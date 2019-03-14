#include "../interface/FiveTracksProducer.h"
#include <tuple>
#include <map>

float FiveTracksProducer::DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
   float p1 = t1.phi();
   float p2 = t2.phi();
   float e1 = t1.eta();
   float e2 = t2.eta();
   auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

bool FiveTracksProducer::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<maxDPtRel &&
	DeltaR(t1,t2) < maxDeltaR);
}

bool
FiveTracksProducer::isAbHadron(int pdgID) {

  if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
  return false;

}

bool
FiveTracksProducer::isAMixedbHadron(int pdgID, int momPdgID) {

  if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) ||
  (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0))
  return true;
  return false;

}

bool
FiveTracksProducer::isTheCandidate(reco::GenParticleRef genY) {

  bool goToJPsi = false;
  bool goToPhi = false;

  for(size_t j = 0; j < genY->numberOfDaughters(); ++j)
  {

    const reco::Candidate * daughter = genY->daughter(j);
    if(daughter->pdgId() == 443)
      goToJPsi=true;
    if(daughter->pdgId() == 333)
    {
      bool kP = false, kN = false;
      for(size_t k = 0; k <daughter->numberOfDaughters(); ++k)
      {
        const reco::Candidate * gdaughter = daughter->daughter(k);

        if(goToPhi && goToJPsi)
        {
          if(gdaughter->pdgId()==321)
            kP=true;
          if(gdaughter->pdgId()==-321)
            kN=true;
        }

      }
      goToPhi = kP && kN;
    }

  }

  return (goToJPsi && goToPhi);

}

bool FiveTracksProducer::isSameTrack(reco::Track t1, reco::Track t2)
{

  float p1 = t1.phi();
  float p2 = t2.phi();
  float e1 = t1.eta();
  float e2 = t2.eta();
  auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

  float deltaR = sqrt((e1-e2)*(e1-e2) + dp*dp);
  float deltaPt = ((t1.pt() - t2.pt())/t1.pt());

  return (deltaR <= 0.001) &&( deltaPt <= 0.01);

}

FiveTracksProducer::FiveTracksProducer(const edm::ParameterSet& iConfig):
  DiMuonDiTrakCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiMuoDiTrak"))),
  TrakCollection_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
  trakPtCut_(iConfig.existsAs<double>("TrakPtCut") ? iConfig.getParameter<double>("TrakPtCut") : 0.8),
  thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
  thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
  TriggerCollection_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  FiveTrakMassCuts_(iConfig.getParameter<std::vector<double>>("FiveTrakCuts")),
  numMasses_(iConfig.getParameter<uint32_t>("NumMasses"))
{
  produces<pat::CompositeCandidateCollection>("FiveTracks");

  nevents = 0;

  kaonmass = 0.493677;
  pionmass = 0.13957061;
  protonmass = 0.93827208;
  psi2smass = 3.686093;

  maxDeltaR = 0.01;
  maxDPtRel = 2.0;
  trackmass = kaonmass;

  ncomboneg = 0;
  ncombo = 0;
  ncomboneu = 0;

  pdgToMass[211] = pionmass;
  pdgToMass[2212] = protonmass;
  pdgToMass[321] = kaonmass;

}

void FiveTracksProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::unique_ptr<pat::CompositeCandidateCollection> fiveCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> dimuonditrak;
  iEvent.getByToken(DiMuonDiTrakCollection_,dimuonditrak);

  edm::Handle<std::vector<pat::PackedCandidate> > trak;
  iEvent.getByToken(TrakCollection_,trak);

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;
  reco::Vertex theBeamSpotV = reco::Vertex(bs.position(), bs.covariance3D());

  edm::Handle<reco::VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trig;
  iEvent.getByToken(TriggerCollection_,trig);

  // const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  // edm::ESHandle<MagneticField> magneticField;
  // iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  // const MagneticField* field = magneticField.product();

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);


  float FiveTrakMassMax = FiveTrakMassCuts_[1];
  float FiveTrakMassMin = FiveTrakMassCuts_[0];

  //Sorting new masses - Kaon is the baseline
  // std::sort( ExtraMasses_.begin(), ExtraMasses_.end() );
  // ExtraMasses_.erase( std::unique( ExtraMasses_.begin(), vec.end() ), ExtraMasses_.end());

  // std::vector<float> extraMasses;
  //
  // for (size_t i = 0; i < ExtraMasses_.size(); i++)
  //   if(pdgToMass.find(ExtraMasses_[i])!=pdgToMass.end())
  //     extraMasses.push_back(pdgToMass[ExtraMasses_[i]]);
  //
  // float deltaMassMin = 0.0, deltaMassMax = 0.0;
  // for(const auto& m : extraMasses)
  // {
  //   deltaMassMin = std::max(deltaMassMin,kaonmass - m);
  //   deltaMassMax = std::max(deltaMassMin,m - kaonmass);
  // }

  FiveTrakMassMax = FiveTrakMassMax + 3*pionmass;
  FiveTrakMassMin = FiveTrakMassMax - 3*pionmass;

  // reco::TrackCollection allTheTracks;
  // for (size_t i = 0; i < trak->size(); i++)
  // {
  //   auto t = trak->at(i);
  //   if(t.pt()<0.5) continue;
  //   if(!(t.hasTrackDetails())) continue;
  //   allTheTracks.push_back(*(t.bestTrack()));
  //
  // }

  KinematicParticleFactoryFromTransientTrack pFactory;

  std::map< std::tuple <int,int,int> ,int> doneFlag;
  std::map<size_t,float> bestVertexPos, bestVertexNeg, bestVertexNeu;
  std::map<size_t,pat::CompositeCandidate> posCollection,negCollection,neuCollection;

  for (size_t d = 0; d < dimuonditrak->size(); d++) {

       auto dimuonditrakCand = dimuonditrak->at(d);

       if(dimuonditrakCand.userFloat("vProb")<0.01)
         continue;

       // const reco::Vertex thePrimaryV = *(dimuonditrakCand.userData<reco::Vertex>("bestPV"));
       // const reco::Vertex thePrimaryV = *dimuonditrakCand.userData<reco::Vertex>("PVwithmuons");
       const pat::CompositeCandidate * dimuon_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrakCand.daughter("dimuon"));

       const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(dimuonditrakCand.daughter("dimuon")->daughter("highMuon"));
       const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(dimuonditrakCand.daughter("dimuon")->daughter("lowMuon"));
       // const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(pmu1->originalObject());
       // const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(pmu2->originalObject());


       //I want the positive and negative track to build psi2S with charge == 0
       const pat::PackedCandidate *tp = dynamic_cast <pat::PackedCandidate *>(dimuonditrakCand.daughter("ditrak")->daughter("highTrak"));
       const pat::PackedCandidate *tm = dynamic_cast <pat::PackedCandidate *>(dimuonditrakCand.daughter("ditrak")->daughter("lowTrak"));
       int tpId = dimuonditrakCand.userInt("pId");
       int tmId = dimuonditrakCand.userInt("mId");

       std::vector<double> oneMasses,twoMasses,threeMasses;
       oneMasses.push_back(pionmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(kaonmass); // p p k
       oneMasses.push_back(kaonmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(pionmass); // k p p
       oneMasses.push_back(pionmass);  twoMasses.push_back(kaonmass);  threeMasses.push_back(pionmass); // p k p
       oneMasses.push_back(pionmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(pionmass); // p p p


       //Adding the fifth track
       //Possibilities:
       // B+  -> Psi' K+ -> JPsi π+ π- K+
       // B-  -> Psi' K- -> JPsi π+ π- K-
       // B+  -> J/Psi π+ π- K+
       // B-  -> J/Psi K+ K- K-
       // B+  -> J/Psi K+ K- K+ (Phi->KK or incoherent)
       // B-  -> J/Psi K+ K- K- (Phi->KK or incoherent)
       // [ 6 tracks: B0s -> Psi' Phi -> J/Psi π π K K
       // B0 -> J/Psi Phi K0 -> J/Psi K K K0
       //
       
       const unsigned int numMasses = 4; //(int) numMasses_;

       for (size_t i = 0; i < trak->size(); i++) {

         std::array<float,numMasses> fiveTracksMass;

         float minDR_third = 10000.0;
         float minDP_third = 10000.0;
         float minDE_third = 10000.0;

         auto fifthTrack = trak->at(i);


         if(fifthTrack.pt()<trakPtCut_) continue;
         if(fifthTrack.charge() == 0) continue;
	       //if(!isMC_ and fabs(fifthTrack.pdgId())!=211) continue;
	       if(!(fifthTrack.trackHighPurity())) continue;
         if(!(fifthTrack.hasTrackDetails())) continue;

         if (IsTheSame(fifthTrack,*tp) || int(i) == tpId) continue;
         if (IsTheSame(fifthTrack,*tm) || int(i) == tmId) continue;
         if (IsTheSame(fifthTrack,*pmu1) || IsTheSame(fifthTrack,*pmu2) ) continue;

         std::vector<int> ids;
         ids.push_back(i);ids.push_back(tpId);ids.push_back(tmId);
         std::sort(ids.begin(),ids.end());
         if(doneFlag.find(std::tuple<int,int,int>(ids[0],ids[1],ids[2]))!=doneFlag.end())
          continue;
         else
            doneFlag[std::tuple<int,int,int>(ids[0],ids[1],ids[2])] = 1.0;

         pat::CompositeCandidate fiveCand = makeFiveCandidate(dimuonditrakCand, fifthTrack);

         double five_ma_fit = 14000.;
         double five_vp_fit = -9999.;
         double five_x2_fit = 10000.;
         double five_nd_fit = 10000.;

         std::vector<pat::CompositeCandidate> fiveCands, ref_fiveCands;

         bool atLeastOne = false;
         std::vector< bool > insideMass;

         // int numMasses = (int) extraMasses.size();

         std::vector<reco::TransientTrack> fiveTracks;
         std::vector<RefCountedKinematicParticle> kaonParticles, pionParticles;

         float kinChi = 0.;
         float kinNdf = 0.;

         fiveTracks.push_back((*theB).build(*(pmu1->innerTrack()))); // µ
         fiveTracks.push_back((*theB).build(*(pmu2->innerTrack()))); // µ
         fiveTracks.push_back((*theB).build(*(tp->bestTrack()))); // K/π
         fiveTracks.push_back((*theB).build(*(tm->bestTrack()))); // K/π
         fiveTracks.push_back((*theB).build(*(fifthTrack.bestTrack()))); // K/π

         kaonParticles.push_back(pFactory.particle(fiveTracks[0],muonMass,kinChi,kinNdf,muonSigma));
         kaonParticles.push_back(pFactory.particle(fiveTracks[1],muonMass,kinChi,kinNdf,muonSigma));
         kaonParticles.push_back(pFactory.particle(fiveTracks[2],trakMassP,kinChi,kinNdf,trakSigmaP));
         kaonParticles.push_back(pFactory.particle(fiveTracks[3],trakMassM,kinChi,kinNdf,trakSigmaM));
         kaonParticles.push_back(pFactory.particle(fiveTracks[4],fifthMass,kinChi,kinNdf,fifthSigma));

         KinematicConstrainedVertexFitter fiveFitter;
         MultiTrackKinematicConstraint *jpsi_mtc = new  TwoTrackMassKinematicConstraint(jpsiMass);
         RefCountedKinematicTree fiveVertexFitTree;
         fiveVertexFitTree = fiveFitter.fit(kaonParticles,jpsi_mtc);

         if (fiveVertexFitTree->isEmpty()) continue;

         fiveVertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle fitF = fiveVertexFitTree->currentParticle();
         RefCountedKinematicVertex fitFVertex = fiveVertexFitTree->currentDecayVertex();

         if (!(fitF->currentState().isValid())) continue;

         hasRefit[j] = 1.0;

         five_ma_fit = fitF->currentState().mass();
         five_x2_fit = fitFVertex->chiSquared();
         five_nd_fit = (double)(fitFVertex->degreesOfFreedom());
         five_vp_fit = ChiSquaredProbability(five_x2_fit,five_nd_fit);


         if(five_vp_fit < 0.001) continue;

         float minDR_pos = 10000.0, minDR_neg = 10000.0;
         float minDP_pos = 10000.0, minDP_neg = 10000.0;
         float minDE_pos = 10000.0, minDE_neg = 10000.0;

         for (size_t ii = 0; ii < muons->size(); ii++)
         {
            auto thisMuon = muons->at(ii);

            float DeltaEta = fabs(thisMuon.eta()-fifthTrack.eta());
            float DeltaP   = fabs(thisMuon.p()-fifthTrack.p());
            float DeltaPt = ((fifthTrack.pt() - thisMuon.pt())/fifthTrack.pt());

            minDR_pos = -std::max(minDR_pos,DeltaEta);
            minDP_pos = -std::max(minDP_pos,DeltaP);
            minDE_pos = -std::max(minDE_pos,DeltaPt);
         }


         int    five_ch_fit = fiveCand.charge();
         double five_px_fit = fitF->currentState().kinematicParameters().momentum().x();
         double five_py_fit = fitF->currentState().kinematicParameters().momentum().y();
         double five_pz_fit = fitF->currentState().kinematicParameters().momentum().z();
         double five_en_fit = sqrt(five_ma_fit*five_ma_fit+five_px_fit*five_px_fit+five_py_fit*five_py_fit+five_pz_fit*five_pz_fit);
         double five_vx_fit = fitFVertex->position().x();
         double five_vy_fit = fitFVertex->position().y();
         double five_vz_fit = fitFVertex->position().z();

         reco::CompositeCandidate recoFive_rf(five_ch_fit,math::XYZTLorentzVector(five_px_fit,five_py_fit,five_pz_fit,five_en_fit),
                                           math::XYZPoint(five_vx_fit,five_vy_fit,five_vz_fit),0);

         pat::CompositeCandidate fiveCand_rf(recoFive_rf);

         //////////////////////////////////////////////////////////////////////////////
         //PV Selection(s)
         int numVertex = 4;
         std::array <float,numVertex> sumPTPV, cosAlpha, ctauPV, ctauErrPV, fromPV;
         // std::cout << "debug    9 "<< std::endl;
         TVector3 vtx, vdiff, pvtx;
         VertexDistanceXY vdistXY;
         reco::Vertex thePrimaryV,thePrimaryVDZ, thePrimaryZero, thePrimaryVCA;
         TwoTrackMinimumDistance ttmd;

         double x_px_fit = fitX->currentState().kinematicParameters().momentum().x();
         double x_py_fit = fitX->currentState().kinematicParameters().momentum().y();
         double x_pz_fit = fitX->currentState().kinematicParameters().momentum().z();
         // double x_en_fit = sqrt(x_ma_fit*x_ma_fit+x_px_fit*x_px_fit+x_py_fit*x_py_fit+x_pz_fit*x_pz_fit);
         double x_vx_fit = fitXVertex->position().x();
         double x_vy_fit = fitXVertex->position().y();
         double x_vz_fit = fitXVertex->position().z();
         vtx.SetXYZ(x_vx_fit,x_vy_fit,0);

         bool status = ttmd.calculate( GlobalTrajectoryParameters(
           GlobalPoint(x_vx_fit,x_vy_fit,x_vz_fit),
           GlobalVector(x_px_fit,x_py_fit,x_pz_fit),TrackCharge(0),&(*magneticField)),
           GlobalTrajectoryParameters(
             GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
             GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(0),&(*magneticField)));
         float extrapZ=-9E20;

         if (status) extrapZ=ttmd.points().first.z();

         TVector3 pperp(x_px_fit, x_py_fit, 0);
         AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);

         reco::VertexCollection verteces;
         std::vector<int> vKeys;
         verteces.push_back(theBeamSpotV);

         thePrimaryZero = reco::Vertex(*(priVtxs->begin()));
         verteces.push_back(thePrimaryV);
         vKeys.push_back(0);

         float minDz = 999999.;
         double maxCosAlpha = -1.0;
         if ( (priVtxs->begin() == priVtxs->end()) )
         {
           // std::cout << "debug    10 "<< std::endl;
           thePrimaryVCA = reco::Vertex(*(priVtxs->begin()));
           thePrimaryVDZ = reco::Vertex(*(priVtxs->begin()));
           verteces.push_back(thePrimaryV);
           verteces.push_back(thePrimaryVDZ);
           vKeys.push_back(0);
           vKeys.push_back(0);
         }else
         {

           reco::Vertex p,pz;
           int thisp,thispz;
           for(size_t pV = 0; pV<priVtxs->size();++pV)
           {
             auto thisPV = priVtxs->at(pV);

             pvtx.SetXYZ(thePrimaryVCA.position().x(),thePrimaryVCA.position().y(),0);
             vdiff = vtx - pvtx;
             double thisCosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
             if(thisCosAlpha>maxCosAlpha)
             {
               thePrimaryVCA = reco::Vertex(thisPV);
               maxCosAlpha = thisCosAlpha;
               p = reco::Vertex(thisPV);
               thisp = pV;
             }

             float deltaZ = fabs(extrapZ - thisPV.position().z()) ;
             if ( deltaZ < minDz ) {
               minDz = deltaZ;
               pz = reco::Vertex(thisPV);
               thispz =pV;
             }
           }
           verteces.push_back(p);
           vKeys.push_back(thisp);
           verteces.push_back(pz);
           vKeys.push_back(thispz);
         }


         vtx.SetXYZ(five_vx_fit,five_vy_fit,0);
         TVector3 pperp(five_px_fit, five_py_fit, 0);
         AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);

         //////////////////////////////////////////////////
         //Refit PVs (not BS)

         for(size_t i = 0; i < numVertex; i++)
         {

           ctauPV[i] = (-1000.0);
           ctauErrPV[i] = (-1000.0);
           cosAlpha[i] = (-1000.0);
           fromPV[i] = (-1000.0);
         }

         // std::cout << "debug    13 "<< std::endl;
         for(size_t i = 0; i < verteces.size(); i++)
         {
           pvtx.SetXYZ(verteces[i].position().x(),verteces[i].position().y(),0);
           vdiff = vtx - pvtx;
           cosAlpha.push_back(vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp()));
           Measurement1D distXY = vdistXY.distance(reco::Vertex(*fitXVertex), verteces[i]);
           ctauPV.push_back(distXY.value()*cosAlpha[i] * x_ma_fit/pperp.Perp());
           GlobalError v1e = (reco::Vertex(*fitXVertex)).error();
           GlobalError v2e = verteces[i].error();
           AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
           ctauErrPV.push_back(sqrt(ROOT::Math::Similarity(vpperp,vXYe))*x_ma_fit/(pperp.Perp2()));
         }


         fiveCand.addUserFloat("mass_ref_0",five_ma_fit);
         fiveCand.addDaughter(fiveCand_rf,"ref_cand");

         fiveCand.addUserFloat("thirdTrackMuonDR",minDR_pos);
         fiveCand.addUserFloat("thirdTrackMuonDP",minDP_pos);
         fiveCand.addUserFloat("thirdTrackMuonDE",minDE_pos);

         pat::CompositeCandidate thisFive;

         for(size_t j = 0; j<numMasses_;j++)
          fiveTracksMass[j] = makeFiveCandidateMixed(*dimuon_cand, *tp, *tm, fifthTrack,oneMasses[j] ,twoMasses[j] ,threeMasses[j]).mass();

             fiveCand.addUserInt("dimuontt_index",int(d));
             fiveCand.addUserInt("pId",tpId);
             fiveCand.addUserInt("mId",tmId);
             fiveCand.addUserInt("tId",i);

             fiveCand.addUserFloat(name,fiveTracksMass[j]);

             std::string name;
             for(size_t j = 1; j<numMasses_+1;j++)
             {
              name = "mass_ref_" + std::to_string(j);
              fiveCand.addUserFloat(name,fiveTracksMass[j-1]);
             }

             // fiveCand.addDaughter(dimuonditrakCand,"dimuonditrak");
             fiveCand.addDaughter(*tp,"trakOne");
             fiveCand.addDaughter(*tm,"trakTwo");
             fiveCand.addDaughter(fifthTrack,"trakThree");


             fiveCand.addUserData("bestPV",reco::Vertex(thePrimaryZero));
             fiveCand.addUserData("cosPV",reco::Vertex(thePrimaryVCA));
             fiveCand.addUserData("zPV",reco::Vertex(thePrimaryVDZ));
             fiveCand.addUserData("bs",reco::Vertex(thePrimaryV));

             fiveCand.addUserFloat("vtxX",x_vx_fit);
             fiveCand.addUserFloat("vtxY",x_vy_fit);
             fiveCand.addUserFloat("vtxZ",x_vz_fit);


             fiveCand.addUserFloat("cosAlphaBS",cosAlpha[0]);
             fiveCand.addUserFloat("ctauPVBS",ctauPV[0]);
             fiveCand.addUserFloat("ctauErrPVBS",ctauErrPV[0]);

             fiveCand.addUserFloat("ttFromPVBS",float(fromPV[0]));


             fiveCand.addUserFloat("cosAlpha",cosAlpha[1]);
             fiveCand.addUserFloat("ctauPV",ctauPV[1]);
             fiveCand.addUserFloat("ctauErrPV",ctauErrPV[1]);

             fiveCand.addUserFloat("ttFromPV",float(fromPV[1]));

             fiveCand.addUserFloat("cosAlpha_alpha",cosAlpha[2]);
             fiveCand.addUserFloat("ctauPV_alpha",ctauPV[2]);
             fiveCand.addUserFloat("ctauErrPV_alpha",ctauErrPV[2]);

             fiveCand.addUserFloat("ttFromPV_alpha",float(fromPV[2]));

             fiveCand.addUserFloat("cosAlphaDZ",cosAlpha[3]);
             fiveCand.addUserFloat("ctauPVDZ",ctauPV[3]);
             fiveCand.addUserFloat("ctauErrPVDZ",ctauErrPV[3]);

             fiveCand.addUserFloat("ttFromPVDZ",float(fromPV[3]));
             ///DCA
             std::vector<float> DCAs;

             for(size_t j = 0; j < fiveTracks.size() - 1;++j)
             {
               TrajectoryStateClosestToPoint TS1 = fiveTracks[fiveTracks.size()-1].impactPointTSCP();
               TrajectoryStateClosestToPoint TS2 = fiveTracks[j].impactPointTSCP();
               float dca = 1E20;
               if (TS1.isValid() && TS2.isValid()) {
                 ClosestApproachInRPhi cApp;
                 cApp.calculate(TS1.theState(), TS2.theState());
                 if (cApp.status() ) dca = cApp.distance();
               }
               DCAs.push_back(dca);
             }


             fiveCand.addUserFloat("vProb",fiveTracksVProb[0]);
             fiveCand.addUserFloat("ctauPV",fiveTracksCTau[0]);
             fiveCand.addUserFloat("ctauErrPV",fiveTracksCTauErr[0]);
             fiveCand.addUserFloat("cosAlpha",fiveTracksCosAlpha[0]);
             fiveCand.addUserFloat("nDof",fiveTracksVNDof[0]);
             fiveCand.addUserFloat("vChi2",fiveTracksVChi2[0]);

             fiveCand.addUserFloat("dca_m1t3",DCAs[0]);
             fiveCand.addUserFloat("dca_m2t3",DCAs[1]);
             fiveCand.addUserFloat("dca_t1t3",DCAs[2]);
             fiveCand.addUserFloat("dca_t2t3",DCAs[3]);

             // std::cout << std::endl;
             if(atLeastOne)
              fiveCandColl->push_back(fiveCand);

           }

        }


     ncombo += fiveCandColl->size();

  iEvent.put(std::move(fiveCandColl),"FiveTracks");

  nevents++;
}

void FiveTracksProducer::endJob(){
  std::cout << "#########################################" << std::endl;
  std::cout << "FiveTracks Candidate producer report:" << std::endl;
  std::cout << "#########################################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "No. dimtt + trk candidates " << ncombo << std::endl;
  std::cout << "#########################################" << std::endl;
}

bool FiveTracksProducer::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.01 && DeltaP < 0.01) return true;
  return false;
}

bool FiveTracksProducer::IsTheSame(const pat::PackedCandidate& t1, const pat::PackedCandidate& t2){
  double DeltaEta = fabs(t1.eta()-t2.eta());
  double DeltaP   = fabs(t1.p()-t2.p());
  if (DeltaEta < 0.01 && DeltaP < 0.01) return true;
  return false;
}


pat::CompositeCandidate FiveTracksProducer::makeFiveCandidate(
                                          const pat::CompositeCandidate& dimuonditrak,
                                          const pat::PackedCandidate& trak
                                         ){

  pat::CompositeCandidate fiveCand, dimuontrak;
  fiveCand.addDaughter(dimuonditrak,"dimuonditrak");
  fiveCand.addDaughter(trak,"fifth");
  fiveCand.setCharge(dimuonditrak.charge()+trak.charge());

  double m_trak = trackmass;
  math::XYZVector mom_trak = trak.momentum();
  double e_trak = sqrt(m_trak*m_trak + mom_trak.Mag2());
  math::XYZTLorentzVector p4_trak = math::XYZTLorentzVector(mom_trak.X(),mom_trak.Y(),mom_trak.Z(),e_trak);

  reco::Candidate::LorentzVector v = p4_trak + dimuonditrak.p4();

  fiveCand.setP4(v);

  return fiveCand;
}

pat::CompositeCandidate FiveTracksProducer::makeFiveCandidateMixed(
                                          const pat::CompositeCandidate& dimuon,
                                          const pat::PackedCandidate& trakP,
                                          const pat::PackedCandidate& trakN,
                                          const pat::PackedCandidate& trak3,
                                          double massOne,
                                          double massTwo,
                                          double massThree
                                         ){

  pat::CompositeCandidate fiveCand, trakOne, trakTwo, trakThree;
  pat::CompositeCandidate dimuonDiTrakOne, dimuonDiTrakTwo, dimuonDiTrakThree;
  pat::CompositeCandidate triTrak;

  fiveCand.addDaughter(dimuon,"dimuon");

  triTrak.setCharge(trakP.charge()+trakN.charge()+trak3.charge());
  fiveCand.setCharge(dimuon.charge()+trakP.charge()+trakN.charge()+trak3.charge());

  math::XYZVector mom_trakP = trakP.momentum();
  double e_trakP = sqrt(massOne*massOne + mom_trakP.Mag2());
  math::XYZTLorentzVector p4_trakP = math::XYZTLorentzVector(mom_trakP.X(),mom_trakP.Y(),mom_trakP.Z(),e_trakP);
  trakOne.setCharge(trakP.charge());
  trakOne.setP4(p4_trakP);

  math::XYZVector mom_trakN = trakN.momentum();
  double e_trakN = sqrt(massTwo*massTwo + mom_trakN.Mag2());
  math::XYZTLorentzVector p4_trakN = math::XYZTLorentzVector(mom_trakN.X(),mom_trakN.Y(),mom_trakN.Z(),e_trakN);
  trakTwo.setCharge(trakN.charge());
  trakTwo.setP4(p4_trakN);

  math::XYZVector mom_trak3 = trak3.momentum();
  double e_trak3 = sqrt(massThree*massThree + mom_trak3.Mag2());
  math::XYZTLorentzVector p4_trak3 = math::XYZTLorentzVector(mom_trak3.X(),mom_trak3.Y(),mom_trak3.Z(),e_trak3);
  trakThree.setCharge(trak3.charge());
  trakThree.setP4(p4_trak3);

  fiveCand.addDaughter(trakOne,"trakOne");
  fiveCand.addDaughter(trakTwo,"trakTwo");
  fiveCand.addDaughter(trakThree,"trakThree");

  dimuonDiTrakOne     = makePsi2SCandidate(dimuon,trakOne,trakTwo);
  dimuonDiTrakTwo     = makePsi2SCandidate(dimuon,trakOne,trakThree);
  dimuonDiTrakThree   = makePsi2SCandidate(dimuon,trakTwo,trakThree);

  fiveCand.addDaughter(dimuonDiTrakOne,"dimuonDiTrakOne");
  fiveCand.addDaughter(dimuonDiTrakTwo,"dimuonDiTrakTwo");
  fiveCand.addDaughter(dimuonDiTrakThree,"dimuonDiTrakThree");

  reco::Candidate::LorentzVector v = p4_trakP + p4_trak3 + p4_trak3 + dimuon.p4();
  reco::Candidate::LorentzVector vT = trakP.p4() + trakN.p4() + trak3.p4();

  triTrak.setP4(vT);

  fiveCand.addDaughter(triTrak,"triTrak");

  fiveCand.setP4(v);

  return fiveCand;
}

pat::CompositeCandidate FiveTracksProducer::makeFiveCandidateMixed(
                                          const pat::CompositeCandidate& dimuon,
                                          const pat::CompositeCandidate& trakP,
                                          const pat::CompositeCandidate& trakN,
                                          const pat::CompositeCandidate& trak3
                                         ){

  pat::CompositeCandidate fiveCand, trakOne, trakTwo, trakThree;
  pat::CompositeCandidate dimuonDiTrakOne, dimuonDiTrakTwo, dimuonDiTrakThree;
  pat::CompositeCandidate triTrak;

  fiveCand.addDaughter(dimuon,"dimuon");

  fiveCand.setCharge(dimuon.charge()+trakP.charge()+trakN.charge()+trak3.charge());
  triTrak.setCharge(trakP.charge()+trakN.charge()+trak3.charge());

  fiveCand.addDaughter(trakP,"trakOne");
  fiveCand.addDaughter(trakN,"trakTwo");
  fiveCand.addDaughter(trak3,"trakThree");

  dimuonDiTrakOne     = makePsi2SCandidate(dimuon,trakP,trakN);
  dimuonDiTrakTwo     = makePsi2SCandidate(dimuon,trakP,trak3);
  dimuonDiTrakThree   = makePsi2SCandidate(dimuon,trakN,trak3);

  fiveCand.addDaughter(dimuonDiTrakOne,"dimuonDiTrakOne");
  fiveCand.addDaughter(dimuonDiTrakTwo,"dimuonDiTrakTwo");
  fiveCand.addDaughter(dimuonDiTrakThree,"dimuonDiTrakThree");

  reco::Candidate::LorentzVector v  = trakP.p4() + trakN.p4() + trak3.p4() + dimuon.p4();
  reco::Candidate::LorentzVector vT = trakP.p4() + trakN.p4() + trak3.p4();

  triTrak.setP4(vT);

  fiveCand.addDaughter(triTrak,"triTrak");

  fiveCand.setP4(v);

  return fiveCand;
}

pat::CompositeCandidate FiveTracksProducer::makePsi2SCandidate(
                                          const pat::CompositeCandidate& dimuon,
                                          const pat::CompositeCandidate& t1,
                                          const pat::CompositeCandidate& t2
                                         ){

  pat::CompositeCandidate psi2sCand, ditrak;
  psi2sCand.setCharge(dimuon.charge()+t1.charge()+t2.charge());
  ditrak.setCharge(t1.charge()+t2.charge());
  psi2sCand.addDaughter(dimuon,"dimuon");
  psi2sCand.addDaughter(t1,"trakOne");
  psi2sCand.addDaughter(t2,"trakTwo");

  reco::Candidate::LorentzVector v  = t1.p4() + t2.p4() + dimuon.p4();
  reco::Candidate::LorentzVector vT = t1.p4() + t2.p4();

  ditrak.setP4(vT);

  psi2sCand.addDaughter(ditrak,"ditrak");

  psi2sCand.setP4(v);

  return psi2sCand;
}

reco::Candidate::LorentzVector FiveTracksProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(FiveTracksProducer);
