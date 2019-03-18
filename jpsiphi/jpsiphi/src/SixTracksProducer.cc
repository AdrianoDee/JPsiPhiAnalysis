#include "../interface/SixTracksProducer.h"
#include <tuple>
#include <map>

float SixTracksProducer::DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
   float p1 = t1.phi();
   float p2 = t2.phi();
   float e1 = t1.eta();
   float e2 = t2.eta();
   auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

bool SixTracksProducer::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<maxDPtRel &&
	DeltaR(t1,t2) < maxDeltaR);
}

bool
SixTracksProducer::isAbHadron(int pdgID) {

  if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
  return false;

}

bool
SixTracksProducer::isAMixedbHadron(int pdgID, int momPdgID) {

  if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) ||
  (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0))
  return true;
  return false;

}

bool
SixTracksProducer::isTheCandidate(reco::GenParticleRef genY) {

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

bool SixTracksProducer::isSameTrack(reco::Track t1, reco::Track t2)
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

SixTracksProducer::SixTracksProducer(const edm::ParameterSet& iConfig):
  FiveTrakCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("FiveCollection"))),
  TrakCollection_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
  trakPtCut_(iConfig.existsAs<double>("TrakPtCut") ? iConfig.getParameter<double>("TrakPtCut") : 0.8),
  thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
  thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
  TriggerCollection_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  SixTrakMassCuts_(iConfig.getParameter<std::vector<double>>("SixTrakCuts")),
  numMasses_(iConfig.getParameter<uint32_t>("NumMasses"))
{
  produces<pat::CompositeCandidateCollection>("SixTracks");

  nevents = 0;

  kaonmass = 0.493677;
  pionmass = 0.13957061;
  psi2smass = 3.686093;

  maxDeltaR = 0.01;
  maxDPtRel = 2.0;
  trackmass = kaonmass;

  ncomboneg = 0;
  ncombo = 0;
  ncomboneu = 0;

  allMuons_ = consumes<pat::MuonCollection>((edm::InputTag)"slimmedMuons");

}

void SixTracksProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::unique_ptr<pat::CompositeCandidateCollection> sixCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(allMuons_,muons);

  edm::Handle<pat::CompositeCandidateCollection> fivetrak;
  iEvent.getByToken(FiveTrakCollection_,fivetrak);

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

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  const MagneticField* field = magneticField.product();

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);


  float SixTrakMassMax = SixTrakMassCuts_[1];
  float SixTrakMassMin = SixTrakMassCuts_[0];


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

  std::map< std::tuple <int,int,int,int> ,bool> doneFlag;
  std::map<size_t,float> bestVertexPos, bestVertexNeg, bestVertexNeu;
  std::map<size_t,pat::CompositeCandidate> posCollection,negCollection,neuCollection;

  const int numMasses = 6;//numMasses_;

  for (size_t d = 0; d < fivetrak->size(); d++) {

       auto fivetrakCand = fivetrak->at(d);

       if(fivetrakCand.userFloat("vProb")<0.001)
         continue;
       // const reco::Vertex thePrimaryV = *(fivetrakCand.userData<reco::Vertex>("bestPV"));
       // const reco::Vertex thePrimaryV = *fivetrakCand.userData<reco::Vertex>("PVwithmuons");
       const pat::CompositeCandidate *dimuonditrak_cand = dynamic_cast <const pat::CompositeCandidate *>(fivetrakCand.daughter("dimuonditrak"));
       const pat::CompositeCandidate *dimuon_cand = dynamic_cast <const pat::CompositeCandidate *>(dimuonditrak_cand->daughter("dimuon"));

       const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(dimuonditrak_cand->daughter("dimuon")->daughter("highMuon"));
       const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(dimuonditrak_cand->daughter("dimuon")->daughter("lowMuon"));
       // const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(pmu1->originalObject());
       // const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(pmu2->originalObject());


       //I want the positive and negative track to build psi2S with charge == 0
       const pat::PackedCandidate *tp = dynamic_cast <const pat::PackedCandidate *>(dimuonditrak_cand->daughter("ditrak")->daughter("highTrak"));
       const pat::PackedCandidate *tm = dynamic_cast <const pat::PackedCandidate *>(dimuonditrak_cand->daughter("ditrak")->daughter("lowTrak"));
       const pat::PackedCandidate *tt = dynamic_cast <const pat::PackedCandidate *>(fivetrakCand.daughter("fifth"));

       int tpId = fivetrakCand.userInt("pId");
       int tmId = fivetrakCand.userInt("mId");
       int ttId = fivetrakCand.userInt("tId");

       std::array< std::array <float,4>, numMasses-1> trackMasses;
       trackMasses[0][0] = kaonmass; trackMasses[0][1] = pionmass; trackMasses[0][2] = kaonmass; trackMasses[0][3] = pionmass; // k p k p
       trackMasses[1][0] = kaonmass; trackMasses[1][1] = pionmass; trackMasses[1][2] = pionmass; trackMasses[1][3] = kaonmass; // k p k p
       trackMasses[2][0] = pionmass; trackMasses[2][1] = kaonmass; trackMasses[2][2] = kaonmass; trackMasses[2][3] = pionmass; // k p k p
       trackMasses[3][0] = pionmass; trackMasses[3][1] = kaonmass; trackMasses[3][2] = pionmass; trackMasses[3][3] = kaonmass; // k p k p
       trackMasses[4][0] = pionmass; trackMasses[4][1] = pionmass; trackMasses[4][2] = kaonmass; trackMasses[4][3] = kaonmass; // k p k p


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

       for (size_t i = 0; i < trak->size(); i++) {

         std::vector<float> sixTracksMass;
         std::vector<float> psi2sOne, psi2sTwo;
         std::vector<float> fiveTracksCTau, fiveTracksCTauErr, fiveTracksCosAlpha;

         float minDR_fourth = 10000.0;
         float minDP_fourth = 10000.0;
         float minDPt_fourth = 10000.0;

         auto sixthTrack = trak->at(i);

         if(sixthTrack.pt()<trakPtCut_) continue;
         //if(sixthTrack.charge() == 0) continue;
	       //if(!isMC_ and fabs(sixthTrack.pdgId())!=211) continue;
	       if(!(sixthTrack.trackHighPurity())) continue;
         if(!(sixthTrack.hasTrackDetails())) continue;

         if (IsTheSame(sixthTrack,*tp) || int(i) == tpId) continue;
         if (IsTheSame(sixthTrack,*tm) || int(i) == tmId) continue;
         if (IsTheSame(sixthTrack,*tm) || int(i) == ttId) continue;
         if (IsTheSame(sixthTrack,*pmu1) || IsTheSame(sixthTrack,*pmu2) ) continue;

         std::vector<int> ids;
         ids.push_back(i);ids.push_back(tpId);ids.push_back(tmId);ids.push_back(ttId);
         std::sort(ids.begin(),ids.end());
         if(doneFlag.find(std::tuple<int,int,int,int>(ids[0],ids[1],ids[2],ids[3]))!=doneFlag.end())
          continue;
         else
            doneFlag[std::tuple<int,int,int,int>(ids[0],ids[1],ids[2],ids[3])] = true;

         pat::CompositeCandidate sixCand = makeSixCandidate(fivetrakCand, sixthTrack);

         if (sixCand.mass() > SixTrakMassMin || sixCand.mass() < SixTrakMassMax)
         continue;

         sixCand.addUserFloat("sixCandMass",sixCand.mass());

         double six_ma_fit = 14000.;
         double six_vp_fit = -9999.;
         double six_x2_fit = 10000.;
         double six_nd_fit = 10000.;

         std::vector<pat::CompositeCandidate> fiveCands, ref_fiveCands;

         bool atLeastOne = false;
         std::vector< bool > insideMass;

         const ParticleMass muonMass(0.1056583);
         float muonSigma = muonMass*1E-6;
         const ParticleMass kaonMass(kaonmass);
         float kaonSigma = kaonMass*1E-6;
         const ParticleMass pionMass(pionmass);
         float pionSigma = pionMass*1E-6;

         std::vector<reco::TransientTrack> sixTracks;
         std::vector<RefCountedKinematicParticle> dParticles;

         float kinChi = 0.;
         float kinNdf = 0.;

         sixTracks.push_back((*theB).build(*(pmu1->innerTrack()))); // µ
         sixTracks.push_back((*theB).build(*(pmu2->innerTrack()))); // µ
         sixTracks.push_back((*theB).build(*(tp->bestTrack()))); // K/π
         sixTracks.push_back((*theB).build(*(tm->bestTrack()))); // K/π
         sixTracks.push_back((*theB).build(*(tt->bestTrack()))); // K/π
         sixTracks.push_back((*theB).build(*(sixthTrack.bestTrack()))); // K/π

         dParticles.push_back(pFactory.particle(sixTracks[0],muonMass,kinChi,kinNdf,muonSigma));
         dParticles.push_back(pFactory.particle(sixTracks[1],muonMass,kinChi,kinNdf,muonSigma));
         dParticles.push_back(pFactory.particle(sixTracks[2],kaonmass,kinChi,kinNdf,kaonSigma));
         dParticles.push_back(pFactory.particle(sixTracks[3],kaonmass,kinChi,kinNdf,kaonSigma));
         dParticles.push_back(pFactory.particle(sixTracks[4],pionmass,kinChi,kinNdf,pionSigma));
         dParticles.push_back(pFactory.particle(sixTracks[5],pionmass,kinChi,kinNdf,pionSigma));


         //Mass Constrained fit
         double jpsiMass = 3.096916;

         KinematicConstrainedVertexFitter sixFitter;
         MultiTrackKinematicConstraint *jpsi_mtc = new  TwoTrackMassKinematicConstraint(jpsiMass);
         RefCountedKinematicTree sixVertexFitTree;

         sixVertexFitTree = sixFitter.fit(dParticles,jpsi_mtc);

         if (sixVertexFitTree->isEmpty()) continue;

         sixVertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle fitS = sixVertexFitTree->currentParticle();
         RefCountedKinematicVertex fitSVertex = sixVertexFitTree->currentDecayVertex();

         if (!(fitS->currentState().isValid())) continue;

         six_ma_fit = fitS->currentState().mass();
         six_x2_fit = fitSVertex->chiSquared();
         six_nd_fit = (double)(fitSVertex->degreesOfFreedom());
         six_vp_fit = ChiSquaredProbability(six_x2_fit,six_nd_fit);

         if(six_vp_fit < 0.0005) continue;

         for (size_t ii = 0; ii < muons->size(); ii++)
         {
            auto thisMuon = muons->at(ii);

            float DeltaEta = fabs(thisMuon.eta()-sixthTrack.eta());
            float DeltaP   = fabs(thisMuon.p()-sixthTrack.p());
            float DeltaPt = ((sixthTrack.pt() - thisMuon.pt())/sixthTrack.pt());

            minDR_fourth = -std::max(-minDR_fourth,-DeltaEta);
            minDP_fourth = -std::max(-minDP_fourth,-DeltaP);
            minDPt_fourth = -std::max(-minDPt_fourth,-DeltaPt);
         }


         int    six_ch_fit = sixCand.charge();
         double six_px_fit = fitS->currentState().kinematicParameters().momentum().x();
         double six_py_fit = fitS->currentState().kinematicParameters().momentum().y();
         double six_pz_fit = fitS->currentState().kinematicParameters().momentum().z();
         double six_en_fit = sqrt(six_ma_fit*six_ma_fit+six_px_fit*six_px_fit+six_py_fit*six_py_fit+six_pz_fit*six_pz_fit);
         double six_vx_fit = fitSVertex->position().x();
         double six_vy_fit = fitSVertex->position().y();
         double six_vz_fit = fitSVertex->position().z();

         reco::CompositeCandidate recoSIx_rf(six_ch_fit,math::XYZTLorentzVector(six_px_fit,six_py_fit,six_pz_fit,six_en_fit),
                                           math::XYZPoint(six_vx_fit,six_vy_fit,six_vz_fit),0);

         pat::CompositeCandidate sixCand_rf(recoSIx_rf);

         ////////////////
         //Vertexing
         TVector3 vtx, vdiff, pvtx;
         VertexDistanceXY vdistXY;
         vtx.SetXYZ(six_vx_fit,six_vy_fit,0);
         TVector3 pperp(six_px_fit, six_py_fit, 0);
         AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
         reco::Vertex thePrimaryV,thePrimaryVDZ, thePrimaryZero, thePrimaryVCA;
         TwoTrackMinimumDistance ttmd;

         reco::VertexCollection verteces;
         std::vector<int> vKeys;
         verteces.push_back(theBeamSpotV);

         bool status = ttmd.calculate( GlobalTrajectoryParameters(
           GlobalPoint(six_vx_fit,six_vy_fit,six_vz_fit),
           GlobalVector(six_px_fit,six_py_fit,six_pz_fit),TrackCharge(0),&(*magneticField)),
           GlobalTrajectoryParameters(
             GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
             GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(0),&(*magneticField)));
         float extrapZ=-9E20;

         if (status) extrapZ=ttmd.points().first.z();

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


         std::vector <double> sumPTPV,cosAlpha,ctauPV,ctauErrPV;




         for(size_t i = 0; i < verteces.size(); i++)
         {
           pvtx.SetXYZ(verteces[i].position().x(),verteces[i].position().y(),0);
           vdiff = vtx - pvtx;
           cosAlpha[i] = (vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp()));
           Measurement1D distXY = vdistXY.distance(reco::Vertex(*fitSVertex), verteces[i]);
           ctauPV[i] = (distXY.value()*cosAlpha[i] * six_ma_fit/pperp.Perp());
           GlobalError v1e = (reco::Vertex(*fitSVertex)).error();
           GlobalError v2e = verteces[i].error();
           AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
           ctauErrPV[i] = (sqrt(ROOT::Math::Similarity(vpperp,vXYe))*six_ma_fit/(pperp.Perp2()));
         }

         std::vector<double> oneMasses, twoMasses, threeMasses, fourMasses;
         oneMasses.push_back(pionmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(kaonmass); fourMasses.push_back(kaonmass); // p p k k

         oneMasses.push_back(pionmass);  twoMasses.push_back(kaonmass);  threeMasses.push_back(pionmass); fourMasses.push_back(kaonmass); // p k p k
         oneMasses.push_back(pionmass);  twoMasses.push_back(kaonmass);  threeMasses.push_back(kaonmass); fourMasses.push_back(pionmass); // p k k p

         oneMasses.push_back(kaonmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(kaonmass); fourMasses.push_back(pionmass); // k p k p
         oneMasses.push_back(kaonmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(pionmass); fourMasses.push_back(kaonmass); // k p p k

         oneMasses.push_back(kaonmass);  twoMasses.push_back(kaonmass);  threeMasses.push_back(pionmass); fourMasses.push_back(pionmass); // k k p p


         for(size_t j = 0; j<numMasses_;j++)
          sixTracksMass[j] = makeSixCandidateMixed(*dimuon_cand, *tp, *tm, *tt,sixthTrack,oneMasses[j] ,twoMasses[j] ,threeMasses[j],fourMasses[j]).mass();

          std::string name;
          for(size_t j = 1; j<numMasses_+1;j++)
          {
           name = "mass_ref_" + std::to_string(j);
           sixCand.addUserFloat(name,sixTracksMass[j-1]);
          }

          sixCand.addUserInt("five_id",int(d));

          sixCand.addDaughter(*tp,"trakOne");
          sixCand.addDaughter(*tm,"trakTwo");
          sixCand.addDaughter(*tm,"trakThree");
          sixCand.addDaughter(sixthTrack,"trakFour");

          sixCand.addDaughter(sixCand_rf,"ref_cand");

          sixCand.addUserFloat("mass_ref_0",six_ma_fit);

          sixCand.addUserFloat("fourthTrackMuonDR",minDR_fourth);
          sixCand.addUserFloat("fourthTrackMuonDP",minDP_fourth);
          sixCand.addUserFloat("fourthTrackMuonDPt",minDPt_fourth);

           if(atLeastOne)
            sixCandColl->push_back(sixCand);

           }

        }


     ncombo += sixCandColl->size();

  iEvent.put(std::move(sixCandColl),"SixTracks");

  nevents++;
}

void SixTracksProducer::endJob(){
  std::cout << "#########################################" << std::endl;
  std::cout << "SixTracks Candidate producer report:" << std::endl;
  std::cout << "#########################################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "No. dim + 4 trk candidates " << ncombo << std::endl;
  std::cout << "#########################################" << std::endl;
}

bool SixTracksProducer::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.01 && DeltaP < 0.01) return true;
  return false;
}

bool SixTracksProducer::IsTheSame(const pat::PackedCandidate& t1, const pat::PackedCandidate& t2){
  double DeltaEta = fabs(t1.eta()-t2.eta());
  double DeltaP   = fabs(t1.p()-t2.p());
  if (DeltaEta < 0.01 && DeltaP < 0.01) return true;
  return false;
}


pat::CompositeCandidate SixTracksProducer::makeSixCandidate(
                                          const pat::CompositeCandidate& fivetrak,
                                          const pat::PackedCandidate& trak
                                         ){

  pat::CompositeCandidate fiveCand;
  fiveCand.addDaughter(fivetrak,"fivetrak");
  fiveCand.addDaughter(trak,"sixth");
  fiveCand.setCharge(fivetrak.charge()+trak.charge());

  double m_trak = trackmass;
  math::XYZVector mom_trak = trak.momentum();
  double e_trak = sqrt(m_trak*m_trak + mom_trak.Mag2());
  math::XYZTLorentzVector p4_trak = math::XYZTLorentzVector(mom_trak.X(),mom_trak.Y(),mom_trak.Z(),e_trak);

  reco::Candidate::LorentzVector v = p4_trak + fiveCand.p4();

  fiveCand.setP4(v);

  return fiveCand;
}

pat::CompositeCandidate SixTracksProducer::makeSixCandidateMixed(
                                          const pat::CompositeCandidate& dimuon,
                                          const pat::PackedCandidate& trakP,
                                          const pat::PackedCandidate& trakN,
                                          const pat::PackedCandidate& trak3,
                                          const pat::PackedCandidate& trak4,
                                          double massOne,
                                          double massTwo,
                                          double massThree,
                                          double massFour
                                         ){

  pat::CompositeCandidate sixCand, trakOne, trakTwo, trakThree, trakFour;

  sixCand.setCharge(dimuon.charge()+trakP.charge()+trakN.charge()+trak3.charge());

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

  math::XYZVector mom_trak4 = trak4.momentum();
  double e_trak4 = sqrt(massThree*massThree + mom_trak4.Mag2());
  math::XYZTLorentzVector p4_trak4 = math::XYZTLorentzVector(mom_trak4.X(),mom_trak4.Y(),mom_trak4.Z(),e_trak4);
  trakFour.setCharge(trak4.charge());
  trakFour.setP4(p4_trak4);

  reco::Candidate::LorentzVector v = p4_trakP + p4_trakN + p4_trak3 + p4_trak4 + dimuon.p4();

  sixCand.setP4(v);

  return sixCand;
}


reco::Candidate::LorentzVector SixTracksProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(SixTracksProducer);
