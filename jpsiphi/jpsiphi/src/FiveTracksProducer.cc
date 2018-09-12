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
  FiveTrakMassCuts_(iConfig.getParameter<std::vector<double>>("FiveTrakCuts"))
{
  produces<pat::CompositeCandidateCollection>("FiveTracks");

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

       const reco::Vertex thePrimaryV = *(dimuonditrakCand.userData<reco::Vertex>("bestPV"));
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

       std::vector<double> oneMasses,twoMasses,threeMasses, hasRefit;
       oneMasses.push_back(kaonmass); oneMasses.push_back(pionmass);oneMasses.push_back(kaonmass);oneMasses.push_back(pionmass);oneMasses.push_back(pionmass);
       twoMasses.push_back(kaonmass); twoMasses.push_back(pionmass);twoMasses.push_back(pionmass);twoMasses.push_back(kaonmass);twoMasses.push_back(pionmass);
       threeMasses.push_back(kaonmass); threeMasses.push_back(kaonmass);threeMasses.push_back(pionmass);threeMasses.push_back(pionmass);threeMasses.push_back(pionmass);


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
       for (size_t i = 0; i < trak->size(); i++) {

         std::vector<float> fiveTracksMass, fiveTracksVProb, fiveTracksVChi2, fiveTracksVNDof;
         std::vector<float> psi2sOne, psi2sTwo;
         std::vector<float> fiveTracksCTau, fiveTracksCTauErr, fiveTracksCosAlpha;


         auto fifthTrack = trak->at(i);


         if(fifthTrack.pt()<trakPtCut_) continue;
         //if(fifthTrack.charge() == 0) continue;
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
         trackmass = kaonmass;
         pat::CompositeCandidate fiveCandKaon = makeFiveCandidate(dimuonditrakCand, fifthTrack);
         trackmass = pionmass;
         pat::CompositeCandidate fiveCandPion = makeFiveCandidate(dimuonditrakCand, fifthTrack);


         double five_ma_fit = 14000.;
         double five_vp_fit = -9999.;
         double five_x2_fit = 10000.;
         double five_nd_fit = 10000.;

         std::vector<pat::CompositeCandidate> fiveCands, ref_fiveCands;

         bool atLeastOne = false;
         std::vector< bool > insideMass;
         for(size_t j = 0; j<oneMasses.size();j++)
         {

             fiveTracksMass.push_back(FiveTrakMassMin-0.2);
             fiveTracksVNDof.push_back(-1.0);
             fiveTracksVChi2.push_back(-1.0);
             fiveTracksVProb.push_back(-0.1);
             fiveTracksCTau.push_back(-1000.0);
             fiveTracksCTauErr.push_back(-1000.0);
             fiveTracksCosAlpha.push_back(-1.1);
             insideMass.push_back(false);
             hasRefit.push_back(0.0);

             //if(i!=5 && fifthTrack.charge() == 0) continue;
             pat::CompositeCandidate thisFive = makeFiveCandidateMixed(*dimuon_cand, *tp, *tm, fifthTrack,oneMasses[j] ,twoMasses[j] ,threeMasses[j]);
             fiveCands.push_back(thisFive);

             //Kinematic Fit
             const ParticleMass muonMass(0.1056583);
             float muonSigma = muonMass*1E-6;
             const ParticleMass trakMassP(oneMasses[j]);
             float trakSigmaP = trakMassP*1E-6;
             const ParticleMass trakMassM(twoMasses[j]);
             float trakSigmaM = trakMassM*1E-6;
             const ParticleMass fifthMass(threeMasses[j]);
             float fifthSigma = fifthMass*1E-6;


             if (thisFive.mass() < FiveTrakMassMax && thisFive.mass() > FiveTrakMassMin)
             atLeastOne = true;

             if (thisFive.mass() < FiveTrakMassMax && thisFive.mass() > FiveTrakMassMin)
             insideMass[j] = true;

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

             KinematicParticleVertexFitter fiveFitter;
             RefCountedKinematicTree fiveVertexFitTree;
             fiveVertexFitTree = fiveFitter.fit(kaonParticles);

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

             //five_ndof_fit = (double)(fitFVertex->degreesOfFreedom());

             int    five_ch_fit = fiveCandKaon.charge();
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

             TVector3 vtx;
             TVector3 pvtx;
             VertexDistanceXY vdistXY;

             vtx.SetXYZ(five_vx_fit,five_vy_fit,0);
             TVector3 pperp(five_px_fit, five_py_fit, 0);
             AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
             pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
             TVector3 vdiff = vtx - pvtx;
             double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
             Measurement1D distXY = vdistXY.distance(reco::Vertex(*fitFVertex), thePrimaryV);

             double ctauPV = distXY.value()*cosAlpha * five_ma_fit/pperp.Perp();
             GlobalError v1e = (reco::Vertex(*fitFVertex)).error();
             GlobalError v2e = thePrimaryV.error();
             AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
             double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*five_ma_fit/(pperp.Perp2());

             fiveTracksMass[j]      = five_ma_fit;
             fiveTracksVProb[j]     = five_vp_fit;
             fiveTracksCTau[j]      = ctauPV;
             fiveTracksCTauErr[j]   = ctauErrPV;
             fiveTracksCosAlpha[j]  = cosAlpha;
             fiveTracksVNDof[j]     = five_nd_fit;
             fiveTracksVChi2[j]     = five_x2_fit;

 // get first muon
             bool child = fiveVertexFitTree->movePointerToTheFirstChild();
             RefCountedKinematicParticle fitPart = fiveVertexFitTree->currentParticle();
             if (!child) break;
             float m1_ma_fit = fitPart->currentState().mass();
             int   m1_ch_fit = fitPart->currentState().particleCharge();
             float m1_px_fit = fitPart->currentState().kinematicParameters().momentum().x();
             float m1_py_fit = fitPart->currentState().kinematicParameters().momentum().y();
             float m1_pz_fit = fitPart->currentState().kinematicParameters().momentum().z();
             float m1_en_fit = sqrt(m1_ma_fit*m1_ma_fit+m1_px_fit*m1_px_fit+m1_py_fit*m1_py_fit+m1_pz_fit*m1_pz_fit);
             reco::CompositeCandidate recoMu1(m1_ch_fit,math::XYZTLorentzVector(m1_px_fit,m1_py_fit,m1_pz_fit,m1_en_fit),
                                              math::XYZPoint(five_vx_fit,five_vy_fit,five_vz_fit),13);
             pat::CompositeCandidate patMu1(recoMu1);
 // get second muon
             child = fiveVertexFitTree->movePointerToTheNextChild();
             fitPart = fiveVertexFitTree->currentParticle();
             if (!child) break;
             float m2_ma_fit = fitPart->currentState().mass();
             int   m2_ch_fit = fitPart->currentState().particleCharge();
             float m2_px_fit = fitPart->currentState().kinematicParameters().momentum().x();
             float m2_py_fit = fitPart->currentState().kinematicParameters().momentum().y();
             float m2_pz_fit = fitPart->currentState().kinematicParameters().momentum().z();
             float m2_en_fit = sqrt(m2_ma_fit*m2_ma_fit+m2_px_fit*m2_px_fit+m2_py_fit*m2_py_fit+m2_pz_fit*m2_pz_fit);
             reco::CompositeCandidate recoMu2(m2_ch_fit,math::XYZTLorentzVector(m2_px_fit,m2_py_fit,m2_pz_fit,m2_en_fit),
                                              math::XYZPoint(five_vx_fit,five_vy_fit,five_vz_fit),13);
             pat::CompositeCandidate patMu2(recoMu2);

 // Define psi from two muons
 	           pat::CompositeCandidate psi;
 	           psi.addDaughter(patMu1,"highMuon");
             psi.addDaughter(patMu2,"lowMuon");
             psi.setP4(patMu1.p4()+patMu2.p4());
 // get tn
             child = fiveVertexFitTree->movePointerToTheNextChild();
             fitPart = fiveVertexFitTree->currentParticle();
             if (!child) break;
             float tk1_ma_fit = fitPart->currentState().mass();
             int   tk1_ch_fit = fitPart->currentState().particleCharge();
             float tk1_px_fit = fitPart->currentState().kinematicParameters().momentum().x();
             float tk1_py_fit = fitPart->currentState().kinematicParameters().momentum().y();
             float tk1_pz_fit = fitPart->currentState().kinematicParameters().momentum().z();
             float tk1_en_fit = sqrt(tk1_ma_fit*tk1_ma_fit+tk1_px_fit*tk1_px_fit+tk1_py_fit*tk1_py_fit+tk1_pz_fit*tk1_pz_fit);
             reco::CompositeCandidate recoTk1(tk1_ch_fit,math::XYZTLorentzVector(tk1_px_fit,tk1_py_fit,tk1_pz_fit,tk1_en_fit),
                                              math::XYZPoint(five_vx_fit,five_vy_fit,five_vz_fit),321);
             pat::CompositeCandidate patTk1(recoTk1);

 // get tn
             child = fiveVertexFitTree->movePointerToTheNextChild();
             fitPart = fiveVertexFitTree->currentParticle();
             if (!child) break;
             float tk2_ma_fit = fitPart->currentState().mass();
             int   tk2_ch_fit = fitPart->currentState().particleCharge();
             float tk2_px_fit = fitPart->currentState().kinematicParameters().momentum().x();
             float tk2_py_fit = fitPart->currentState().kinematicParameters().momentum().y();
             float tk2_pz_fit = fitPart->currentState().kinematicParameters().momentum().z();
             float tk2_en_fit = sqrt(tk2_ma_fit*tk2_ma_fit+tk2_px_fit*tk2_px_fit+tk2_py_fit*tk2_py_fit+tk2_pz_fit*tk2_pz_fit);
             reco::CompositeCandidate recoTk2(tk2_ch_fit,math::XYZTLorentzVector(tk2_px_fit,tk2_py_fit,tk2_pz_fit,tk2_en_fit),
                                              math::XYZPoint(five_vx_fit,five_vy_fit,five_vz_fit),321);
             pat::CompositeCandidate patTk2(recoTk2);
 // get tn
             child = fiveVertexFitTree->movePointerToTheNextChild();
             fitPart = fiveVertexFitTree->currentParticle();
             if (!child) break;
             float tk3_ma_fit = fitPart->currentState().mass();
             int   tk3_ch_fit = fitPart->currentState().particleCharge();
             float tk3_px_fit = fitPart->currentState().kinematicParameters().momentum().x();
             float tk3_py_fit = fitPart->currentState().kinematicParameters().momentum().y();
             float tk3_pz_fit = fitPart->currentState().kinematicParameters().momentum().z();
             float tk3_en_fit = sqrt(tk3_ma_fit*tk3_ma_fit+tk3_px_fit*tk3_px_fit+tk3_py_fit*tk3_py_fit+tk3_pz_fit*tk3_pz_fit);

             reco::CompositeCandidate recoTk3(tk3_ch_fit,math::XYZTLorentzVector(tk3_px_fit,tk3_py_fit,tk3_pz_fit,tk3_en_fit),
                                              math::XYZPoint(five_vx_fit,five_vy_fit,five_vz_fit),321);
             pat::CompositeCandidate patTk3(recoTk3);

               fiveCands[j].addDaughter(makeFiveCandidateMixed(psi,patTk1,patTk2,patTk3),"fiveRef");

            }

             fiveCandKaon.addUserInt("dimuontt_index",int(d));

             // fiveCandKaon.addDaughter(dimuonditrakCand,"dimuonditrak");
             fiveCandKaon.addDaughter(*tp,"trakOne");
             fiveCandKaon.addDaughter(*tm,"trakTwo");
             fiveCandKaon.addDaughter(fifthTrack,"trakThree");

             std::string name;
             for(size_t j = 0; j<fiveCands.size();j++)
             {

               fiveCands[j].addUserFloat("vProb",fiveTracksVProb[j]);
               fiveCands[j].addUserFloat("ctauPV",fiveTracksCTau[j]);
               fiveCands[j].addUserFloat("ctauErrPV",fiveTracksCTauErr[j]);
               fiveCands[j].addUserFloat("cosAlpha",fiveTracksCosAlpha[j]);
               fiveCands[j].addUserFloat("nDof",fiveTracksVNDof[j]);
               fiveCands[j].addUserFloat("vChi2",fiveTracksVChi2[j]);
               fiveCands[j].addUserFloat("mass_ref",fiveTracksMass[j]);
               fiveCands[j].addUserFloat("has_ref",hasRefit[j]);

               name = "fiveCand_" + std::to_string(j);

               fiveCandKaon.addDaughter(fiveCands[j],name);
              // std::cout << fiveCands[j].mass() << " - ";

             }

             // std::cout << std::endl;
             if(atLeastOne)
              fiveCandColl->push_back(fiveCandKaon);

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

  pat::CompositeCandidate fiveCandKaon, dimuontrak;
  fiveCandKaon.addDaughter(dimuonditrak,"dimuonditrak");
  fiveCandKaon.addDaughter(trak,"fifth");
  fiveCandKaon.setCharge(dimuonditrak.charge()+trak.charge());

  double m_trak = trackmass;
  math::XYZVector mom_trak = trak.momentum();
  double e_trak = sqrt(m_trak*m_trak + mom_trak.Mag2());
  math::XYZTLorentzVector p4_trak = math::XYZTLorentzVector(mom_trak.X(),mom_trak.Y(),mom_trak.Z(),e_trak);

  reco::Candidate::LorentzVector v = p4_trak + dimuonditrak.p4();

  fiveCandKaon.setP4(v);

  return fiveCandKaon;
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
