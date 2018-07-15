#include "../interface/FiveTracksProducerFit.h"

float FiveTracksProducerFit::DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
   float p1 = t1.phi();
   float p2 = t2.phi();
   float e1 = t1.eta();
   float e2 = t2.eta();
   auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

bool FiveTracksProducerFit::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<maxDPtRel &&
	DeltaR(t1,t2) < maxDeltaR);
}

bool
FiveTracksProducerFit::isAbHadron(int pdgID) {

  if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
  return false;

}

bool
FiveTracksProducerFit::isAMixedbHadron(int pdgID, int momPdgID) {

  if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) ||
  (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0))
  return true;
  return false;

}

bool
FiveTracksProducerFit::isTheCandidate(reco::GenParticleRef genY) {

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

bool FiveTracksProducerFit::isSameTrack(reco::Track t1, reco::Track t2)
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

FiveTracksProducerFit::FiveTracksProducerFit(const edm::ParameterSet& iConfig):
  DiMuonDiTrakCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiMuoDiTrak"))),
  TrakCollection_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
  thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
  thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
  TriggerCollection_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  FiveTrakMassCuts_(iConfig.getParameter<std::vector<double>>("FiveTrakCuts")),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters")),
  isMC_(iConfig.getParameter<bool>("IsMC")),
  addMCTruth_(iConfig.getParameter<bool>("AddMCTruth")),
  addSameSig_(iConfig.getParameter<bool>("AddSS"))
{
  produces<pat::CompositeCandidateCollection>("FiveTracksKaon");

  nevents = 0;

  kaonmass = 0.493677;
  pionmass = 0.13957061;

  maxDeltaR = 0.01;
  maxDPtRel = 2.0;
  trackmass = kaonmass;

  ncombo = 0;

}

void FiveTracksProducerFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::unique_ptr<pat::CompositeCandidateCollection> fiveCandKaonColl(new pat::CompositeCandidateCollection);

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

  for (size_t d = 0; d < dimuonditrak->size(); d++) {

       auto dimuonditrakCand = dimuonditrak->at(d);

       if(dimuonditrakCand.userFloat("vProb")<0.0)
         continue;

       const reco::Vertex thePrimaryV = *dimuonditrakCand.userData<reco::Vertex>("bestPV");
       // const reco::Vertex thePrimaryV = *dimuonditrakCand.userData<reco::Vertex>("PVwithmuons");

       const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(dimuonditrakCand.daughter("dimuon")->daughter("highMuon"));
       const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(dimuonditrakCand.daughter("dimuon")->daughter("lowMuon"));
       // const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(pmu1->originalObject());
       // const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(pmu2->originalObject());

       const pat::PackedCandidate *tp = dynamic_cast <pat::PackedCandidate *>(dimuonditrakCand.daughter("ditrak")->daughter("highTrak"));
       const pat::PackedCandidate *tm = dynamic_cast <pat::PackedCandidate *>(dimuonditrakCand.daughter("ditrak")->daughter("lowTrak"));
       int tpId = dimuonditrakCand.userInt("pId");
       int tmId = dimuonditrakCand.userInt("mId");

//Adding a kaon
       for (size_t i = 0; i < trak->size(); i++) {
         auto fifthTrack = trak->at(i);

         if(fifthTrack.pt()<0.7) continue;
         if(fifthTrack.charge() == 0) continue;
	       //if(!isMC_ and fabs(fifthTrack.pdgId())!=211) continue;
	       if(!(fifthTrack.trackHighPurity())) continue;
         if(!(fifthTrack.hasTrackDetails())) continue;

         if (IsTheSame(fifthTrack,*tp) || int(i) == tpId) continue;
         if (IsTheSame(fifthTrack,*tm) || int(i) == tmId) continue;
         if ( IsTheSame(fifthTrack,*pmu1) || IsTheSame(fifthTrack,*pmu2) ) continue;

         trackmass = kaonmass;
         pat::CompositeCandidate fiveCandKaon = makeFiveCandidate(dimuonditrakCand, fifthTrack);
         trackmass = pionmass;
         pat::CompositeCandidate fiveCandPion = makeFiveCandidate(dimuonditrakCand, fifthTrack);

         if (fiveCandKaon.mass() > FiveTrakMassMax || fiveCandKaon.mass() < FiveTrakMassMin)
         if (fiveCandPion.mass() > FiveTrakMassMax || fiveCandPion.mass() < FiveTrakMassMin)
          continue;

         //KaonRefit
         const ParticleMass muonMass(0.1056583);
         float muonSigma = muonMass*1E-6;
         const ParticleMass trakMassP(tp->mass());
         float trakSigmaP = trakMassP*1E-6;
         const ParticleMass trakMassM(tm->mass());
         float trakSigmaM = trakMassM*1E-6;
         const ParticleMass fifthKaonMass(kaonmass);
         float fifthSigma = fifthKaonMass*1E-6;
         const ParticleMass fifthPionMass(pionmass);

         std::vector<reco::TransientTrack> fiveTracks;
         std::vector<RefCountedKinematicParticle> kaonParticles, pionParticles;

         float kinChi = 0.;
         float kinNdf = 0.;

         fiveTracks.push_back((*theB).build(*(pmu1->innerTrack()))); // K+
         fiveTracks.push_back((*theB).build(*(pmu2->innerTrack()))); // K+
         fiveTracks.push_back((*theB).build(*(tp->bestTrack()))); // K+
         fiveTracks.push_back((*theB).build(*(tm->bestTrack()))); // K+
         fiveTracks.push_back((*theB).build(*(fifthTrack.bestTrack()))); // K+

         kaonParticles.push_back(pFactory.particle(fiveTracks[0],muonMass,kinChi,kinNdf,muonSigma));
         kaonParticles.push_back(pFactory.particle(fiveTracks[1],muonMass,kinChi,kinNdf,muonSigma));
         kaonParticles.push_back(pFactory.particle(fiveTracks[2],trakMassP,kinChi,kinNdf,trakSigmaP));
         kaonParticles.push_back(pFactory.particle(fiveTracks[3],trakMassM,kinChi,kinNdf,trakSigmaM));
         kaonParticles.push_back(pFactory.particle(fiveTracks[4],fifthKaonMass,kinChi,kinNdf,fifthSigma));

         KinematicParticleVertexFitter kaonFitter;
         RefCountedKinematicTree kaonVertexFitTree;
         kaonVertexFitTree = kaonFitter.fit(kaonParticles);

         if (kaonVertexFitTree->isEmpty()) continue;

         kaonVertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle fitF = kaonVertexFitTree->currentParticle();
         RefCountedKinematicVertex fitFVertex = kaonVertexFitTree->currentDecayVertex();

         double kaon_ma_fit = 14000.;
         double kaon_vp_fit = -9999.;
         double kaon_x2_fit = 10000.;
         double kaon_ndof_fit = 10000.;

         if (!(fitF->currentState().isValid())) continue;

         kaon_ma_fit = fitF->currentState().mass();
         kaon_x2_fit = fitFVertex->chiSquared();
         kaon_vp_fit = ChiSquaredProbability(kaon_x2_fit,
                                              (double)(fitFVertex->degreesOfFreedom()));
         kaon_ndof_fit = (double)(fitFVertex->degreesOfFreedom());

         kinChi = 0.;
         kinNdf = 0.;
         fiveTracks.clear();
         fiveTracks.push_back((*theB).build(*(pmu1->innerTrack()))); // K+
         fiveTracks.push_back((*theB).build(*(pmu2->innerTrack()))); // K+
         fiveTracks.push_back((*theB).build(*(tp->bestTrack()))); // K+
         fiveTracks.push_back((*theB).build(*(tm->bestTrack()))); // K+
         fiveTracks.push_back((*theB).build(*(fifthTrack.bestTrack()))); // K+

         pionParticles.push_back(pFactory.particle(fiveTracks[0],muonMass,kinChi,kinNdf,muonSigma));
         pionParticles.push_back(pFactory.particle(fiveTracks[1],muonMass,kinChi,kinNdf,muonSigma));
         pionParticles.push_back(pFactory.particle(fiveTracks[2],trakMassP,kinChi,kinNdf,trakSigmaP));
         pionParticles.push_back(pFactory.particle(fiveTracks[3],trakMassM,kinChi,kinNdf,trakSigmaM));
         pionParticles.push_back(pFactory.particle(fiveTracks[4],fifthPionMass,kinChi,kinNdf,fifthSigma));

         KinematicParticleVertexFitter pionFitter;
         RefCountedKinematicTree pionVertexFitTree;
         pionVertexFitTree = pionFitter.fit(pionParticles);

         if (pionVertexFitTree->isEmpty()) continue;

         pionVertexFitTree->movePointerToTheTop();
         fitF = pionVertexFitTree->currentParticle();
         fitFVertex = pionVertexFitTree->currentDecayVertex();

         double pion_ma_fit = 14000.;

         if (!(fitF->currentState().isValid())) continue;

         pion_ma_fit = fitF->currentState().mass();

         TVector3 vtx;
         TVector3 pvtx;
         VertexDistanceXY vdistXY;
         // int   kaon_ch_fit = fiveCandKaon.charge();
         double kaon_px_fit = fitF->currentState().kinematicParameters().momentum().x();
         double kaon_py_fit = fitF->currentState().kinematicParameters().momentum().y();
         // double kaon_pz_fit = fitF->currentState().kinematicParameters().momentum().z();
         // double kaon_en_fit = sqrt(kaon_ma_fit*kaon_ma_fit+kaon_pkaon_fit*kaon_pkaon_fit+kaon_py_fit*kaon_py_fit+kaon_pz_fit*kaon_pz_fit);
         double kaon_vx_fit = fitFVertex->position().x();
         double kaon_vy_fit = fitFVertex->position().y();
         // double kaon_vz_fit = fitFVertex->position().z();

         vtx.SetXYZ(kaon_vx_fit,kaon_vy_fit,0);
         TVector3 pperp(kaon_px_fit, kaon_py_fit, 0);
         AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
         pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
         TVector3 vdiff = vtx - pvtx;
         double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
         Measurement1D distXY = vdistXY.distance(reco::Vertex(*fitFVertex), thePrimaryV);
         double ctauPV = distXY.value()*cosAlpha * kaon_ma_fit/pperp.Perp();
         GlobalError v1e = (reco::Vertex(*fitFVertex)).error();
         GlobalError v2e = thePrimaryV.error();
         AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
         double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*kaon_ma_fit/(pperp.Perp2());

         fiveCandKaon.addUserFloat("mass_kaon_rf",kaon_ma_fit);
         fiveCandKaon.addUserFloat("mass_pion_rf",pion_ma_fit);
         fiveCandKaon.addUserFloat("mass_pion",fiveCandPion.mass());
         fiveCandKaon.addUserFloat("vProb",kaon_vp_fit);
         fiveCandKaon.addUserFloat("vChi2",kaon_x2_fit);
         fiveCandKaon.addUserFloat("nDof",kaon_ndof_fit);
         fiveCandKaon.addUserFloat("cosAlpha",cosAlpha);
         fiveCandKaon.addUserFloat("ctauPV",ctauPV);
         fiveCandKaon.addUserFloat("ctauErrPV",ctauErrPV);

         fiveCandKaon.addUserInt("index",d);

         fiveCandKaonColl->push_back(fiveCandKaon);

         ++ncombo;
       }
     }

  iEvent.put(std::move(fiveCandKaonColl),"FiveTracksKaon");
  nevents++;
}

void FiveTracksProducerFit::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "FiveTracks Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "No. dimuonditrak + trk candidates " << ncombo << std::endl;
  std::cout << "###########################" << std::endl;
}

bool FiveTracksProducerFit::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.01 && DeltaP < 0.01) return true;
  return false;
}

bool FiveTracksProducerFit::IsTheSame(const pat::PackedCandidate& t1, const pat::PackedCandidate& t2){
  double DeltaEta = fabs(t1.eta()-t2.eta());
  double DeltaP   = fabs(t1.p()-t2.p());
  if (DeltaEta < 0.01 && DeltaP < 0.01) return true;
  return false;
}


pat::CompositeCandidate FiveTracksProducerFit::makeFiveCandidate(
                                          const pat::CompositeCandidate& dimuonditrak,
                                          const pat::PackedCandidate& trak
                                         ){

  pat::CompositeCandidate fiveCandKaon, dimuon, dimuontrak;
  fiveCandKaon.addDaughter(dimuonditrak,"dimuonditrak");
  fiveCandKaon.addDaughter(trak,"fifth");
  fiveCandKaon.setCharge(dimuonditrak.charge()+trak.charge());
  dimuon = dynamic_cast <pat::CompositeCandidate *>(dimuonditrak.daughter("dimuon"));

  double m_trak = trackmass;
  math::XYZVector mom_trak = trak.momentum();
  double e_trak = sqrt(m_trak*m_trak + mom_trak.Mag2());
  math::XYZTLorentzVector p4_trak = math::XYZTLorentzVector(mom_trak.X(),mom_trak.Y(),mom_trak.Z(),e_trak);

  reco::Candidate::LorentzVector v = p4_trak + dimuonditrak.p4();
  reco::Candidate::LorentzVector j = p4_trak + dimuontrak.p4();

  fiveCandKaon.setP4(v);
  dimuontrak.setP4(j);
  fiveCandKaon.addDaughter(trak,"dimuontrak");

  return fiveCandKaon;
}


reco::Candidate::LorentzVector FiveTracksProducerFit::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(FiveTracksProducerFit);
