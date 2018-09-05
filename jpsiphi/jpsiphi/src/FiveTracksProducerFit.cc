#include "../interface/FiveTracksProducerFit.h"
#include <tuple>
#include <map>

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
  FiveTrakMassCuts_(iConfig.getParameter<std::vector<double>>("FiveTrakCuts"))
{
  produces<pat::CompositeCandidateCollection>("FiveTracksPos");
  produces<pat::CompositeCandidateCollection>("FiveTracksNeg");
  produces<pat::CompositeCandidateCollection>("FiveTracksNeu");

  nevents = 0;

  kaonmass = 0.493677;
  pionmass = 0.13957061;

  maxDeltaR = 0.01;
  maxDPtRel = 2.0;
  trackmass = kaonmass;

  ncomboneg = 0;
  ncombo = 0;
  ncomboneu = 0;
}

void FiveTracksProducerFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

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

       const pat::PackedCandidate *tp = dynamic_cast <pat::PackedCandidate *>(dimuonditrakCand.daughter("ditrak")->daughter("highTrak"));
       const pat::PackedCandidate *tm = dynamic_cast <pat::PackedCandidate *>(dimuonditrakCand.daughter("ditrak")->daughter("lowTrak"));
       int tpId = dimuonditrakCand.userInt("pId");
       int tmId = dimuonditrakCand.userInt("mId");

       std::vector<float> oneMasses,twoMasses,threeMasses;
       oneMasses.push_back(kaonmass);oneMasses.push_back(pionmass);oneMasses.push_back(kaonmass);oneMasses.push_back(pionmass);oneMasses.push_back(pionmass);
       twoMasses.push_back(kaonmass);twoMasses.push_back(pionmass);twoMasses.push_back(pionmass);twoMasses.push_back(kaonmass);twoMasses.push_back(pionmass);
       threeMasses.push_back(kaonmass);threeMasses.push_back(kaonmass);threeMasses.push_back(pionmass);threeMasses.push_back(pionmass);threeMasses.push_back(pionmass);


//Adding a kaon
       for (size_t i = 0; i < trak->size(); i++) {

         std::vector<float> fiveTracksMass, fiveTracksVProb;// fiveTracksPt, fiveTracksCharge;
         std::vector<float> psi2sOne, psi2sTwo;
         std::vector<float> fiveTracksCTau, fiveTracksCTauErr, fiveTracksCosAlpha;
         auto fifthTrack = trak->at(i);

         if(fifthTrack.pt()<0.7) continue;
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
         trackmass = kaonmass;
         pat::CompositeCandidate fiveCandKaon = makeFiveCandidate(dimuonditrakCand, fifthTrack);
         trackmass = pionmass;
         pat::CompositeCandidate fiveCandPion = makeFiveCandidate(dimuonditrakCand, fifthTrack);

         if (fiveCandKaon.mass() > FiveTrakMassMax || fiveCandKaon.mass() < FiveTrakMassMin)
         if (fiveCandPion.mass() > FiveTrakMassMax || fiveCandPion.mass() < FiveTrakMassMin)
          continue;

         double kaon_ma_fit = 14000.;
         double kaon_vp_fit = -9999.;
         double kaon_x2_fit = 10000.;
         // double kaon_ndof_fit = 10000.;

         for(size_t i = 0; i<oneMasses.size();i++)
         {

             fiveTracksMass.push_back(-1.0);
             // fiveTracksPt.push_back(-1.0);
             fiveTracksVProb.push_back(-1.0);
             psi2sOne.push_back(-1.0);
             psi2sTwo.push_back(-1.0);
             fiveTracksCTau.push_back(-1.0);
             fiveTracksCTauErr.push_back(-1.0);
             fiveTracksCosAlpha.push_back(-1.0);
             //fiveTracksMassRef.push_back(-1.0);
             // fiveTracksCharge.push_back(-5.0);

             //KaonRefit
             const ParticleMass muonMass(0.1056583);
             float muonSigma = muonMass*1E-6;
             const ParticleMass trakMassP(oneMasses[i]);
             float trakSigmaP = trakMassP*1E-6;
             const ParticleMass trakMassM(twoMasses[i]);
             float trakSigmaM = trakMassM*1E-6;
             const ParticleMass fifthMass(threeMasses[i]);
             float fifthSigma = fifthMass*1E-6;


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
             kaonParticles.push_back(pFactory.particle(fiveTracks[4],fifthMass,kinChi,kinNdf,fifthSigma));

             KinematicParticleVertexFitter kaonFitter;
             RefCountedKinematicTree kaonVertexFitTree;
             kaonVertexFitTree = kaonFitter.fit(kaonParticles);

             if (kaonVertexFitTree->isEmpty()) continue;

             kaonVertexFitTree->movePointerToTheTop();
             RefCountedKinematicParticle fitF = kaonVertexFitTree->currentParticle();
             RefCountedKinematicVertex fitFVertex = kaonVertexFitTree->currentDecayVertex();

             if (!(fitF->currentState().isValid())) continue;

             kaon_ma_fit = fitF->currentState().mass();
             kaon_x2_fit = fitFVertex->chiSquared();
             kaon_vp_fit = ChiSquaredProbability(kaon_x2_fit,
                                                  (double)(fitFVertex->degreesOfFreedom()));
             //kaon_ndof_fit = (double)(fitFVertex->degreesOfFreedom());


             double kaon_px_fit = fitF->currentState().kinematicParameters().momentum().x();
             double kaon_py_fit = fitF->currentState().kinematicParameters().momentum().y();
             // double kaon_pz_fit = fitF->currentState().kinematicParameters().momentum().z();
             // double kaon_en_fit = sqrt(kaon_ma_fit*kaon_ma_fit+kaon_pkaon_fit*kaon_pkaon_fit+kaon_py_fit*kaon_py_fit+kaon_pz_fit*kaon_pz_fit);
             double kaon_vx_fit = fitFVertex->position().x();
             double kaon_vy_fit = fitFVertex->position().y();
         // double kaon_vz_fit = fitFVertex->position().z();

             TVector3 vtx;
             TVector3 pvtx;
             VertexDistanceXY vdistXY;

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

             fiveTracksMass[i] = kaon_ma_fit;
//             fiveTracksPt[i] = fiveCandKaon.pt();
             fiveTracksVProb[i] = kaon_vp_fit;
//             fiveTracksCharge[i] = fiveCandKaon.charge();
             fiveTracksCTau[i] = ctauPV;
             fiveTracksCTauErr[i] = ctauErrPV;
             fiveTracksCosAlpha[i] = cosAlpha;

             pat::CompositeCandidate onePsi2S,twoPsi2S;
             onePsi2S.setMass(-1.0); twoPsi2S.setMass(-1.0);
             if(i==4)
             {
               onePsi2S = makePsi2SCandidate(*dimuon_cand,*tp,*tm);

               if(fifthTrack.charge()<0)
               {
                 twoPsi2S = makePsi2SCandidate(*dimuon_cand,*tp,fifthTrack);
               }else
               {
                 twoPsi2S = makePsi2SCandidate(*dimuon_cand,fifthTrack,*tm);
               }
             }
             if(i==1)
               onePsi2S = makePsi2SCandidate(*dimuon_cand,*tp,*tm);
             if(i==2 && fifthTrack.charge()>0)
               onePsi2S = makePsi2SCandidate(*dimuon_cand,fifthTrack,*tm);
             if(i==3 && fifthTrack.charge()<0)
               onePsi2S = makePsi2SCandidate(*dimuon_cand,*tp,fifthTrack);

             psi2sOne[i] = onePsi2S.mass();
             psi2sTwo[i] = onePsi2S.mass();

             }

             fiveCandKaon.addUserInt("index",d);
             for(size_t i = 0; i<oneMasses.size();i++)
             {
               std::string name = "mass_" + std::to_string(i);
               fiveCandKaon.addUserFloat(name,fiveTracksMass[i]);
               name = "vProb_" + std::to_string(i);
               fiveCandKaon.addUserFloat(name,fiveTracksVProb[i]);
               name = "ctau_" + std::to_string(i);
               fiveCandKaon.addUserFloat(name,fiveTracksCTau[i]);
               name = "ctauErr_" + std::to_string(i);
               fiveCandKaon.addUserFloat(name,fiveTracksCTauErr[i]);
               name = "cosAlpha_" + std::to_string(i);
               fiveCandKaon.addUserFloat(name,fiveTracksCosAlpha[i]);
               name = "onePsi2S_" + std::to_string(i);
               fiveCandKaon.addUserFloat(name,psi2sOne[i]);
               name = "twoPsi2S_" + std::to_string(i);
               fiveCandKaon.addUserFloat(name,psi2sTwo[i]);

             }

             fiveCandColl->push_back(fiveCandKaon);
           }

        }


     ncombo += fiveCandColl->size();

  iEvent.put(std::move(fiveCandColl),"FiveTracksPos");

  nevents++;
}

void FiveTracksProducerFit::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "FiveTracks Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "No. dimtt + trk pos candidates " << ncombo << std::endl;
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

pat::CompositeCandidate FiveTracksProducerFit::makePsi2SCandidate(
                                          const pat::CompositeCandidate& dimuonditrak,
                                          const pat::PackedCandidate& t1,
                                          const pat::PackedCandidate& t2
                                         ){

  pat::CompositeCandidate psi2sCand;
  psi2sCand.setCharge(dimuonditrak.charge()+t1.charge()+t2.charge());

  double m_trak = pionmass;
  math::XYZVector mom_trak_1 = t1.momentum();
  math::XYZVector mom_trak_2 = t2.momentum();
  double e_trak_1 = sqrt(m_trak*m_trak + mom_trak_1.Mag2());
  double e_trak_2 = sqrt(m_trak*m_trak + mom_trak_2.Mag2());

  math::XYZTLorentzVector p4_trak_1 = math::XYZTLorentzVector(mom_trak_1.X(),mom_trak_1.Y(),mom_trak_1.Z(),e_trak_1);
  math::XYZTLorentzVector p4_trak_2 = math::XYZTLorentzVector(mom_trak_2.X(),mom_trak_2.Y(),mom_trak_2.Z(),e_trak_2);

  reco::Candidate::LorentzVector v = p4_trak_1 + p4_trak_2 + dimuonditrak.p4();

  psi2sCand.setP4(v);

  return psi2sCand;
}

reco::Candidate::LorentzVector FiveTracksProducerFit::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(FiveTracksProducerFit);
