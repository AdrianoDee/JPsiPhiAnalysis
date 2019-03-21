#include "../interface/SixTracksProducer.h"
#include <tuple>
#include <map>


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

float SixTracksProducer::DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
   float p1 = t1.phi();
   float p2 = t2.phi();
   float e1 = t1.eta();
   float e2 = t2.eta();
   auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

float SixTracksProducer::DeltaPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
   return (fabs(t1.pt()-t2.pt())/t2.pt());
}

bool SixTracksProducer::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<MaxDPtRel_ &&
	DeltaR(t1,t2) < MaxDeltaR_);
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
  FiveTrackCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("FiveCollection"))),
  TrackCollection_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
  TriggerCollection_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
  trackPtCut_(iConfig.existsAs<double>("TrackPtCut") ? iConfig.getParameter<double>("TrackPtCut") : 0.8),
  MaxDeltaR_(iConfig.existsAs<double>("DRCut") ? iConfig.getParameter<double>("DRCut") : 0.01),
  MaxDPtRel_(iConfig.existsAs<double>("DPtCut") ? iConfig.getParameter<double>("DPtCut") : 2.0),
  TrackGenMap_(consumes<edm::Association<reco::GenParticleCollection>>(iConfig.getParameter<edm::InputTag>("TrackMatcher"))),
  thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
  thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
  TriggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  SixTrackMassCuts_(iConfig.getParameter<std::vector<double>>("SixTrackCuts")),
  HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters")),
  IsMC_(iConfig.getParameter<bool>("isMC"))
{
  produces<pat::CompositeCandidateCollection>("SixTracks");

  nevents = 0;

  kaonmass = 0.493677;
  pionmass = 0.13957061;
  psi2smass = 3.686093;

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

  edm::Handle<pat::CompositeCandidateCollection> fivetrack;
  iEvent.getByToken(FiveTrackCollection_,fivetrack);

  edm::Handle<edm::View<pat::PackedCandidate> > track;
  iEvent.getByToken(TrackCollection_,track);

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;
  reco::Vertex theBeamSpotV = reco::Vertex(bs.position(), bs.covariance3D());

  edm::Handle<reco::VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trig;
  iEvent.getByToken(TriggerCollection_,trig);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( TriggerResults_ , triggerResults_handle);

  // const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);


  float SixTrackMassMax = SixTrackMassCuts_[1];
  float SixTrackMassMin = SixTrackMassCuts_[0];


  // reco::TrackCollection allTheTracks;
  // for (size_t i = 0; i < track->size(); i++)
  // {
  //   auto t = track->at(i);
  //   if(t.pt()<0.5) continue;
  //   if(!(t.hasTrackDetails())) continue;
  //   allTheTracks.push_back(*(t.bestTrack()));
  //
  // }

  pat::TriggerObjectStandAloneCollection filteredColl;
  std::map<int,pat::TriggerObjectStandAlone> matchedColl;
  std::map<size_t,double> trackDeltaR,trackDeltaPt;
  std::vector < UInt_t > filterResults;
  std::map<int,UInt_t> filters;

  //Trigger Collections
  for ( size_t iTrigObj = 0; iTrigObj < trig->size(); ++iTrigObj ) {

    pat::TriggerObjectStandAlone unPackedTrigger( trig->at( iTrigObj ) );

    unPackedTrigger.unpackPathNames( names );
    unPackedTrigger.unpackFilterLabels(iEvent,*triggerResults_handle);

    bool filtered = false;
    UInt_t thisFilter = 0;

    for (size_t i = 0; i < HLTFilters_.size(); i++)
      if(unPackedTrigger.hasFilterLabel(HLTFilters_[i]))
        {
          thisFilter += (1<<i);
          filtered = true;
        }

    if(filtered)
    {
      filteredColl.push_back(unPackedTrigger);
      filterResults.push_back(thisFilter);
    }
  }

  //Tracks Collections Trigger Matching
  for (size_t i = 0; i < track->size(); i++) {

    auto t = track->at(i);

    bool matched = false;
    for (std::vector<pat::TriggerObjectStandAlone>::const_iterator trigger = filteredColl.begin(), triggerEnd=filteredColl.end(); trigger!= triggerEnd; ++trigger)
  for ( size_t iTrigObj = 0; iTrigObj < filteredColl.size(); ++iTrigObj )
    {
      auto thisTrig = filteredColl.at(iTrigObj);
      if(MatchByDRDPt(t,filteredColl[iTrigObj]))
      {
        if(matched)
        {
          if(trackDeltaR[i] > DeltaR(t,thisTrig))
          {
            filters[i] = filterResults[iTrigObj];
            matchedColl[i] = thisTrig;
            trackDeltaR[i] = fabs(DeltaR(t,thisTrig));
            trackDeltaPt[i] = fabs(DeltaPt(t,thisTrig));
          }
        }else
        {
          filters[i] = filterResults[iTrigObj];
          matchedColl[i] = thisTrig;
          trackDeltaR[i] = fabs(DeltaR(t,thisTrig));
          trackDeltaPt[i] = fabs(DeltaPt(t,thisTrig));
        }

        matched = true;
      }
    }
    if(!matched)
    {
      filters[i] = 0;
      trackDeltaR[i] = -1.0;
      trackDeltaPt[i] = -1.0;
    }

  }

  edm::Handle<edm::Association<reco::GenParticleCollection>> theGenMap;
  iEvent.getByToken(TrackGenMap_,theGenMap);

  KinematicParticleFactoryFromTransientTrack pFactory;

  std::map< std::tuple <int,int,int,int> ,bool> doneFlag;
  std::map<size_t,float> bestVertexPos, bestVertexNeg, bestVertexNeu;
  std::map<size_t,pat::CompositeCandidate> posCollection,negCollection,neuCollection;

  const int numMasses = 6;//numMasses;

  for (size_t d = 0; d < fivetrack->size(); d++) {

       auto fivetrackCand = fivetrack->at(d);

       if(fivetrackCand.userFloat("vProb")<0.001)
         continue;
       // const reco::Vertex thePrimaryV = *(fivetrackCand.userData<reco::Vertex>("bestPV"));
       // const reco::Vertex thePrimaryV = *fivetrackCand.userData<reco::Vertex>("PVwithmuons");
       const pat::CompositeCandidate *dimuonditrack_cand = dynamic_cast <const pat::CompositeCandidate *>(fivetrackCand.daughter("dimuonditrack"));
       const pat::CompositeCandidate *dimuon_cand = dynamic_cast <const pat::CompositeCandidate *>(dimuonditrack_cand->daughter("dimuon"));

       const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(dimuonditrack_cand->daughter("dimuon")->daughter("highMuon"));
       const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(dimuonditrack_cand->daughter("dimuon")->daughter("lowMuon"));
       // const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(pmu1->originalObject());
       // const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(pmu2->originalObject());


       //I want the positive and negative track to build psi2S with charge == 0
       const pat::PackedCandidate *tp = dynamic_cast <const pat::PackedCandidate *>(dimuonditrack_cand->daughter("ditrack")->daughter("highTrack"));
       const pat::PackedCandidate *tm = dynamic_cast <const pat::PackedCandidate *>(dimuonditrack_cand->daughter("ditrack")->daughter("lowTrack"));
       const pat::PackedCandidate *tt = dynamic_cast <const pat::PackedCandidate *>(fivetrackCand.daughter("fifth"));

       int tpId = fivetrackCand.userInt("pId");
       int tmId = fivetrackCand.userInt("mId");
       int ttId = fivetrackCand.userInt("tId");

       std::array< std::array <float,4>, numMasses-1> trackMasses;
       trackMasses[0][0] = kaonmass; trackMasses[0][1] = pionmass; trackMasses[0][2] = kaonmass; trackMasses[0][3] = pionmass; // k p k p
       trackMasses[1][0] = kaonmass; trackMasses[1][1] = pionmass; trackMasses[1][2] = pionmass; trackMasses[1][3] = kaonmass; // k p k p
       trackMasses[2][0] = pionmass; trackMasses[2][1] = kaonmass; trackMasses[2][2] = kaonmass; trackMasses[2][3] = pionmass; // k p k p
       trackMasses[3][0] = pionmass; trackMasses[3][1] = kaonmass; trackMasses[3][2] = pionmass; trackMasses[3][3] = kaonmass; // k p k p
       trackMasses[4][0] = pionmass; trackMasses[4][1] = pionmass; trackMasses[4][2] = kaonmass; trackMasses[4][3] = kaonmass; // k p k p

       int mmtt_id = fivetrackCand.userInt("dimuontt_index");
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

       for (size_t i = 0; i < track->size(); i++) {

         std::vector<float> sixTracksMass;
         std::vector<float> psi2sOne, psi2sTwo;
         std::vector<float> fiveTracksCTau, fiveTracksCTauErr, fiveTracksCosAlpha;

         float minDR_fourth = 10000.0;
         float minDP_fourth = 10000.0;
         float minDPt_fourth = 10000.0;

         auto sixthTrack = track->at(i);

         if(sixthTrack.pt()<trackPtCut_) continue;
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

         pat::CompositeCandidate sixCand = makeSixCandidate(fivetrackCand, sixthTrack);

         if (sixCand.mass() > SixTrackMassMin || sixCand.mass() < SixTrackMassMax)
         continue;

         sixCand.addUserFloat("sixCandMass",sixCand.mass());

         double six_ma_fit = 14000.;
         double six_vp_fit = -9999.;
         double six_x2_fit = 10000.;
         double six_nd_fit = 10000.;

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


         int numVertex = 6;
         std::array <float,6> cosAlpha, ctauPV, ctauErrPV, fromPV;

         for(int i = 0; i < numVertex; i++)
         {

           ctauPV[i] = (-1000.0);
           ctauErrPV[i] = (-1000.0);
           cosAlpha[i] = (-1000.0);
           fromPV[i] = (-1000.0);
         }

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

         ///DCA
         std::vector<float> DCAs;

         TrajectoryStateClosestToPoint TS1 = sixTracks[sixTracks.size()-1].impactPointTSCP();

         if (TS1.isValid())
         {
             for(size_t j = 0; j < sixTracks.size() - 1;++j)
             {
               TrajectoryStateClosestToPoint TS2 = sixTracks[j].impactPointTSCP();
               float dca = 1E20;
               if (TS1.isValid() && TS2.isValid()) {
                 ClosestApproachInRPhi cApp;
                 cApp.calculate(TS1.theState(), TS2.theState());
                 if (cApp.status() ) dca = cApp.distance();
               }
               DCAs.push_back(dca);
             }
         }

         std::vector<double> oneMasses, twoMasses, threeMasses, fourMasses;
         oneMasses.push_back(pionmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(kaonmass); fourMasses.push_back(kaonmass); // p p k k

         oneMasses.push_back(pionmass);  twoMasses.push_back(kaonmass);  threeMasses.push_back(pionmass); fourMasses.push_back(kaonmass); // p k p k
         oneMasses.push_back(pionmass);  twoMasses.push_back(kaonmass);  threeMasses.push_back(kaonmass); fourMasses.push_back(pionmass); // p k k p

         oneMasses.push_back(kaonmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(kaonmass); fourMasses.push_back(pionmass); // k p k p
         oneMasses.push_back(kaonmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(pionmass); fourMasses.push_back(kaonmass); // k p p k

         oneMasses.push_back(kaonmass);  twoMasses.push_back(kaonmass);  threeMasses.push_back(pionmass); fourMasses.push_back(pionmass); // k k p p

         auto thisSix = makeSixCandidateMixed(*dimuon_cand, *tp, *tm, *tt,sixthTrack,kaonmass,kaonmass,kaonmass,kaonmass);

         for(size_t j = 0; j<numMasses;j++)
          sixTracksMass[j] = makeSixCandidateMixed(*dimuon_cand, *tp, *tm, *tt,sixthTrack,oneMasses[j] ,twoMasses[j] ,threeMasses[j],fourMasses[j]).mass();

          sixCand.addUserInt("five_index",int(d));
          sixCand.addUserInt("pId",tpId);
          sixCand.addUserInt("mId",tmId);
          sixCand.addUserInt("tId",ttId);
          sixCand.addUserInt("fId",mmtt_id);

          std::string name;
          for(size_t j = 1; j<numMasses+1;j++)
          {
           name = "mass_ref_" + std::to_string(j);
           sixCand.addUserFloat(name,sixTracksMass[j-1]);
          }

          sixCand.addUserInt("five_id",int(d));

          sixCand.addDaughter(*tp,"trackOne");
          sixCand.addDaughter(*tm,"trackTwo");
          sixCand.addDaughter(*tm,"trackThree");
          sixCand.addDaughter(sixthTrack,"trackFour");

          sixCand.addDaughter(sixCand_rf,"ref_cand");
          sixCand.addDaughter(thisSix,"first_six_ref");

          sixCand.addUserFloat("mass_ref_0",six_ma_fit);

          sixCand.addUserFloat("fourthTrackMuonDR",minDR_fourth);
          sixCand.addUserFloat("fourthTrackMuonDP",minDP_fourth);
          sixCand.addUserFloat("fourthTrackMuonDPt",minDPt_fourth);

          sixCand.addUserData("bestPV",reco::Vertex(thePrimaryZero));
          sixCand.addUserData("cosPV",reco::Vertex(thePrimaryVCA));
          sixCand.addUserData("zPV",reco::Vertex(thePrimaryVDZ));
          sixCand.addUserData("bs",reco::Vertex(thePrimaryV));

          sixCand.addUserFloat("vtxX",six_vx_fit);
          sixCand.addUserFloat("vtxY",six_vy_fit);
          sixCand.addUserFloat("vtxZ",six_vz_fit);

          sixCand.addUserFloat("cosAlphaBS",cosAlpha[0]);
          sixCand.addUserFloat("ctauPVBS",ctauPV[0]);
          sixCand.addUserFloat("ctauErrPVBS",ctauErrPV[0]);
          sixCand.addUserFloat("tFFromPVBS",float(fromPV[0]));

          sixCand.addUserFloat("cosAlpha",cosAlpha[1]);
          sixCand.addUserFloat("ctauPV",ctauPV[1]);
          sixCand.addUserFloat("ctauErrPV",ctauErrPV[1]);
          sixCand.addUserFloat("tFFromPV",float(fromPV[1]));

          sixCand.addUserFloat("cosAlpha_alpha",cosAlpha[2]);
          sixCand.addUserFloat("ctauPV_alpha",ctauPV[2]);
          sixCand.addUserFloat("ctauErrPV_alpha",ctauErrPV[2]);
          sixCand.addUserFloat("tFFromPV_alpha",float(fromPV[2]));

          sixCand.addUserFloat("cosAlphaDZ",cosAlpha[3]);
          sixCand.addUserFloat("ctauPVDZ",ctauPV[3]);
          sixCand.addUserFloat("ctauErrPVDZ",ctauErrPV[3]);
          sixCand.addUserFloat("tFFromPVDZ",float(fromPV[3]));

          sixCand.addUserInt("fourthKaonMatch",filters[i]);

          sixCand.addUserFloat("dca_m1t4",DCAs[0]);
          sixCand.addUserFloat("dca_m2t4",DCAs[1]);
          sixCand.addUserFloat("dca_t1t4",DCAs[2]);
          sixCand.addUserFloat("dca_t2t4",DCAs[3]);
          sixCand.addUserFloat("dca_t3t4",DCAs[4]);

          sixCand.addUserFloat("vProb",six_vp_fit);
          sixCand.addUserFloat("nDof",six_nd_fit);
          sixCand.addUserFloat("vChi2",six_x2_fit);

          if(IsMC_)
           {
             float hasGen = -1.0;

             if(theGenMap.isValid())
             {
               //
               auto ref = track->refAt(i);

               if(theGenMap->contains(ref.id()))
               {
                if(((*theGenMap)[edm::Ref<edm::View<pat::PackedCandidate>>(track, i)]).isNonnull())
                {
                  auto genP = ((*theGenMap)[edm::Ref<edm::View<pat::PackedCandidate>>(track, i)]);
                  hasGen = 1.0;
                  sixCand.addDaughter(*genP,"fourthTrackGen");

                 }
               }

             }
             sixCand.addUserFloat("hasFourthGen",hasGen);
           }

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
                                          const pat::CompositeCandidate& fivetrack,
                                          const pat::PackedCandidate& track
                                         ){

  pat::CompositeCandidate sixCand;
  sixCand.addDaughter(fivetrack,"fiveCand");
  sixCand.addDaughter(track,"fourthTrack");
  sixCand.setCharge(fivetrack.charge()+track.charge());

  double m_track = trackmass;
  math::XYZVector mom_track = track.momentum();
  double e_track = sqrt(m_track*m_track + mom_track.Mag2());
  math::XYZTLorentzVector p4_track = math::XYZTLorentzVector(mom_track.X(),mom_track.Y(),mom_track.Z(),e_track);

  reco::Candidate::LorentzVector v = p4_track + fivetrack.p4();

  sixCand.setP4(v);

  return sixCand;
}

pat::CompositeCandidate SixTracksProducer::makeSixCandidateMixed(
                                          const pat::CompositeCandidate& dimuon,
                                          const pat::PackedCandidate& trackP,
                                          const pat::PackedCandidate& trackN,
                                          const pat::PackedCandidate& track3,
                                          const pat::PackedCandidate& track4,
                                          double massOne,
                                          double massTwo,
                                          double massThree,
                                          double massFour
                                         ){

  pat::CompositeCandidate sixCand, trackOne, trackTwo, trackThree, trackFour;

  pat::CompositeCandidate fiveTrackOne, fiveTrackTwo, fiveTrackThree, fiveTrackFour;

  sixCand.setCharge(dimuon.charge()+trackP.charge()+trackN.charge()+track3.charge());

  math::XYZVector mom_trackP = trackP.momentum();
  double e_trackP = sqrt(massOne*massOne + mom_trackP.Mag2());
  math::XYZTLorentzVector p4_trackP = math::XYZTLorentzVector(mom_trackP.X(),mom_trackP.Y(),mom_trackP.Z(),e_trackP);
  trackOne.setCharge(trackP.charge());
  trackOne.setP4(p4_trackP);

  math::XYZVector mom_trackN = trackN.momentum();
  double e_trackN = sqrt(massTwo*massTwo + mom_trackN.Mag2());
  math::XYZTLorentzVector p4_trackN = math::XYZTLorentzVector(mom_trackN.X(),mom_trackN.Y(),mom_trackN.Z(),e_trackN);
  trackTwo.setCharge(trackN.charge());
  trackTwo.setP4(p4_trackN);

  math::XYZVector mom_track3 = track3.momentum();
  double e_track3 = sqrt(massThree*massThree + mom_track3.Mag2());
  math::XYZTLorentzVector p4_track3 = math::XYZTLorentzVector(mom_track3.X(),mom_track3.Y(),mom_track3.Z(),e_track3);
  trackThree.setCharge(track3.charge());
  trackThree.setP4(p4_track3);

  math::XYZVector mom_track4 = track4.momentum();
  double e_track4 = sqrt(massThree*massThree + mom_track4.Mag2());
  math::XYZTLorentzVector p4_track4 = math::XYZTLorentzVector(mom_track4.X(),mom_track4.Y(),mom_track4.Z(),e_track4);
  trackFour.setCharge(track4.charge());
  trackFour.setP4(p4_track4);

  fiveTrackOne     = makeFiveCandidateMixed(dimuon,trackP,trackN,track3);
  fiveTrackTwo     = makeFiveCandidateMixed(dimuon,trackP,trackN,track4);
  fiveTrackThree     = makeFiveCandidateMixed(dimuon,trackP,track3,track4);
  fiveTrackFour     = makeFiveCandidateMixed(dimuon,trackN,track3,track4);

  // P N 3
  // P N 4
  // P 3 4
  // N 3 4
  //sixCand.addDaughter(fiveTrackOne,"fiveTrackOne"); is already there
  sixCand.addDaughter(fiveTrackTwo,"fiveTrackTwo");
  sixCand.addDaughter(fiveTrackThree,"fiveTrackThree");
  sixCand.addDaughter(fiveTrackFour,"fiveTrackFour");

  reco::Candidate::LorentzVector v = p4_trackP + p4_trackN + p4_track3 + p4_track4 + dimuon.p4();

  sixCand.setP4(v);

  return sixCand;
}


reco::Candidate::LorentzVector SixTracksProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}

pat::CompositeCandidate SixTracksProducer::makeFiveCandidateMixed(
                                          const pat::CompositeCandidate& dimuon,
                                          const pat::PackedCandidate& trackP,
                                          const pat::PackedCandidate& trackN,
                                          const pat::PackedCandidate& track3
                                         ){

  pat::CompositeCandidate fiveCand, trackOne, trackTwo, trackThree;
  pat::CompositeCandidate dimuonDiTrackOne, dimuonDiTrackTwo, dimuonDiTrackThree;
  pat::CompositeCandidate triTrack;

  fiveCand.addDaughter(dimuon,"dimuon");

  fiveCand.setCharge(dimuon.charge()+trackP.charge()+trackN.charge()+track3.charge());
  triTrack.setCharge(trackP.charge()+trackN.charge()+track3.charge());

  fiveCand.addDaughter(trackP,"trackOne");
  fiveCand.addDaughter(trackN,"trackTwo");
  fiveCand.addDaughter(track3,"trackThree");

  dimuonDiTrackOne     = makeDimuonDiTrackCandidate(dimuon,trackP,trackN);
  dimuonDiTrackTwo     = makeDimuonDiTrackCandidate(dimuon,trackP,track3);
  dimuonDiTrackThree   = makeDimuonDiTrackCandidate(dimuon,trackN,track3);

  fiveCand.addDaughter(dimuonDiTrackOne,"dimuonDiTrackOne");
  fiveCand.addDaughter(dimuonDiTrackTwo,"dimuonDiTrackTwo");
  fiveCand.addDaughter(dimuonDiTrackThree,"dimuonDiTrackThree");

  reco::Candidate::LorentzVector v  = trackP.p4() + trackN.p4() + track3.p4() + dimuon.p4();
  reco::Candidate::LorentzVector vT = trackP.p4() + trackN.p4() + track3.p4();

  triTrack.setP4(vT);

  fiveCand.addDaughter(triTrack,"triTrack");

  fiveCand.setP4(v);

  return fiveCand;
}

pat::CompositeCandidate SixTracksProducer::makeDimuonDiTrackCandidate(
                                          const pat::CompositeCandidate& dimuon,
                                          const pat::PackedCandidate& t1,
                                          const pat::PackedCandidate& t2
                                         ){

  pat::CompositeCandidate psi2sCand, ditrack;
  psi2sCand.setCharge(dimuon.charge()+t1.charge()+t2.charge());
  ditrack.setCharge(t1.charge()+t2.charge());
  psi2sCand.addDaughter(dimuon,"dimuon");
  psi2sCand.addDaughter(t1,"trackOne");
  psi2sCand.addDaughter(t2,"trackTwo");

  reco::Candidate::LorentzVector v  = t1.p4() + t2.p4() + dimuon.p4();
  reco::Candidate::LorentzVector vT = t1.p4() + t2.p4();

  ditrack.setP4(vT);

  psi2sCand.addDaughter(ditrack,"ditrack");

  psi2sCand.setP4(v);

  return psi2sCand;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SixTracksProducer);
