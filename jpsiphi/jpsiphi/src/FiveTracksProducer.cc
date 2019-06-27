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

float FiveTracksProducer::DeltaPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
   return (fabs(t1.pt()-t2.pt())/t2.pt());
}

bool FiveTracksProducer::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<MaxDPtRel_ &&
	DeltaR(t1,t2) < MaxDeltaR_);
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
  DiMuonDiTrackCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiMuoDiTrack"))),
  TrackCollection_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
  TriggerCollection_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
  TrackPtCut_(iConfig.existsAs<double>("TrackPtCut") ? iConfig.getParameter<double>("TrackPtCut") : 0.8),
  MaxDeltaR_(iConfig.existsAs<double>("DRCut") ? iConfig.getParameter<double>("DRCut") : 0.01),
  MaxDPtRel_(iConfig.existsAs<double>("DPtCut") ? iConfig.getParameter<double>("DPtCut") : 2.0),
  TrackGenMap_(consumes<edm::Association<reco::GenParticleCollection>>(iConfig.getParameter<edm::InputTag>("TrackMatcher"))),
  thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
  thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
  TriggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  FiveTrackMassCuts_(iConfig.getParameter<std::vector<double>>("FiveTrackCuts")),
  HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters")),
  IsMC_(iConfig.getParameter<bool>("IsMC"))
{
  produces<pat::CompositeCandidateCollection>("FiveTracks");

  nevents = 0;

  kaonmass = 0.493677;
  pionmass = 0.13957061;
  protonmass = 0.93827208;
  psi2smass = 3.686093;
  jpsiMass = 3.096916;

  trackmass = kaonmass;

  ncomboneg = 0;
  ncombo = 0;
  ncomboneu = 0;

  pdgToMass[211] = pionmass;
  pdgToMass[2212] = protonmass;
  pdgToMass[321] = kaonmass;

  allMuons_ = consumes<pat::MuonCollection>((edm::InputTag)"slimmedMuons");

}

void FiveTracksProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::unique_ptr<pat::CompositeCandidateCollection> fiveCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(allMuons_,muons);

  edm::Handle<pat::CompositeCandidateCollection> dimuonditrack;
  iEvent.getByToken(DiMuonDiTrackCollection_,dimuonditrack);

  edm::Handle<edm::View<pat::PackedCandidate> > track;
  iEvent.getByToken(TrackCollection_,track);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trig;
  iEvent.getByToken(TriggerCollection_,trig);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( TriggerResults_ , triggerResults_handle);

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;
  reco::Vertex theBeamSpotV = reco::Vertex(bs.position(), bs.covariance3D());

  edm::Handle<reco::VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);


  // const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);


  float FiveTrackMassMax = FiveTrackMassCuts_[1];
  float FiveTrackMassMin = FiveTrackMassCuts_[0];

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

  FiveTrackMassMax = FiveTrackMassMax + 3*pionmass;
  FiveTrackMassMin = FiveTrackMassMin - 3*pionmass;


  edm::Handle<edm::Association<reco::GenParticleCollection>> theGenMap;
  iEvent.getByToken(TrackGenMap_,theGenMap);

  KinematicParticleFactoryFromTransientTrack pFactory;

  std::map< std::tuple <int,int,int> ,int> doneFlag;
  std::map<size_t,float> bestVertexPos, bestVertexNeg, bestVertexNeu;
  std::map<size_t,pat::CompositeCandidate> posCollection,negCollection,neuCollection;

  for (size_t d = 0; d < dimuonditrack->size(); d++) {

       auto dimuonditrackCand = dimuonditrack->at(d);

       if(dimuonditrackCand.userFloat("vProb")<0.01)
         continue;

       // const reco::Vertex thePrimaryV = *(dimuonditrackCand.userData<reco::Vertex>("bestPV"));
       // const reco::Vertex thePrimaryV = *dimuonditrackCand.userData<reco::Vertex>("PVwithmuons");
       const pat::CompositeCandidate * dimuon_cand = dynamic_cast <pat::CompositeCandidate *>(dimuonditrackCand.daughter("dimuon"));

       const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(dimuonditrackCand.daughter("dimuon")->daughter("highMuon"));
       const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(dimuonditrackCand.daughter("dimuon")->daughter("lowMuon"));
       // const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(pmu1->originalObject());
       // const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(pmu2->originalObject());


       //I want the positive and negative track to build psi2S with charge == 0
       const pat::PackedCandidate *tp = dynamic_cast <pat::PackedCandidate *>(dimuonditrackCand.daughter("ditrack")->daughter("highTrack"));
       const pat::PackedCandidate *tm = dynamic_cast <pat::PackedCandidate *>(dimuonditrackCand.daughter("ditrack")->daughter("lowTrack"));
       int tpId = dimuonditrackCand.userInt("pId");
       int tmId = dimuonditrackCand.userInt("mId");

       std::vector<double> oneMasses,twoMasses,threeMasses;
       oneMasses.push_back(pionmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(kaonmass); // p p k
       oneMasses.push_back(kaonmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(pionmass); // k p p
       oneMasses.push_back(pionmass);  twoMasses.push_back(kaonmass);  threeMasses.push_back(pionmass); // p k p
       oneMasses.push_back(pionmass);  twoMasses.push_back(pionmass);  threeMasses.push_back(pionmass); // p p p

       const ParticleMass muonMass(0.1056583);
       float muonSigma = muonMass*1E-6;
       const ParticleMass kaonMass(kaonmass);
       float kaonSigma = kaonMass*1E-6;

       //Adding the third track
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

       const unsigned int numMasses = 4; //(int) numMasses;

       for (size_t i = 0; i < track->size(); i++) {

         std::array<float,numMasses> fiveTracksMass;

         float minDR_third = 10000.0;
         float minDP_third = 10000.0;
         float minDPt_third = 10000.0;

         auto thirdTrack = track->at(i);


         if(thirdTrack.pt()<TrackPtCut_) continue;
         if(thirdTrack.charge() == 0) continue;
	       //if(!IsMC_ and fabs(thirdTrack.pdgId())!=211) continue;
	       if(!(thirdTrack.trackHighPurity())) continue;
         if(!(thirdTrack.hasTrackDetails())) continue;

         if (IsTheSame(thirdTrack,*tp) || int(i) == tpId) continue;
         if (IsTheSame(thirdTrack,*tm) || int(i) == tmId) continue;
         if (IsTheSame(thirdTrack,*pmu1) || IsTheSame(thirdTrack,*pmu2) ) continue;

         std::vector<int> ids;
         ids.push_back(i);ids.push_back(tpId);ids.push_back(tmId);
         std::sort(ids.begin(),ids.end());
         if(doneFlag.find(std::tuple<int,int,int>(ids[0],ids[1],ids[2]))!=doneFlag.end())
          continue;
         else
            doneFlag[std::tuple<int,int,int>(ids[0],ids[1],ids[2])] = 1.0;

         pat::CompositeCandidate fiveCand = makeFiveCandidate(dimuonditrackCand, thirdTrack);

         if (fiveCand.mass() > FiveTrackMassMax || fiveCand.mass() < FiveTrackMassMin) continue;

         double five_ma_fit = 14000.;
         double five_vp_fit = -9999.;
         double five_x2_fit = 10000.;
         double five_nd_fit = 10000.;

         // int numMasses = (int) extraMasses.size();

         std::vector<reco::TransientTrack> fiveTracks;
         std::vector<RefCountedKinematicParticle> kaonParticles, pionParticles;

         float kinChi = 0.;
         float kinNdf = 0.;

         fiveTracks.push_back((*theB).build(*(pmu1->innerTrack()))); // µ
         fiveTracks.push_back((*theB).build(*(pmu2->innerTrack()))); // µ
         fiveTracks.push_back((*theB).build(*(tp->bestTrack()))); // K/π
         fiveTracks.push_back((*theB).build(*(tm->bestTrack()))); // K/π
         fiveTracks.push_back((*theB).build(*(thirdTrack.bestTrack()))); // K/π

         kaonParticles.push_back(pFactory.particle(fiveTracks[0],muonMass,kinChi,kinNdf,muonSigma));
         kaonParticles.push_back(pFactory.particle(fiveTracks[1],muonMass,kinChi,kinNdf,muonSigma));
         kaonParticles.push_back(pFactory.particle(fiveTracks[2],kaonMass,kinChi,kinNdf,kaonSigma));
         kaonParticles.push_back(pFactory.particle(fiveTracks[3],kaonMass,kinChi,kinNdf,kaonSigma));
         kaonParticles.push_back(pFactory.particle(fiveTracks[4],kaonMass,kinChi,kinNdf,kaonSigma));

         KinematicConstrainedVertexFitter fiveFitter;
         MultiTrackKinematicConstraint *jpsi_mtc = new  TwoTrackMassKinematicConstraint(jpsiMass);
         RefCountedKinematicTree fiveVertexFitTree;
         fiveVertexFitTree = fiveFitter.fit(kaonParticles,jpsi_mtc);

         if (fiveVertexFitTree->isEmpty()) continue;

         fiveVertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle fitF = fiveVertexFitTree->currentParticle();
         RefCountedKinematicVertex fitFVertex = fiveVertexFitTree->currentDecayVertex();

         if (!(fitF->currentState().isValid())) continue;

         five_ma_fit = fitF->currentState().mass();
         five_x2_fit = fitFVertex->chiSquared();
         five_nd_fit = (double)(fitFVertex->degreesOfFreedom());
         five_vp_fit = ChiSquaredProbability(five_x2_fit,five_nd_fit);


         if(five_vp_fit < 0.001) continue;


         for (size_t ii = 0; ii < muons->size(); ii++)
         {
            auto thisMuon = muons->at(ii);

            float DeltaEta = fabs(thisMuon.eta()-thirdTrack.eta());
            float DeltaP   = fabs(thisMuon.p()-thirdTrack.p());
            float DeltaPt = ((thirdTrack.pt() - thisMuon.pt())/thirdTrack.pt());

            minDR_third = -std::max(-minDR_third,-DeltaEta);
            minDP_third = -std::max(-minDP_third,-DeltaP);
            minDPt_third = -std::max(-minDPt_third,-DeltaPt);
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
         std::array <float,4> cosAlpha, ctauPV, ctauErrPV, fromPV;
         // std::cout << "debug    9 "<< std::endl;
         TVector3 vtx, vdiff, pvtx;
         VertexDistanceXY vdistXY;
         reco::Vertex thePrimaryVDZ, thePrimaryZero, thePrimaryVCA;
         TwoTrackMinimumDistance ttmd;

         bool status = ttmd.calculate( GlobalTrajectoryParameters(
           GlobalPoint(five_vx_fit,five_vy_fit,five_vz_fit),
           GlobalVector(five_px_fit,five_py_fit,five_pz_fit),TrackCharge(0),&(*magneticField)),
           GlobalTrajectoryParameters(
             GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
             GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(0),&(*magneticField)));
         float extrapZ=-9E20;

         if (status) extrapZ=ttmd.points().first.z();

         vtx.SetXYZ(five_vx_fit,five_vy_fit,0);
         TVector3 pperp(five_px_fit, five_py_fit, 0);
         AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);

         reco::VertexCollection verteces;
         std::vector<int> vKeys;
         verteces.push_back(theBeamSpotV);
         vKeys.push_back(0);

         thePrimaryZero = reco::Vertex(*(priVtxs->begin()));
         verteces.push_back(thePrimaryZero);
         vKeys.push_back(0);

         float minDz = 999999.;
         double maxCosAlpha = -1.0;
         if ( (priVtxs->begin() == priVtxs->end()) )
         {
           // std::cout << "debug    10 "<< std::endl;
           thePrimaryVCA = reco::Vertex(*(priVtxs->begin()));
           thePrimaryVDZ = reco::Vertex(*(priVtxs->begin()));
           verteces.push_back(thePrimaryVCA);
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




         //////////////////////////////////////////////////
         //Refit PVs (not BS)

         for(int i = 0; i < numVertex; i++)
         {

           ctauPV[i] = (-1000.0);
           ctauErrPV[i] = (-1000.0);
           cosAlpha[i] = (-1000.0);
           fromPV[i] = (-1000.0);
         }

         for(size_t i = 1; i < verteces.size(); i++)
         {
           fromPV[i] = thirdTrack.fromPV(vKeys[i]);
         }

         // std::cout << "debug    13 "<< std::endl;
         for(size_t i = 0; i < verteces.size(); i++)
         {
           pvtx.SetXYZ(verteces[i].position().x(),verteces[i].position().y(),0);
           vdiff = vtx - pvtx;
           cosAlpha[i] = (vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp()));
           Measurement1D distXY = vdistXY.distance(reco::Vertex(*fitFVertex), verteces[i]);
           ctauPV[i] = (distXY.value()*cosAlpha[i] * five_ma_fit/pperp.Perp());
           GlobalError v1e = (reco::Vertex(*fitFVertex)).error();
           GlobalError v2e = verteces[i].error();
           AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
           ctauErrPV[i] = (sqrt(ROOT::Math::Similarity(vpperp,vXYe))*five_ma_fit/(pperp.Perp2()));
         }


         fiveCand.addUserFloat("mass_ref_0",five_ma_fit);
         fiveCand.addDaughter(fiveCand_rf,"ref_cand");

         fiveCand.addUserFloat("thirdTrackMuonDR",minDR_third);
         fiveCand.addUserFloat("thirdTrackMuonDP",minDP_third);
         fiveCand.addUserFloat("thirdTrackMuonDPt",minDPt_third);

         auto thisFive = makeFiveCandidateMixed(*dimuon_cand, *tp, *tm, thirdTrack,kaonmass,kaonmass,kaonmass);

         // for(size_t j = 0; j<numMasses;j++)
         //  fiveTracksMass[j] = makeFiveCandidateMixed(*dimuon_cand, *tp, *tm, thirdTrack,oneMasses[j] ,twoMasses[j] ,threeMasses[j]).mass();

             fiveCand.addUserInt("dimuontt_index",int(d));
             fiveCand.addUserInt("pId",tpId);
             fiveCand.addUserInt("mId",tmId);
             fiveCand.addUserInt("tId",i);

             std::string name;
             for(size_t j = 1; j<numMasses+1;j++)
             {
              name = "mass_ref_" + std::to_string(j);
              fiveCand.addUserFloat(name,fiveTracksMass[j-1]);
             }

             fiveCand.addDaughter(dimuonditrackCand,"dimuonditrack");
             fiveCand.addDaughter(*dimuon_cand,"dimuon");

             fiveCand.addDaughter(*tp,"trackOne");
             fiveCand.addDaughter(*tm,"trackTwo");
             fiveCand.addDaughter(thirdTrack,"trackThree");


             fiveCand.addUserData("bestPV",reco::Vertex(thePrimaryZero));
             fiveCand.addUserData("cosPV",reco::Vertex(thePrimaryVCA));
             fiveCand.addUserData("zPV",reco::Vertex(thePrimaryVDZ));
             fiveCand.addUserData("bS",reco::Vertex(theBeamSpotV));

             fiveCand.addUserFloat("vtxX",five_vx_fit);
             fiveCand.addUserFloat("vtxY",five_vy_fit);
             fiveCand.addUserFloat("vtxZ",five_vz_fit);


             fiveCand.addUserFloat("cosAlphaBS",cosAlpha[0]);
             fiveCand.addUserFloat("ctauPVBS",ctauPV[0]);
             fiveCand.addUserFloat("ctauErrPVBS",ctauErrPV[0]);

             // fiveCand.addUserFloat("tTFromPVBS",float(fromPV[0]));


             fiveCand.addUserFloat("cosAlpha",cosAlpha[1]);
             fiveCand.addUserFloat("ctauPV",ctauPV[1]);
             fiveCand.addUserFloat("ctauErrPV",ctauErrPV[1]);

             fiveCand.addUserFloat("tTFromPV",float(fromPV[1]));

             fiveCand.addUserFloat("cosAlphaCA",cosAlpha[2]);
             fiveCand.addUserFloat("ctauPVCA",ctauPV[2]);
             fiveCand.addUserFloat("ctauErrPVCA",ctauErrPV[2]);

             fiveCand.addUserFloat("tTFromPVCA",float(fromPV[2]));

             fiveCand.addUserFloat("cosAlphaDZ",cosAlpha[3]);
             fiveCand.addUserFloat("ctauPVDZ",ctauPV[3]);
             fiveCand.addUserFloat("ctauErrPVDZ",ctauErrPV[3]);

             fiveCand.addUserFloat("tTFromPVDZ",float(fromPV[3]));
             ///DCA
             std::vector<float> DCAs;
             TrajectoryStateClosestToPoint TS1 = fiveTracks[fiveTracks.size()-1].impactPointTSCP();

             if (TS1.isValid())
             {
                 for(size_t j = 0; j < fiveTracks.size() - 1;++j)
                 {
                   TrajectoryStateClosestToPoint TS2 = fiveTracks[j].impactPointTSCP();
                   float dca = 1E20;
                   if (TS1.isValid() && TS2.isValid()) {
                     ClosestApproachInRPhi cApp;
                     cApp.calculate(TS1.theState(), TS2.theState());
                     if (cApp.status() ) dca = cApp.distance();
                   }
                   DCAs.push_back(dca);
                 }
             }



             fiveCand.addUserFloat("vProb",five_vp_fit);
             fiveCand.addUserFloat("nDof",five_nd_fit);
             fiveCand.addUserFloat("vChi2",five_x2_fit);

             fiveCand.addUserFloat("dca_m1t3",DCAs[0]);
             fiveCand.addUserFloat("dca_m2t3",DCAs[1]);
             fiveCand.addUserFloat("dca_t1t3",DCAs[2]);
             fiveCand.addUserFloat("dca_t2t3",DCAs[3]);

             fiveCand.addDaughter(thisFive,"first_five_ref");

             fiveCand.addUserInt("thirdKaonMatch",filters[i]);


             if(IsMC_)
              {
                float hasGen = -1.0;

                if(theGenMap.isValid())
                {
                  //posTrack
                  auto ref = track->refAt(i);

                  if(theGenMap->contains(ref.id()))
                  {
                   if(((*theGenMap)[edm::Ref<edm::View<pat::PackedCandidate>>(track, i)]).isNonnull())
                   {
                     auto genP = ((*theGenMap)[edm::Ref<edm::View<pat::PackedCandidate>>(track, i)]);
                     hasGen = 1.0;
                     fiveCand.addDaughter(*genP,"thirdTrackGen");

                    }
                  }

                }
                // if(hasHighGen * hasLowGen >= 0.0)
                //   std::cout << "Has some gen ref " << std::endl;
                fiveCand.addUserFloat("hasThirdGen",hasGen);


              }

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
                                          const pat::CompositeCandidate& dimuonditrack,
                                          const pat::PackedCandidate& track
                                         ){

  pat::CompositeCandidate fiveCand;
  // fiveCand.addDaughter(dimuonditrack,"dimuonditrack");
  // fiveCand.addDaughter(track,"third");
  fiveCand.setCharge(dimuonditrack.charge()+track.charge());

  double m_track = trackmass;
  math::XYZVector mom_track = track.momentum();
  double e_track = sqrt(m_track*m_track + mom_track.Mag2());
  math::XYZTLorentzVector p4_track = math::XYZTLorentzVector(mom_track.X(),mom_track.Y(),mom_track.Z(),e_track);

  reco::Candidate::LorentzVector v = p4_track + dimuonditrack.p4();

  fiveCand.setP4(v);

  return fiveCand;
}

pat::CompositeCandidate FiveTracksProducer::makeFiveCandidateMixed(
                                          const pat::CompositeCandidate& dimuon,
                                          const pat::PackedCandidate& trackP,
                                          const pat::PackedCandidate& trackN,
                                          const pat::PackedCandidate& track3,
                                          double massOne,
                                          double massTwo,
                                          double massThree
                                         ){

  pat::CompositeCandidate fiveCand, trackOne, trackTwo, trackThree;
  pat::CompositeCandidate dimuonDiTrackOne, dimuonDiTrackTwo, dimuonDiTrackThree;
  pat::CompositeCandidate triTrack;

  fiveCand.addDaughter(dimuon,"dimuon");

  triTrack.setCharge(trackP.charge()+trackN.charge()+track3.charge());
  fiveCand.setCharge(dimuon.charge()+trackP.charge()+trackN.charge()+track3.charge());

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

  fiveCand.addDaughter(trackOne,"trackOne");
  fiveCand.addDaughter(trackTwo,"trackTwo");
  fiveCand.addDaughter(trackThree,"trackThree");

  dimuonDiTrackOne     = makeDimuonDiTrackCandidate(dimuon,trackOne,trackTwo);
  dimuonDiTrackTwo     = makeDimuonDiTrackCandidate(dimuon,trackOne,trackThree);
  dimuonDiTrackThree   = makeDimuonDiTrackCandidate(dimuon,trackTwo,trackThree);

  fiveCand.addDaughter(dimuonDiTrackOne,"dimuonDiTrackOne");
  fiveCand.addDaughter(dimuonDiTrackTwo,"dimuonDiTrackTwo");
  fiveCand.addDaughter(dimuonDiTrackThree,"dimuonDiTrackThree");

  reco::Candidate::LorentzVector v = p4_trackP + p4_trackN + p4_track3 + dimuon.p4();
  reco::Candidate::LorentzVector vT = trackP.p4() + trackN.p4() + track3.p4();

  triTrack.setP4(vT);

  fiveCand.addDaughter(triTrack,"triTrack");

  fiveCand.setP4(v);

  return fiveCand;
}

pat::CompositeCandidate FiveTracksProducer::makeFiveCandidateMixed(
                                          const pat::CompositeCandidate& dimuon,
                                          const pat::CompositeCandidate& trackP,
                                          const pat::CompositeCandidate& trackN,
                                          const pat::CompositeCandidate& track3
                                         ){

  pat::CompositeCandidate fiveCand, trackOne, trackTwo, trackThree;
  pat::CompositeCandidate dimuonDiTrackOne, dimuonDiTrackTwo, dimuonDiTrackThree;
  pat::CompositeCandidate triTrack;

  // fiveCand.addDaughter(dimuon,"dimuon");

  fiveCand.setCharge(dimuon.charge()+trackP.charge()+trackN.charge()+track3.charge());
  // triTrack.setCharge(trackP.charge()+trackN.charge()+track3.charge());
  //
  // fiveCand.addDaughter(trackP,"trackOne");
  // fiveCand.addDaughter(trackN,"trackTwo");
  // fiveCand.addDaughter(track3,"trackThree");
  //
  // dimuonDiTrackOne     = makeDimuonDiTrackCandidate(dimuon,trackP,trackN);
  // dimuonDiTrackTwo     = makeDimuonDiTrackCandidate(dimuon,trackP,track3);
  // dimuonDiTrackThree   = makeDimuonDiTrackCandidate(dimuon,trackN,track3);
  //
  // fiveCand.addDaughter(dimuonDiTrackOne,"dimuonDiTrackOne");
  // fiveCand.addDaughter(dimuonDiTrackTwo,"dimuonDiTrackTwo");
  // fiveCand.addDaughter(dimuonDiTrackThree,"dimuonDiTrackThree");

  reco::Candidate::LorentzVector v  = trackP.p4() + trackN.p4() + track3.p4() + dimuon.p4();
  // reco::Candidate::LorentzVector vT = trackP.p4() + trackN.p4() + track3.p4();

  // triTrack.setP4(vT);

  // fiveCand.addDaughter(triTrack,"triTrack");

  fiveCand.setP4(v);

  return fiveCand;
}

pat::CompositeCandidate FiveTracksProducer::makeDimuonDiTrackCandidate(
                                          const pat::CompositeCandidate& dimuon,
                                          const pat::CompositeCandidate& t1,
                                          const pat::CompositeCandidate& t2
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

reco::Candidate::LorentzVector FiveTracksProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(FiveTracksProducer);
