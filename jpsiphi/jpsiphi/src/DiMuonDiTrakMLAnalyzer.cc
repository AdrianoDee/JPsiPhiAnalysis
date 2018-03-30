// -*- C++ -*-
//
// Package:    DiMuonDiTrakMLAnalyzer
// Class:      DiMuonDiTrakMLAnalyzer
//
// Author:  Adriano Di Florio
//

#include "../interface/DiMuonDiTrakMLAnalyzer.h"
#include "../interface/DiMuonVtxReProducer.h"


//
// constructors and destructor
//

bool DiMuonDiTrakMLAnalyzer::IsTheSame(const reco::Muon& mu1, const reco::Muon& mu2){
  double DeltaEta = fabs(mu1.eta()-mu2.eta());
  double DeltaP   = fabs(mu1.p()-mu2.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;

  return false;
}

bool DiMuonDiTrakMLAnalyzer::IsTheSame(const reco::Muon& mu, const reco::Track& tk){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;

  return false;
}

bool DiMuonDiTrakMLAnalyzer::IsTheSame( const reco::Track& tk1,  const reco::Track& tk2)
{
double DeltaEta = fabs(tk1.eta()-tk2.eta());
double DeltaP   = fabs(tk1.p()-tk2.p());
if (DeltaEta < 0.02 && DeltaP < 0.02) return true;

return false;
}

float DiMuonDiTrakMLAnalyzer::DeltaR(const reco::Track t1, const pat::TriggerObject t2)
{
   float p1 = t1.phi();
   float p2 = t2.phi();
   float e1 = t1.eta();
   float e2 = t2.eta();
   auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

bool DiMuonDiTrakMLAnalyzer::MatchByDRDPt(const reco::Track t1, const pat::TriggerObject t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<maxDPtRel &&
	DeltaR(t1,t2) < maxDeltaR);
}

DiMuonDiTrakMLAnalyzer::DiMuonDiTrakMLAnalyzer(const edm::ParameterSet& iConfig):
muons_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("Muons"))),
traks_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("Tracks"))),
triggerEvent_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerEvent"))),
thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
triggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
DiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonCuts")),
DiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiTrakCuts")),
DiMuonMass_(iConfig.getParameter<double>("DiMuonMass")),
HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs")),
HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters"))
// thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
// thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
// DiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonCuts")),
// DiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiTrakCuts")),
// DiMuonDiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonDiTrakCuts")),
// massCands_(iConfig.getParameter<std::vector<double>>("CandsMasses"))
{
  // produces<pat::CompositeCandidateCollection>();
  muon_mass = 0.1056583715;
  cands = 0;
  dimuoncands = 0;
  trigger = 0;
  maxDeltaR = 0.01;
  maxDPtRel = 2.0;

  edm::Service < TFileService > fs;
  ml_tree = fs->make < TTree > ("DiMuonDiTrackML", "Tree of DiTrakDiMuon");

  ml_tree->Branch("run",      &run,      "run/i");
  ml_tree->Branch("event",    &event,    "event/l");
  ml_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  ml_tree->Branch("trigger",  &trigger,  "trigger/i");

  ml_tree->Branch("dimuonditrak_p4", "TLorentzVector", &dimuonditrak_p4);
  ml_tree->Branch("ditrak_p4", "TLorentzVector", &ditrak_p4);
  ml_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
  //
  // ml_tree->Branch("nditrak",  &nditrak,    "nditrak/i");
  // ml_tree->Branch("trigger",  &trigger,  "trigger/i");
  // ml_tree->Branch("charge",   &charge,   "charge/I");
  //
  // ml_tree->Branch("isBest",   &isBest,   "isBest/O");
  //
  // if(addTrigger_)
  // {
  //   ml_tree->Branch("trigs_pt",   &trigs_pt);
  //   ml_tree->Branch("trigs_eta",   &trigs_eta);
  //   ml_tree->Branch("trigs_phi",   &trigs_phi);
  //   ml_tree->Branch("trigs_m",   &trigs_m);
  //   ml_tree->Branch("trigs_filters", &trigs_filters);
  // }
  // ml_tree->Branch("ditrak_p4", "TLorentzVector", &ditrak_p4);
  // ml_tree->Branch("trakP_p4",  "TLorentzVector", &trakP_p4);
  // ml_tree->Branch("trakN_p4",  "TLorentzVector", &trakN_p4);
  //
  // ml_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
  // ml_tree->Branch("vProb",     &vProb,      "vProb/F");
  // ml_tree->Branch("DCA",       &DCA,        "DCA/F");
  // ml_tree->Branch("ctauPV",    &ctauPV,     "ctauPV/F");
  // ml_tree->Branch("ctauErrPV", &ctauErrPV,  "ctauErrPV/F");
  // ml_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");
  // ml_tree->Branch("lxy",       &lxyPV,      "lxy/F");
  // ml_tree->Branch("lxyErrPV",    &lxyErrPV,      "lxyErr/F");
  //
  // ml_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");


}

DiMuonDiTrakMLAnalyzer::~DiMuonDiTrakMLAnalyzer() {}

//
// member functions
//

const reco::Candidate* DiMuonDiTrakMLAnalyzer::GetAncestor(const reco::Candidate* p) {
  if (p->numberOfMothers()) {
    if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
    else return p->mother(0);
  }
  return p;
}


/* Grab Trigger information. Save it in variable trigger, trigger is an uint between 0 and 256, in binary it is:
(pass 2)(pass 1)(pass 0)
ex. 7 = pass 0, 1 and 2
ex. 6 = pass 1, 2
ex. 1 = pass 0
*/


UInt_t DiMuonDiTrakMLAnalyzer::getTriggerBits(const edm::Event& iEvent ) {

  UInt_t trigger = 0;

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_ , triggerResults_handle);

  if (triggerResults_handle.isValid()) {
    const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
    unsigned int NTRIGGERS = HLTs_.size();

    for (unsigned int i = 0; i < NTRIGGERS; i++) {
      for (int version = 1; version < 20; version++) {
        std::stringstream ss;
        ss << HLTs_[i] << "_v" << version;
        unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
        if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
          trigger += (1<<i);
          break;
        }
      }
    }
  } else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;

  return trigger;
}

// ------------ method called for each event  ------------
void DiMuonDiTrakMLAnalyzer::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  using namespace reco;

  std::vector<int> pixelDets{0,1,2,3,14,15,16,29,30,31};

  std::vector<double> mmMasses,kMasses;
  mmMasses.push_back( 0.1056583715 );
  mmMasses.push_back( 0.1056583715 );
  kMasses.push_back( 0.493677 );
  kMasses.push_back( 0.493677 );

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(traks_,tracks);

  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muons_,muons);

  edm::Handle<trigger::TriggerEvent> triggerEvent;
  iEvent.getByToken(triggerEvent_,triggerEvent);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_ , triggerResults_handle);

  trigger = 0;

  if (triggerResults_handle.isValid())
    trigger = getTriggerBits(iEvent);//,triggerResults_handle);
  else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;

  std::vector < UInt_t > filterResults;
  trigger::TriggerObjectCollection filteredColl;
  reco::MuonCollection filteredMuons;
  reco::TrackCollection filteredTracks;
  std::vector<unsigned int> muonTrigs,trackTrigs;
  const trigger::size_type nFilters(triggerEvent->sizeFilters());

  for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter)
  {
    //get the filter name
    std::string filterTag = triggerEvent->filterTag(iFilter).encode();
    //search for this filter in the one we want
    if(std::find(HLTFilters_.begin(),HLTFilters_.end(),filterTag)==HLTFilters_.end())
      continue;
    trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
    const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());

    for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey)
    {
      trigger::size_type objKey = objectKeys.at(iKey);
      filteredColl.push_back(triggerObjects[objKey]);
    }
  }


  for(reco::MuonCollection::const_iterator muon = muons->begin();muon != muons->end(); ++muon )
  {
    bool matched = false;
    // for (TriggerObjectCollection::const_iterator trigger = filteredColl.begin(), triggerEnd=filteredColl.end(); trigger!= triggerEnd; ++trigger)
    // {
    for (size_t i = 0; i < filteredColl.size(); i++)
    {
      if(MatchByDRDPt(*muon,filteredColl[i]))
      {
        if(matched)
        {
          if(DeltaR(*muon,filteredColl[matchedColl.back()]) > DeltaR(*muon,filteredColl[i]))
          {
            muonTrigs.pop_back();
            muonTrigs.push_back(i);

          }
        }

        if(!matched)
          {
            filteredTracks.push_back(*muon);
            muonTrigs.push_back(i);
          }

        matched = true;
      }
    }
  }

  for(reco::TrackCollection::const_iterator trak = tracks->begin();trak != tracks->end(); ++trak )
  {
    bool matched = false;
    // for (TriggerObjectCollection::const_iterator trigger = filteredColl.begin(), triggerEnd=filteredColl.end(); trigger!= triggerEnd; ++trigger)
    // {
    for (size_t i = 0; i < filteredColl.size(); i++)
    {
      if(MatchByDRDPt(*trak,filteredColl[i]))
      {
        if(matched)
        {
          if(DeltaR(*trak,filteredColl[matchedColl.back()]) > DeltaR(*trak,filteredColl[i]))
          {
            trackTrigs.pop_back();
            trackTrigs.push_back(i);

          }
        }

        if(!matched)
          {
            filteredMuons.push_back(*trak);
            trackTrigs.push_back(i);
          }

        matched = true;
      }
    }
  }


  reco::Vertex thePrimaryV;

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  const MagneticField* field = magneticField.product();

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;

  edm::Handle<reco::VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);
  if ( priVtxs->begin() != priVtxs->end() ) {
    thePrimaryV = reco::Vertex(*(priVtxs->begin()));
  }
  else {
    thePrimaryV = reco::Vertex(bs.position(), bs.covariance3D());
  }


  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  int max = 0;

  float vProb, vNDF, vChi2, minDz = 999999.;
  float vProb_mumu, vNDF_mumu, vChi2_mumu;
  float cosAlpha, ctauPV, ctauErrPV, dca;
  float l_xy, lErr_xy;

  for(reco::MuonCollection::const_iterator mPos = filteredMuons.begin();mPos != filteredMuons.end(); ++mPos )
  {
    if(mPos->charge()<=0.0) continue;
    // if (!(mPos->bestTrackRef().isNonnull())) continue;
    if (!(mPos->innerTrack().isNonnull())) continue;

    for(reco::MuonCollection::const_iterator mNeg = filteredMuons.begin();mNeg != filteredMuons.end(); ++mNeg )
    {
      if(mNeg->charge()>=0.0) continue;
      // if (!(mNeg->bestTrackRef().isNonnull())) continue;
      if (!(mNeg->innerTrack().isNonnull())) continue;

      if(mNeg==mPos) continue;

      if(IsTheSame(*mPos,*mNeg)) continue;

      std::vector<TransientVertex> vDiMuon;

      // Candidate::LorentzVector mumu = mNeg->p4() + mPos->p4();

      TLorentzVector mu1, mu2,mumuP4;

      mu1.SetXYZM(mNeg->track()->px(),mNeg->track()->py(),mNeg->track()->pz(),muon_mass);
      mu2.SetXYZM(mPos->track()->px(),mPos->track()->py(),mPos->track()->pz(),muon_mass);

      mumuP4=mu1+mu2;
      // mumucand.setP4(mumu);
      // mumucand.setCharge(mNeg->charge()+mPos->charge());

      if(mumuP4.M() > DiMuonMassCuts_[1]) continue;
      if(mumuP4.M() < DiMuonMassCuts_[0]) continue;

      std::vector<reco::TransientTrack> mm_ttks;

      mm_ttks.push_back(theTTBuilder->build(mNeg->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
      mm_ttks.push_back(theTTBuilder->build(mPos->track()));

      TransientVertex mumuVertex = vtxFitter.vertex(mm_ttks);
      CachingVertex<5> VtxForInvMass = vtxFitter.vertex( mm_ttks );

      Measurement1D MassWErr(mPos->mass(),-9999.);
      if ( field->nominalValue() > 0 )
      MassWErr = massCalculator.invariantMass( VtxForInvMass, mmMasses );
      else
      mumuVertex = TransientVertex();                      // with no arguments it is invalid

      if (!(mumuVertex.isValid()))
      continue;

      vChi2_mumu = mumuVertex.totalChiSquared();
      vNDF_mumu  = mumuVertex.degreesOfFreedom();
      vProb_mumu = TMath::Prob(vChi2_mumu,(int)vNDF_mumu);

      if(vProb_mumu < 0.0) continue;
      dimuoncands++;
      // int pv_index = -1;

      for(reco::TrackCollection::const_iterator posTrack = filteredTracks->begin();posTrack != filteredTracks->end(); ++posTrack )
      {
        if(!(posTrack->extra())) continue;
        if(posTrack->charge()<=0.0) continue;
        if(posTrack->pt() < 0.5) continue;

        if(!(posTrack->extra().isNonnull())) continue;

        if(IsTheSame(*mPos,*posTrack)) continue;
        if(IsTheSame(*mNeg,*posTrack)) continue;

        // if(!(posTrack->isNonnull())) continue;


        for(reco::TrackCollection::const_iterator negTrack = filteredTracks->begin();negTrack != filteredTracks->end(); ++negTrack )
        {
          if(!(negTrack->extra())) continue;
          if(negTrack->charge()>=0.0) continue;
          if(negTrack->pt() < 0.5) continue;

          if(!(negTrack->extra().isNonnull())) continue;
          if(negTrack==posTrack) continue;
          // if(!(negTrack->isNonnull())) continue;

          if(IsTheSame(*negTrack,*posTrack)) continue;

          if(IsTheSame(*mPos,*negTrack)) continue;
          if(IsTheSame(*mNeg,*negTrack)) continue;

          TLorentzVector kNeg, kPos, kkP4;

          kNeg.SetXYZM(negTrack->px(),negTrack->py(),negTrack->pz(),muon_mass);
          kPos.SetXYZM(posTrack->px(),posTrack->py(),posTrack->pz(),muon_mass);

          kkP4=kNeg+kPos;

          float deltaphi  = kPos.Phi() - kNeg.Phi();
          while (deltaphi < -M_PI) deltaphi += 2*M_PI;
          while (deltaphi >  M_PI) deltaphi -= 2*M_PI;
          float deltaeta  = kPos.Eta() - kNeg.Eta();
          float deltar    = sqrt(pow(deltaphi,2) + pow(deltaeta,2));

          if(kkP4.M() > DiTrakMassCuts_[1]) continue;
          if(kkP4.M() < DiTrakMassCuts_[0]) continue;


          //
          //build the dikaon secondary vertex
          std::vector<reco::TransientTrack> t_tks;
          t_tks.push_back(theTTBuilder->build(*posTrack));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
          t_tks.push_back(theTTBuilder->build(*negTrack)); // otherwise the vertex will have transient refs inside.
          TransientVertex myVertex = vtxFitter.vertex(t_tks);

          CachingVertex<5> VtxForInvMass = vtxFitter.vertex( t_tks );
          //
          Measurement1D MassWErr(kkP4.M(),-9999.);
          if ( field->nominalValue() > 0 ) MassWErr = massCalculator.invariantMass( VtxForInvMass, kMasses );
          else myVertex = TransientVertex();                   // this is an invalid vertex by definition
          //
          //
          if (myVertex.isValid()) {
            float vChi2 = myVertex.totalChiSquared();
            float vNDF  = myVertex.degreesOfFreedom();
            float vProb(TMath::Prob(vChi2,(int)vNDF));

            if(vProb < 0.01) continue;

            TVector3 vtx;
            TVector3 pvtx;
            VertexDistanceXY vdistXY;

            vtx.SetXYZ(myVertex.position().x(),myVertex.position().y(),0);
            TVector3 pperp(kkP4.Px(), kkP4.Py(), 0);
            AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);

            float minDz = 99999.;
            float extrapZ=-9E20;
            bool resolveAmbiguity_ = true;

            if (resolveAmbiguity_) {
              TwoTrackMinimumDistance ttmd;
              bool status = ttmd.calculate( GlobalTrajectoryParameters(
                GlobalPoint(myVertex.position().x(), myVertex.position().y(), myVertex.position().z()),
                GlobalVector(kkP4.Px(),kkP4.Py(),kkP4.Pz()),TrackCharge(0),&(*magneticField)),
                GlobalTrajectoryParameters(
                  GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
                  GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(0),&(*magneticField)));
                  if (status) extrapZ=ttmd.points().first.z();

                  // int ii_pv = -1;
                  for (VertexCollection::const_iterator itv = priVtxs->begin(), itvend = priVtxs->end(); itv != itvend; ++itv) {
                    // ii_pv++;
                    float deltaZ = fabs(extrapZ - itv->position().z()) ;
                    if ( deltaZ < minDz ) {
                      minDz = deltaZ;
                      thePrimaryV = Vertex(*itv);
                      // pv_index = ii_pv;
                    }
                  }
                } else {
                  minDz = -1;
                  // pv_index = 0;
                  thePrimaryV = (*priVtxs)[0];
                  extrapZ = thePrimaryV.position().z();
                }
                //
                //         myPhi.addUserInt("oniaPV",which_vertex);
                //   myPhi.addUserInt("iPV",pv_index);
                //   myPhi.addUserFloat("dzPV",minDz);
                //   myPhi.addUserFloat("extrapZPV",extrapZ);
                //
                //   // count the number of high Purity tracks with pT > 500 MeV attached to the chosen vertex
                //   double vertexWeight = 0., sumPTPV = 0.;
                //   int countTksOfPV = 0;
                //         for (size_t kk=1; kk<(size_t)ntracks_pv; kk++) {
                //             const pat::GenericParticle *it3 = &(kaons->at(kk));
                //             if (!it3->track().isNonnull())                  continue;
                //             reco::Track track = *it3->track();
                //             if(track.pt() < 0.5)                            continue;
                //             if(!track.quality(reco::TrackBase::highPurity)) continue;
                //             TransientTrack tt = theTTBuilder->build(track);
                //             pair<bool,Measurement1D> tkPVdist = IPTools::absoluteImpactParameter3D(tt,thePrimaryV);
                //             if (!tkPVdist.first)                  continue;
                //             if (tkPVdist.second.significance()>3) continue;
                //             if (track.ptError()/track.pt()>0.1)   continue;
                //             if (it3 == it2 || it3 == it)          continue;
                //             countTksOfPV++;
                //             sumPTPV += track.pt();
                //             vertexWeight += thePrimaryV.trackWeight(it3->track());
                //         }
                //
                //   myPhi.addUserInt("countTksOfPV", countTksOfPV);
                //   myPhi.addUserFloat("vertexWeight", (float) vertexWeight);
                //   myPhi.addUserFloat("sumPTPV", (float) sumPTPV);
                //
                //   ///DCA
                TrajectoryStateClosestToPoint kNegTS = t_tks[0].impactPointTSCP();
                TrajectoryStateClosestToPoint kPosTS = t_tks[1].impactPointTSCP();
                float dca = 1E20;
                if (kNegTS.isValid() && kPosTS.isValid()) {
                  ClosestApproachInRPhi cApp;
                  cApp.calculate(kNegTS.theState(), kPosTS.theState());
                  if (cApp.status() ) dca = cApp.distance();
                }

                if(dca > 1000000000000000000) continue;

                //   myPhi.addUserFloat("DCA", dca );
                //   ///end DCA
                //
                //   myPhi.addUserData("PVwithkaons",thePrimaryV);
                //   npvtracks = thePrimaryV.nTracks();
                //
                //   // lifetime using PV
                pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
                TVector3 vdiff = vtx - pvtx;
                double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
                Measurement1D distXY = vdistXY.distance(Vertex(myVertex), thePrimaryV);
                double ctauPV = distXY.value()*cosAlpha * kkP4.M()/pperp.Perp();
                GlobalError v1e = (Vertex(myVertex)).error();
                GlobalError v2e = thePrimaryV.error();
                AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
                double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*kkP4.M()/(pperp.Perp2());

                std::vector<reco::TransientTrack> MuMuKK;
                // reco::TrackRef The4Tks[4]={,mNeg->innerTrack(),(*posTrack).extra(),(*negTrack).extra()};
                reco::TrackRef  the2MuTks[2] = {mPos->innerTrack(),mNeg->innerTrack()};

                MuMuKK.push_back(theTTBuilder->build(&the2MuTks[0])); // Mu+
                MuMuKK.push_back(theTTBuilder->build(&the2MuTks[1])); // Mu-
                MuMuKK.push_back(theTTBuilder->build(*posTrack)); // K+
                MuMuKK.push_back(theTTBuilder->build(*negTrack)); // K+

                KinematicParticleFactoryFromTransientTrack pFactory;

                const ParticleMass mMass(0.1056583745);
                float mSigma = mMass*1E-6;
                const ParticleMass kMass(0.493677);
                float kSigma = kMass*1E-6;

                std::vector<RefCountedKinematicParticle> allDaughters;
                allDaughters.push_back(pFactory.particle (MuMuKK[0], mMass, float(0), float(0), mSigma));
                allDaughters.push_back(pFactory.particle (MuMuKK[1], mMass, float(0), float(0), mSigma));
                allDaughters.push_back(pFactory.particle (MuMuKK[2], kMass, float(0), float(0), kSigma));
                allDaughters.push_back(pFactory.particle (MuMuKK[3], kMass, float(0), float(0), kSigma));

                KinematicConstrainedVertexFitter constVertexFitter;
                MultiTrackKinematicConstraint *mumu_c = new  TwoTrackMassKinematicConstraint(DiMuonMass_);
                RefCountedKinematicTree allTree = constVertexFitter.fit(allDaughters,mumu_c);

                if (allTree->isEmpty()) continue;

                allTree->movePointerToTheTop();
                RefCountedKinematicParticle TheParticle = allTree->currentParticle();
                RefCountedKinematicVertex TheVertex = allTree->currentDecayVertex();
                float mmkk_ma_fit = 14000.;
                float mmkk_vp_fit = -9999.;
                float mmkk_x2_fit = 10000.;
                if (TheParticle->currentState().isValid()) {
                    mmkk_ma_fit = TheParticle->currentState().mass();
                    mmkk_x2_fit = TheVertex->chiSquared();
                    mmkk_vp_fit = ChiSquaredProbability(mmkk_x2_fit,TheVertex->degreesOfFreedom());
                }
                if ( mmkk_ma_fit < 5.15 || mmkk_ma_fit > 5.55 || mmkk_vp_fit < 0.005 ) continue;
                cands++;
                // VertexDistanceXY vdistXY;
                float mmkk_px_fit = TheParticle->currentState().kinematicParameters().momentum().x();
                float mmkk_py_fit = TheParticle->currentState().kinematicParameters().momentum().y();
                float mmkk_pz_fit = TheParticle->currentState().kinematicParameters().momentum().z();
                float mmkk_en_fit = sqrt(mmkk_ma_fit*mmkk_ma_fit+mmkk_px_fit*mmkk_px_fit+mmkk_py_fit*mmkk_py_fit+mmkk_pz_fit*mmkk_pz_fit);
                float mmkk_vx_fit = TheVertex->position().x();
                float mmkk_vy_fit = TheVertex->position().y();
                float mmkk_vz_fit = TheVertex->position().z();

                reco::CompositeCandidate recoMMKK(0.,math::XYZTLorentzVector(mmkk_px_fit,mmkk_py_fit,mmkk_pz_fit,mmkk_en_fit),
                                                     math::XYZPoint(mmkk_vx_fit,mmkk_vy_fit,mmkk_vz_fit));
                pat::CompositeCandidate patMMKK(recoMMKK);

                dimuonditrak_p4.SetPtEtaPhiM(recoMMKK.pt(),recoMMKK.eta(),recoMMKK.phi(),recoMMKK.mass());
                dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
                ditrak_p4.SetPtEtaPhiM(0.,0.,0.,0.);

                //
                //   myPhi.addUserFloat("ppdlPV",ctauPV);
                //         myPhi.addUserFloat("ppdlErrPV",ctauErrPV);
                //   myPhi.addUserFloat("cosAlpha",cosAlpha);
                //
                //   // lifetime using BS
                //   pvtx.SetXYZ(theBeamSpotV.position().x(),theBeamSpotV.position().y(),0);
                //   vdiff = vtx - pvtx;
                //   cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
                //   distXY = vdistXY.distance(Vertex(myVertex), theBeamSpotV);
                //   double ctauBS = distXY.value()*cosAlpha*kkP4.M()/pperp.Perp();
                //   GlobalError v1eB = (Vertex(myVertex)).error();
                //   GlobalError v2eB = theBeamSpotV.error();
                //         AlgebraicSymMatrix33 vXYeB = v1eB.matrix()+ v2eB.matrix();
                //   double ctauErrBS = sqrt(ROOT::Math::Similarity(vpperp,vXYeB))*kkP4.M()/(pperp.Perp2());
                //
                //   myPhi.addUserFloat("ppdlBS",ctauBS);
                //         myPhi.addUserFloat("ppdlErrBS",ctauErrBS);
                //
                //   if (addCommonVertex_) myPhi.addUserData("commonVertex",Vertex(myVertex));
                //         myPhi.addUserInt("npvtracks", npvtracks);
                //   myPhi.addUserInt("ntracks_pv",ntracks_pv );
                //
                //         // ---- If here push back to output ----
                //         phiOutput->push_back(myPhi);
                //
              }


            }


          }


        }

      }

      // 	int padHalfSize = 8;
      // 	int padSize = padHalfSize*2;
      //   int maxpix = 0;
      //
      //   for(reco::TrackCollection::const_iterator itTrack = tracks->begin();itTrack != tracks->end(); ++itTrack )
      //     {
      // 	int noPixels= 0,noStripOne = 0, noStripTwo = 0;
      // 	int counter = 0;
      // 	float clusterSize = 0.0;
      // //	std::cout<<"On "<< itTrack->found() ;
      // 	for ( trackingRecHit_iterator recHit = (*itTrack).recHitsBegin();recHit != (*itTrack).recHitsEnd(); ++recHit )
      // 	{
      // 		counter++;
      // /*		 if(!(*recHit))
      //                  continue;
      //
      //                  if (!((*recHit)->isValid()))
      //                  continue;
      //
      //                  if(!((*recHit)->hasPositionAndError()))
      //                  continue;
      // */
      // 		TrackerSingleRecHit const * hit= dynamic_cast<TrackerSingleRecHit const *>(*recHit);
      // //    BaseTrackerRecHit const * bhit = dynamic_cast<BaseTrackerRecHit const *>(recHit);
      //
      // 		DetId detid = (*recHit)->geographicalId();
      // 		unsigned int subdetid = detid.subdetId();
      //
      //
      // 	        //if(!(siPix))
      // 	        //continue;
      // 	        //
      //
      //     	 if(detid.det() != DetId::Tracker) continue;
      // //	 if (!((subdetid==1) || (subdetid==2))) continue;
      // // 	 if()
      // 		if (dynamic_cast<SiPixelRecHit const *>(hit))
      // 		{		noPixels++;
      // 				clusterSize += float(dynamic_cast<SiPixelRecHit const *>(hit)->cluster()->size());
      // 				clusterSize /= float(counter);
      //
      // 		auto clust = dynamic_cast<SiPixelRecHit const *>(hit)->cluster();
      // 		TH2F hClust("hClust","hClust",
      //               padSize,
      //               clust->x()-padHalfSize,
      //               clust->x()+padHalfSize,
      //               padSize,
      //               clust->y()-padHalfSize,
      //               clust->y()+padHalfSize);
      //
      // 		for (int nx = 0; nx < padSize; ++nx)
      //               for (int ny = 0; ny < padSize; ++ny)
      //               hClust.SetBinContent(nx,ny,0.0);
      //
      //               for (int k = 0; k < clust->size(); ++k)
      //               hClust.SetBinContent(hClust.FindBin((float)clust->pixel(k).x, (float)clust->pixel(k).y),(float)clust->pixel(k).adc);
      //
      // 		for (int ny = padSize; ny>0; --ny)
      //               {
      //                 for(int nx = 0; nx<padSize; nx++)
      //                 {
      //                   int n = (ny+2)*(padSize + 2) - 2 -2 - nx - padSize; //see TH2 reference for clarification
      //
      //        //          std::cout << hClust.GetBinContent(n) << " ";
      // 		}
      //               }
      // 	//	std::cout << std::endl;
      //
      // 		}
      // 		if (dynamic_cast<SiStripRecHit1D const *>(hit))
      // 		noStripOne++;
      //
      // 	if (dynamic_cast<SiStripRecHit2D const *>(hit))
      // 			noStripTwo++;
      //
      // 	}
      // //	std::cout << " n. pixels = " << noPixels<< " 1DStrips = " << noStripOne << " 2DStrips = " << noStripTwo<< " clustsize : "<< clusterSize <<std::endl;
      //
      // 	     max = std::max(max,int(itTrack->found()));
      // 		maxpix = std::max(maxpix,noPixels);
      // }
      //   std::cout<<"Max = " << max<< " Max pixels " << maxpix << std::endl;


    }

    // ------------ method called once each job just before starting event loop  ------------
    void DiMuonDiTrakMLAnalyzer::beginJob() {}

    // ------------ method called once each job just after ending the event loop  ------------
    void DiMuonDiTrakMLAnalyzer::endJob() {

      std::cout << "No. candidates : "<<cands<<std::endl;
      std::cout << "No. dimuoncand : "<<dimuoncands<<std::endl;
    }

    // ------------ method called when starting to processes a run  ------------
    void DiMuonDiTrakMLAnalyzer::beginRun(edm::Run const &, edm::EventSetup const &) {}

    // ------------ method called when ending the processing of a run  ------------
    void DiMuonDiTrakMLAnalyzer::endRun(edm::Run const &, edm::EventSetup const &) {}

    // ------------ method called when starting to processes a luminosity block  ------------
    void DiMuonDiTrakMLAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

    // ------------ method called when ending the processing of a luminosity block  ------------
    void DiMuonDiTrakMLAnalyzer::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

    // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
    void DiMuonDiTrakMLAnalyzer::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
      //The following says we do not know what parameters are allowed so do no validation
      // Please change this to state exactly what you do use, even if it is no parameters
      edm::ParameterSetDescription desc;
      desc.setUnknown();
      descriptions.addDefault(desc);
    }

    //define this as a plug-in
    DEFINE_FWK_MODULE(DiMuonDiTrakMLAnalyzer);
