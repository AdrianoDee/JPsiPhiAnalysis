#include "../interface/DiMuonDiTrakML.h"
#include "../interface/DiMuonVtxReProducer.h"

DiMuonDiTrakML::DiMuonDiTrakML(const edm::ParameterSet& iConfig):
muons_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("Muons"))),
traks_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("Tracks"))),
thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
DiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonCuts")),
DiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiTrakCuts"))
// thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
// thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
// DiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonCuts")),
// DiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiTrakCuts")),
// DiMuonDiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonDiTrakCuts")),
// massCands_(iConfig.getParameter<std::vector<double>>("CandsMasses"))
{
  // produces<pat::CompositeCandidateCollection>();
  muon_mass = 0.1056583715;
}


DiMuonDiTrakML::~DiMuonDiTrakML()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


void
DiMuonDiTrakML::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::vector<double> mmMasses;
  mmMasses.push_back( 0.1056583715 );
  mmMasses.push_back( 0.1056583715 );

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(traks_,tracks);

  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muons_,muons);

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

  for(reco::MuonCollection::const_iterator mPos = muons->begin();mPos != muons->end(); ++mPos )
  {
    if(mPos->charge()<=0.0) continue;
    // if (!(mPos->bestTrack()->isNonnull())) continue;
    // if (!(mPos->innerTrack()->isNonnull())) continue;

    for(reco::MuonCollection::const_iterator mNeg = muons->begin();mNeg != muons->end(); ++mNeg )
    {
      if(mNeg->charge()>=0.0) continue;
      // if (!(mNeg->bestTrack()->isNonnull())) continue;
      // if (!(mNeg->innerTrack()->isNonnull())) continue;

      std::vector<TransientVertex> vDiMuon;

      // Candidate::LorentzVector mumu = mNeg->p4() + mPos->p4();

      TLorentzVector mu1, mu2,mumuP4;

      mu1.SetXYZM(mNeg->track()->px(),mNeg->track()->py(),mNeg->track()->pz(),muon_mass);
      mu2.SetXYZM(mPos->track()->px(),mPos->track()->py(),mPos->track()->pz(),muon_mass);

      mumuP4=mu1+mu2;
      // mumucand.setP4(mumu);
      // mumucand.setCharge(mNeg->charge()+mPos->charge());

      if(mumuP4.M() < DiMuonMassCuts_[1]) continue;
      if(mumuP4.M() > DiMuonMassCuts_[0]) continue;

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

      for(reco::TrackCollection::const_iterator posTrack = tracks->begin();posTrack != tracks->end(); ++posTrack )
      {
        if(posTrack->charge()<=0.0) continue;
        // if(!(posTrack->isNonnull())) continue;

        for(reco::TrackCollection::const_iterator negTrack = tracks->begin();negTrack != tracks->end(); ++negTrack )
        {
          if(negTrack->charge()>=0.0) continue;
          // if(!(negTrack->isNonnull())) continue;


        }


      }


    }

  }
    for(reco::TrackCollection::const_iterator itTrack = tracks->begin();itTrack != tracks->end(); ++itTrack )
          max = std::max(max,int(itTrack->numberOfValidHits()));

  std::cout<<max<<std::endl;
  //loop on
  // std::sort(mmttCollection->begin(),mmttCollection->end(),vPComparator_);
  // iEvent.put(std::move(mmttCollection));

}


// ------------ method called once each job just before starting event loop  ------------
void
DiMuonDiTrakML::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DiMuonDiTrakML::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakML);
