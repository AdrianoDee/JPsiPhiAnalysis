#include "../interface/DiMuon.h"
#include "../interface/DiMuonVtxReProducer.h"

DiMuonPAT::DiMuonPAT(const edm::ParameterSet& iConfig):
muons_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("Muons"))),
thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
dimuonSelection_(iConfig.existsAs<std::string>("DiMuonCuts") ? iConfig.getParameter<std::string>("DiMuonCuts") : "")
{
  revtxtrks_ = consumes<reco::TrackCollection>((edm::InputTag)"generalTracks"); //if that is not true, we will raise an exception
  revtxbs_ = consumes<reco::BeamSpot>((edm::InputTag)"offlineBeamSpot");

  produces<pat::CompositeCandidateCollection>();
  muon_mass = 0.1056583715;
}


DiMuonPAT::~DiMuonPAT()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


void
DiMuonPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  vector<double> mmMasses;
  mmMasses.push_back( 0.1056583715 );
  mmMasses.push_back( 0.1056583715 );

  typedef Candidate::LorentzVector LorentzVector;

  std::unique_ptr<pat::CompositeCandidateCollection> mumuCollection(new pat::CompositeCandidateCollection);

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muons_,muons);

  Vertex thePrimaryV;

  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  const MagneticField* field = magneticField.product();

  edm::Handle<BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  BeamSpot bs = *theBeamSpot;

  edm::Handle<VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);
  if ( priVtxs->begin() != priVtxs->end() ) {
    thePrimaryV = Vertex(*(priVtxs->begin()));
  }
  else {
    thePrimaryV = Vertex(bs.position(), bs.covariance3D());
  }

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  TrackCollection muonLess;

  float DiMuonMassMax_ = dimuonMassCuts_[1];
  float DiMuonMassMin_ = dimuonMassCuts_[0];

  float vProb, vNDF, vChi2, minDz = 999999.;
  float cosAlpha, ctauPV, ctauErrPV, dca;
  float l_xy, lErr_xy;

  for(View<pat::Muon>::const_iterator mNeg = muons->begin(), itend = muons->end(); mNeg != itend; ++mNeg){

    if(mNeg->charge()>=0.0) continue;

    for(View<pat::Muon>::const_iterator mPos = muons->begin(), itend = muons->end(); mPos != itend; ++mPos)
    {
      if(mNeg == mPos) continue;

      if(mPos->charge()<=0.0) continue;

      if (!(mNeg->track().isNonnull() && mPos->track().isNonnull())) continue;

      vProb = -1.0; vNDF = -1.0; vChi2 = -1.0;
      cosAlpha = -1.0; ctauPV = -1.0; ctauErrPV = -1.0;
      minDz = 999999.; dca = 1E20;

      pat::CompositeCandidate mumucand;

      vector<TransientVertex> pvs;

      mumucand.addDaughter(*mNeg,"muonN");
      mumucand.addDaughter(*mPos,"muonP");

      LorentzVector mumu = mNeg->p4() + mPos->p4();
      TLorentzVector mu1, mu2,mumuP4;

      mu1.SetXYZM(mNeg->track()->px(),mNeg->track()->py(),mNeg->track()->pz(),muon_mass);
      mu2.SetXYZM(mPos->track()->px(),mPos->track()->py(),mPos->track()->pz(),muon_mass);

      mumuP4=mu1+mu2;
      mumucand.setP4(mumu);
      mumucand.setCharge(mNeg->charge()+mPos->charge());

      if(!dimuonSelection_(mumucand)) continue;

      std::vector<TransientTrack> mm_ttks;
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

      vChi2 = mumuVertex.totalChiSquared();
      vNDF  = mumuVertex.degreesOfFreedom();
      vProb = TMath::Prob(vChi2,(int)vNDF);

      //Vertex parameters
      TVector3 vtx,vtx3D;
      TVector3 pvtx,pvtx3D;
      VertexDistanceXY vdistXY;

      vtx.SetXYZ(mumuVertex.position().x(),mumuVertex.position().y(),0);
      vtx3D.SetXYZ(mumuVertex.position().x(),mumuVertex.position().y(),mumuVertex.position().z());
      TVector3 pperp(mumu.px(), mumu.py(), 0);
      TVector3 pperp3D(mumu.px(), mumu.py(), mumu.pz());
      AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
      AlgebraicVector3 vpperp3D(pperp.x(),pperp.y(),pperp.z());

      //Resolving pileup ambiguity with two trak min distance
      TwoTrackMinimumDistance ttmd;
      bool status = ttmd.calculate( GlobalTrajectoryParameters(
        GlobalPoint(mumuVertex.position().x(), mumuVertex.position().y(), mumuVertex.position().z()),
        GlobalVector(mumucand.px(),mumucand.py(),mumucand.pz()),TrackCharge(0),&(*magneticField)),
        GlobalTrajectoryParameters(
          GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
          GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(0),&(*magneticField)));

      float extrapZ=-9E20;

      if (status) extrapZ=ttmd.points().first.z();

      for(VertexCollection::const_iterator itv = priVtxs->begin(), itvend = priVtxs->end(); itv != itvend; ++itv)
      {
          float deltaZ = fabs(extrapZ - itv->position().z()) ;
          if ( deltaZ < minDz ) {
              minDz = deltaZ;
              thePrimaryV = Vertex(*itv);
            }
        }

      //Distance of Closest Approach
      TrajectoryStateClosestToPoint mu1TS = mm_ttks[0].impactPointTSCP();
      TrajectoryStateClosestToPoint mu2TS = mm_ttks[1].impactPointTSCP();

      if (mu1TS.isValid() && mu2TS.isValid()) {
        ClosestApproachInRPhi cApp;
        cApp.calculate(mu1TS.theState(), mu2TS.theState());
        if (cApp.status() ) dca = cApp.distance();
      }

      //Lifetime calculations
      pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
      TVector3 vdiff = vtx - pvtx;
      cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());

      Measurement1D distXY = vdistXY.distance(Vertex(mumuVertex), thePrimaryV);
      ctauPV = distXY.value()*cosAlpha * mumucand.mass()/pperp.Perp();

      GlobalError v1e = (Vertex(mumuVertex)).error();
      GlobalError v2e = thePrimaryV.error();
      AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
      ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*mumucand.mass()/(pperp.Perp2());

      AlgebraicVector3 vDiff;
      vDiff[0] = vdiff.x(); vDiff[1] = vdiff.y(); vDiff[2] = 0 ;
      l_xy = vdiff.Perp();
      lErr_xy = sqrt(ROOT::Math::Similarity(vDiff,vXYe)) / vdiff.Perp();

      //MuonLess PV

      Vertex muonLessPV = Vertex();
      muonLess.clear();
      muonLess.reserve(thePrimaryV.tracksSize());

      if( thePrimaryV.tracksSize()>2) {
        // Primary vertex matched to the dimuon, now refit it removing the two muons
        DiMuonVtxReProducer revertex(priVtxs, iEvent);
        edm::Handle<reco::TrackCollection> pvtracks;
        iEvent.getByToken(revtxtrks_,   pvtracks);
        if( !pvtracks.isValid()) { std::cout << "pvtracks NOT valid " << std::endl; }
        else {
          edm::Handle<reco::BeamSpot> pvbeamspot;
          iEvent.getByToken(revtxbs_, pvbeamspot);
          if (pvbeamspot.id() != theBeamSpot.id()) edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";

          const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(mNeg->originalObject());
          const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(mPos->originalObject());

          if (rmu1 != nullptr && rmu2 != nullptr && rmu1->track().id() == pvtracks.id() && rmu2->track().id() == pvtracks.id()) {

            if( thePrimaryV.hasRefittedTracks() ) {

              std::vector<reco::Track>::const_iterator itRefittedTrack = thePrimaryV.refittedTracks().begin();
              std::vector<reco::Track>::const_iterator refittedTracksEnd = thePrimaryV.refittedTracks().end();
              for( ; itRefittedTrack != refittedTracksEnd; ++itRefittedTrack )
              {
                if( thePrimaryV.originalTrack(*itRefittedTrack).key() == rmu1->track().key() ) continue;
                if( thePrimaryV.originalTrack(*itRefittedTrack).key() == rmu2->track().key() ) continue;

                muonLess.push_back(*(thePrimaryV.originalTrack(*itRefittedTrack)));
              }
            }
            else {
              std::vector<reco::TrackBaseRef>::const_iterator itPVtrack = thePrimaryV.tracks_begin();
              for( ; itPVtrack != thePrimaryV.tracks_end(); ++itPVtrack ) if (itPVtrack->isNonnull()) {
                if( itPVtrack->key() == rmu1->track().key() ) continue;
                if( itPVtrack->key() == rmu2->track().key() ) continue;

                muonLess.push_back(**itPVtrack);
              }
            }
            if (muonLess.size()>1 && muonLess.size() < thePrimaryV.tracksSize()){
              pvs = revertex.makeVertices(muonLess, *pvbeamspot, iSetup) ;
              if (!pvs.empty()) {
                muonLessPV = Vertex(pvs.front());
              }
            }
          }
        }
      }

      mumucand.addUserFloat("vNChi2",vChi2/vNDF);
      mumucand.addUserFloat("vProb",vProb);
      mumucand.addUserFloat("DCA", dca );
      mumucand.addUserFloat("MassErr",MassWErr.error());
      mumucand.addUserFloat("ctauPV",ctauPV);
      mumucand.addUserFloat("ctauErrPV",ctauErrPV);
      mumucand.addUserFloat("lxy",l_xy);
      mumucand.addUserFloat("lErrxy",lErr_xy);
      mumucand.addUserFloat("cosAlpha",cosAlpha);
      mumucand.addUserData("thePV",Vertex(thePrimaryV));
      mumucand.addUserData("theMuLessPV",Vertex(muonLessPV));
      mumucand.addUserData("theVertex",Vertex(mumuVertex));


      mumuCollection->push_back(mumucand);

    }
  }

  iEvent.put(std::move(mumuCollection));

}


// ------------ method called once each job just before starting event loop  ------------
void
DiMuonPAT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DiMuonPAT::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonPAT);
