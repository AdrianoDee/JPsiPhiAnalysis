#include "../interface/DiTrak.h"

DiTrakPAT::DiTrakPAT(const edm::ParameterSet& iConfig):
traks_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Traks"))),
thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
ditrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiTrakCuts")),
massTraks_(iConfig.getParameter<std::vector<double>>("TraksMasses"))
{
  produces<pat::CompositeCandidateCollection>();
}


DiTrakPAT::~DiTrakPAT()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

const pat::CompositeCandidate DiTrakPAT::makeTTCandidate(
                                          const pat::PackedCandidate& trakP,
                                          const pat::PackedCandidate& trakN
                                         ){

  pat::CompositeCandidate TTCand;
  TTCand.addDaughter(trakP,"trakP");
  TTCand.addDaughter(trakN,"trakN");
  TTCand.setCharge(trakP.charge()+trakN.charge());

  double m_trakP = massTraks_[0];
  math::XYZVector mom_trakP = trakP.momentum();
  double e_trakP = sqrt(m_trakP*m_trakP + mom_trakP.Mag2());
  math::XYZTLorentzVector p4_trakP = math::XYZTLorentzVector(mom_trakP.X(),mom_trakP.Y(),mom_trakP.Z(),e_trakP);
  double m_trakN = massTraks_[1];
  math::XYZVector mom_trakN = trakN.momentum();
  double e_trakN = sqrt(m_trakN*m_trakN + mom_trakN.Mag2());
  math::XYZTLorentzVector p4_trakN = math::XYZTLorentzVector(mom_trakN.X(),mom_trakN.Y(),mom_trakN.Z(),e_trakN);
  reco::Candidate::LorentzVector vTT = p4_trakP + p4_trakN;
  TTCand.setP4(vTT);

  return TTCand;
}


void
DiTrakPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  typedef Candidate::LorentzVector LorentzVector;

  std::unique_ptr<pat::CompositeCandidateCollection> trakCollection(new pat::CompositeCandidateCollection);

  edm::Handle<std::vector<pat::PackedCandidate> > traks;
  iEvent.getByToken(traks_,traks);

  Vertex thePrimaryV;
  Vertex theBeamSpotV;

  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  const MagneticField* field = magneticField.product();

  edm::Handle<BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  BeamSpot bs = *theBeamSpot;
  theBeamSpotV = Vertex(bs.position(), bs.covariance3D());

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

  float TrakTrakMassMax_ = ditrakMassCuts_[1];
  float TrakTrakMassMin_ = ditrakMassCuts_[0];

  ParticleMass trakP_mass = massTraks_[0];
  ParticleMass trakN_mass = massTraks_[0];

  float trakP_sigma = trakP_mass*1.e-6;
  float trakN_sigma = trakN_mass*1.e-6;

  float vProb, vNDF, vChi2, minDz = 999999.;
  float cosAlpha, ctauPV, ctauErrPV, dca;

  Vertex thePrimaryV;

  for (size_t i = 0; i < traks->size(); i++)
  {
    auto posTrack = traks->at(i);

    if(posTrack.charge() <= 0 ) continue;
    if(posTrack.pt()<0.5) continue;
    if(fabs(posTrack.pdgId())!=211) continue;


    for (size_t j = 0; j < traks->size(); j++){

      vProb = -1.0; vNDF = -1.0; vChi2 = -1.0;
      cosAlpha = -1.0; ctauPV = -1.0; ctauErrPV = -1.0;
      dca = -1.0; minDz = 999999.; dca = 1E20;

      if (i == j) continue;

      auto negTrack = traks->at(j);

      if(negTrack.charge() >= 0 ) continue;
      if(negTrack.pt()<0.5) continue;
      if(fabs(negTrack.pdgId())!=211) continue;

      pat::CompositeCandidate TTCand = makeTTCandidate(posTrack,negTrack);
      vector<TransientVertex> pvs;

      if ( !(TTCand.mass() < TrakTrakMassMax_ && TTCand.mass() > TrakTrakMassMin_) )
        continue;

      vector<TransientTrack> muon_ttks;
      muon_ttks.push_back(theTTBuilder->build(mNeg.track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
      muon_ttks.push_back(theTTBuilder->build(mPos.track()));

      TransientVertex mumuVertex = vtxFitter.vertex(muon_ttks);
      CachingVertex<5> VtxForInvMass = vtxFitter.vertex( muon_ttks );

      Measurement1D MassWErr(mumu.M(),-9999.);
      if ( field->nominalValue() > 0 )
          MassWErr = massCalculator.invariantMass( VtxForInvMass, muMasses );
      else
          mumuVertex = TransientVertex();                      // with no arguments it is invalid

      if (!(mumuVertex.isValid()))
          continue;

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
      TrajectoryStateClosestToPoint mu1TS = muon_ttks[0].impactPointTSCP();
      TrajectoryStateClosestToPoint mu2TS = muon_ttks[1].impactPointTSCP();

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
      ppdlPV = distXY.value()*cosAlpha * mumucand.mass()/pperp.Perp();

      GlobalError v1e = (Vertex(mumuVertex)).error();
      GlobalError v2e = thePrimaryV.error();
      AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
      ppdlErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*mumucand.mass()/(pperp.Perp2());

      AlgebraicVector3 vDiff;
      vDiff[0] = vdiff.x(); vDiff[1] = vdiff.y(); vDiff[2] = 0 ;
      ctauPV = vdiff.Perp();
      ctauErrPV = sqrt(ROOT::Math::Similarity(vDiff,vXYe)) / vdiff.Perp();


      mumucand.addUserFloat("vNChi2",vChi2/vNDF);
      mumucand.addUserFloat("vProb",vProb);
      mumucand.addUserFloat("DCA", dca );
      mumucand.addUserFloat("MassErr",MassWErr.error());
      mumucand.addUserFloat("ctauPV",ctauPV);
      mumucand.addUserFloat("ctauErrPV",ctauErrPV);
      mumucand.addUserFloat("cosAlpha",cosAlpha);
      mumucand.addUserData("thePV",Vertex(thePrimaryV));
      mumucand.addUserData("theVertex",Vertex(mumuVertex));

      trakCollection->push_back(TTCand);


    } // loop over second track
  }

  iEvent.put(std::move(trakCollection));

}


// ------------ method called once each job just before starting event loop  ------------
void
DiTrakPAT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DiTrakPAT::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiTrakPAT);
