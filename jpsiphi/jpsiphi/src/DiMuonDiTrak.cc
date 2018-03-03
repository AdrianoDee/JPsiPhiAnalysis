#include "../interface/DiMuonDiTrak.h"

DiMuonDiTrakPAT::DiMuonDiTrakPAT(const edm::ParameterSet& iConfig):
dimuons_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiMuons"))),
ditraks_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiTraks"))),
thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
DiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonCuts")),
DiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiTrakCuts")),
DiMuonDiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonDiTrakCuts"))
{
  produces<pat::CompositeCandidateCollection>();
}


DiMuonDiTrakPAT::~DiMuonDiTrakPAT()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

const pat::CompositeCandidate DiMuonDiTrakPAT::makeDiMuonTTCandidate(
                                          const pat::CompositeCandidate& dimuon,
				          const pat::CompositeCandidate& tt
                                         ){

  pat::CompositeCandidate DiMuonTCand;
  DiMuonTCand.addDaughter(dimuon,"dimuon");
  DiMuonTCand.addDaughter(tt,"ditrak");
  DiMuonTCand.setVertex(dimuon.vertex());
  DiMuonTCand.setCharge(tt.charge());

  reco::Candidate::LorentzVector vDiMuonT = dimuon.p4() + tt.p4();
  DiMuonTCand.setP4(vDiMuonT);

  return DiMuonTCand;

}


void
DiMuonDiTrakPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  typedef Candidate::LorentzVector LorentzVector;

  edm::Handle<pat::CompositeCandidateCollection> dimuon;
  iEvent.getByToken(dimuons_,dimuon);

  edm::Handle<pat::CompositeCandidateCollection> ditrak;
  iEvent.getByToken(ditraks_,ditrak);

  edm::Handle<BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  BeamSpot bs = *theBeamSpot;

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  std::unique_ptr<pat::CompositeCandidateCollection> mmttCollection(new pat::CompositeCandidateCollection);

  float DiMuonDiTrakMassMax_ = DiMuonDiTrakMassCuts_[1];
  float DiMuonDiTrakMassMin_ = DiMuonDiTrakMassCuts_[0];

  float DiMuonMassMax_ = DiMuonMassCuts_[1];
  float DiMuonMassMin_ = DiMuonMassCuts_[0];

  float DiTrakMassMax_ = DiTrakMassCuts_[1];
  float DiTrakMassMin_ = DiTrakMassCuts_[0];

  float vProb, vNDF, vChi2, minDz = 999999.;
  float cosAlpha, ctauPV, ctauErrPV;
  float l_xy, lErr_xy;

  for (pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuon->begin(); dimuonCand != dimuon->end(); ++dimuonCand)
  {
    if ( !(dimuonCand->mass() < DiMuonMassMax_  && dimuonCand->mass() > DiMuonMassMin_) )
      continue;

    for (pat::CompositeCandidateCollection::const_iterator ditrakCand = ditrak->begin(); ditrakCand != ditrak->end(); ++ditrakCand)
    {
      if ( !(ditrakCand->mass() < DiTrakMassMax_  && ditrakCand->mass() > DiTrakMassMin_) )
        continue;

        vProb = -1.0; vNDF = -1.0; vChi2 = -1.0;
        cosAlpha = -1.0; ctauPV = -1.0; ctauErrPV = -1.0;
        minDz = 999999.;

        pat::CompositeCandidate mmttCand = makeDiMuonTTCandidate(*dimuonCand, *&TTCand);

        if ( !(mmttCand->mass() < DiMuonDiTrakMassMax_  && mmttCand->mass() > DiMuonDiTrakMassMin_) )
          continue;

        const pat::PackedCandidate *trakP = dynamic_cast<const pat::PackedCandidate*>(ditrakCand->daughter("trakP"));
        const pat::PackedCandidate *trakN = dynamic_cast<const pat::PackedCandidate*>(ditrakCand->daughter("trakN"));

        const pat::Muon *muonP = (dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muonP") ) )->innerTrack();
        const pat::Muon *muonN = (dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muonN") ) )->innerTrack();

        std::vector<reco::TransientTrack> MuMuTT;
        MuMuTT.push_back((*theB).build(muonP));
        MuMuTT.push_back((*theB).build(muonN));

        if(!trakP->hasTrackDetails())
          continue;
        else if(trakP->bestTrack())
          MuMuTT.push_back((*theB).build(*(trakP->bestTrack())));
        else
          MuMuTT.push_back((*theB).build((trakP->pseudoTrack())));


        if(!trakN->hasTrackDetails())
          continue;
        else if(trakN->bestTrack())
          MuMuTT.push_back((*theB).build(*(trakN->bestTrack())));
        else
          MuMuTT.push_back((*theB).build((trakN->pseudoTrack())));

        TransientVertex mmttVertex = vtxFitter.vertex(MuMuTT);
        CachingVertex<5> VtxForInvMass = vtxFitter.vertex( MuMuTT );

        Measurement1D MassWErr(mPos->mass(),-9999.);
        if ( field->nominalValue() > 0 )
            MassWErr = massCalculator.invariantMass( VtxForInvMass, mmMasses );
        else
            mmttVertex = TransientVertex();                      // with no arguments it is invalid

        if (!(mmttVertex.isValid()))
            continue;

        vChi2 = mmttVertex.totalChiSquared();
        vNDF  = mmttVertex.degreesOfFreedom();
        vProb(TMath::Prob(vChi2,(int)vNDF));

        //Vertex parameters
        TVector3 vtx,vtx3D;
        TVector3 pvtx,pvtx3D;
        VertexDistanceXY vdistXY;

        vtx.SetXYZ(mmttVertex.position().x(),mmttVertex.position().y(),0);
        vtx3D.SetXYZ(mmttVertex.position().x(),mmttVertex.position().y(),mmttVertex.position().z());
        TVector3 pperp(mumu.px(), mumu.py(), 0);
        TVector3 pperp3D(mumu.px(), mumu.py(), mumu.pz());
        AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
        AlgebraicVector3 vpperp3D(pperp.x(),pperp.y(),pperp.z());

        thePrimaryV = dimuontt->daughter("dimuon").userData<Vertex>("thePV");

        //Lifetime calculations
        pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
        TVector3 vdiff = vtx - pvtx;
        cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());

        Measurement1D distXY = vdistXY.distance(Vertex(mmttVertex), thePrimaryV);
        ctauPV = distXY.value()*cosAlpha * mmttCand.mass()/pperp.Perp();

        GlobalError v1e = (Vertex(mmttVertex)).error();
        GlobalError v2e = thePrimaryV.error();
        AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
        ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*mmttCand.mass()/(pperp.Perp2());

        AlgebraicVector3 vDiff;
        vDiff[0] = vdiff.x(); vDiff[1] = vdiff.y(); vDiff[2] = 0 ;
        l_xy = vdiff.Perp();
        lErr_xy = sqrt(ROOT::Math::Similarity(vDiff,vXYe)) / vdiff.Perp();


        mmttCand.addUserFloat("vNChi2",vChi2/vNDF);
        mmttCand.addUserFloat("vProb",vProb);
        mmttCand.addUserFloat("MassErr",MassWErr.error());
        mmttCand.addUserFloat("ctauPV",ctauPV);
        mmttCand.addUserFloat("ctauErrPV",ctauErrPV);
        mmttCand.addUserFloat("lxy",l_xy);
        mmttCand.addUserFloat("lErrxy",lErr_xy);
        mmttCand.addUserFloat("cosAlpha",cosAlpha);
        mmttCand.addUserData("thePV",Vertex(thePrimaryV));
        mmttCand.addUserData("theVertex",Vertex(mmttVertex));

        mmttCollection->push_back(mmttCand);

    }
  }

  iEvent.put(std::move(mmttCollection));

}


// ------------ method called once each job just before starting event loop  ------------
void
DiMuonDiTrakPAT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DiMuonDiTrakPAT::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakPAT);
