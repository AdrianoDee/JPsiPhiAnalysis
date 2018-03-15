// -*- C++ -*-
//
// Package:    DiMuonDiTrakKinematicFit
// Class:      DiMuonDiTrakKinematicFit
//
/**\class DiMuonDiTrakKinematicFit jpsiphi/jpsiphi/src/DiMuonDiTrakKinematicFit.cc

 Description: performs vertex kinematical fit for DiMuon-DiTrack

**/
//
// Original Author:  Adriano Di Florio
//         Created:  Based on Alberto Sanchez Hernandez PsiTrkTrk Code

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

///For kinematic fit:
#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include <boost/foreach.hpp>
#include <string>

//
// class declaration
//

class DiMuonDiTrakKinematicFit : public edm::EDProducer {
   public:
      explicit DiMuonDiTrakKinematicFit(const edm::ParameterSet&);
      ~DiMuonDiTrakKinematicFit() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginJob() override ;
      void produce(edm::Event&, const edm::EventSetup&) override;
      void endJob() override ;

      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

// ----------member data ---------------------------
  edm::EDGetTokenT<pat::CompositeCandidateCollection> DiMuonTTCand_;
  double DiMuonMass_;
  std::vector<double> DiMuonTTMassCuts_;
  std::vector<double> MassTraks_;
  std::string Product_name_;
  std::vector<double> massCands_;

  InvariantMassFromVertex massCalculator;

  template<typename T>
  struct GreaterByVProb {
     typedef T first_argument_type;
     typedef T second_argument_type;
     bool operator()( const T & t1, const T & t2 ) const {
        return t1.userFloat("vProb") > t2.userFloat("vProb");
     }
  };
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DiMuonDiTrakKinematicFit::DiMuonDiTrakKinematicFit(const edm::ParameterSet& iConfig) {
  DiMuonTTCand_       = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiMuonDiTrak"));
  DiMuonMass_         = iConfig.getParameter<double>("DiMuonMass");
  DiMuonTTMassCuts_   = iConfig.getParameter<std::vector<double>>("DiMuonTrakTrakMassCuts");
  MassTraks_          = iConfig.getParameter<std::vector<double>>("MassTraks");
  Product_name_       = iConfig.getParameter<std::string>("Product");
  massCands_          = iConfig.getParameter<std::vector<double>>("CandsMasses");

// kinematic refit collections
  produces<pat::CompositeCandidateCollection>(Product_name_);

// now do what ever other initialization is needed
}

DiMuonDiTrakKinematicFit::~DiMuonDiTrakKinematicFit() {
// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void DiMuonDiTrakKinematicFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  vector<double> fourMasses;
  fourMasses.push_back( massCands_[0] );
  fourMasses.push_back( massCands_[1] );
  fourMasses.push_back( massCands_[2] );
  fourMasses.push_back( massCands_[3] );

  // Grab paramenters
  edm::Handle<pat::CompositeCandidateCollection> PsiTCandHandle;
  iEvent.getByToken(DiMuonTTCand_, PsiTCandHandle);

// Kinematic refit collection
  std::unique_ptr< pat::CompositeCandidateCollection > mmttCollectionVertex(new pat::CompositeCandidateCollection);

// Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  const MagneticField* field = magneticField.product();

  int indexPsiT=-1;

  float vProb, vNDF, vChi2;
  float cosAlpha, ctauPV, ctauErrPV;
  float l_xy, lErr_xy;

  for (pat::CompositeCandidateCollection::iterator dimuontt=PsiTCandHandle->begin(); dimuontt!=PsiTCandHandle->end(); ++dimuontt) {

    const pat::CompositeCandidate *dimuonC = dynamic_cast<const pat::CompositeCandidate *>(dimuontt->daughter("dimuon"));
    const pat::CompositeCandidate *ditrakC = dynamic_cast<const pat::CompositeCandidate*>(dimuontt->daughter("ditrak"));

    if(dimuonC->userFloat("vProb")<0.0)
      continue;

    indexPsiT++;

    std::vector <reco::Track> JpsiTk;

    const pat::PackedCandidate *trakP = dynamic_cast<const pat::PackedCandidate*>(ditrakC->daughter("trakP"));
    const pat::PackedCandidate *trakN = dynamic_cast<const pat::PackedCandidate*>(ditrakC->daughter("trakN"));

    JpsiTk.push_back(*( dynamic_cast<const pat::Muon*>(dimuontt->daughter("dimuon")->daughter("muonP") ) )->innerTrack());
    JpsiTk.push_back(*( dynamic_cast<const pat::Muon*>(dimuontt->daughter("dimuon")->daughter("muonN") ) )->innerTrack());

    std::vector<reco::TransientTrack> MuMuTT;
    MuMuTT.push_back((*theB).build(&JpsiTk[0]));
    MuMuTT.push_back((*theB).build(&JpsiTk[1]));

    if(!trakP->hasTrackDetails())
      continue;
    else if(trakP->bestTrack())
      MuMuTT.push_back((*theB).build(*(trakP->bestTrack()))); // K+
    else
      MuMuTT.push_back((*theB).build((trakP->pseudoTrack()))); // K+


    if(!trakN->hasTrackDetails())
      continue;
    else if(trakN->bestTrack())
      MuMuTT.push_back((*theB).build(*(trakN->bestTrack()))); // K+
    else
      MuMuTT.push_back((*theB).build((trakN->pseudoTrack()))); // K+

    const reco::Vertex thePrimaryV = *dimuonC->userData<reco::Vertex>("PVwithmuons");

    TransientVertex mmttVertex = vtxFitter.vertex(MuMuTT);
    CachingVertex<5> VtxForInvMass = vtxFitter.vertex( MuMuTT );

    Measurement1D MassWErr(dimuontt->mass(),-9999.);
    if ( field->nominalValue() > 0 )
        MassWErr = massCalculator.invariantMass( VtxForInvMass, fourMasses );
    else
        mmttVertex = TransientVertex();                      // with no arguments it is invalid

    if (!(mmttVertex.isValid()))
        continue;

    LorentzVector mumutrktrk = trakP->p4() + trakN->p4() + muonP->p4() + muonN->p4();

    vChi2 = mmttVertex.totalChiSquared();
    vNDF  = mmttVertex.degreesOfFreedom();
    vProb = TMath::Prob(vChi2,(int)vNDF);

    //Vertex parameters
    TVector3 vtx,vtx3D;
    TVector3 pvtx,pvtx3D;
    VertexDistanceXY vdistXY;

    vtx.SetXYZ(mmttVertex.position().x(),mmttVertex.position().y(),0);
    vtx3D.SetXYZ(mmttVertex.position().x(),mmttVertex.position().y(),mmttVertex.position().z());
    TVector3 pperp(mumutrktrk.px(), mumutrktrk.py(), 0);
    TVector3 pperp3D(mumutrktrk.px(), mumutrktrk.py(), mumutrktrk.pz());
    AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
    AlgebraicVector3 vpperp3D(pperp.x(),pperp.y(),pperp.z());

    //Lifetime calculations
    pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
    TVector3 vdiff = vtx - pvtx;
    cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());

    Measurement1D distXY = vdistXY.distance(Vertex(mmttVertex), thePrimaryV);
    ctauPV = distXY.value()*cosAlpha * dimuontt->mass()/pperp.Perp();

    GlobalError v1e = (Vertex(mmttVertex)).error();
    GlobalError v2e = thePrimaryV.error();
    AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
    ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*dimuontt->mass()/(pperp.Perp2());

    AlgebraicVector3 vDiff;
    vDiff[0] = vdiff.x(); vDiff[1] = vdiff.y(); vDiff[2] = 0 ;
    l_xy = vdiff.Perp();
    lErr_xy = sqrt(ROOT::Math::Similarity(vDiff,vXYe)) / vdiff.Perp();

    dimuontt->addUserFloat("vNChi2",vChi2/vNDF);
    dimuontt->addUserFloat("vProb",vProb);
    dimuontt->addUserFloat("MassErr",MassWErr.error());
    dimuontt->addUserFloat("ctauPV",ctauPV);
    dimuontt->addUserFloat("ctauErrPV",ctauErrPV);
    dimuontt->addUserFloat("lxy",l_xy);
    dimuontt->addUserFloat("lErrxy",lErr_xy);
    dimuontt->addUserFloat("cosAlpha",cosAlpha);
    dimuontt->addUserData("thePV",Vertex(thePrimaryV));
    dimuontt->addUserData("theVertex",Vertex(mmttVertex));

    mmttCollectionVertex->push_back(dimuontt);

  }
// End kinematic fit

// now sort by vProb
  DiMuonDiTrakKinematicFit::GreaterByVProb<pat::CompositeCandidate> vPComparator;
  std::sort(mmttCollectionVertex->begin(),mmttCollectionVertex->end(),vPComparator);

  iEvent.put(std::move(mmttCollectionVertex),Product_name_);
  // std::cout << "mmttCollectionVertex size: "<< mmttCollectionVertex->size()<<std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
void DiMuonDiTrakKinematicFit::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiMuonDiTrakKinematicFit::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiMuonDiTrakKinematicFit::beginRun(edm::Run&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void DiMuonDiTrakKinematicFit::endRun(edm::Run&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiMuonDiTrakKinematicFit::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiMuonDiTrakKinematicFit::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMuonDiTrakKinematicFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakKinematicFit);
