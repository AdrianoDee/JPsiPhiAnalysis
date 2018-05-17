// -*- C++ -*-
//
// Package:    MuMuProducerPAT
// Class:      MuMuProducerPAT
//
/**\class MuMuProducerPAT MuMuProducerPAT.cc myProducers/MuMuProducerPAT/src/MuMuProducerPAT.cc

Description: <one line class summary>
Make rootTuple for JPsiKK reconstruction

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Keith Ulmer   keith.ulmer@colorado.edu
//
//


/// system include files
#include <memory>

/// user include files
#include "../interface/MuMuProducerPAT.h"
#include "../interface/VertexReProducer.h"
//#include "DataFormats/Candidate/interface/OverlapChecker.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

//#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

/// for 53x
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include <vector>
#include <utility>

#include "TMath.h"
#include "Math/VectorUtil.h"

/// useless so far
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
//#include "HepMC/GenVertex.h"
//#include <HepMC/GenVertex.h>
//#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

///
/// constants, enums and typedefs
///

typedef math::Error<3>::type CovarianceMatrix;

const ParticleMass muon_mass = 0.10565837; //pdg mass
const ParticleMass kaon_mass = 0.493667; //pdg mass
ParticleMass JPsi_mass = 3.096916;
const ParticleMass Phi_mass = 1.0194;


/// Setting insignificant mass sigma to avoid singularities in the covariance matrix.
float small_sigma = muon_mass*1.e-6;
//float small_sigma = kaon_mass*1.e-6; /// SEMRA

///
/// static data member definitions
///

///
/// constructors and destructor
///

DiMuonProducerPAT::DiMuonProducerPAT(const edm::ParameterSet& iConfig):

hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults",edm::InputTag("TriggerResults::HLT")) ),
inputGEN_(iConfig.getUntrackedParameter<edm::InputTag>("inputGEN",edm::InputTag("genParticles"))),
vtxSample__(iConfig.getUntrackedParameter<std::string>("vtxSample_",std::string("offlinePrimaryVertices"))),

jspiMassCuts_(iConfig.getParameter<std::vector<float>>("JPsiMassCuts")),
psiMassCuts_(iConfig.getParameter<std::vector<float>>("PsiMassCuts")),

doData_( iConfig.getUntrackedParameter<bool>("DoDataAnalysis", true) ),
doMC__( iConfig.getUntrackedParameter<bool>("DoMonteCarloTree", true) ),

addCommonVertex_(iConfig.getParameter<bool>("addCommonVertex")),
resolveAmbiguity_(iConfig.getParameter<bool>("resolvePileUpAmbiguity")),

addMCTruth_(iConfig.getParameter<bool>("addMCTruth")),
MCParticle_( iConfig.getUntrackedParameter<int>("MonteCarloParticleId", 20443) ), /// 20443 X, 100443 Psi(2S), 9120443 X from B / decide later for X(4140)
MCExclusiveDecay__( iConfig.getUntrackedParameter<bool>("MonteCarloExclusiveDecay", true) ),
MCMother_( iConfig.getUntrackedParameter<int>("MonteCarloMotherId", 511) ), /// 511 B0 (=anti-B0), 531 B0 / decide later MCMotherId for X(4140)
MCDaughtersN_( iConfig.getUntrackedParameter<int>(" MonteCarloDaughtersN", 3) ), /// will be same

MuMinPixHits_(iConfig.getUntrackedParameter<int>("MinNumMuPixHits", 0)),
MuMinSiHits_(iConfig.getUntrackedParameter<int>("MinNumMuSiHits", 0)),
MuMaxNormChi_(iConfig.getUntrackedParameter<float>("MaxMuNormChi2", 1000)),
MuMaxD0_(iConfig.getUntrackedParameter<float>("MaxMuD0", 1000)),

TriggerCut_(iConfig.getUntrackedParameter<bool>("TriggerCut",true)),
HLTFileters_(iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching")),
FiltersForMatching_(iConfig.getUntrackedParameter<std::vector<std::string> >("FiltersForMatching")),
Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output",true)),
{
  // revtxtrks_ = "generalTracks"; //if that is not true, we will raise an exception
  // revtxbs_ = "offlineBeamSpot";
  // genCands_ = "genParticles";

  produces<pat::CompositeCandidateCollection>();

  /// now do what ever initialization is needed


}

MuMuProducerPAT::~MuMuProducerPAT()
{
  /// do anything here that needs to be done at desctruction time
  /// (e.g. close files, deallocate resources etc.)

}

UInt_t MuMuProducerPAT::isTriggerMatched(pat::Muon* posMuon, pat::Muon* negMuon) {

  UInt_t matched = 0;  // if no list is given, is not matched

  // if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<FiltersForMatching_.size(); iTr++ ) {
    // std::cout << HLTFilters_[iTr] << std::endl;
    const pat::TriggerObjectStandAloneCollection mu1HLTMatches = posMuon->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    const pat::TriggerObjectStandAloneCollection mu2HLTMatches = negMuon->triggerObjectMatchesByFilter(HLTFilters_[iTr]);
    if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr);
    // if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) std::cout << std::endl << HLTFilters_[iTr] << std::endl;
  }

  return matched;
}

// UInt_t MuMuProducerPAT::getTriggerBits(const edm::Event& iEvent ) {
//
//   UInt_t trigger = 0;
//
//   edm::Handle< edm::TriggerResults > triggerResults_handle;
//   iEvent.getByLabel( triggerResults_Label , triggerResults_handle);
//
//   if (triggerResults_handle.isValid()) {
//      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
//      unsigned int NTRIGGERS = HLTs_.size();
//
//      for (unsigned int i = 0; i < NTRIGGERS; i++) {
//         for (int version = 1; version < 20; version++) {
//            std::stringstream ss;
//            ss << HLTs_[i] << "_v" << version;
//            unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
//            if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
//               trigger += (1<<i);
//               break;
//            }
//         }
//      }
//    } else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
//
//    return trigger;
// }


///
/// member functions
///

/// ------------ method called to for each event  ------------
void MuMuProducerPAT::produce(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;
  using namespace reco;

  /// get event content information
  int debug = 0;

  if(Debug_) std::cout <<"Debug : " << ++debug << std::endl;

  std::auto_ptr<pat::CompositeCandidateCollection> oniaOutput(new pat::CompositeCandidateCollection);

  bool decayChainOK = false;
  int runNum = iEvent.id().run();
  int evtNum = iEvent.id().event();
  int lumiNum = iEvent.id().luminosityBlock();


  bool hasRequestedTrigger = false;
  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  // /// first get HLT results
  // unsigned int triggerBit = -1;
  //
  // std::map<string,int> HLTPreScaleMap;
  // edm::Handle<edm::TriggerResults> hltresults;
  // try {
  //   iEvent.getByLabel(hlTriggerResults_, hltresults);
  // }
  // catch ( ... ) {
  //   std::cout << "Couldn't get handle on HLT Trigger!" << std::endl;
  // }
  // if (!hltresults.isValid()) {
  //   std::cout << "No Trigger Results!" << std::endl;
  // }
  // else {
  //
  //   int ntrigs = hltresults->size();
  //
  //   if (ntrigs==0)
  //     std::cout << "No trigger name given in TriggerResults of the input " << std::endl;
  //
  //   triggerBit = getTriggerBits(iEvent);
  //
  // } /// end valid trigger


  Vertex thePrimaryVtx, theBeamSpotVtx;
  math::XYZPoint RefVtx;

  Int_t thePrimaryVtx_multiplicity = -1 ;

  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
  if ( beamSpotHandle.isValid() ) {
    beamSpot = *beamSpotHandle;
    theBeamSpotVtx = Vertex(beamSpot.position(), beamSpot.covariance3D());
  }
  else std::cout << "No beam spot available from EventSetup" << std::endl;

  Handle<VertexCollection> recVtxs;
  iEvent.getByLabel(vtxSample_, recVtxs);
  unsigned int nVtxTrks = 0;
  if ( recVtxs->begin() != recVtxs->end() ) {
    thePrimaryVtx_multiplicity = recVtxs->size() ;

    if (resolveAmbiguity_) {
      //thePrimaryVtx = Vertex(*(recVtxs->begin()));
      //std::cout <<"here" <<std::endl;
      thePrimaryVtx = *(recVtxs->begin());
    }
    else {
      for ( reco::VertexCollection::const_iterator vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx) {
        if (nVtxTrks < vtx->tracksSize() ) {
          nVtxTrks = vtx->tracksSize();
          thePrimaryVtx = Vertex(*vtx);
        }
      }
    }
  } else {
    thePrimaryVtx = Vertex(beamSpot.position(), beamSpot.covariance3D());
    thePrimaryVtx_multiplicity = 1 ;
  }

  std::cout <<"Debug : " << ++debug << std::endl;

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  RefVtx = thePrimaryVtx.position(); /// reference primary vertex choosen
  reco::Vertex pV = thePrimaryVtx;
  int n_pV = thePrimaryVtx_multiplicity ;
  float priVtx_VProb = ChiSquaredProbability( (float)(thePrimaryVtx.chi2()), (float)(thePrimaryVtx.ndof())) ;

  // VertexHigherPtSquared vertexHigherPtSquared ;
  // tracksPtSq_pV = vertexHigherPtSquared.sumPtSquared(thePrimaryVtx) ;

  /// /// /// /// /// /// /// /// /// /// /// /// /// ///
  /// MUONS
  /// /// /// /// /// /// /// /// /// /// /// /// /// ///

  Handle< std::vector<pat::Muon> > thePATMuonHandle;
  iEvent.getByLabel("patMuonsWithTrigger", thePATMuonHandle);


  /// /// /// /// /// /// /// /// /// /// /// /// /// ///
  /// MC Truth
  /// /// /// /// /// /// /// /// /// /// /// /// /// ///
  //
  //   if (doMC_) {
  //     /*
  //     // Get generated event
  //     //Handle<edm::HepMCProduct> hepEv;
  //     //iEvent.getByLabel("generator", hepEv);
  //     Handle<GenEventInfoProduct> genEvtInfo;
  //     iEvent.getByLabel("generator", genEvtInfo);
  //
  //     //const HepMC::GenEvent *myGenEvent = hepEv->GetEvent();
  //     const HepMC::GenEvent *myGenEvent = genEvtInfo->GetEvent();
  //     n_genEvtVtx = myGenEvent->vertices_size() ;
  //
  //     HepMC::GenVertex* primaryGenVtx = *(myGenEvent->vertices_begin()) ;
  //
  //     genEvtVtx_X->push_back( primaryGenVtx->point3d().x() );
  //     genEvtVtx_Y->push_back( primaryGenVtx->point3d().y() );
  //     genEvtVtx_Z->push_back( primaryGenVtx->point3d().z() );
  //     //genEvtVtx_XE = (primaryGenVtx->xError()) ;
  //     //genEvtVtx_YE = (primaryGenVtx->yError()) ;
  //     //genEvtVtx_ZE = (primaryGenVtx->zError()) ;
  //     //genEvtVtx_NormChi2 = (primaryGenVtx->normalizedChi2()) ;
  //     //genEvtVtx_Chi2 = primaryGenVtx->chi2() ;
  //     //genEvtVtx_VProb = ChiSquaredProbability( (float)(primaryGenVtx.chi2()), (float)(primaryGenVtx.ndof())) ;
  //     genEvtVtx_particles->push_back( primaryGenVtx->particles_out_size() );
  //     */
  //
  //     Handle< std::vector< PileupSummaryInfo > >  PupInfo;
  //     iEvent.getByLabel("addPileupInfo", PupInfo);
  //     std::vector<PileupSummaryInfo>::const_iterator PVI;
  //     if (Debug_) std::cout <<"\nBunchXing multiplicity = " <<PupInfo->size() <<std::endl ;
  //     for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
  //     if (Debug_) std::cout <<"Pileup Information: bunchXing, nvtx: " <<PVI->getBunchCrossing() <<" " <<PVI->getPU_NumInteractions() <<std::endl;
  //
  //     Handle<GenParticleCollection> genParticles;
  //     iEvent.getByLabel("genParticles", genParticles);
  //     if (Debug_) std::cout << "############### GenParticles Analysis ###############" << std::endl;
  //     float jpsiPx=0., jpsiPy=0., jpsiPz=0.;
  //     float  mupPx=0., mupPy=0., mupPz=0., mumPx=0., mumPy=0., mumPz=0.;
  //     float phiPx=0., phiPy=0., phiPz=0.;
  //     float  kpPx=0., kpPy=0., kpPz=0., kmPx=0., kmPy=0., kmPz=0.;
  //     //float pionPx=0., pionPy=0., pionPz=0., kaonPx=0., kaonPy=0., kaonPz=0.;
  //     //int pionCh=0, kaonCh=0 ;
  //
  //     for (size_t i = 0; i < genParticles->size(); ++ i) {
  //       nMCAll++;
  //       const reco::GenParticle &p = (*genParticles)[i];
  //       int pdgid = p.pdgId() ;
  //       int dauNum = p.numberOfDaughters();
  //       MCPdgIdAll->push_back( pdgid );
  //       MCDanNumAll->push_back( dauNum );
  //
  //       if ( MCExclusiveDecay_ ) {
  //         /// check if there is a MCMother which has MCDaughtersN daughters
  //         if ( abs(pdgid) == MCMother  &&  dauNum == MCDaughtersN ) {
  //           bool mumuOK = false;
  //           bool kkOK = false;
  //           //bool pionOK = false, kaonOK = false;
  //
  //           for (int j=0; j<dauNum; ++j) {
  //             const Candidate *dau = p.daughter(j);
  //             if (Debug_) std::cout << "dauPdgId = " << dau->pdgId() << std::endl;
  //
  //             /// check if one of B0 daughters is a psi(nS) whitch has 2 muons as daughters /// SEMRA ask again !!!
  //             int mumuId = 0 ;
  //             if (skipJPsi) /// SEMRA cleaned skipPsi2S
  //             if (Debug_) std::cout <<"Skipping J/psi!" <<std::endl ; /// SEMRA cleaned skipPsi2S
  //             //else if (skipPsi2S) /// SEMRA
  //             //  mumuId = 443 ; /// SEMRA (JPsi ID)
  //
  //             if ( ((skipJPsi) && (dau->pdgId() == mumuId)) ||
  //             ((!skipJPsi) && (dau->pdgId()%1000 == 443)) ) {
  //               jpsiPx = dau->px(); jpsiPy = dau->py(); jpsiPz = dau->pz();
  //               int jpsiDauNum = dau->numberOfDaughters();
  //               if (Debug_) std::cout << "jpsiDauNum = " << jpsiDauNum << std::endl;
  //               int muNum = 0;
  //               for (int k=0; k<jpsiDauNum; ++k) {
  //                 const Candidate *grandDau = dau->daughter(k);
  //                 if (Debug_)  std::cout << "grandDauPdgId = " << grandDau->pdgId() << std::endl;
  //                 if ( abs(grandDau->pdgId()) == 13 ) {
  //                   muNum++;
  //                   if (grandDau->pdgId() < 0) {
  //                     mupPx = grandDau->px(); mupPy = grandDau->py(); mupPz = grandDau->pz();
  //                   } else {
  //                     mumPx = grandDau->px(); mumPy = grandDau->py(); mumPz = grandDau->pz();
  //                   }
  //                 }
  //               }
  //               if ( muNum == 2 ) mumuOK = true ;
  //
  //             } /// end check if one of the MCMother daughters is a J/Psi or psi'
  //
  //             /// for Phi
  //             phiPx = dau->px(); phiPy = dau->py(); phiPz = dau->pz();
  //             int phiDauNum = dau->numberOfDaughters();
  //             if (Debug_) std::cout << "phiDauNum = " << phiDauNum << std::endl;
  //             int kNum = 0;
  //             for (int n=0; n<phiDauNum; ++n) {
  //               const Candidate *grandDau = dau->daughter(n);
  //               if (Debug_)  std::cout << "grandDauPdgId = " << grandDau->pdgId() << std::endl;
  //               if ( abs(grandDau->pdgId()) == 321 ) {
  //                 kNum++;
  //                 if (grandDau->pdgId() < 0) {
  //                   kpPx = grandDau->px(); kpPy = grandDau->py(); kpPz = grandDau->pz();
  //                 } else {
  //                   kmPx = grandDau->px(); kmPy = grandDau->py(); kmPz = grandDau->pz();
  //                 }
  //               }
  //             }
  //             if ( kNum == 2 ) kkOK = true ;
  //
  //
  //             /*else if ( abs(dau->pdgId()) == 211 ) { // check if one of B0 daughters is a pion /// SEMRA ask again !!!
  //             pionPx = dau->px(); pionPy = dau->py(); pionPz = dau->pz();
  //             pionCh = (dau->pdgId() == 211)? 1 : -1;
  //             pionOK = true; /// SEMRA pions change with kaons for B0 ?
  //           } else if ( abs(dau->pdgId()) == 321 ) { // check if one of B0 daughters is a kaon /// SEMRA ask again !!!
  //           kaonPx = dau->px(); kaonPy=dau->py(); kaonPz=dau->pz();
  //           kaonCh = (dau->pdgId() == 321)? 1 : -1;
  //           kaonOK = true;
  //         }*/
  //
  //       } /// end loop on MCMother daughters
  //
  //       if (Debug_) std::cout << "mumuOK = " << mumuOK << ", kkOK = " << kkOK << std::endl;
  //       if ( mumuOK && kkOK ) {
  //         if (Debug_) {
  //           std::cout <<"\nnumber of X mothers = " <<p.numberOfMothers() <<std::endl ;
  //           std::cout <<"X mother pdgID = " <<p.mother(0)->pdgId() <<std::endl ;
  //         }
  //         ++nMCX ;
  //         PriVtxGen_X->push_back( p.vx() ) ;
  //         PriVtxGen_Y->push_back( p.vy() ) ;
  //         PriVtxGen_Z->push_back( p.vz() ) ;
  //         PriVtxGen_VProb->push_back( p.vertexNormalizedChi2() ) ;
  //         PriVtxGen_Chi2->push_back( p.vertexChi2() ) ;
  //         PriVtxGen_Ndof->push_back( p.vertexNdof() ) ;
  //
  //         Bool_t status = kTRUE ;
  //         const Candidate *x_ancestor = p.mother(0) ; /// a particle can have several mothers
  //         Int_t n_ancestors = 1 ;
  //         while ( status ) {
  //           if ( abs(x_ancestor->pdgId()) <= 8 || x_ancestor->pdgId() == 21 || x_ancestor->status() == 3 ) {
  //             status = kFALSE ;
  //             if (Debug_) std::cout <<"X ancestor ID = " <<x_ancestor->pdgId() <<std::endl ;
  //             genEvtVtx_X->push_back( x_ancestor->daughter(0)->vx() ) ;
  //             genEvtVtx_Y->push_back( x_ancestor->daughter(0)->vy() ) ;
  //             genEvtVtx_Z->push_back( x_ancestor->daughter(0)->vz() ) ;
  //             genEvtVtx_particles->push_back( x_ancestor->numberOfDaughters() ) ;
  //             n_XAncestors->push_back( n_ancestors ) ;
  //           }
  //           else {
  //             x_ancestor = x_ancestor->mother(0) ;
  //             n_ancestors++ ;
  //           }
  //         }
  //
  //         MCJPsiPx->push_back(jpsiPx); MCJPsiPy->push_back(jpsiPy); MCJPsiPz->push_back(jpsiPz);
  //         MCmupPx->push_back(mupPx); MCmupPy->push_back(mupPy); MCmupPz->push_back(mupPz);
  //         MCmumPx->push_back(mumPx); MCmumPy->push_back(mumPy); MCmumPz->push_back(mumPz);
  //         MCPhiPx->push_back(phiPx); MCPhiPy->push_back(phiPy); MCPhiPz->push_back(phiPz);
  //         MCkpPx->push_back(kpPx); MCkpPy->push_back(kpPy); MCkpPz->push_back(kpPz);
  //         MCkmPx->push_back(kmPx); MCkmPy->push_back(kmPy); MCkmPz->push_back(kmPz);
  //         //MCpionPx->push_back(pionPx); MCpionPy->push_back(pionPy); MCpionPz->push_back(pionPz);
  //         //MCkaonPx->push_back(kaonPx); MCkaonPy->push_back(kaonPy); MCkaonPz->push_back(kaonPz);
  //         //MCpionCh->push_back(pionCh) ; MCkaonCh->push_back(kaonCh) ;
  //         decayChainOK = true;
  //         MCPx->push_back( p.px() );
  //         MCPy->push_back( p.py() );
  //         MCPz->push_back( p.pz() );
  //       }
  //       if (Debug_) std::cout << "decayChainOK = " << decayChainOK << std::endl;
  //     } // if ( abs(pdgid) == MCMother  &&  dauNum == 3 )
  //   } // if ( !MCExclusiveDecay_ )
  //
  // } // for (size_t i = 0; i < genParticles->size(); ++ i)
  // } // if (doMC_)
  //

  /// reconstruction only for events with B decaying in psi(nS)+Pi+K /// SEMRA JPsiPhi !!!

  if ( (doMC_ && !MCExclusiveDecay_) || (doMC_ && (MCExclusiveDecay_ && decayChainOK)) || doData_ ) {

    bool isEventWithInvalidMu = false;

    if (Debug_) std::cout << "Starting event with " << thePATMuonHandle->size() << " muons" << std::endl;

    if ((thePATMuonHandle->size()) > 10000)
      std::cout << "Too many Muons: " << thePATMuonHandle->size() << std::endl;
    else //if (thePATMuonHandle->size() >= 2) { // check
      if (thePATMuonHandle->size() >= 2  && (hasRequestedTrigger || !TriggerCut_)) {

        if (Debug_) std::cout <<"============================  Evt: " <<evtNum <<" accept event with 2 mu and trigger ==============================================" <<std::endl;

        int posMuonType, negMuonType, posMuonTrackType, MuonTrackType;
        float posMuonDzVtx, posMuonDxyVtx;

        int nMatchedStationsPos, nMatchedStationsNeg,nOverlapMusPos, nOverlapMusNeg, nSharingSegWithPos, nSharingSegWithNeg;


        float mumuVProb;

        /// get MuMu cands
        for ( std::vector<pat::Muon>::const_iterator posMuon = thePATMuonHandle->begin(); posMuon != thePATMuonHandle->end(); ++posMuon ) {


          /// push back all muon information
          const reco::Muon* recoPosMuon = dynamic_cast<const reco::Muon * >(posMuon->originalObject());
          // muPx->push_back(recoPosMuon->px());
          // muPy->push_back(recoPosMuon->py());
          // muPz->push_back(recoPosMuon->pz());
          // muCharge->push_back(recoPosMuon->charge());

          if (recoPosMuon->charge()<=0.0) continue;

          if (recoPosMuon->pt()<=0.5) continue;

          if (recoPosMuon->track().isNull()) continue;

          if (recoPosMuon->bestTrackRef().isNull()) continue;

          if (recoPosMuon->innerTrack().isNull() && recoPosMuon->outerTrack().isNull() && recoPosMuon->globalTrack().isNull()) continue;

          if ( !(recoPosMuon->track()->hitPattern().trackerLayersWithMeasurement()) ) {
            isEventWithInvalidMu = true;
            if (Debug_) std::cout <<"evt:" <<evtNum <<" problem with trackerLayersWithMeasurement" <<std::endl;
            continue ;
          }
          if ( !(recoPosMuon->track()->hitPattern().pixelLayersWithMeasurement()) ) {
            if (Debug_) std::cout <<"evt:" <<evtNum <<" problem with pixelLayersWithMeasurement" <<std::endl;
            isEventWithInvalidMu = true;
            continue ;
          }

          /// Pixel Cuts On
          if (recoPosMuon->track()->hitPattern().numberOfValidPixelHits() < MuMinPixHits_
          || recoPosMuon->track()->hitPattern().numberOfValidStripHits() < MuMinSiHits_
          || recoPosMuon->track()->chi2()/recoPosMuon->track()->ndof() > MuMaxNormChi_
          || fabs(recoPosMuon->track()->dxy(RefVtx)) > MuMaxD0_) {
            continue ;
          }

          TrackRef muPosTrack = posMuon->track() ; //Needed for Refitting

          if ( muPosTrack.isNull() )
          continue ;

          const reco::Track recoPosMuonTrack = *(recoPosMuon->bestTrack());

          posMuonTrackType = recoPosMuon->MuonTrackType();
          posMuonType = (int)(recoPosMuon->type());

          const reco::Track* recoPosMuonInTrack = 0, *recoPosMuonOutTrack = 0, *recoPosMuonGlobTrack = 0;

          if (!(recoPosMuon->innerTrack().isNull()))
          recoPosMuonInTrack = (recoPosMuon->track()).get();
          if (!(recoPosMuon->outerTrack().isNull()))
          recoPosMuonOutTrack = (recoPosMuon->outerTrack()).get();
          if (!(recoPosMuon->globalTrack().isNull()))
          recoPosMuonGlobTrack = (recoPosMuon->globalTrack()).get();

          posMuonDzVtx = recoPosMuon->bestTrackRef()->dz(RefVtx);
          posMuonDxyVtx = recoPosMuon->bestTrackRef()->dxy(RefVtx);

          nMatchedStationsPos = recoPosMuon->numberOfMatchedStations();

          ////////////////// Muons Overlap Checks //////////////////

          int nOverlapMus = 0, nSharingSegWith = -1;
          int nSegments1 = recoPosMuon->numberOfMatches(reco::Muon::SegmentArbitration);

          for ( std::vector<pat::Muon>::const_iterator otherMuon = posMuon+1; otherMuon != thePATMuonHandle->end(); ++otherMuon) {

            const reco::Muon* recoOtherMuon = dynamic_cast<const reco::Muon*>(otherMuon->originalObject());
            if ( isSameMuon(*recoPosMuon, *recoOtherMuon)) continue;
            if ( !muon::isGoodMuon(*recoOtherMuon, muon::TMOneStationTight) ) continue;
            /// geometric overlap
            if ( muon::overlap( *recoPosMuon, *recoOtherMuon ) )
            nOverlapMus++ ;
            /// shared segments
            int nSegments2 = recoOtherMuon->numberOfMatches(reco::Muon::SegmentArbitration);

            if (nSegments2 == 0 || nSegments1 == 0) continue;

            float sf = muon::sharedSegments(*recoPosMuon, *recoOtherMuon) / std::min<float>(nSegments1, nSegments2);

            if (sf > 0.5) {
              nSharingSegWith = 0;
              if ( !isBetterMuon(*recoPosMuon, *recoOtherMuon) )
              nSharingSegWith++ ;
            }
          }

          nOverlapMusPos = nOverlapMus;
          nSharingSegWithPos = nSharingSegWith;

          // muNOverlap->push_back( nOverlapMus ) ;
          // muNSharingSegWith->push_back( nSharingSegWith ) ;

          ////////////////// check for muon2 //////////////////
          for ( std::vector<pat::Muon>::const_iterator negMuon = posMuon+1; negMuon != thePATMuonHandle->end(); ++negMuon) {

            if(negMuon->charge() >= 0) continue ;

            if(negMuon->pt() <= 0.5) continue ;

            const reco::Muon* recoNegMuon = dynamic_cast<const reco::Muon *>(negMuon->originalObject()) ;

            if (muon::overlap(*recoPosMuon, *recoNegMuon) )
            continue ;

            if (recoNegMuon->track().isNull()) continue;

            if (recoNegMuon->bestTrackRef().isNull()) continue;

            if (recoNegMuon->innerTrack().isNull() && recoNegMuon->outerTrack().isNull() && recoNegMuon->globalTrack().isNull()) continue;

            TrackRef muNegTrack = negMuon->track() ;
            if ( muNegTrack.isNull() )
            continue ;

            /// cuts on muon2
            if (recoNegMuon->track()->hitPattern().numberOfValidPixelHits() < MuMinPixHits_
            || recoNegMuon->track()->hitPattern().numberOfValidStripHits() < MuMinSiHits_
            || recoNegMuon->track()->chi2()/recoPosMuon->track()->ndof() > MuMaxNormChi_
            || fabs(recoNegMuon->track()->dxy(RefVtx)) > MuMaxD0_)

            continue ;

            nMatchedStationsNeg = recoNegMuon->numberOfMatchedStations();

            ////////////////// Muons Overlap Checks //////////////////

            int nOverlapMus = 0, nSharingSegWith = -1;
            int nSegments1 = recoNegMuon->numberOfMatches(reco::Muon::SegmentArbitration);

            for ( std::vector<pat::Muon>::const_iterator otherMuon = posMuon+1; otherMuon != thePATMuonHandle->end(); ++otherMuon) {

              const reco::Muon* recoOtherMuon = dynamic_cast<const reco::Muon*>(otherMuon->originalObject());
              if ( isSameMuon(*recoNegMuon, *recoOtherMuon)) continue;
              if ( !muon::isGoodMuon(*recoOtherMuon, muon::TMOneStationTight) ) continue;
              /// geometric overlap
              if ( muon::overlap( *recoNegMuon, *recoOtherMuon ) )
              nOverlapMus++ ;
              /// shared segments
              int nSegments2 = recoOtherMuon->numberOfMatches(reco::Muon::SegmentArbitration);

              if (nSegments2 == 0 || nSegments1 == 0) continue;

              float sf = muon::sharedSegments(*recoNegMuon, *recoOtherMuon) / std::min<float>(nSegments1, nSegments2);

              if (sf > 0.5) {
                nSharingSegWith = 0;
                if ( !isBetterMuon(*recoNegMuon, *recoOtherMuon) )
                nSharingSegWith++ ;
              }
            }

            nOverlapMusNeg = nOverlapMus;
            nSharingSegWithNeg = nSharingSegWith;

            pat::CompositeCandidate mumuCandidateidate;

            // ---- define and set candidate's 4momentum  ----
            math::XYZTLorentzVector muP, muN,mumuP4;

            muP.SetXYZM(recoPosMuon->bestTrackRef()->px(),recoPosMuon->bestTrackRef()->py(),recoPosMuon->bestTrackRef()->pz(),muon_mass);
            muN.SetXYZM(recoNegMuon->bestTrackRef()->px(),recoNegMuon->bestTrackRef()->py(),recoNegMuon->bestTrackRef()->pz(),muon_mass);
            // LorentzVector mumu;

            mumuP4 = muP + muN;

            mumuCandidate.setP4(mumuP4);
            mumuCandidate.setCharge(recoPosMuon->charge() + recoNegMuon->charge());

            float deltaRMuMu = reco::deltaR2(recoPosMuon->eta(),recoPosMuon->phi(),recoNegMuon->eta(),recoNegMuon->phi());


            ////////////////// get the MuMu information //////////////////
            TransientTrack muonPosTT( muPosTrack, &(*bFieldHandle) );
            TransientTrack muonNegTT( muNegTrack, &(*bFieldHandle) );
            KinematicParticleFactoryFromTransientTrack pFactory;

            /// initial chi2 and ndf before kinematic fits
            float chi = 0., ndf = 0.;
            std::vector<RefCountedKinematicParticle> muons; /// the final state muons produced by the KinematicParticleFactory
            muons.push_back( pFactory.particle( muonPosTT, muon_mass, chi, ndf, small_sigma));
            muons.push_back( pFactory.particle( muonNegTT, muon_mass, chi, ndf, small_sigma));
            KinematicParticleVertexFitter MuMuFitter; /// creating the vertex fitter for JPsi
            RefCountedKinematicTree MuMuVertexFitTree;
            MuMuVertexFitTree = MuMuFitter.fit(muons);

            if (!MuMuVertexFitTree->isValid())
            continue ;

            MuMuVertexFitTree->movePointerToTheTop();
            RefCountedKinematicParticle mumuCandidate_fromFit = MuMuVertexFitTree->currentParticle();
            RefCountedKinematicVertex mumuCandidate_vertex_fromFit = MuMuVertexFitTree->currentDecayVertex();
            MuMuVertexFitTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle MuPosCand_fromFit = MuMuVertexFitTree->currentParticle();
            MuMuVertexFitTree->movePointerToTheNextChild();
            RefCountedKinematicParticle MuNegCand_fromFit = MuMuVertexFitTree->currentParticle();

            float dimuon_vx_fit = mumuCandidate_vertex_fromFit->position().x();
            float dimuon_vy_fit = mumuCandidate_vertex_fromFit->position().y();
            float dimuon_vz_fit = mumuCandidate_vertex_fromFit->position().z();

            mumuVProb = ChiSquaredProbability((float)( mumuCandidate_vertex_fromFit->chiSquared()),(float)( mumuCandidate_vertex_fromFit->degreesOfFreedom()));

            if (mumuVProb < 0.001)
            continue;

            mumuChi2 = mumuCandidate_vertex_fromFit->chiSquared();

            ////////////////// fill the MuMu vectors //////////////////
            if (mumuCandidate_fromFit->currentState().mass() < JPsiMinMass  ||  mumuCandidate_fromFit->currentState().mass() > JPsiMaxMass)
            continue ;

            float dimuon_ma_fit = mumuCandidate_fromFit->currentState().mass();
            int   dimuon_ch_fit = mumuCandidate_fromFit->currentState().particleCharge();
            float dimuon_px_fit = mumuCandidate_fromFit->currentState().kinematicParameters().momentum().x();
            float dimuon_py_fit = mumuCandidate_fromFit->currentState().kinematicParameters().momentum().y();
            float dimuon_pz_fit = mumuCandidate_fromFit->currentState().kinematicParameters().momentum().z();
            float dimuon_en_fit = sqrt(dimuon_ma_fit*dimuon_ma_fit+dimuon_px_fit*dimuon_px_fit+dimuon_py_fit*dimuon_py_fit+dimuon_pz_fit*dimuon_pz_fit);

            reco::CompositeCandidate reco_ref_JPsi(dimuon_ch_fit,math::XYZTLorentzVector(dimuon_px_fit,dimuon_py_fit,dimuon_pz_fit,dimuon_en_fit),
            math::XYZPoint(dimuon_vx_fit,dimuon_vy_fit,dimuon_vz_fit),443);
            pat::CompositeCandidate pat_ref_JPsi(reco_ref_JPsi);

            //////////////////// For Lifetimes Calculations ////////////////////
            TVector3 MuMu_vtx((*mumuCandidate_vertex_fromFit).position().x(), (*mumuCandidate_vertex_fromFit).position().y(), 0) ;
            TVector3 MuMu_pperp(mumuCandidate_fromFit->currentState().globalMomentum().x(), mumuCandidate_fromFit->currentState().globalMomentum().y(), 0);
            TVector3 MuMu_vtx3D((*mumuCandidate_vertex_fromFit).position().x(), (*mumuCandidate_vertex_fromFit).position().y(), (*mumuCandidate_vertex_fromFit).position().z()) ;
            TVector3 MuMu_pperp3D(mumuCandidate_fromFit->currentState().globalMomentum().x(),mumuCandidate_fromFit->currentState().globalMomentum().y(), mumuCandidate_fromFit->currentState().globalMomentum().z());


            float muonPos_ma_fit = MuPosCand_fromFit->currentState().mass();
            int   muonPos_ch_fit = MuPosCand_fromFit->currentState().particleCharge();
            float muonPos_px_fit = MuPosCand_fromFit->currentState().kinematicParameters().momentum().x();
            float muonPos_py_fit = MuPosCand_fromFit->currentState().kinematicParameters().momentum().y();
            float muonPos_pz_fit = MuPosCand_fromFit->currentState().kinematicParameters().momentum().z();
            float muonPos_en_fit = sqrt(muonPos_ma_fit*muonPos_ma_fit+muonPos_px_fit*muonPos_px_fit+muonPos_py_fit*muonPos_py_fit+muonPos_pz_fit*muonPos_pz_fit);

            reco::CompositeCandidate reco_ref_PM(muonPos_ch_fit,math::XYZTLorentzVector(muonPos_px_fit,muonPos_py_fit,muonPos_pz_fit,muonPos_en_fit),math::XYZPoint(dimuon_vx_fit,dimuon_vy_fit,dimuon_vz_fit),-13);
            pat::CompositeCandidate pat_ref_PM(reco_ref_PM);

            float muonNeg_ma_fit = MuNegCand_fromFit->currentState().mass();
            int   muonNeg_ch_fit = MuNegCand_fromFit->currentState().particleCharge();
            float muonNeg_px_fit = MuNegCand_fromFit->currentState().kinematicParameters().momentum().x();
            float muonNeg_py_fit = MuNegCand_fromFit->currentState().kinematicParameters().momentum().y();
            float muonNeg_pz_fit = MuNegCand_fromFit->currentState().kinematicParameters().momentum().z();
            float muonNeg_en_fit = sqrt(muonNeg_ma_fit*muonNeg_ma_fit+muonNeg_px_fit*muonNeg_px_fit+muonNeg_py_fit*muonNeg_py_fit+muonNeg_pz_fit*muonNeg_pz_fit);

            reco::CompositeCandidate reco_ref_NM(muonNeg_ch_fit,math::XYZTLorentzVector(muonNeg_px_fit,muonNeg_py_fit,muonNeg_pz_fit,muonNeg_en_fit),math::XYZPoint(dimuon_vx_fit,dimuon_vy_fit,dimuon_vz_fit),13);
            pat::CompositeCandidate pat_ref_NM(reco_ref_NM);


            mumuCandidate.addDaughter(*muonPos, "muonPos");
            mumuCandidate.addDaughter(*muonNeg,"muonNeg");

            pat_ref_JPsi.addUserFloat("deltaR",deltaRMuMu);
            pat_ref_JPsi.addUserFloat("mumuP4",mumuP4.M());

            pat_ref_JPsi.addDaughter(*pat_ref_PM, "ref_muonPos");
            pat_ref_JPsi.addDaughter(*pat_ref_NM, "ref_muonNeg");
            pat_ref_JPsi.addDaughter(*mumuCandidate, "mumuCandidate");
            // pat_ref_JPsi.addDaughter(*pat_ref_PM, "muonPos");
            // pat_ref_JPsi.addDaughter(*pat_ref_NM, "muonNeg");

            pat_ref_JPsi.addUserInt("nMatchedStationsPos",    nMatchedStationsPos);
            pat_ref_JPsi.addUserInt("nOverlapMusPos",         nOverlapMusPos);
            pat_ref_JPsi.addUserInt("nSharingSegWithPos",     nSharingSegWithPos);

            pat_ref_JPsi.addUserInt("nMatchedStationsNeg",    nMatchedStationsNeg);
            pat_ref_JPsi.addUserInt("nOverlapMusNeg",         nOverlapMusNeg);
            pat_ref_JPsi.addUserInt("nSharingSegWithNeg",     nSharingSegWithNeg);

            // ref_Jpsi.push_back(pat_ref_JPsi);
            // ref_mupos.push_back(pat_ref_PM);
            // ref_muneg.push_back(pat_ref_NM);
            //
            // Jpsi_p4.push_back(recoNegMuon->p4() + recoPosMuon->p4());
            // mupos_p4.push_back(recoPosMuon->p4());
            // muneg_p4.push_back(recoNegMuon->p4());

            pat_ref_JPsi.addUserFloat("muPos_Chi2", MuPosCand_fromFit->chiSquared());
            pat_ref_JPsi.addUserFloat("muPos_NDF",  MuPosCand_fromFit->degreesOfFreedom());
            pat_ref_JPsi.addUserFloat("muNeg_Chi2", MuNegCand_fromFit->chiSquared());
            pat_ref_JPsi.addUserFloat("muNeg_NDF",  MuNegCand_fromFit->degreesOfFreedom());

            pat_ref_JPsi.addUserFloat("VProb", mumuVProb);
            pat_ref_JPsi.addUserFloat("Chi2",  mumuChi2);

            AlgebraicVector3 MuMu_v3pperp ;
            MuMu_v3pperp[0] = MuMu_pperp.x(); MuMu_v3pperp[1] = MuMu_pperp.y(); MuMu_v3pperp[2] = 0.;

            TVector3 MuMu_pvtx, MuMu_pvtx3D, MuMu_vdiff, MuMu_vdiff3D ;
            float MuMu_cosAlpha, MuMu_cosAlpha3D, MuMu_ctau ;
            VertexDistanceXY MuMu_vdistXY ;
            Measurement1D MuMu_distXY ;
            GlobalError MuMu_v1e = (Vertex(*mumuCandidate_vertex_fromFit)).error(), MuMu_v2e;
            AlgebraicSymMatrix33 MuMu_vXYe ;
            float MuMu_ctauErr ;
            float MuMu_lxy, MuMu_lxyErr, MuMu_lxyz, MuMu_lxyzErr ;
            ROOT::Math::SVector<float, 3> MuMu_vDiff, MuMu_vDiff3D ; // needed by Similarity method

            ////////////////// Lifetime wrt PV for MuMu //////////////////
            MuMu_v2e = thePrimaryVtx.error();
            MuMu_vXYe = MuMu_v1e.matrix() + MuMu_v2e.matrix() ;

            /// 2D
            MuMu_pvtx.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), 0) ;
            MuMu_vdiff = MuMu_vtx - MuMu_pvtx ;
            MuMu_cosAlpha = MuMu_vdiff.Dot(MuMu_pperp) / (MuMu_vdiff.Perp()*MuMu_pperp.Perp()) ;
            MuMu_lxy = MuMu_vdiff.Perp();
            MuMu_vDiff[0] = MuMu_vdiff.x(); MuMu_vDiff[1] = MuMu_vdiff.y(); MuMu_vDiff[2] = 0 ; // needed by Similarity method
            MuMu_lxyErr = sqrt(ROOT::Math::Similarity(MuMu_vDiff,MuMu_vXYe)) / MuMu_vdiff.Perp();
            MuMu_distXY = MuMu_vdistXY.distance(Vertex(*mumuCandidate_vertex_fromFit), Vertex(thePrimaryVtx));
            MuMu_ctau = MuMu_distXY.value() * MuMu_cosAlpha * mumuCandidate_fromFit->currentState().mass() / MuMu_pperp.Perp();
            MuMu_ctauErr = sqrt(ROOT::Math::Similarity(MuMu_v3pperp,MuMu_vXYe)) * mumuCandidate_fromFit->currentState().mass() / (MuMu_pperp.Perp2()) ;

            /// 3D
            MuMu_pvtx3D.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), thePrimaryVtx.position().z());
            MuMu_vdiff3D = MuMu_vtx3D - MuMu_pvtx3D;
            MuMu_cosAlpha3D = MuMu_vdiff3D.Dot(MuMu_pperp3D)/(MuMu_vdiff3D.Mag()*MuMu_pperp3D.Mag());
            MuMu_lxyz = MuMu_vdiff3D.Mag();
            MuMu_vDiff3D[0] = MuMu_vdiff3D.x(); MuMu_vDiff3D[1] = MuMu_vdiff3D.y(); MuMu_vDiff3D[2] = MuMu_vdiff3D.z() ;
            MuMu_lxyzErr = sqrt(ROOT::Math::Similarity(MuMu_vDiff3D,MuMu_vXYe)) / MuMu_vdiff3D.Mag();

            pat_ref_JPsi.addUserFloat("cosAlpha",    MuMu_cosAlpha);
            pat_ref_JPsi.addUserFloat("cosAlpha3D",  MuMu_cosAlpha3D);
            pat_ref_JPsi.addUserFloat("ctau",        MuMu_ctau);
            pat_ref_JPsi.addUserFloat("ctauErr",     MuMu_ctauErr);
            pat_ref_JPsi.addUserFloat("lxy",         MuMu_lxy);
            pat_ref_JPsi.addUserFloat("lxyErr",      MuMu_lxyErr);
            pat_ref_JPsi.addUserFloat("lxyz",        MuMu_lxyz);
            pat_ref_JPsi.addUserFloat("lxyzErr",     MuMu_lxyzErr);

            ////////////////// Lifetime wrt BS for MuMu //////////////////
            MuMu_v2e = theBeamSpotVtx.error();
            MuMu_vXYe = MuMu_v1e.matrix() + MuMu_v2e.matrix();

            /// 2D
            MuMu_pvtx.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), 0);
            MuMu_vdiff = MuMu_vtx - MuMu_pvtx;
            MuMu_cosAlpha = MuMu_vdiff.Dot(MuMu_pperp)/(MuMu_vdiff.Perp()*MuMu_pperp.Perp());
            MuMu_lxy = MuMu_vdiff.Perp();
            MuMu_vDiff[0] = MuMu_vdiff.x(); MuMu_vDiff[1] = MuMu_vdiff.y(); MuMu_vDiff[2] = 0 ; // needed by Similarity method
            MuMu_lxyErr = sqrt(ROOT::Math::Similarity(MuMu_vDiff,MuMu_vXYe)) / MuMu_vdiff.Perp();
            MuMu_distXY = MuMu_vdistXY.distance(Vertex(*XCand_vertex_fromMCFit), Vertex(theBeamSpotVtx));
            MuMu_ctau = MuMu_distXY.value() * MuMu_cosAlpha * (XCand_fromMCFit->currentState().mass() / MuMu_pperp.Perp()) ;
            MuMu_ctauErr = sqrt(ROOT::Math::Similarity(MuMu_v3pperp,MuMu_vXYe)) * XCand_fromMCFit->currentState().mass()/MuMu_pperp.Perp2();

            /// 3D
            MuMu_pvtx3D.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), theBeamSpotVtx.position().z());
            MuMu_vdiff3D = MuMu_vtx3D - MuMu_pvtx3D;
            MuMu_cosAlpha3D = MuMu_vdiff3D.Dot(MuMu_pperp3D)/(MuMu_vdiff3D.Mag()*MuMu_pperp3D.Mag());
            MuMu_lxyz = MuMu_vdiff3D.Mag();
            MuMu_vDiff3D[0] = MuMu_vdiff3D.x(); MuMu_vDiff3D[1] = MuMu_vdiff3D.y(); MuMu_vDiff3D[2] = MuMu_vdiff3D.z() ;
            MuMu_lxyzErr = sqrt(ROOT::Math::Similarity(MuMu_vDiff3D,MuMu_vXYe)) / MuMu_vdiff3D.Mag();

            pat_ref_JPsi.addUserFloat("BS_cosAlpha",    MuMu_cosAlpha);
            pat_ref_JPsi.addUserFloat("BS_cosAlpha3D",  MuMu_cosAlpha3D);
            pat_ref_JPsi.addUserFloat("BS_ctau",        MuMu_ctau);
            pat_ref_JPsi.addUserFloat("BS_ctauErr",     MuMu_ctauErr);
            pat_ref_JPsi.addUserFloat("BS_lxy",         MuMu_lxy);
            pat_ref_JPsi.addUserFloat("BS_lxyErr",      MuMu_lxyErr);
            pat_ref_JPsi.addUserFloat("BS_lxyz",        MuMu_lxyz);
            pat_ref_JPsi.addUserFloat("BS_lxyzErr",     MuMu_lxyzErr);

            // MuMuDecayVtx_XE->push_back( sqrt( mumuCandidate_vertex_fromFit->error().cxx()) );
            // MuMuDecayVtx_YE->push_back( sqrt( mumuCandidate_vertex_fromFit->error().cyy()) );
            // MuMuDecayVtx_ZE->push_back( sqrt( mumuCandidate_vertex_fromFit->error().czz()) );


            Int_t dimuonType = 0;   //0 nothing,  1 J/psi  , 2 psi(2S)

            if ( mumuCandidate_fromFit->currentState().mass() > jspiMassCuts_[0]  &&  mumuCandidate_fromFit->currentState().mass() < jspiMassCuts_[1] )
            dimuonType = 1 ;
            else
            if ( mumuCandidate_fromFit->currentState().mass() > psiMassCuts_[0]  &&  mumuCandidate_fromFit->currentState().mass() < psiMassCuts_[1] )
            dimuonType = 2;

            pat_ref_JPsi.addUserInt("dimuonType",    dimuonType);

            if (Debug_) std::cout <<dimuonType <<std::endl;

            if (Debug_) std::cout <<"evt:" <<evtNum <<" MuMu with diMuonType = " <<dimuonType <<std::endl;
            //if (Debug_) std::cout << "POINT 0" << std::endl;
            // MuMuType->push_back(dimuonType);
            //if (Debug_) std::cout << "POINT  2" << std::endl;

            pat_ref_JPsi.addUserInt("isTriggerMatched",isTriggerMatched(posMuon,negMuon));

            int isTriggerMatched = 0;

            int ntriggers = HLTFileters_.size();
            if (Debug_) std::cout << "ntriggers: " << ntriggers << std::endl;

            for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++)
            {
              if (Debug_) std::cout << "MatchingTriggerResult[" << MatchTrig << "]: " << MatchingTriggerResult[MatchTrig] << std::endl;
              if ( MatchingTriggerResult[MatchTrig]!=0 )
              {
                if (Debug_) std::cout << "CHECKING FiltersForMatching_[" << MatchTrig << "]: " << FiltersForMatching_[MatchTrig] << std::endl;
                pat::TriggerObjectStandAloneCollection mu1HLTMatches = posMuon->triggerObjectMatchesByFilter( FiltersForMatching_[MatchTrig] );
                pat::TriggerObjectStandAloneCollection mu2HLTMatches = negMuon->triggerObjectMatchesByFilter( FiltersForMatching_[MatchTrig] );
                bool pass1 = mu1HLTMatches.size() > 0;
                bool pass2 = mu2HLTMatches.size() > 0;

                if ((pass1) && (pass2))
                {
                  isTriggerMatched = 1;
                  MuMuMuonTrigMatch->push_back(true);
                  if (Debug_) std::cout <<"Matched MuMu" <<std::endl ;
                } else
                //if (Debug_) std::cout << "POINT 9" << std::endl;
                MuMuMuonTrigMatch->push_back(false);
              }
              else
              //if (Debug_) std::cout << "POINT 10" << std::endl;
              MuMuMuonTrigMatch->push_back(false);
            }


            /// vertex without matched muons
            std::vector<TransientVertex> pvs ;
            Vertex mumuLessPV = thePrimaryVtx ;

            if (addMuMulessPrimaryVertex_)
            {
              VertexReProducer revertex(recVtxs, iEvent);
              Handle<TrackCollection> pvtracks;
              iEvent.getByLabel(revertex.inputTracks(),   pvtracks);
              Handle<BeamSpot>        pvbeamspot;
              iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);

              if ( pvbeamspot.isValid() < 0 )
              continue ;
              if (pvbeamspot.id() != beamSpotHandle.id()) {
                edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this Producer.";
              }
              const reco::Muon *rmu_1 = dynamic_cast<const reco::Muon*>( posMuon->originalObject() ) ;
              const reco::Muon *rmu_2 = dynamic_cast<const reco::Muon*>( negMuon->originalObject() ) ;

              if (rmu_1 != 0  &&  rmu_2 != 0  &&  rmu_1->track().id() == pvtracks.id()  &&  rmu_2->track().id() == pvtracks.id() ) {
                TrackCollection MuMuLess;
                MuMuLess.reserve(pvtracks->size());
                for (size_t i = 0, n = pvtracks->size(); i < n; ++i) {
                  if (i == rmu_1->track().key()) continue;
                  if (i == rmu_2->track().key()) continue;
                  MuMuLess.push_back((*pvtracks)[i]);
                }
                if (Debug_) std::cout <<"pvbeamspot.isValid() = " <<pvbeamspot.isValid() <<std::endl ;
                pvs = revertex.makeVertices(MuMuLess, *pvbeamspot, iSetup) ;
                if (!pvs.empty()) {
                  mumuLessPV = Vertex(pvs.front());
                }
              }
            }
            mumuLessPvs_n.push_back( pvs.size() );
            mumuLessPVs.push_back( mumuLessPV);

            muons.clear();


          } // 2nd loop over muons (look for mu-)
        } //first loop over muons (look for mu+)
      } // if (thePATMuonHandle->size() >= 2  && hasRequestedTrigger) {
      } // if (doMC_ || doData_)
      // AT THE END OF THE EVENT fill the tree and clear the vectors
      // ===========================================================


    }
    //}/// produce
    /// ------------ method called once each job just before starting event loop  ------------
    void MuMuProducerPAT::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
    {
    }
    void MuMuProducerPAT::beginJob()
    {


    }/// begin Job

    /// ------------ method called once each job just after ending the event loop  ------------
    void MuMuProducerPAT::endJob() {
      mumukktree->GetDirectory()->cd();
      mumukktree->Write();
    }/// endjob


    bool MuMuProducerPAT::isAbHadron(int pdgID) {

      if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
      return false;

    }

    bool MuMuProducerPAT::isAMixedbHadron(int pdgID, int momPdgID) {

      if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) ||
      (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0))
      return true;
      return false;

    }

    std::pair<int, float> MuMuProducerPAT::findCandMCInfo(reco::GenParticleRef genCand) {

      int momJpsiID = 0;
      float trueLife = -99.;
      //std::cout <<"externalmodule"<<std::endl;

      if (genCand->numberOfMothers()>0) {

        TVector3 trueVtx(0.0,0.0,0.0);
        TVector3 trueP(0.0,0.0,0.0);
        TVector3 trueVtxMom(0.0,0.0,0.0);

        trueVtx.SetXYZ(genCand->vertex().x(),genCand->vertex().y(),genCand->vertex().z());
        trueP.SetXYZ(genCand->momentum().x(),genCand->momentum().y(),genCand->momentum().z());

        bool aBhadron = false;
        reco::GenParticleRef Candmom = genCand->motherRef();       // find mothers
        if (Candmom.isNull()) {
          std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
          return result;
        } else {
          reco::GenParticleRef CandGrandMom = Candmom->motherRef();
          if (isAbHadron(Candmom->pdgId())) {
            if (CandGrandMom.isNonnull() && isAMixedbHadron(Candmom->pdgId(),CandGrandMom->pdgId())) {
              momJpsiID = CandGrandMom->pdgId();
              trueVtxMom.SetXYZ(CandGrandMom->vertex().x(),CandGrandMom->vertex().y(),CandGrandMom->vertex().z());
            } else {
              momJpsiID = Candmom->pdgId();
              trueVtxMom.SetXYZ(Candmom->vertex().x(),Candmom->vertex().y(),Candmom->vertex().z());
            }
            aBhadron = true;
          } else {
            if (CandGrandMom.isNonnull() && isAbHadron(CandGrandMom->pdgId())) {
              reco::GenParticleRef JpsiGrandgrandmom = CandGrandMom->motherRef();
              if (JpsiGrandgrandmom.isNonnull() && isAMixedbHadron(CandGrandMom->pdgId(),JpsiGrandgrandmom->pdgId())) {
                momJpsiID = JpsiGrandgrandmom->pdgId();
                trueVtxMom.SetXYZ(JpsiGrandgrandmom->vertex().x(),JpsiGrandgrandmom->vertex().y(),JpsiGrandgrandmom->vertex().z());
              } else {
                momJpsiID = CandGrandMom->pdgId();
                trueVtxMom.SetXYZ(CandGrandMom->vertex().x(),CandGrandMom->vertex().y(),CandGrandMom->vertex().z());
              }
              aBhadron = true;
            }
          }
          if (!aBhadron) {
            momJpsiID = Candmom->pdgId();
            trueVtxMom.SetXYZ(Candmom->vertex().x(),Candmom->vertex().y(),Candmom->vertex().z());
          }
        }

        TVector3 vdiff = trueVtx - trueVtxMom;
        trueLife = vdiff.Perp()*genCand->mass()/trueP.Perp();
      }
      std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
      return result;

    }

    float MuMuProducerPAT::getSigmaOfLogdEdx(float logde)
    {
      return 0.3;
    }

    float MuMuProducerPAT::getEnergyLoss(const reco::TrackRef & track)
    {
      if (iexception_dedx==1) return 9999.;
      const reco::DeDxDataValueMap & eloss = *energyLoss;
      return eloss[track].dEdx();
    }

    float MuMuProducerPAT::nsigmaofdedx(const reco::TrackRef & track, float & theo, float & sigma)
    {

      // no usable dE/dx if p > 2
      float nsigma = 99 ;
      if (iexception_dedx==1) return nsigma ;

      float m  = 0.13957;
      float bg = track->p() / m;

      theo = getLogdEdx(bg);


      int nhitr = track->numberOfValidHits();
      float meas = log(getEnergyLoss(track));
      sigma = getSigmaOfLogdEdx(theo) * pow(nhitr,-0.65);
      if (sigma>0)
      nsigma = (meas-theo) / sigma ;
      return nsigma;
    }


    float MuMuProducerPAT::getLogdEdx(float bg)
    {
      const float a =  3.25 ;
      const float b =  0.288;
      const float c = -0.852;

      float beta = bg/sqrt(bg*bg + 1);
      float dedx = log( a/(beta*beta) + b * log(bg) + c );

      return dedx;

    }


    float MuMuProducerPAT::GetMass(const reco::TrackRef & track){
      float P = track->p();
      float C = 2.625;
      float K = 2.495;
      float I = getEnergyLoss(track);
      return sqrt((I-C)/K)*P;
    }


    template<typename T>
    bool MuMuProducerPAT::isBetterMuon(const T &mu1, const T &mu2) const {
      if (mu2.track().isNull()) return true;
      if (mu1.track().isNull()) return false;
      if (mu1.isPFMuon() != mu2.isPFMuon()) return mu1.isPFMuon();
      if (mu1.isGlobalMuon() != mu2.isGlobalMuon()) return mu1.isGlobalMuon();
      if (mu1.charge() == mu2.charge() && deltaR2(mu1,mu2) < 0.0009) {
        return mu1.track()->ptError()/mu1.track()->pt() < mu2.track()->ptError()/mu2.track()->pt();
      } else {
        int nm1 = mu1.numberOfMatches(reco::Muon::SegmentArbitration);
        int nm2 = mu2.numberOfMatches(reco::Muon::SegmentArbitration);
        return (nm1 != nm2 ? nm1 > nm2 : mu1.pt() > mu2.pt());
      }
    }

    bool MuMuProducerPAT::isSameMuon(const reco::Muon &mu1, const reco::Muon &mu2) const {
      return (& mu1 == & mu2) ||
      //(mu1.originalObjectRef() == mu2.originalObjectRef()) ||
      (mu1.reco::Muon::innerTrack().isNonnull() ?
      mu1.reco::Muon::innerTrack() == mu2.reco::Muon::innerTrack() :
      mu1.reco::Muon::outerTrack() == mu2.reco::Muon::outerTrack());
    }


    /// define this as a plug-in
    DEFINE_FWK_MODULE(MuMuProducerPAT);

    // rsync -vut --existing src/MuMuProducerPAT.cc semrat@lxplus.cern.ch:/afs/cern.ch/user/s/semrat/scratch0/CMSSW_5_3_22/src/X4140/MuMuProducerPAT/src/MuMuProducerPAT.cc
