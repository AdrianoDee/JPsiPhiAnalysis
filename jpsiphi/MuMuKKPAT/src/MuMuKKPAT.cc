// -*- C++ -*-
//
// Package:    MuMuKKPAT
// Class:      MuMuKKPAT
//
/**\class MuMuKKPAT MuMuKKPAT.cc myAnalyzers/MuMuKKPAT/src/MuMuKKPAT.cc

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
#include "../interface/MuMuKKPAT.h"
#include "../interface/VertexReProducer.h"
//#include "DataFormats/Candidate/interface/OverlapChecker.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
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
MuMuKKPAT::MuMuKKPAT(const edm::ParameterSet& iConfig) :
hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults",edm::InputTag("TriggerResults::HLT")) ),
inputGEN_(iConfig.getUntrackedParameter<edm::InputTag>("inputGEN",edm::InputTag("genParticles"))),
vtxSample(iConfig.getUntrackedParameter<std::string>("VtxSample",std::string("offlinePrimaryVertices"))),

doData( iConfig.getUntrackedParameter<bool>("DoDataAnalysis", true) ),
doMC( iConfig.getUntrackedParameter<bool>("DoMonteCarloTree", true) ),
MCParticle( iConfig.getUntrackedParameter<int>("MonteCarloParticleId", 20443) ), /// 20443 X, 100443 Psi(2S), 9120443 X from B / decide later for X(4140)
MCExclusiveDecay( iConfig.getUntrackedParameter<bool>("MonteCarloExclusiveDecay", true) ),
MCMother( iConfig.getUntrackedParameter<int>("MonteCarloMotherId", 511) ), /// 511 B0 (=anti-B0), 531 B0 / decide later MCMotherId for X(4140)
MCDaughtersN( iConfig.getUntrackedParameter<int>(" MonteCarloDaughtersN", 3) ), /// will be same
doMuMuMassConst( iConfig.getUntrackedParameter<bool>("DoMuMuMassConstraint", true) ),
skipJPsi(iConfig.getUntrackedParameter<bool>("SkipJPsi", false)),

MuMinPixHits(iConfig.getUntrackedParameter<int>("MinNumMuPixHits", 0)),
MuMinSiHits(iConfig.getUntrackedParameter<int>("MinNumMuSiHits", 0)),
MuMaxNormChi(iConfig.getUntrackedParameter<double>("MaxMuNormChi2", 1000)),
MuMaxD0(iConfig.getUntrackedParameter<double>("MaxMuD0", 1000)),
sharedFraction(iConfig.getUntrackedParameter<double>("sharedFraction", 0.5)),

TrMinSiHits(iConfig.getUntrackedParameter<int>("MinNumTrSiHits", 0)),
TrMinPt(iConfig.getUntrackedParameter<double>("MinTrPt", 0)),
TrMaxNormChi2(iConfig.getUntrackedParameter<double>("MaxTrChi2NDF", 10)),
TriggersForMatching_(iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching")),
FiltersForMatching_(iConfig.getUntrackedParameter<std::vector<std::string> >("FiltersForMatching")),
resolveAmbiguity_(iConfig.getUntrackedParameter<bool>("resolvePileUpAmbiguity",true)),
addMuMulessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addMuMulessPrimaryVertex", true)),

JPsiMinMass(iConfig.getUntrackedParameter<double>("MinJPsiMass", 2.8)),
JPsiMaxMass(iConfig.getUntrackedParameter<double>("MaxJPsiMass", 3.4)),
PhiMinMass(iConfig.getUntrackedParameter<double>("MinPhiMass", 0.97)),
PhiMaxMass(iConfig.getUntrackedParameter<double>("MaxPhiMass", 1.07)),
XMassMin(iConfig.getUntrackedParameter<double>("XMassMin", 4.0)),
XMassMax(iConfig.getUntrackedParameter<double>("XMassMax", 6.0)),
MuMuTrackMaxDR(iConfig.getUntrackedParameter<double>("MaxMuMuTrackDR", 1)),

XTrackMaxDR(iConfig.getUntrackedParameter<double>("MaxXCandTrackDR", 1.1)),
UseXDR(iConfig.getUntrackedParameter<bool>("UseXDR", false)),
MuMuKKMinB0Mass(iConfig.getUntrackedParameter<double>("MinMuMuKKB0Mass", 0)),
MuMuKKMaxB0Mass(iConfig.getUntrackedParameter<double>("MaxMuMuKKB0Mass", 10)),
// MuMuKKMaxXMass(iConfig.getUntrackedParameter<double>("MaxMuMuKKXMass", 10)),
addXlessPrimaryVertex_(iConfig.getUntrackedParameter<bool>("addXlessPrimaryVertex", true)),

Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output",true)),
DeDxEstimator_(iConfig.getUntrackedParameter<std::string>("DeDxEstimator", std::string("dedxHarmonic2"))),
m_dEdxDiscrimTag(iConfig.getUntrackedParameter<std::string>("DeDxDiscriminator", std::string("dedxHarmonic2"))),
m_dEdxDiscrimTag_kaon(iConfig.getUntrackedParameter<std::string>("DeDxDiscriminator", std::string("dedxHarmonic2"))),

mumukktree(0),
runNum(0), evtNum(0), lumiNum(0),
trigRes(0), trigNames(0), L1TT(0), MatchTriggerNames(0),
/// counters for X
nMu(0), nMuMu(0), nX(0), nKK(0),
nX_pre0(0), nX_pre1(0), nX_pre2(0), nX_pre3(0), nX_pre4(0), nX_pre5(0), nX_pre6(0), nX_pre7(0), nX_pre8(0), nX_pre9(0), nX_pre10(0), nX_pre11(0), nX_pre12(0), nX_pre13(0), nX_pre14(0), nX_pre15(0),


n_pV(0),
tracks(0),
tracksPtSq_pV(0),
/// indices

mu1Idx(0), mu2Idx(0), MuMuType(0), ka1Idx(0), ka2Idx(0),
X_MuMuIdx(0), X_ka1Idx(0), X_ka2Idx(0),
/// MC Analysis /// n_B0Ancestors & no for X
n_genEvtVtx(0), genEvtVtx_X(0), genEvtVtx_Y(0), genEvtVtx_Z(0), genEvtVtx_particles(0), n_XAncestors(0),
nMCAll(0), nMCX(0), /*nMCXVtx(0),*/ MCPdgIdAll(0), MCDanNumAll(0),
// Gen Primary Vertex
PriVtxGen_X(0), PriVtxGen_Y(0), PriVtxGen_Z(0), PriVtxGen_EX(0), PriVtxGen_EY(0), PriVtxGen_EZ(0),
PriVtxGen_Chi2(0), PriVtxGen_CL(0), PriVtxGen_Ndof(0), PriVtxGen_tracks(0),
MCJPsiPx(0), MCJPsiPy(0), MCJPsiPz(0),
MCmupPx(0), MCmupPy(0), MCmupPz(0),
MCmumPx(0), MCmumPy(0), MCmumPz(0),
MCPhiPx(0), MCPhiPy(0), MCPhiPz(0),
MCkpPx(0), MCkpPy(0), MCkpPz(0),
MCkmPx(0), MCkmPy(0), MCkmPz(0),
MCpionPx(0), MCpionPy(0), MCpionPz(0),
MCkaonPx(0), MCkaonPy(0), MCkaonPz(0),
MCpionCh(0), MCkaonCh(0),
MCPx(0), MCPy(0), MCPz(0),
/// generic muons
muPx(0), muPy(0), muPz(0), muCharge(0),
muPhits(0), muShits(0), muLayersTr(0), muLayersPix(0),
muD0(0),  muD0E(0), muDz(0), muChi2(0), muNDF(0),
mufHits(0), muFirstBarrel(0), muFirstEndCap(0),
muDzVtx(0), muDxyVtx(0), muDzVtxErr(0), muKey(0),
muIsGlobal(0), muIsPF(0),
muGlMuHits(0), muGlChi2(0), muGlNDF(0), muGlMatchedStation(0),
muGlDzVtx(0), muGlDxyVtx(0),
nMatchedStations(0),
muType(0), muQual(0), muTrack(0), muNOverlap(0), muNSharingSegWith(0),
/// generic tracks
trNotRef(0), trRef(0),
trPx(0), trPy(0), trPz(0), trE(0),
trNDF(0), trPhits(0), trShits(0), trChi2(0),
trD0(0), trD0E(0), trCharge(0),
trfHits(0), trFirstBarrel(0), trFirstEndCap(0),
trDzVtx(0), trDxyVtx(0),
trQualityHighPurity(0), trQualityTight(0),
tr_nsigdedx(0), tr_dedx(0), tr_dedxMass(0), tr_theo(0), tr_sigma(0),
tr_dedx_byHits(0), tr_dedxErr_byHits(0), tr_saturMeas_byHits(0), tr_Meas_byHits(0),
/// MuMu cand & KK cand
MuMuMass(0), MuMuPx(0), MuMuPy(0), MuMuPz(0),
MuMuVtx_CL(0), MuMuVtx_Chi2(0),
MuMuDecayVtx_X(0), MuMuDecayVtx_Y(0), MuMuDecayVtx_Z(0),
MuMuDecayVtx_XE(0), MuMuDecayVtx_YE(0), MuMuDecayVtx_ZE(0),
MuMuMuonTrigMatch(0),
KKMass(0), KKPx(0), KKPy(0), KKPz(0),
KKVtx_CL(0), KKVtx_Chi2(0),
KKDecayVtx_X(0), KKDecayVtx_Y(0), KKDecayVtx_Z(0),
KKDecayVtx_XE(0), KKDecayVtx_YE(0), KKDecayVtx_ZE(0),
/// muons after JPsi (MuMu) fit &kaons after Phi (KK) fit
muPos_MuMu_Px(0), muPos_MuMu_Py(0), muPos_MuMu_Pz(0), muPos_MuMu_Chi2(0), muPos_MuMu_NDF(0),
muNeg_MuMu_Px(0), muNeg_MuMu_Py(0), muNeg_MuMu_Pz(0), muNeg_MuMu_Chi2(0), muNeg_MuMu_NDF(0),
kaonPos_KK_Px(0), kaonPos_KK_Py(0), kaonPos_KK_Pz(0), kaonPos_KK_Chi2(0), kaonPos_KK_NDF(0),
kaonNeg_KK_Px(0), kaonNeg_KK_Py(0), kaonNeg_KK_Pz(0), kaonNeg_KK_Chi2(0), kaonNeg_KK_NDF(0),
// DR_MuMu_K1(0), DR_MuMu_K2(0), DR_MuMuKK_K1(0), DR_MuMuKK_K2(0),
/// Primary Vertex with "MuMu correction"
mumuLessPvs_n(0),
PriVtxMuMuCorr_X(0), PriVtxMuMuCorr_Y(0), PriVtxMuMuCorr_Z(0), PriVtxMuMuCorr_EX(0), PriVtxMuMuCorr_EY(0), PriVtxMuMuCorr_EZ(0),
PriVtxMuMuCorr_Chi2(0), PriVtxMuMuCorr_CL(0), PriVtxMuMuCorr_tracks(0),
nTrk(0),
/// X candidates
xMass(0), xVtx_CL(0), xVtx_Chi2(0),
xPx(0), xPy(0), xPz(0), xPxE(0), xPyE(0), xPzE(0),
xDecayVtx_X(0), xDecayVtx_Y(0), xDecayVtx_Z(0), xDecayVtx_XE(0), xDecayVtx_YE(0), xDecayVtx_ZE(0),
/// Muons and tracks after X candidates fit
mu1Px_MuMuKK(0), mu1Py_MuMuKK(0), mu1Pz_MuMuKK(0), mu1E_MuMuKK(0),
mu2Px_MuMuKK(0), mu2Py_MuMuKK(0), mu2Pz_MuMuKK(0), mu2E_MuMuKK(0),
k1Px_MuMuKK(0), k1Py_MuMuKK(0), k1Pz_MuMuKK(0), k1E_MuMuKK(0),
kaonPos_nsigdedx(0), kaonPos_dedx(0), kaonPos_dedxMass(0), kaonPos_theo(0), kaonPos_sigma(0),
kaonPos_dedx_byHits(0), kaonPos_dedxErr_byHits(0), kaonPos_saturMeas_byHits(0), kaonPos_Meas_byHits(0),
k2Px_MuMuKK(0), k2Py_MuMuKK(0), k2Pz_MuMuKK(0), k2E_MuMuKK(0),
kaonNeg_nsigdedx(0), kaonNeg_dedx(0), kaonNeg_dedxMass(0), kaonNeg_theo(0), kaonNeg_sigma(0),
kaonNeg_dedx_byHits(0), kaonNeg_dedxErr_byHits(0), kaonNeg_saturMeas_byHits(0), kaonNeg_Meas_byHits(0),
/// Primary Vertex with largest B0_cos(alpha) no less values for X
PriVtx_XCosAlpha_n(0),
PriVtx_XCosAlpha_X(0), PriVtx_XCosAlpha_Y(0), PriVtx_XCosAlpha_Z(0), PriVtx_XCosAlpha_EX(0), PriVtx_XCosAlpha_EY(0), PriVtx_XCosAlpha_EZ(0),
PriVtx_XCosAlpha_Chi2(0), PriVtx_XCosAlpha_CL(0), PriVtx_XCosAlpha_tracks(0),
PriVtx_XCosAlpha3D_n(0),
PriVtx_XCosAlpha3D_X(0), PriVtx_XCosAlpha3D_Y(0), PriVtx_XCosAlpha3D_Z(0), PriVtx_XCosAlpha3D_EX(0), PriVtx_XCosAlpha3D_EY(0), PriVtx_XCosAlpha3D_EZ(0),
PriVtx_XCosAlpha3D_Chi2(0), PriVtx_XCosAlpha3D_CL(0), PriVtx_XCosAlpha3D_tracks(0),
XLessPV_tracksPtSq(0), XLessPV_4tracksPtSq(0),
PriVtxXLess_n(0),
PriVtxXLess_X(0), PriVtxXLess_Y(0), PriVtxXLess_Z(0), PriVtxXLess_EX(0), PriVtxXLess_EY(0), PriVtxXLess_EZ(0),
PriVtxXLess_Chi2(0), PriVtxXLess_CL(0), PriVtxXLess_tracks(0),
PriVtxXLess_XCosAlpha_n(0),
PriVtxXLess_XCosAlpha_X(0), PriVtxXLess_XCosAlpha_Y(0), PriVtxXLess_XCosAlpha_Z(0), PriVtxXLess_XCosAlpha_EX(0), PriVtxXLess_XCosAlpha_EY(0), PriVtxXLess_XCosAlpha_EZ(0),
PriVtxXLess_XCosAlpha_Chi2(0), PriVtxXLess_XCosAlpha_CL(0), PriVtxXLess_XCosAlpha_tracks(0),
PriVtxXLess_XCosAlpha3D_n(0),
PriVtxXLess_XCosAlpha3D_X(0), PriVtxXLess_XCosAlpha3D_Y(0), PriVtxXLess_XCosAlpha3D_Z(0), PriVtxXLess_XCosAlpha3D_EX(0), PriVtxXLess_XCosAlpha3D_EY(0), PriVtxXLess_XCosAlpha3D_EZ(0),
PriVtxXLess_XCosAlpha3D_Chi2(0), PriVtxXLess_XCosAlpha3D_CL(0), PriVtxXLess_XCosAlpha3D_tracks(0),
/// Primary Vertex with "B0 correction"
PriVtxXCorr_n(0),
PriVtxXCorr_X(0), PriVtxXCorr_Y(0), PriVtxXCorr_Z(0), PriVtxXCorr_EX(0), PriVtxXCorr_EY(0), PriVtxXCorr_EZ(0),
PriVtxXCorr_Chi2(0), PriVtxXCorr_CL(0), PriVtxXCorr_tracks(0),
/// Lifetime variables for B0
xCosAlphaBS(0), xCosAlpha3DBS(0), xCTauBS(0), xCTauBSE(0), xLxyBS(0), xLxyBSE(0), xLxyzBS(0), xLxyzBSE(0),
xCosAlphaPV(0), xCosAlpha3DPV(0), xCTauPV(0), xCTauPVE(0), xLxyPV(0), xLxyPVE(0), xLxyzPV(0), xLxyzPVE(0),
xCosAlphaPVCosAlpha(0), xCosAlpha3DPVCosAlpha(0), xCTauPVCosAlpha(0), xCTauPVCosAlphaE(0), xLxyPVCosAlpha(0), xLxyPVCosAlphaE(0), xLxyzPVCosAlpha(0), xLxyzPVCosAlphaE(0),
xCosAlphaPVCosAlpha3D(0), xCosAlpha3DPVCosAlpha3D(0), xCTauPVCosAlpha3D(0), xCTauPVCosAlpha3DE(0), xLxyPVCosAlpha3D(0), xLxyPVCosAlpha3DE(0), xLxyzPVCosAlpha3D(0), xLxyzPVCosAlpha3DE(0),
xCosAlphaXLessPV(0), xCosAlpha3DXLessPV(0),xCTauXLessPV(0),
xCTauXLessPVE(0), xLxyXLessPV(0), xLxyXLessPVE(0), xLxyzXLessPV(0), xLxyzXLessPVE(0),
xCosAlphaXLessPVCosAlpha(0), xCosAlpha3DXLessPVCosAlpha(0), xCTauXLessPVCosAlpha(0), xCTauXLessPVCosAlphaE(0), xLxyXLessPVCosAlpha(0), xLxyXLessPVCosAlphaE(0), xLxyzXLessPVCosAlpha(0), xLxyzXLessPVCosAlphaE(0),
xCosAlphaXLessPVCosAlpha3D(0), xCosAlpha3DXLessPVCosAlpha3D(0), xCTauXLessPVCosAlpha3D(0), xCTauXLessPVCosAlpha3DE(0), xLxyXLessPVCosAlpha3D(0), xLxyXLessPVCosAlpha3DE(0), xLxyzXLessPVCosAlpha3D(0), xLxyzXLessPVCosAlpha3DE(0),
xCosAlphaPVX(0), xCTauPVX(0), xCTauPVXE(0), xLxyPVX(0), xLxyzPVX(0), xLxyzPVXE(0),
xCTauPVX_3D(0), xCTauPVX_3D_err(0),
/// dxy, dz, dxyE, dzE for kaons from PV, BS, B0LessPV
kaonPos_dxy_PV(0), kaonPos_dz_PV(0), kaonNeg_dxy_PV(0), kaonNeg_dz_PV(0),
kaonPos_dxy_BS(0), kaonPos_dz_BS(0), kaonNeg_dxy_BS(0), kaonNeg_dz_BS(0),
kaonPos_dxy_XLessPV(0), kaonPos_dz_XLessPV(0), kaonNeg_dxy_XLessPV(0), kaonNeg_dz_XLessPV(0),
kaonPos_dxyE(0), kaonPos_dzE(0), kaonNeg_dxyE(0), kaonNeg_dzE(0),

kaonPosFromPV(0), kaonNegFromPV(0)

{
  /// now do what ever initialization is needed
  MuMuMinMass = JPsiMinMass;
  MuMuMaxMass = JPsiMaxMass;
  KKMinMass = PhiMinMass;
  KKMaxMass = PhiMaxMass;
  MinXMass = XMassMin;
  MaxXMass = XMassMax;

}

MuMuKKPAT::~MuMuKKPAT()
{
  /// do anything here that needs to be done at desctruction time
  /// (e.g. close files, deallocate resources etc.)

}


///
/// member functions
///

/// ------------ method called to for each event  ------------
void MuMuKKPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /// get event content information
  int debug = 0;

  if(Debug_) std::cout <<"Debug : " << ++debug << std::endl;

  bool decayChainOK = false;
  runNum = iEvent.id().run();
  evtNum = iEvent.id().event();
  lumiNum = iEvent.id().luminosityBlock();


  bool hasRequestedTrigger = false;
  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  /// first get HLT results
  std::map<string,int> HLTPreScaleMap;
  edm::Handle<edm::TriggerResults> hltresults;
  try {
    iEvent.getByLabel(hlTriggerResults_, hltresults);
  }
  catch ( ... ) {
    std::cout << "Couldn't get handle on HLT Trigger!" << std::endl;
  }
  if (!hltresults.isValid()) {
    std::cout << "No Trigger Results!" << std::endl;
  }
  else {
    int ntrigs = hltresults->size();
    if (ntrigs==0){
      std::cout << "No trigger name given in TriggerResults of the input " << std::endl;
    }

    /// get hold of trigger names - based on TriggerResults object!
    edm::TriggerNames triggerNames_;
    triggerNames_ = iEvent.triggerNames(*hltresults);
    int ntriggers = TriggersForMatching_.size();
    for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++) { // initialize MatchingTriggerResult array
      MatchingTriggerResult[MatchTrig] = 0;
    }

    for (int itrig = 0; itrig < ntrigs; itrig++) {
      string trigName = triggerNames_.triggerName(itrig);
      int hltflag = (*hltresults)[itrig].accept();
      if (Debug_) if (hltflag) std::cout << trigName << " " <<hltflag <<std::endl;
      trigRes->push_back(hltflag);
      trigNames->push_back(trigName);

      int ntriggers = TriggersForMatching_.size();
      for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++) {
        if (TriggersForMatching_[MatchTrig] == triggerNames_.triggerName(itrig)){
          MatchingTriggerResult[MatchTrig] = hltflag;
          if (hltflag==1) hasRequestedTrigger = true;
          break;
        }
      }
    }
    for (int MatchTrig = 0; MatchTrig<ntriggers; MatchTrig++){
      if (Debug_) std::cout << TriggersForMatching_[MatchTrig]<<std::endl;
      MatchTriggerNames->push_back(TriggersForMatching_[MatchTrig]);
    }

    ///
    /// Get HLT map : triggername associated with its prescale, saved only for accepted trigger
    ///
    for (unsigned int i=0; i<triggerNames_.size(); i++){
      if ( hltresults->accept(i) ) { //  save trigger info only for accepted paths
        /// get the prescale from the HLTConfiguration, initialized at beginRun
        int prescale = hltConfig_.prescaleValue(iEvent,iSetup,triggerNames_.triggerNames().at(i));
        if (Debug_) std::cout<<" HLT===> "<<triggerNames_.triggerNames().at(i)<<" prescale ="<<prescale<<std::endl;
        HLTPreScaleMap[triggerNames_.triggerNames().at(i)] = prescale;
      }
    }
    HLTTrig = &HLTPreScaleMap; // store in the branch

  } /// end valid trigger

  std::cout <<"Debug : " << ++debug << std::endl;

  /// get L1 trigger info
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
  iEvent.getByLabel( edm::InputTag("gtDigis"), gtRecord);
  const DecisionWord dWord = gtRecord->decisionWord();
  const TechnicalTriggerWord ttWord = gtRecord->technicalTriggerWord();
  for(unsigned int l1i = 0; l1i != ttWord.size(); ++l1i){
    L1TT->push_back(ttWord.at(l1i));
  }

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
  iEvent.getByLabel(vtxSample, recVtxs);
  unsigned int nVtxTrks = 0;
  if ( recVtxs->begin() != recVtxs->end() ) {
    thePrimaryVtx_multiplicity = recVtxs->size() ;

    if (addMuMulessPrimaryVertex_ || addXlessPrimaryVertex_ || resolveAmbiguity_) {
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
  pV = thePrimaryVtx;
  n_pV = thePrimaryVtx_multiplicity ;
  // priVtx_X = (thePrimaryVtx.position().x()) ;
  // priVtx_Y = (thePrimaryVtx.position().y()) ;
  // priVtx_Z = (thePrimaryVtx.position().z()) ;
  // priVtx_XE = (thePrimaryVtx.xError()) ;
  // priVtx_YE = (thePrimaryVtx.yError()) ;
  // priVtx_ZE = (thePrimaryVtx.zError()) ;
  // priVtx_NormChi2 = (thePrimaryVtx.normalizedChi2()) ;
  // priVtx_Chi2 = thePrimaryVtx.chi2() ;
  // priVtx_CL = ChiSquaredProbability( (double)(thePrimaryVtx.chi2()), (double)(thePrimaryVtx.ndof())) ;
  // priVtx_tracks = thePrimaryVtx.tracksSize() ;
  VertexHigherPtSquared vertexHigherPtSquared ;
  tracksPtSq_pV = vertexHigherPtSquared.sumPtSquared(thePrimaryVtx) ;

  /// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// try reconstruction without fitting modules
  /// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  Handle< std::vector<pat::GenericParticle> > thePATTrackHandle;
  iEvent.getByLabel("cleanPatTrackCands", thePATTrackHandle); /// container of tracks with pion mass hypothesis
  Handle< std::vector<pat::GenericParticle> > theKaonRefittedPATTrackHandle;
  iEvent.getByLabel("cleanPatTrackKaonCands", theKaonRefittedPATTrackHandle); /// container of tracks with kaon mass hypothesis

  // for ( std::vector<pat::GenericParticle>::const_iterator TrackNotRefitted = thePATTrackHandle->begin(); TrackNotRefitted != thePATTrackHandle->end(); ++TrackNotRefitted ) {
  //   for ( std::vector<pat::GenericParticle>::const_iterator TrackRefitted = theKaonRefittedPATTrackHandle->begin(); TrackRefitted != theKaonRefittedPATTrackHandle->end(); ++TrackRefitted ) {
  //     if ( TrackNotRefitted->track().key() == TrackRefitted->track().key() ) {
  //       trNotRef->push_back( TrackNotRefitted->p() ) ;
  //       trRef->push_back( TrackRefitted->p() ) ;
  //       break ;
  //     }
  //   }
  //   break ;
  // }

  Handle< std::vector<pat::Muon> > thePATMuonHandle;
  iEvent.getByLabel("patMuonsWithTrigger", thePATMuonHandle);

  Handle<reco::DeDxDataValueMap> elossCollection;
  energyLoss = 0;
  iexception_dedx = 0;
  try {
    iEvent.getByLabel(DeDxEstimator_, elossCollection);
    energyLoss = elossCollection.product();
  } catch ( cms::Exception& ex ) {
    if (evtNum < 100) edm::LogError("Analyzer") <<"Warning can't get collection with label: elossCollection";
    iexception_dedx = 1;
  }

  std::cout <<"Debug : " << ++debug << std::endl;

  /// dE/dx hits
  Handle<edm::ValueMap<reco::DeDxData> > dEdxTrackHandle;
  try {
    iEvent.getByLabel(m_dEdxDiscrimTag, dEdxTrackHandle);
    dEdxTrack = *dEdxTrackHandle.product();
  } catch ( cms::Exception& ex ) {
    if (evtNum < 100) edm::LogError("Analyzer") <<"Warning can't get collection with label: dEdxTrackHandle";
    iexception_dedx = 1;
  }

  Handle<edm::ValueMap<reco::DeDxData> > dEdxTrackHandle_Kaon;
  try {
    iEvent.getByLabel(m_dEdxDiscrimTag_kaon, dEdxTrackHandle_Kaon);
    dEdxTrack_Kaon = *dEdxTrackHandle_Kaon.product();
  } catch ( cms::Exception& ex ) {
    if (evtNum < 100) edm::LogError("Analyzer") <<"Warning can't get collection with label: dEdxTrackHandle_Kaon";
    iexception_dedx = 1;
  }


  std::cout <<"Debug : " << ++debug << std::endl;

  ////////////////// check MC truth //////////////////
  if (doMC) {
    /*
    // Get generated event
    //Handle<edm::HepMCProduct> hepEv;
    //iEvent.getByLabel("generator", hepEv);
    Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByLabel("generator", genEvtInfo);

    //const HepMC::GenEvent *myGenEvent = hepEv->GetEvent();
    const HepMC::GenEvent *myGenEvent = genEvtInfo->GetEvent();
    n_genEvtVtx = myGenEvent->vertices_size() ;

    HepMC::GenVertex* primaryGenVtx = *(myGenEvent->vertices_begin()) ;

    genEvtVtx_X->push_back( primaryGenVtx->point3d().x() );
    genEvtVtx_Y->push_back( primaryGenVtx->point3d().y() );
    genEvtVtx_Z->push_back( primaryGenVtx->point3d().z() );
    //genEvtVtx_XE = (primaryGenVtx->xError()) ;
    //genEvtVtx_YE = (primaryGenVtx->yError()) ;
    //genEvtVtx_ZE = (primaryGenVtx->zError()) ;
    //genEvtVtx_NormChi2 = (primaryGenVtx->normalizedChi2()) ;
    //genEvtVtx_Chi2 = primaryGenVtx->chi2() ;
    //genEvtVtx_CL = ChiSquaredProbability( (double)(primaryGenVtx.chi2()), (double)(primaryGenVtx.ndof())) ;
    genEvtVtx_particles->push_back( primaryGenVtx->particles_out_size() );
    */

    Handle< std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    if (Debug_) std::cout <<"\nBunchXing multiplicity = " <<PupInfo->size() <<std::endl ;
    for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
    if (Debug_) std::cout <<"Pileup Information: bunchXing, nvtx: " <<PVI->getBunchCrossing() <<" " <<PVI->getPU_NumInteractions() <<std::endl;

    Handle<GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
    if (Debug_) std::cout << "############### GenParticles Analysis ###############" << std::endl;
    float jpsiPx=0., jpsiPy=0., jpsiPz=0.;
    float  mupPx=0., mupPy=0., mupPz=0., mumPx=0., mumPy=0., mumPz=0.;
    float phiPx=0., phiPy=0., phiPz=0.;
    float  kpPx=0., kpPy=0., kpPz=0., kmPx=0., kmPy=0., kmPz=0.;
    //float pionPx=0., pionPy=0., pionPz=0., kaonPx=0., kaonPy=0., kaonPz=0.;
    //int pionCh=0, kaonCh=0 ;

    for (size_t i = 0; i < genParticles->size(); ++ i) {
      nMCAll++;
      const reco::GenParticle &p = (*genParticles)[i];
      int pdgid = p.pdgId() ;
      int dauNum = p.numberOfDaughters();
      MCPdgIdAll->push_back( pdgid );
      MCDanNumAll->push_back( dauNum );

      if ( MCExclusiveDecay ) {
        /// check if there is a MCMother which has MCDaughtersN daughters
        if ( abs(pdgid) == MCMother  &&  dauNum == MCDaughtersN ) {
          bool mumuOK = false;
          bool kkOK = false;
          //bool pionOK = false, kaonOK = false;

          for (int j=0; j<dauNum; ++j) {
            const Candidate *dau = p.daughter(j);
            if (Debug_) std::cout << "dauPdgId = " << dau->pdgId() << std::endl;

            /// check if one of B0 daughters is a psi(nS) whitch has 2 muons as daughters /// SEMRA ask again !!!
            int mumuId = 0 ;
            if (skipJPsi) /// SEMRA cleaned skipPsi2S
            if (Debug_) std::cout <<"Skipping J/psi!" <<std::endl ; /// SEMRA cleaned skipPsi2S
            //else if (skipPsi2S) /// SEMRA
            //  mumuId = 443 ; /// SEMRA (JPsi ID)

            if ( ((skipJPsi) && (dau->pdgId() == mumuId)) ||
            ((!skipJPsi) && (dau->pdgId()%1000 == 443)) ) {
              jpsiPx = dau->px(); jpsiPy = dau->py(); jpsiPz = dau->pz();
              int jpsiDauNum = dau->numberOfDaughters();
              if (Debug_) std::cout << "jpsiDauNum = " << jpsiDauNum << std::endl;
              int muNum = 0;
              for (int k=0; k<jpsiDauNum; ++k) {
                const Candidate *grandDau = dau->daughter(k);
                if (Debug_)  std::cout << "grandDauPdgId = " << grandDau->pdgId() << std::endl;
                if ( abs(grandDau->pdgId()) == 13 ) {
                  muNum++;
                  if (grandDau->pdgId() < 0) {
                    mupPx = grandDau->px(); mupPy = grandDau->py(); mupPz = grandDau->pz();
                  } else {
                    mumPx = grandDau->px(); mumPy = grandDau->py(); mumPz = grandDau->pz();
                  }
                }
              }
              if ( muNum == 2 ) mumuOK = true ;

            } /// end check if one of the MCMother daughters is a J/Psi or psi'

            /// for Phi
            phiPx = dau->px(); phiPy = dau->py(); phiPz = dau->pz();
            int phiDauNum = dau->numberOfDaughters();
            if (Debug_) std::cout << "phiDauNum = " << phiDauNum << std::endl;
            int kNum = 0;
            for (int n=0; n<phiDauNum; ++n) {
              const Candidate *grandDau = dau->daughter(n);
              if (Debug_)  std::cout << "grandDauPdgId = " << grandDau->pdgId() << std::endl;
              if ( abs(grandDau->pdgId()) == 321 ) {
                kNum++;
                if (grandDau->pdgId() < 0) {
                  kpPx = grandDau->px(); kpPy = grandDau->py(); kpPz = grandDau->pz();
                } else {
                  kmPx = grandDau->px(); kmPy = grandDau->py(); kmPz = grandDau->pz();
                }
              }
            }
            if ( kNum == 2 ) kkOK = true ;


            /*else if ( abs(dau->pdgId()) == 211 ) { // check if one of B0 daughters is a pion /// SEMRA ask again !!!
            pionPx = dau->px(); pionPy = dau->py(); pionPz = dau->pz();
            pionCh = (dau->pdgId() == 211)? 1 : -1;
            pionOK = true; /// SEMRA pions change with kaons for B0 ?
          } else if ( abs(dau->pdgId()) == 321 ) { // check if one of B0 daughters is a kaon /// SEMRA ask again !!!
          kaonPx = dau->px(); kaonPy=dau->py(); kaonPz=dau->pz();
          kaonCh = (dau->pdgId() == 321)? 1 : -1;
          kaonOK = true;
        }*/

      } /// end loop on MCMother daughters

      if (Debug_) std::cout << "mumuOK = " << mumuOK << ", kkOK = " << kkOK << std::endl;
      if ( mumuOK && kkOK ) {
        if (Debug_) {
          std::cout <<"\nnumber of X mothers = " <<p.numberOfMothers() <<std::endl ;
          std::cout <<"X mother pdgID = " <<p.mother(0)->pdgId() <<std::endl ;
        }
        ++nMCX ;
        PriVtxGen_X->push_back( p.vx() ) ;
        PriVtxGen_Y->push_back( p.vy() ) ;
        PriVtxGen_Z->push_back( p.vz() ) ;
        PriVtxGen_CL->push_back( p.vertexNormalizedChi2() ) ;
        PriVtxGen_Chi2->push_back( p.vertexChi2() ) ;
        PriVtxGen_Ndof->push_back( p.vertexNdof() ) ;

        Bool_t status = kTRUE ;
        const Candidate *x_ancestor = p.mother(0) ; /// a particle can have several mothers
        Int_t n_ancestors = 1 ;
        while ( status ) {
          if ( abs(x_ancestor->pdgId()) <= 8 || x_ancestor->pdgId() == 21 || x_ancestor->status() == 3 ) {
            status = kFALSE ;
            if (Debug_) std::cout <<"X ancestor ID = " <<x_ancestor->pdgId() <<std::endl ;
            genEvtVtx_X->push_back( x_ancestor->daughter(0)->vx() ) ;
            genEvtVtx_Y->push_back( x_ancestor->daughter(0)->vy() ) ;
            genEvtVtx_Z->push_back( x_ancestor->daughter(0)->vz() ) ;
            genEvtVtx_particles->push_back( x_ancestor->numberOfDaughters() ) ;
            n_XAncestors->push_back( n_ancestors ) ;
          }
          else {
            x_ancestor = x_ancestor->mother(0) ;
            n_ancestors++ ;
          }
        }

        MCJPsiPx->push_back(jpsiPx); MCJPsiPy->push_back(jpsiPy); MCJPsiPz->push_back(jpsiPz);
        MCmupPx->push_back(mupPx); MCmupPy->push_back(mupPy); MCmupPz->push_back(mupPz);
        MCmumPx->push_back(mumPx); MCmumPy->push_back(mumPy); MCmumPz->push_back(mumPz);
        MCPhiPx->push_back(phiPx); MCPhiPy->push_back(phiPy); MCPhiPz->push_back(phiPz);
        MCkpPx->push_back(kpPx); MCkpPy->push_back(kpPy); MCkpPz->push_back(kpPz);
        MCkmPx->push_back(kmPx); MCkmPy->push_back(kmPy); MCkmPz->push_back(kmPz);
        //MCpionPx->push_back(pionPx); MCpionPy->push_back(pionPy); MCpionPz->push_back(pionPz);
        //MCkaonPx->push_back(kaonPx); MCkaonPy->push_back(kaonPy); MCkaonPz->push_back(kaonPz);
        //MCpionCh->push_back(pionCh) ; MCkaonCh->push_back(kaonCh) ;
        decayChainOK = true;
        MCPx->push_back( p.px() );
        MCPy->push_back( p.py() );
        MCPz->push_back( p.pz() );
      }
      if (Debug_) std::cout << "decayChainOK = " << decayChainOK << std::endl;
    } // if ( abs(pdgid) == MCMother  &&  dauNum == 3 )
  } // if ( !MCExclusiveDecay )

} // for (size_t i = 0; i < genParticles->size(); ++ i)
} // if (doMC)


/// reconstruction only for events with B decaying in psi(nS)+Pi+K /// SEMRA JPsiPhi !!!
if ( (doMC && !MCExclusiveDecay) || (doMC && (MCExclusiveDecay && decayChainOK)) || doData ) {

  bool isEventWithInvalidMu = false;

  if (Debug_) std::cout << "Starting event with " << thePATTrackHandle->size() << " tracks, and " << thePATMuonHandle->size() << " muons" << std::endl;

  if ((thePATMuonHandle->size()) * (thePATTrackHandle->size()) > 20000) {
    std::cout << "Too many Muons: " << thePATMuonHandle->size() << ", and Tracks: " << thePATTrackHandle->size() << std::endl;
  } else //if (thePATMuonHandle->size() >= 2) { // check
    if (thePATMuonHandle->size() >= 2  && hasRequestedTrigger) {

      if (Debug_) std::cout <<"============================  Evt: " <<evtNum <<" accept event with 2 mu and trigger ==============================================" <<std::endl;

      bool CollectTracks_ = false;

      ////////////////// filling track tree //////////////////
      if(CollectTracks_)
        {
        for ( std::vector<pat::GenericParticle>::const_iterator iTr = thePATTrackHandle->begin(); iTr != thePATTrackHandle->end(); ++iTr ) {
          pat::GenericParticle tr = *iTr;
          tracks.push_back(*tr.track());

          // trPx->push_back(tr.px());
          // trPy->push_back(tr.py());
          // trPz->push_back(tr.pz());
          // trE->push_back(tr.energy());
          // trPhits->push_back(tr.track()->hitPattern().numberOfValidPixelHits());
          // trShits->push_back(tr.track()->hitPattern().numberOfValidStripHits());
          // trChi2->push_back(tr.track()->chi2());
          // trNDF->push_back(tr.track()->ndof());
          // trD0->push_back(tr.track()->d0());
          // trD0E->push_back(tr.track()->d0Error());
          // trCharge->push_back(tr.charge());
          // float hits = (1.0*tr.track()->found() )/ (tr.track()->found()+ tr.track()->lost() + tr.track()->trackerExpectedHitsInner().numberOfHits() + tr.track()->trackerExpectedHitsOuter().numberOfHits());
          // trfHits->push_back(hits);
          // trFirstBarrel->push_back(tr.track()->hitPattern().hasValidHitInFirstPixelBarrel());
          // trFirstEndCap->push_back(tr.track()->hitPattern().hasValidHitInFirstPixelEndcap());
          trDzVtx->push_back(tr.track()->dz(RefVtx));
          trDxyVtx->push_back(tr.track()->dxy(RefVtx));
          double theo = 0., sigma = 0. ;
          tr_nsigdedx->push_back(nsigmaofdedx(tr.track(),theo,sigma));
          tr_dedx->push_back(getEnergyLoss(tr.track()));
          tr_dedxMass->push_back(GetMass(tr.track()));
          tr_theo->push_back(theo);
          tr_sigma->push_back(sigma);
          tr_dedx_byHits->push_back( (dEdxTrack)[tr.track()].dEdx() );
          tr_dedxErr_byHits->push_back( (dEdxTrack)[tr.track()].dEdxError() );
          tr_saturMeas_byHits->push_back( (dEdxTrack)[tr.track()].numberOfSaturatedMeasurements() );
          tr_Meas_byHits->push_back( (dEdxTrack)[tr.track()].numberOfMeasurements() );
          /// Track quality:
          /// loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4, looseSetWithPV=5, highPuritySetWithPV=6
          // bool ishighPurity = tr.track()->quality(reco::TrackBase::highPurity);
          // trQualityHighPurity->push_back(ishighPurity);
          // trQualityTight->push_back(tr.track()->quality(reco::TrackBase::tight));
        }
     }

      /// get MuMu cands
      for ( std::vector<pat::Muon>::const_iterator posMuon = thePATMuonHandle->begin(); posMuon != thePATMuonHandle->end(); ++posMuon ) {

        /// push back all muon information
        ++nMu;
        const reco::Muon* recoPosMuon = dynamic_cast<const reco::Muon * >(posMuon->originalObject());
        // muPx->push_back(recoPosMuon->px());
        // muPy->push_back(recoPosMuon->py());
        // muPz->push_back(recoPosMuon->pz());
        // muCharge->push_back(recoPosMuon->charge());

        if (recoPosMuon->track().isNull()) continue;
        if (recoPosMuon->charge()<=0.0) continue;

        // if (recoPosMuon->track().isNull()) { // rmu->track() returns innerTrack();
        //   std::cout << "no track for " << std::distance(thePATMuonHandle->begin(), posMuon) << " filling defaults" << std::endl;
        //   /// AF
        //   muD0->push_back(0);
        //   muDz->push_back(0);
        //   muChi2->push_back(0);
        //   muNDF->push_back(-1);
        //   muPhits->push_back(0);
        //   muShits->push_back(0);
        //   muLayersTr->push_back(0);
        //   muLayersPix->push_back(0);
        //   muDzVtx->push_back(0);
        //   muDxyVtx->push_back(0);
        //   mufHits->push_back(0);
        //   muFirstBarrel->push_back(0);
        //   muFirstEndCap->push_back(0);
        //   muD0E->push_back(0);
        //   muDzVtxErr->push_back(0);
        //   muKey->push_back(0);
        //   muGlChi2->push_back(0);
        //   muGlNDF->push_back(-1);
        //   muGlMuHits->push_back(0);
        //   muGlMatchedStation->push_back(0);
        //   muGlDzVtx->push_back(0);
        //   muGlDxyVtx->push_back(0);
        //   nMatchedStations->push_back(0) ;
        //
        //   if (Debug_) std::cout <<"evt:" <<evtNum << "no track for PAT muon " <<std::distance(thePATMuonHandle->begin(), posMuon) <<" skipping muon... should skip event instead" <<std::endl;
        //   isEventWithInvalidMu = true;
        //   continue;
        // }
        // else {
        //   muD0->push_back(recoPosMuon->track()->d0());
        //   muDz->push_back(recoPosMuon->track()->dz());
        //   muChi2->push_back(recoPosMuon->track()->chi2());
        //   muNDF->push_back(recoPosMuon->track()->ndof());
        //   muPhits->push_back(recoPosMuon->track()->hitPattern().numberOfValidPixelHits());
        //   muShits->push_back(recoPosMuon->track()->hitPattern().numberOfValidStripHits());
        //   if (Debug_) std::cout <<"evt:" <<evtNum <<" trackerLayersWithMeasurement=" <<recoPosMuon->track()->hitPattern().trackerLayersWithMeasurement() <<std::endl;
        //   if ( !(recoPosMuon->track()->hitPattern().trackerLayersWithMeasurement()) ) {
        //     isEventWithInvalidMu = true;
        //     if (Debug_) std::cout <<"evt:" <<evtNum <<" problem with trackerLayersWithMeasurement" <<std::endl;
        //     continue ;
        //   }
        //   if ( !(recoPosMuon->track()->hitPattern().pixelLayersWithMeasurement()) ) {
        //     isEventWithInvalidMu = true;
        //     continue ;
        //   }
        //   muLayersTr->push_back(recoPosMuon->track()->hitPattern().trackerLayersWithMeasurement());
        //   muLayersPix->push_back(recoPosMuon->track()->hitPattern().pixelLayersWithMeasurement());
        //   muDzVtx->push_back(recoPosMuon->track()->dz(RefVtx));
        //   muDxyVtx->push_back(recoPosMuon->track()->dxy(RefVtx));
        //   mufHits->push_back((1.0*recoPosMuon->track()->found())/ (recoPosMuon->track()->found()+ recoPosMuon->track()->lost() + recoPosMuon->track()->trackerExpectedHitsInner().numberOfHits() + recoPosMuon->track()->trackerExpectedHitsOuter().numberOfHits() ) );
        //   if (Debug_) std::cout <<"mu found " <<recoPosMuon->track()->found() <<" fHits=" <<(1.0*recoPosMuon->track()->found())/ (recoPosMuon->track()->found()+ recoPosMuon->track()->lost() + recoPosMuon->track()->trackerExpectedHitsInner().numberOfHits() + recoPosMuon->track()->trackerExpectedHitsOuter().numberOfHits() ) <<std::endl;
        //   muFirstBarrel->push_back(recoPosMuon->track()->hitPattern().hasValidHitInFirstPixelBarrel());
        //   muFirstEndCap->push_back(recoPosMuon->track()->hitPattern().hasValidHitInFirstPixelEndcap());
        //   muD0E->push_back(recoPosMuon->track()->d0Error());
        //   muDzVtxErr->push_back(recoPosMuon->track()->dzError());
        //   muKey->push_back(recoPosMuon->track().key());
        // }
        //
        //
        // muIsGlobal->push_back( recoPosMuon->isGlobalMuon() ) ;
        // muIsPF->push_back( recoPosMuon->isPFMuon() ) ;
        // if ( recoPosMuon->globalTrack().isNull() ) {
        //   muGlMuHits->push_back(0);
        //   muGlChi2->push_back(0);
        //   muGlNDF->push_back(-1);
        //   muGlMatchedStation->push_back(0);
        //   muGlDzVtx->push_back(-1);
        //   muGlDxyVtx->push_back(-1);
        // }
        // else {
        //   muGlMuHits->push_back(recoPosMuon->globalTrack()->hitPattern().numberOfValidMuonHits());
        //   muGlChi2->push_back(recoPosMuon->globalTrack()->chi2());
        //   muGlNDF->push_back(recoPosMuon->globalTrack()->ndof());
        //   muGlMatchedStation->push_back(recoPosMuon->numberOfMatchedStations());
        //   muGlDzVtx->push_back(recoPosMuon->globalTrack()->dz(RefVtx));
        //   muGlDxyVtx->push_back(recoPosMuon->globalTrack()->dxy(RefVtx));
        // }
        // nMatchedStations->push_back(recoPosMuon->numberOfMatchedStations()) ;
        // muType->push_back(recoPosMuon->type());
        // int qm = 0;
        // for (int qi=1; qi!= 24; ++qi) {
        //   if (muon::isGoodMuon(*recoPosMuon, muon::SelectionType(qi)))
        //   qm += 1<<qi;
        // }
        // muQual->push_back(qm);
        // muTrack->push_back(-1);// not implemented yet
        //
        // ////////////////// muon cleaning //////////////////
        // int nOverlapMus = 0, nSharingSegWith = -1;
        // int nSegments1 = recoPosMuon->numberOfMatches(reco::Muon::SegmentArbitration);
        // for ( std::vector<pat::Muon>::const_iterator negMuon = posMuon+1; negMuon != thePATMuonHandle->end(); ++negMuon) {
        //   const reco::Muon* recoNegMuon = dynamic_cast<const reco::Muon*>(negMuon->originalObject());
        //   if ( isSameMuon(*recoPosMuon, *recoNegMuon)) continue;
        //   if ( !muon::isGoodMuon(*recoNegMuon, muon::TMOneStationTight) ) continue;
        //   /// geometric overlap
        //   if ( muon::overlap( *recoPosMuon, *recoNegMuon ) )
        //   nOverlapMus++ ;
        //   /// shared segments
        //   int nSegments2 = recoNegMuon->numberOfMatches(reco::Muon::SegmentArbitration);
        //   if (nSegments2 == 0 || nSegments1 == 0) continue;
        //   double sf = muon::sharedSegments(*recoPosMuon, *recoNegMuon) / std::min<double>(nSegments1, nSegments2);
        //   if (sf > sharedFraction) {
        //     nSharingSegWith = 0;
        //     if ( !isBetterMuon(*recoPosMuon, *recoNegMuon) )
        //     nSharingSegWith++ ;
        //   }
        // }
        // muNOverlap->push_back( nOverlapMus ) ;
        // muNSharingSegWith->push_back( nSharingSegWith ) ;
        //

        // ////////////////// check for muon1 //////////////////
        TrackRef muPosTrack = posMuon->track();
        if ( muPosTrack.isNull() )
        continue;

        /// cuts on muon1
        if (recoPosMuon->track()->hitPattern().numberOfValidPixelHits() < MuMinPixHits
        || recoPosMuon->track()->hitPattern().numberOfValidStripHits() < MuMinSiHits
        || recoPosMuon->track()->chi2()/recoPosMuon->track()->ndof() > MuMaxNormChi
        || fabs(recoPosMuon->track()->dxy(RefVtx)) > MuMaxD0) {
          continue ;
        }

        ////////////////// check for muon2 //////////////////
        for ( std::vector<pat::Muon>::const_iterator negMuon = posMuon+1; negMuon != thePATMuonHandle->end(); ++negMuon) {

          if(negMuon->charge() >= 0) continue ;

          const reco::Muon* recoNegMuon = dynamic_cast<const reco::Muon *>(negMuon->originalObject()) ;

          if (muon::overlap(*recoPosMuon, *recoNegMuon) )
          continue ;

          TrackRef muNegTrack = negMuon->track() ;
          if ( muNegTrack.isNull() )
          continue ;
          /// cuts on muon2
          if (recoNegMuon->track()->hitPattern().numberOfValidPixelHits() < MuMinPixHits
          || recoNegMuon->track()->hitPattern().numberOfValidStripHits() < MuMinSiHits
          || recoNegMuon->track()->chi2()/recoPosMuon->track()->ndof() > MuMaxNormChi
          || fabs(recoNegMuon->track()->dxy(RefVtx)) > MuMaxD0) {
            continue ;
          }



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
          RefCountedKinematicParticle MuMuCand_fromFit = MuMuVertexFitTree->currentParticle();
          RefCountedKinematicVertex MuMuCand_vertex_fromFit = MuMuVertexFitTree->currentDecayVertex();
          MuMuVertexFitTree->movePointerToTheFirstChild();
          RefCountedKinematicParticle MuPosCand_fromFit = MuMuVertexFitTree->currentParticle();
          MuMuVertexFitTree->movePointerToTheNextChild();
          RefCountedKinematicParticle MuNegCand_fromFit = MuMuVertexFitTree->currentParticle();

          double dimuon_vx_fit = MuMuCand_vertex_fromFit->position().x();
          double dimuon_vy_fit = MuMuCand_vertex_fromFit->position().y();
          double dimuon_vz_fit = MuMuCand_vertex_fromFit->position().z();

          ////////////////// fill the MuMu vectors //////////////////
          if (MuMuCand_fromFit->currentState().mass() < JPsiMinMass  ||  MuMuCand_fromFit->currentState().mass() > JPsiMaxMass)
          continue ;

          float dimuon_ma_fit = MuMuCand_fromFit->currentState().mass();
          int   dimuon_ch_fit = MuMuCand_fromFit->currentState().particleCharge();
          float dimuon_px_fit = MuMuCand_fromFit->currentState().kinematicParameters().momentum().x();
          float dimuon_py_fit = MuMuCand_fromFit->currentState().kinematicParameters().momentum().y();
          float dimuon_pz_fit = MuMuCand_fromFit->currentState().kinematicParameters().momentum().z();
          float dimuon_en_fit = sqrt(dimuon_ma_fit*dimuon_ma_fit+dimuon_px_fit*dimuon_px_fit+dimuon_py_fit*dimuon_py_fit+dimuon_pz_fit*dimuon_pz_fit);

          reco::CompositeCandidate reco_ref_JPsi(dimuon_ch_fit,math::XYZTLorentzVector(dimuon_px_fit,dimuon_py_fit,dimuon_pz_fit,dimuon_en_fit),
                                                   math::XYZPoint(dimuon_vx_fit,dimuon_vy_fit,dimuon_vz_fit),443);
          pat::CompositeCandidate pat_ref_JPsi(reco_ref_JPsi);



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

          ref_Jpsi.push_back(pat_ref_JPsi);
          ref_mupos.push_back(pat_ref_PM);
          ref_muneg.push_back(pat_ref_NM);

          Jpsi_p4.push_back(recoNegMuon->p4() + recoPosMuon->p4());
          mupos_p4.push_back(recoPosMuon->p4());
          muneg_p4.push_back(recoNegMuon->p4());

          MuMuVtx_CL->push_back( ChiSquaredProbability((double)( MuMuCand_vertex_fromFit->chiSquared()),(double)( MuMuCand_vertex_fromFit->degreesOfFreedom())) );
          MuMuVtx_Chi2->push_back( MuMuCand_vertex_fromFit->chiSquared() ) ;

          muPos_MuMu_Chi2->push_back( MuPosCand_fromFit->chiSquared());
          muPos_MuMu_NDF->push_back( MuPosCand_fromFit->degreesOfFreedom());
          muNeg_MuMu_Chi2->push_back( MuNegCand_fromFit->chiSquared());
          muNeg_MuMu_NDF->push_back( MuNegCand_fromFit->degreesOfFreedom());


          // mumu_p4->push_back();
          //
          // MuMuMass->push_back( MuMuCand_fromFit->currentState().mass() );
          // MuMuDecayVtx_X->push_back( MuMuCand_vertex_fromFit->position().x() );
          // MuMuDecayVtx_Y->push_back( MuMuCand_vertex_fromFit->position().y() );
          // MuMuDecayVtx_Z->push_back( MuMuCand_vertex_fromFit->position().z() );
          // MuMuDecayVtx_XE->push_back( sqrt( MuMuCand_vertex_fromFit->error().cxx()) );
          // MuMuDecayVtx_YE->push_back( sqrt( MuMuCand_vertex_fromFit->error().cyy()) );
          // MuMuDecayVtx_ZE->push_back( sqrt( MuMuCand_vertex_fromFit->error().czz()) );
          // MuMuVtx_CL->push_back( ChiSquaredProbability((double)( MuMuCand_vertex_fromFit->chiSquared()),(double)( MuMuCand_vertex_fromFit->degreesOfFreedom())) );
          // MuMuVtx_Chi2->push_back( MuMuCand_vertex_fromFit->chiSquared() ) ;
          // MuMuPx->push_back( Mu1Cand_KP.momentum().x() + Mu2Cand_KP.momentum().x() );
          // MuMuPy->push_back( Mu1Cand_KP.momentum().y() + Mu2Cand_KP.momentum().y() );
          // MuMuPz->push_back( Mu1Cand_KP.momentum().z() + Mu2Cand_KP.momentum().z() );
          // mu1Idx->push_back(std::distance(thePATMuonHandle->begin(), posMuon));
          // mu2Idx->push_back(std::distance(thePATMuonHandle->begin(), negMuon));
          //
          // ////////////////// JPsi (MuMu) fit //////////////////
          // muPos_MuMu_Px->push_back( Mu1Cand_KP.momentum().x()); /// SEMRA for JPsi
          // muPos_MuMu_Py->push_back( Mu1Cand_KP.momentum().y());
          // muPos_MuMu_Pz->push_back( Mu1Cand_KP.momentum().z());
          // muPos_MuMu_Chi2->push_back( MuPosCand_fromFit->chiSquared());
          // muPos_MuMu_NDF->push_back( MuPosCand_fromFit->degreesOfFreedom());
          // muNeg_MuMu_Px->push_back( Mu2Cand_KP.momentum().x());
          // muNeg_MuMu_Py->push_back( Mu2Cand_KP.momentum().y());
          // muNeg_MuMu_Pz->push_back( Mu2Cand_KP.momentum().z());
          // muNeg_MuMu_Chi2->push_back( MuNegCand_fromFit->chiSquared());
          // muNeg_MuMu_NDF->push_back( MuNegCand_fromFit->degreesOfFreedom());


          Int_t dimuonType = 0;   //0 nothing,  1 J/psi  , 2 psi(2S)
          if ( MuMuCand_fromFit->currentState().mass() > JPsiMinMass  &&  MuMuCand_fromFit->currentState().mass() < JPsiMaxMass ) {
            dimuonType = 1 ;
          }

          if (Debug_) std::cout <<dimuonType <<std::endl;

          if (Debug_) std::cout <<"evt:" <<evtNum <<" MuMu with diMuonType = " <<dimuonType <<std::endl;
          //if (Debug_) std::cout << "POINT 0" << std::endl;
          // MuMuType->push_back(dimuonType);
          //if (Debug_) std::cout << "POINT  2" << std::endl;

          int ntriggers = TriggersForMatching_.size();
          if (Debug_) std::cout << "ntriggers: " << ntriggers << std::endl;

          for (int MatchTrig = 0; MatchTrig < ntriggers; MatchTrig++)
          {
            if (Debug_) std::cout << "MatchingTriggerResult[" << MatchTrig << "]: " << MatchingTriggerResult[MatchTrig] << std::endl;
            if ( MatchingTriggerResult[MatchTrig]!=0 )
            {
              //if (Debug_) std::cout << "POINT  3" << std::endl;
              if (Debug_) std::cout << "CHECKING FiltersForMatching_[" << MatchTrig << "]: " << FiltersForMatching_[MatchTrig] << std::endl;
              //if (Debug_) std::cout << "POINT  4" << std::endl;
              pat::TriggerObjectStandAloneCollection mu1HLTMatches = posMuon->triggerObjectMatchesByFilter( FiltersForMatching_[MatchTrig] );
              //if (Debug_) std::cout << "POINT  5" << std::endl;
              pat::TriggerObjectStandAloneCollection mu2HLTMatches = negMuon->triggerObjectMatchesByFilter( FiltersForMatching_[MatchTrig] );
              //if (Debug_) std::cout << "POINT  6" << std::endl;
              bool pass1 = mu1HLTMatches.size() > 0;
              bool pass2 = mu2HLTMatches.size() > 0;
              //if (Debug_) std::cout << "POINT  7" << std::endl;
              if ((pass1) && (pass2))
              {
                //if (Debug_) std::cout << "POINT  8" << std::endl;
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
              edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";
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

          ++nMuMu;
          muons.clear();


          //////////////////////////////////////////////////////////////////////
          /// for B0
          if (Debug_) std::cout <<"evt:"<<evtNum<< " is Invalid Muon ?  " <<isEventWithInvalidMu << std::endl;
          //if (skipJPsi && ( dimuonType == 1 ));
          //if (Debug_) std::cout << "POINT 11" <<std::endl;
          nTrk->push_back( thePATTrackHandle->size() ) ;
          // //if (Debug_) std::cout << "POINT 12" <<std::endl;
          // if (thePATTrackHandle->size() < 2) {
          //   //if (Debug_) std::cout << "POINT 13" <<std::endl;
          //   nX_pre0++;
          // }

          if (Debug_) std::cout<<"nmumu : "<<nMuMu<<std::endl;

          ////////////////// check tracks for kaonPos for B0 //////////////////
          for ( std::vector<pat::GenericParticle>::const_iterator trackPos = theKaonRefittedPATTrackHandle->begin(); trackPos != theKaonRefittedPATTrackHandle->end(); ++trackPos ) {
            if (Debug_) if (Debug_) std::cout << "POINT 15" <<std::endl;
            /// check track doesn't overlap with the MuMu candidate tracks
            if (trackPos->charge() <= 0) continue;
            if (trackPos->track().key() == recoPosMuon->track().key()  ||  trackPos->track().key() == recoNegMuon->track().key())
            continue ;

            if (Debug_) if (Debug_) std::cout << "POINT 16" <<std::endl;
            /// cuts on charged tracks
            if (( trackPos->track()->chi2()/trackPos->track()->ndof() > TrMaxNormChi2 )  ||  trackPos->pt() < TrMinPt)
            continue ; nX_pre3++ ;

            if (Debug_) std::cout << "POINT 17" <<std::endl;



            ////////////////// check tracks for kaonNeg for B0 //////////////////
            for ( std::vector<pat::GenericParticle>::const_iterator trackNeg = trackPos+1; trackNeg != theKaonRefittedPATTrackHandle->end(); ++trackNeg ){

              if (Debug_) std::cout << "POINT 18" <<std::endl;
              /// check that this second track doesn't overlap with the the first track candidate
              if (trackNeg->track().key() == trackPos->track().key())
              continue ; nX_pre4++ ;

              /// check track doesn't overlap with the MuMu candidate tracks
              if (trackNeg->track().key() == recoPosMuon->track().key()  ||  trackNeg->track().key() == recoNegMuon->track().key())
              continue ; nX_pre5++ ;
              if (trackNeg->charge() >= 0)
              continue ; nX_pre6++ ;
              /// cuts on charged tracks
              if ((trackNeg->track()->chi2() / trackNeg->track()->ndof() > TrMaxNormChi2)  ||  trackNeg->pt() < TrMinPt)
              continue; nX_pre7++ ;

              if (Debug_) std::cout << "POINT 19" <<std::endl;

              ////////////////// get the KK information //////////////////
              TransientTrack kaonPosTT( trackPos->track(), &(*bFieldHandle) );
              TransientTrack kaonNegTT( trackNeg->track(), &(*bFieldHandle) );
              KinematicParticleFactoryFromTransientTrack pFactory;

              /// initial chi2 and ndf before kinematic fits
              float chi = 0., ndf = 0.;

              if (Debug_) std::cout << "POINT 20" <<std::endl;

              std::vector<RefCountedKinematicParticle> kaons;
              kaons.push_back( pFactory.particle( kaonPosTT, kaon_mass, chi, ndf, small_sigma));
              kaons.push_back( pFactory.particle( kaonNegTT, kaon_mass, chi, ndf, small_sigma));
              KinematicParticleVertexFitter KKFitter;
              RefCountedKinematicTree KKVertexFitTree;
              KKVertexFitTree = KKFitter.fit(kaons);

              if (Debug_) std::cout << "POINT 21" <<std::endl;

              if (!KKVertexFitTree->isValid())
              continue ;

              KKVertexFitTree->movePointerToTheTop();
              RefCountedKinematicParticle KKCand_fromFit = KKVertexFitTree->currentParticle();
              RefCountedKinematicVertex KKCand_vertex_fromFit = KKVertexFitTree->currentDecayVertex();

              KKVertexFitTree->movePointerToTheFirstChild();
              RefCountedKinematicParticle kaonPosCand_fromFit = KKVertexFitTree->currentParticle();
              KKVertexFitTree->movePointerToTheNextChild();
              RefCountedKinematicParticle kaonNegCand_fromFit = KKVertexFitTree->currentParticle();

              double ditrack_vx_fit = KKCand_vertex_fromFit->position().x();
              double ditrack_vy_fit = KKCand_vertex_fromFit->position().y();
              double ditrack_vz_fit = KKCand_vertex_fromFit->position().z();

              if (Debug_) std::cout << "POINT 22" <<std::endl;

              ////////////////// fill the KK vectors //////////////////
              if (KKCand_fromFit->currentState().mass() < KKMinMass  ||  KKCand_fromFit->currentState().mass() > KKMaxMass)
              continue ;

              float ditrack_ma_fit = MuMuCand_fromFit->currentState().mass();
              int   ditrack_ch_fit = MuMuCand_fromFit->currentState().particleCharge();
              float ditrack_px_fit = MuMuCand_fromFit->currentState().kinematicParameters().momentum().x();
              float ditrack_py_fit = MuMuCand_fromFit->currentState().kinematicParameters().momentum().y();
              float ditrack_pz_fit = MuMuCand_fromFit->currentState().kinematicParameters().momentum().z();
              float ditrack_en_fit = sqrt(ditrack_ma_fit*ditrack_ma_fit+ditrack_px_fit*ditrack_px_fit+ditrack_py_fit*ditrack_py_fit+ditrack_pz_fit*ditrack_pz_fit);

              reco::CompositeCandidate reco_ref_Phi(ditrack_ch_fit,math::XYZTLorentzVector(ditrack_px_fit,ditrack_py_fit,ditrack_pz_fit,ditrack_en_fit),
                                                       math::XYZPoint(ditrack_vx_fit,ditrack_vy_fit,ditrack_vz_fit),443);
              pat::CompositeCandidate pat_ref_Phi(reco_ref_Phi);

              if (Debug_) std::cout << "POINT 23" <<std::endl;


              float kaonPos_ma_fit = kaonPosCand_fromFit->currentState().mass();
              int   kaonPos_ch_fit = kaonPosCand_fromFit->currentState().particleCharge();
              float kaonPos_px_fit = kaonPosCand_fromFit->currentState().kinematicParameters().momentum().x();
              float kaonPos_py_fit = kaonPosCand_fromFit->currentState().kinematicParameters().momentum().y();
              float kaonPos_pz_fit = kaonPosCand_fromFit->currentState().kinematicParameters().momentum().z();
              float kaonPos_en_fit = sqrt(kaonPos_ma_fit*kaonPos_ma_fit+kaonPos_px_fit*kaonPos_px_fit+kaonPos_py_fit*kaonPos_py_fit+kaonPos_pz_fit*kaonPos_pz_fit);

              reco::CompositeCandidate reco_ref_PK(kaonPos_ch_fit,math::XYZTLorentzVector(kaonPos_px_fit,kaonPos_py_fit,kaonPos_pz_fit,kaonPos_en_fit),
                                                       math::XYZPoint(ditrack_vx_fit,ditrack_vy_fit,ditrack_vz_fit),-13);
              pat::CompositeCandidate pat_ref_PK(reco_ref_PK);

              if (Debug_) std::cout << "POINT 24" <<std::endl;

              float kaonNeg_ma_fit = MuNegCand_fromFit->currentState().mass();
              int   kaonNeg_ch_fit = MuNegCand_fromFit->currentState().particleCharge();
              float kaonNeg_px_fit = MuNegCand_fromFit->currentState().kinematicParameters().momentum().x();
              float kaonNeg_py_fit = MuNegCand_fromFit->currentState().kinematicParameters().momentum().y();
              float kaonNeg_pz_fit = MuNegCand_fromFit->currentState().kinematicParameters().momentum().z();
              float kaonNeg_en_fit = sqrt(kaonNeg_ma_fit*kaonNeg_ma_fit+kaonNeg_px_fit*kaonNeg_px_fit+kaonNeg_py_fit*kaonNeg_py_fit+kaonNeg_pz_fit*kaonNeg_pz_fit);

              reco::CompositeCandidate reco_ref_NK(kaonNeg_ch_fit,math::XYZTLorentzVector(kaonNeg_px_fit,kaonNeg_py_fit,kaonNeg_pz_fit,kaonNeg_en_fit),
                                                       math::XYZPoint(ditrack_vx_fit,ditrack_vy_fit,ditrack_vz_fit),13);

              if (Debug_) std::cout << "POINT 24" <<std::endl;


              pat::CompositeCandidate pat_ref_NK(reco_ref_NK);

              ref_Phi.push_back(reco_ref_Phi);
              ref_kaonpos.push_back(reco_ref_PK);
              ref_kaonneg.push_back(reco_ref_NK);

              Phi_p4.push_back(trackPos->p4() + trackNeg->p4());
              kpos_p4.push_back(trackPos->p4());
              kneg_p4.push_back(trackNeg->p4());

              KKVtx_CL->push_back( ChiSquaredProbability((double)( KKCand_vertex_fromFit->chiSquared()),(double)( KKCand_vertex_fromFit->degreesOfFreedom())) );
              KKVtx_Chi2->push_back( MuMuCand_vertex_fromFit->chiSquared() ) ;

              if (Debug_) std::cout << "POINT 25" <<std::endl;

              kaonPos_KK_Chi2->push_back( kaonPosCand_fromFit->chiSquared());
              kaonPos_KK_NDF->push_back( kaonPosCand_fromFit->degreesOfFreedom());
              kaonNeg_KK_Chi2->push_back( kaonNegCand_fromFit->chiSquared());
              kaonNeg_KK_NDF->push_back( kaonNegCand_fromFit->degreesOfFreedom());
              //
              // KKMass->push_back( KKCand_fromFit->currentState().mass() );
              // KKDecayVtx_X->push_back( KKCand_vertex_fromFit->position().x() );
              // KKDecayVtx_Y->push_back( KKCand_vertex_fromFit->position().y() );
              // KKDecayVtx_Z->push_back( KKCand_vertex_fromFit->position().z() );
              // KKDecayVtx_XE->push_back( sqrt( KKCand_vertex_fromFit->error().cxx()) );
              // KKDecayVtx_YE->push_back( sqrt( KKCand_vertex_fromFit->error().cyy()) );
              // KKDecayVtx_ZE->push_back( sqrt( KKCand_vertex_fromFit->error().czz()) );
              //
              // KKPx->push_back( Ka1Cand_KP.momentum().x() + Ka2Cand_KP.momentum().x() );
              // KKPy->push_back( Ka1Cand_KP.momentum().y() + Ka2Cand_KP.momentum().y() );
              // KKPz->push_back( Ka1Cand_KP.momentum().z() + Ka2Cand_KP.momentum().z() );
              // ka1Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), trackPos));
              // ka2Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), trackNeg));
              //
              // ////////////////// Phi (KK) fit //////////////////
              // kaonPos_KK_Px->push_back( Ka1Cand_KP.momentum().x());
              // kaonPos_KK_Py->push_back( Ka1Cand_KP.momentum().y());
              // kaonPos_KK_Pz->push_back( Ka1Cand_KP.momentum().z());
              //
              // kaonNeg_KK_Px->push_back( Ka2Cand_KP.momentum().x());
              // kaonNeg_KK_Py->push_back( Ka2Cand_KP.momentum().y());
              // kaonNeg_KK_Pz->push_back( Ka2Cand_KP.momentum().z());
              //
              if (Debug_) std::cout << "POINT 26" <<std::endl;

              ++nKK;
              kaons.clear();

              // ////////////////// cuts on tracks' delta R for B0 //////////////////
              math::XYZTLorentzVector MuMu = (recoPosMuon->p4() + recoNegMuon->p4());
              //math::XYZTLorentzVector MuMuKK = (MuMu + trackPos->p4() + trackNeg->p4());
              // float MuMu_K1_DR = sqrt( pow(MuMu.eta() - trackPos->p4().eta(),2) + pow(MuMu.phi() - trackPos->p4().phi(), 2) );
              // float MuMu_K2_DR = sqrt( pow(MuMu.eta() - trackNeg->p4().eta(),2) + pow(MuMu.phi() - trackNeg->p4().phi(), 2) );
              // float MuMuKK_K1_DR = sqrt( pow(MuMuKK.eta() - trackPos->p4().eta(),2) + pow(MuMuKK.phi() - trackPos->p4().phi(), 2) );
              // float MuMuKK_K2_DR = sqrt( pow(MuMuKK.eta() - trackNeg->p4().eta(),2) + pow(MuMuKK.phi() - trackNeg->p4().phi(), 2) );
              //
              // DR_MuMu_K1->push_back(MuMu_K1_DR);
              // DR_MuMu_K2->push_back(MuMu_K2_DR);
              // DR_MuMuKK_K1->push_back(MuMuKK_K1_DR);
              // DR_MuMuKK_K2->push_back(MuMuKK_K2_DR);
              //
              //
              // if (UseXDR) {
              //   if (MuMuKK_K1_DR > XTrackMaxDR || MuMuKK_K2_DR > XTrackMaxDR)
              //   XTrackMaxDR = 2;
              // } else {
              //   if (MuMu_K1_DR > MuMuTrackMaxDR || MuMu_K2_DR > MuMuTrackMaxDR)
              //   MuMuTrackMaxDR = 3.5;
              // }
              // nX_pre8++ ;

              if (Debug_) std::cout << "POINT 27" <<std::endl;

              math::XYZTLorentzVector xCand = trackPos->p4() + trackNeg->p4() + MuMu;
              ////////////////// cuts on MuMuKK mass window for B0 //////////////////
              if (xCand.M() > MaxXMass  ||  xCand.M() < MinXMass)
              continue ; nX_pre9++ ;

              /// having two oppositely charged muons, and two oppositely charged tracks: try to vertex them
              //TransientTrack kaonPosTT( trackPos->track(), &(*bFieldHandle) );
              //TransientTrack kaonNegTT( trackNeg->track(), &(*bFieldHandle) );

              if (Debug_) std::cout << "POINT 28" <<std::endl;

              TransientTrack kaonPos_notRefit, kaonNeg_notRefit;
              bool notRefPos = false, notRefNeg = false;

              for ( std::vector<pat::GenericParticle>::const_iterator tNotRef = thePATTrackHandle->begin(); tNotRef != thePATTrackHandle->end(); ++tNotRef )
              {
                  if(notRefNeg && notRefPos)
                    break;

                  if ( tNotRef->track().key() == trackNeg->track().key() && ! notRefNeg) {
                    notRefNeg = true;
                    kaonNeg_notRefit = TransientTrack( tNotRef->track(), &(*bFieldHandle) ) ;
                    continue;
                  }

                  if ( tNotRef->track().key() == trackPos->track().key() && ! notRefPos) {
                    notRefPos = true;
                    kaonPos_notRefit = TransientTrack( tNotRef->track(), &(*bFieldHandle) ) ;
                    continue;
                  }


              }

              if (Debug_) std::cout << "POINT 29" <<std::endl;

              bool notRefittedPartner = notRefPos || notRefNeg;
              /// do mass constraint for MuMu cand and do mass constrained vertex fit for B0
              std::vector<RefCountedKinematicParticle> xDaughters,xDaughters_unref;
              xDaughters.push_back(pFactory.particle( muonPosTT, muon_mass, chi, ndf, small_sigma));
              xDaughters.push_back(pFactory.particle( muonNegTT, muon_mass, chi, ndf, small_sigma));
              xDaughters.push_back(pFactory.particle( kaonPosTT, kaon_mass, chi, ndf, small_sigma));
              xDaughters.push_back(pFactory.particle( kaonNegTT, kaon_mass, chi, ndf, small_sigma));

              RefCountedKinematicTree XVertexFitTree, XVertexFitTree_noKrefit ;
              KinematicConstrainedVertexFitter XFitter ;

              if (Debug_) std::cout << "POINT 30" <<std::endl;

              if (doMuMuMassConst) { // MassConst = 'MC' in the following

                if (Debug_) std::cout << "POINT 30.0" <<std::endl;

                MultiTrackKinematicConstraint *MuMu = 0;
                MuMu = new TwoTrackMassKinematicConstraint(JPsi_mass);

                if (Debug_) std::cout << "POINT 30.1" <<std::endl;

                XVertexFitTree = XFitter.fit( xDaughters, MuMu );

                if (Debug_) std::cout << "POINT 30.2" <<std::endl;

                if (notRefittedPartner && notRefNeg && notRefPos) { // use not refitted kaons

                  if (Debug_) std::cout << "POINT 30.3" <<std::endl;

                  xDaughters_unref.push_back(pFactory.particle( muonPosTT, muon_mass, chi, ndf, small_sigma));
                  if (Debug_) std::cout << "POINT 30.3.1" <<std::endl;
                  xDaughters_unref.push_back(pFactory.particle( muonNegTT, muon_mass, chi, ndf, small_sigma));
                  if (Debug_) std::cout << "POINT 30.3.2" <<std::endl;
                  xDaughters_unref.push_back(pFactory.particle( kaonPos_notRefit, kaon_mass, chi, ndf, small_sigma));
                  if (Debug_) std::cout << "POINT 30.3.3" <<std::endl;
                  xDaughters_unref.push_back(pFactory.particle( kaonNeg_notRefit, kaon_mass, chi, ndf, small_sigma));
                  if (Debug_) std::cout << "POINT 30.3.4" <<std::endl;

                  XVertexFitTree_noKrefit = XFitter.fit( xDaughters_unref, MuMu );
                  if (Debug_) std::cout << "POINT 30.4" <<std::endl;
                }
              }
              else {

                if (Debug_) std::cout << "POINT 30.5" <<std::endl;
                XVertexFitTree = XFitter.fit( xDaughters );
                if (Debug_) std::cout << "POINT 30.6" <<std::endl;

                if (notRefittedPartner && notRefNeg && notRefPos) { // use not refitted kaons

                  if (Debug_) std::cout << "POINT 30.7" <<std::endl;
                  xDaughters_unref.push_back(pFactory.particle( muonPosTT, muon_mass, chi, ndf, small_sigma));
                  xDaughters_unref.push_back(pFactory.particle( muonNegTT, muon_mass, chi, ndf, small_sigma));
                  xDaughters_unref.push_back(pFactory.particle( kaonPos_notRefit, kaon_mass, chi, ndf, small_sigma));
                  xDaughters_unref.push_back(pFactory.particle( kaonNeg_notRefit, kaon_mass, chi, ndf, small_sigma));

                  XVertexFitTree_noKrefit = XFitter.fit( xDaughters_unref );
                  if (Debug_) std::cout << "POINT 30.8" <<std::endl;
                }
              }

              if (Debug_) std::cout << "POINT 31" <<std::endl;

              if ( !XVertexFitTree->isValid() ) /// B0 variables started
              continue ; nX_pre10++ ;

              XVertexFitTree->movePointerToTheTop();
              RefCountedKinematicParticle XCand_fromMCFit = XVertexFitTree->currentParticle();
              RefCountedKinematicVertex XCand_vertex_fromMCFit = XVertexFitTree->currentDecayVertex();

              if ( !XCand_vertex_fromMCFit->vertexIsValid() )
              continue ; nX_pre11++ ;

              if ( XCand_vertex_fromMCFit->chiSquared() < 0  ||  XCand_vertex_fromMCFit->chiSquared() > 10000 )
              continue ; nX_pre12++ ;

              if (XCand_vertex_fromMCFit->chiSquared() / XCand_vertex_fromMCFit->degreesOfFreedom() > 10 )
              continue ; nX_pre13++;

              if ( XCand_fromMCFit->currentState().mass() > 100 )
              continue ; nX_pre14++ ;

              double xVtxProb = ChiSquaredProbability((double)(XCand_vertex_fromMCFit->chiSquared()), (double)(XCand_vertex_fromMCFit->degreesOfFreedom()));
              if ( xVtxProb < 0.001 ) //0.0001 )
              continue ; nX_pre15++ ;


              if (Debug_) std::cout << "POINT 32" <<std::endl;

              //////////////////// Lifetimes calculations for B0 ////////////////////
              TVector3 X_vtx((*XCand_vertex_fromMCFit).position().x(), (*XCand_vertex_fromMCFit).position().y(), 0) ;
              TVector3 X_pperp(XCand_fromMCFit->currentState().globalMomentum().x(), XCand_fromMCFit->currentState().globalMomentum().y(), 0);
              TVector3 X_vtx3D((*XCand_vertex_fromMCFit).position().x(), (*XCand_vertex_fromMCFit).position().y(), (*XCand_vertex_fromMCFit).position().z()) ;
              TVector3 X_pperp3D(XCand_fromMCFit->currentState().globalMomentum().x(),XCand_fromMCFit->currentState().globalMomentum().y(), XCand_fromMCFit->currentState().globalMomentum().z());

              AlgebraicVector3 X_v3pperp ;
              X_v3pperp[0] = X_pperp.x(); X_v3pperp[1] = X_pperp.y(); X_v3pperp[2] = 0.;
              TVector3 X_pvtx, X_pvtx3D, X_vdiff, X_vdiff3D ;
              double X_cosAlpha, X_cosAlpha3D, X_ctau ;
              VertexDistanceXY X_vdistXY ;
              Measurement1D X_distXY ;
              GlobalError X_v1e = (Vertex(*XCand_vertex_fromMCFit)).error();
              GlobalError X_v2e ;
              AlgebraicSymMatrix33 X_vXYe ;
              double X_ctauErr ;
              float X_lxy, X_lxyErr, X_lxyz, X_lxyzErr ;
              ROOT::Math::SVector<double, 3> X_vDiff, X_vDiff3D ; // needed by Similarity method

              if (Debug_) std::cout << "POINT 33" <<std::endl;

              ////////////////// Lifetime wrt PV for B0 //////////////////
              X_v2e = thePrimaryVtx.error();
              X_vXYe = X_v1e.matrix() + X_v2e.matrix() ;
              /// 2D
              X_pvtx.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), 0) ;
              X_vdiff = X_vtx - X_pvtx ;
              X_cosAlpha = X_vdiff.Dot(X_pperp) / (X_vdiff.Perp()*X_pperp.Perp()) ;
              X_lxy = X_vdiff.Perp();
              X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
              X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
              X_distXY = X_vdistXY.distance(Vertex(*XCand_vertex_fromMCFit), Vertex(thePrimaryVtx));
              X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
              X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2()) ;
              /// 3D
              X_pvtx3D.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), thePrimaryVtx.position().z());
              X_vdiff3D = X_vtx3D - X_pvtx3D;
              X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/(X_vdiff3D.Mag()*X_pperp3D.Mag());
              X_lxyz = X_vdiff3D.Mag();
              X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
              X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();


              if (Debug_) std::cout << "POINT 34" <<std::endl;


              ////////////////// Last cuts for B0 //////////////////
              //if ( !(X_ctau/X_ctauErr > 2.8) || !(X_cosAlpha > 0.8) ) /// Alexis, BC_CTau_CTauErr plot has a cut which we don't want so we closed it.
              //continue ;


              ////////////////// fill X candidate variables //////////////////

              float xcand_ma_fit = XCand_fromMCFit->currentState().mass();
              int   xcand_ch_fit = XCand_fromMCFit->currentState().particleCharge();
              float xcand_px_fit = XCand_fromMCFit->currentState().kinematicParameters().momentum().x();
              float xcand_py_fit = XCand_fromMCFit->currentState().kinematicParameters().momentum().y();
              float xcand_pz_fit = XCand_fromMCFit->currentState().kinematicParameters().momentum().z();
              float xcand_en_fit = sqrt(xcand_ma_fit*xcand_ma_fit+xcand_px_fit*xcand_px_fit+xcand_py_fit*xcand_py_fit+xcand_pz_fit*xcand_pz_fit);

              double xcand_vx_fit = XCand_vertex_fromMCFit->position().x();
              double xcand_vy_fit = XCand_vertex_fromMCFit->position().y();
              double xcand_vz_fit = XCand_vertex_fromMCFit->position().z();

              reco::CompositeCandidate reco_X(xcand_ch_fit,math::XYZTLorentzVector(xcand_px_fit,xcand_py_fit,xcand_pz_fit,xcand_en_fit),
                                                       math::XYZPoint(xcand_vx_fit,xcand_vy_fit,xcand_vz_fit),443);
              pat::CompositeCandidate pat_X(reco_X);

              xMass->push_back( XCand_fromMCFit->currentState().mass()) ;
              // xPx->push_back( XCand_fromMCFit->currentState().globalMomentum().x()) ;
              // xPy->push_back( XCand_fromMCFit->currentState().globalMomentum().y()) ;
              // xPz->push_back( XCand_fromMCFit->currentState().globalMomentum().z()) ;
              // xPxE->push_back( sqrt( XCand_fromMCFit->currentState().kinematicParametersError().matrix()(3,3) ) ) ;
              // xPyE->push_back( sqrt( XCand_fromMCFit->currentState().kinematicParametersError().matrix()(4,4) ) ) ;
              // xPzE->push_back( sqrt( XCand_fromMCFit->currentState().kinematicParametersError().matrix()(5,5) ) ) ;
              // xVtx_CL->push_back( xVtxProb );
              // xVtx_Chi2->push_back( XCand_vertex_fromMCFit->chiSquared() ) ;
              // xDecayVtx_X->push_back((*XCand_vertex_fromMCFit).position().x());
              // xDecayVtx_Y->push_back((*XCand_vertex_fromMCFit).position().y());
              // xDecayVtx_Z->push_back((*XCand_vertex_fromMCFit).position().z());
              // xDecayVtx_XE->push_back(sqrt((*XCand_vertex_fromMCFit).error().cxx()));
              // xDecayVtx_YE->push_back(sqrt((*XCand_vertex_fromMCFit).error().cyy()));
              // xDecayVtx_ZE->push_back(sqrt((*XCand_vertex_fromMCFit).error().czz()));

              XVertexFitTree->movePointerToTheFirstChild();
              RefCountedKinematicParticle x_muonP_fromMCFit = XVertexFitTree->currentParticle();
              XVertexFitTree->movePointerToTheNextChild();
              RefCountedKinematicParticle x_muonN_fromMCFit = XVertexFitTree->currentParticle();
              XVertexFitTree->movePointerToTheNextChild();
              RefCountedKinematicParticle x_kaonP_fromMCFit = XVertexFitTree->currentParticle();
              XVertexFitTree->movePointerToTheNextChild();
              RefCountedKinematicParticle x_kaonN_fromMCFit = XVertexFitTree->currentParticle();
              /// muon1 & muon2

              float x_muonP_ma_fit = x_muonP_fromMCFit->currentState().mass();
              int   x_muonP_ch_fit = x_muonP_fromMCFit->currentState().particleCharge();
              float x_muonP_px_fit = x_muonP_fromMCFit->currentState().kinematicParameters().momentum().x();
              float x_muonP_py_fit = x_muonP_fromMCFit->currentState().kinematicParameters().momentum().y();
              float x_muonP_pz_fit = x_muonP_fromMCFit->currentState().kinematicParameters().momentum().z();
              float x_muonP_en_fit = sqrt(x_muonP_ma_fit*x_muonP_ma_fit+x_muonP_px_fit*x_muonP_px_fit+x_muonP_py_fit*x_muonP_py_fit+x_muonP_pz_fit*x_muonP_pz_fit);

              reco::CompositeCandidate reco_X_muonP(x_muonP_ch_fit,math::XYZTLorentzVector(x_muonP_px_fit,x_muonP_py_fit,x_muonP_pz_fit,x_muonP_en_fit),
                                                       math::XYZPoint(xcand_vx_fit,xcand_vy_fit,xcand_vz_fit),-13);
              pat::CompositeCandidate pat_X_muonP(reco_X_muonP);

              if (Debug_) std::cout << "POINT 35" <<std::endl;

              float x_muonN_ma_fit = x_muonN_fromMCFit->currentState().mass();
              int   x_muonN_ch_fit = x_muonN_fromMCFit->currentState().particleCharge();
              float x_muonN_px_fit = x_muonN_fromMCFit->currentState().kinematicParameters().momentum().x();
              float x_muonN_py_fit = x_muonN_fromMCFit->currentState().kinematicParameters().momentum().y();
              float x_muonN_pz_fit = x_muonN_fromMCFit->currentState().kinematicParameters().momentum().z();
              float x_muonN_en_fit = sqrt(x_muonN_ma_fit*x_muonN_ma_fit+x_muonN_px_fit*x_muonN_px_fit+x_muonN_py_fit*x_muonN_py_fit+x_muonN_pz_fit*x_muonN_pz_fit);

              reco::CompositeCandidate reco_X_muonN(x_muonN_ch_fit,math::XYZTLorentzVector(x_muonN_px_fit,x_muonN_py_fit,x_muonN_pz_fit,x_muonN_en_fit),
                                                       math::XYZPoint(xcand_vx_fit,xcand_vy_fit,xcand_vz_fit),13);
              pat::CompositeCandidate pat_X_muonN(reco_X_muonN);

              float x_kaonP_ma_fit = x_kaonP_fromMCFit->currentState().mass();
              int   x_kaonP_ch_fit = x_kaonP_fromMCFit->currentState().particleCharge();
              float x_kaonP_px_fit = x_kaonP_fromMCFit->currentState().kinematicParameters().momentum().x();
              float x_kaonP_py_fit = x_kaonP_fromMCFit->currentState().kinematicParameters().momentum().y();
              float x_kaonP_pz_fit = x_kaonP_fromMCFit->currentState().kinematicParameters().momentum().z();
              float x_kaonP_en_fit = sqrt(x_kaonP_ma_fit*x_kaonP_ma_fit+x_kaonP_px_fit*x_kaonP_px_fit+x_kaonP_py_fit*x_kaonP_py_fit+x_kaonP_pz_fit*x_kaonP_pz_fit);

              reco::CompositeCandidate reco_X_kaonP(x_kaonP_ch_fit,math::XYZTLorentzVector(x_kaonP_px_fit,x_kaonP_py_fit,x_kaonP_pz_fit,x_kaonP_en_fit),
                                                       math::XYZPoint(xcand_vx_fit,xcand_vy_fit,xcand_vz_fit),321);
              pat::CompositeCandidate pat_X_kaonP(reco_X_kaonP);

              if (Debug_) std::cout << "POINT 36" <<std::endl;

              float x_kaonN_ma_fit = x_kaonN_fromMCFit->currentState().mass();
              int   x_kaonN_ch_fit = x_kaonN_fromMCFit->currentState().particleCharge();
              float x_kaonN_px_fit = x_kaonN_fromMCFit->currentState().kinematicParameters().momentum().x();
              float x_kaonN_py_fit = x_kaonN_fromMCFit->currentState().kinematicParameters().momentum().y();
              float x_kaonN_pz_fit = x_kaonN_fromMCFit->currentState().kinematicParameters().momentum().z();
              float x_kaonN_en_fit = sqrt(x_kaonN_ma_fit*x_kaonN_ma_fit+x_kaonN_px_fit*x_kaonN_px_fit+x_kaonN_py_fit*x_kaonN_py_fit+x_kaonN_pz_fit*x_kaonN_pz_fit);

              reco::CompositeCandidate reco_X_kaonN(x_kaonN_ch_fit,math::XYZTLorentzVector(x_kaonN_px_fit,x_kaonN_py_fit,x_kaonN_pz_fit,x_kaonN_en_fit),
                                                       math::XYZPoint(xcand_vx_fit,xcand_vy_fit,xcand_vz_fit),-321);
              pat::CompositeCandidate pat_X_kaonN(reco_X_kaonN);

              // mu1Px_MuMuKK->push_back( muPos_MuMuKK->currentState().globalMomentum().x() );
              // mu1Py_MuMuKK->push_back( muPos_MuMuKK->currentState().globalMomentum().y() );
              // mu1Pz_MuMuKK->push_back( muPos_MuMuKK->currentState().globalMomentum().z() );
              // mu1E_MuMuKK->push_back( muPos_MuMuKK->currentState().kinematicParameters().energy() );
              // mu2Px_MuMuKK->push_back( muNeg_MuMuKK->currentState().globalMomentum().x() );
              // mu2Py_MuMuKK->push_back( muNeg_MuMuKK->currentState().globalMomentum().y() );
              // mu2Pz_MuMuKK->push_back( muNeg_MuMuKK->currentState().globalMomentum().z() );
              // mu2E_MuMuKK->push_back( muNeg_MuMuKK->currentState().kinematicParameters().energy() );
              // /// kaonPos & kaonNeg
              // k1Px_MuMuKK->push_back( k1_MuMuKK->currentState().globalMomentum().x() );
              // k1Py_MuMuKK->push_back( k1_MuMuKK->currentState().globalMomentum().y() );
              // k1Pz_MuMuKK->push_back( k1_MuMuKK->currentState().globalMomentum().z() );
              // k1E_MuMuKK->push_back( k1_MuMuKK->currentState().kinematicParameters().energy() );
              Double_t theo = 0., sigma = 0. ;
              kaonPos_nsigdedx->push_back( nsigmaofdedx(trackPos->track(),theo,sigma) );
              kaonPos_dedx->push_back( getEnergyLoss(trackPos->track()) );
              kaonPos_dedxMass->push_back( GetMass(trackPos->track()) );
              kaonPos_theo->push_back( theo );
              kaonPos_sigma->push_back( sigma );
              kaonPos_dedx_byHits->push_back( (dEdxTrack)[trackPos->track()].dEdx() );
              kaonPos_dedxErr_byHits->push_back( (dEdxTrack)[trackPos->track()].dEdxError() );
              kaonPos_saturMeas_byHits->push_back( (dEdxTrack)[trackPos->track()].numberOfSaturatedMeasurements() );
              kaonPos_Meas_byHits->push_back( (dEdxTrack)[trackPos->track()].numberOfMeasurements() );
              // k2Px_MuMuKK->push_back( k2_MuMuKK->currentState().globalMomentum().x() );
              // k2Py_MuMuKK->push_back( k2_MuMuKK->currentState().globalMomentum().y() );
              // k2Pz_MuMuKK->push_back( k2_MuMuKK->currentState().globalMomentum().z() );
              // k2E_MuMuKK->push_back( k2_MuMuKK->currentState().kinematicParameters().energy() );
              theo = 0.; sigma = 0. ;
              kaonNeg_nsigdedx->push_back(nsigmaofdedx(trackNeg->track(),theo,sigma));
              kaonNeg_dedx->push_back(getEnergyLoss(trackNeg->track()));
              kaonNeg_dedxMass->push_back(GetMass(trackNeg->track()));
              kaonNeg_theo->push_back(theo);
              kaonNeg_sigma->push_back(sigma);
              kaonNeg_dedx_byHits->push_back( (dEdxTrack_Kaon)[trackNeg->track()].dEdx() );
              kaonNeg_dedxErr_byHits->push_back( (dEdxTrack_Kaon)[trackNeg->track()].dEdxError() );
              kaonNeg_saturMeas_byHits->push_back( (dEdxTrack_Kaon)[trackNeg->track()].numberOfSaturatedMeasurements() );
              kaonNeg_Meas_byHits->push_back( (dEdxTrack_Kaon)[trackNeg->track()].numberOfMeasurements() );
              /// PV
              xCosAlphaPV->push_back( X_cosAlpha ); xCosAlpha3DPV->push_back( X_cosAlpha3D );
              xCTauPV->push_back( X_ctau ); xCTauPVE->push_back( X_ctauErr );
              xLxyPV->push_back( X_lxy ); xLxyPVE->push_back( X_lxyErr );
              xLxyzPV->push_back( X_lxyz ); xLxyzPVE->push_back( X_lxyzErr );
              /// dxy, dz, dxyE, dzE for kaons from PV
              kaonPos_dxy_PV->push_back( trackPos->track()->dxy(RefVtx) );
              kaonPos_dz_PV->push_back( trackPos->track()->dz(RefVtx) );
              kaonNeg_dxy_PV->push_back( trackNeg->track()->dxy(RefVtx) );
              kaonNeg_dz_PV->push_back( trackNeg->track()->dz(RefVtx) );


              ////////////////// Lifetime wrt BS for B0 //////////////////
              X_v2e = theBeamSpotVtx.error();
              X_vXYe = X_v1e.matrix() + X_v2e.matrix();
              /// 2D
              X_pvtx.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), 0);
              X_vdiff = X_vtx - X_pvtx;
              X_cosAlpha = X_vdiff.Dot(X_pperp)/(X_vdiff.Perp()*X_pperp.Perp());
              X_lxy = X_vdiff.Perp();
              X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
              X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
              X_distXY = X_vdistXY.distance(Vertex(*XCand_vertex_fromMCFit), Vertex(theBeamSpotVtx));
              X_ctau = X_distXY.value() * X_cosAlpha * (XCand_fromMCFit->currentState().mass() / X_pperp.Perp()) ;
              X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass()/X_pperp.Perp2();
              /// 3D
              X_pvtx3D.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), theBeamSpotVtx.position().z());
              X_vdiff3D = X_vtx3D - X_pvtx3D;
              X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/(X_vdiff3D.Mag()*X_pperp3D.Mag());
              X_lxyz = X_vdiff3D.Mag();
              X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
              X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();


              ////////////////// BS (beam spot) for B0 //////////////////
              xCosAlphaBS->push_back( X_cosAlpha ); xCosAlpha3DBS->push_back( X_cosAlpha3D );
              xCTauBS->push_back( X_ctau ); xCTauBSE->push_back( X_ctauErr );
              xLxyBS->push_back( X_lxy ); xLxyBSE->push_back( X_lxyErr );
              xLxyzBS->push_back( X_lxyz ); xLxyzBSE->push_back( X_lxyzErr );

              std::vector<TransientVertex> X_pvs ;
              Vertex XLessPV = thePrimaryVtx ;


              if (Debug_) std::cout << "POINT 37" <<std::endl;


              if (addXlessPrimaryVertex_)
              {
                VertexReProducer revertex(recVtxs, iEvent);
                Handle<TrackCollection> pvtracks;
                iEvent.getByLabel(revertex.inputTracks(), pvtracks);
                Handle<BeamSpot>        pvbeamspot;
                iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);

                if (pvbeamspot.id() != beamSpotHandle.id() )
                edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";

                const reco::Muon *Xrmu_1 = dynamic_cast<const reco::Muon *>(posMuon->originalObject());
                const reco::Muon *Xrmu_2 = dynamic_cast<const reco::Muon *>(negMuon->originalObject());

                if (Xrmu_1 != 0  &&  Xrmu_2 != 0  &&  Xrmu_1->track().id() == pvtracks.id()  &&  Xrmu_2->track().id() == pvtracks.id()
                &&  trackPos->track().id() == pvtracks.id()  &&  trackNeg->track().id() ==  pvtracks.id()) {
                  std::vector<TransientTrack> XLess; // need TransientTrack to keep the TrackRef
                  XLess.reserve( pvtracks->size() );
                  Double_t removedTrksPtSq = 0. ;
                  for (size_t i = 0, n = pvtracks->size(); i < n; ++i) {
                    if (i == Xrmu_1->track().key()) { removedTrksPtSq += (Xrmu_1->track()->pt())*(Xrmu_1->track()->pt()) ;
                      continue; }
                      if (i == Xrmu_2->track().key()) { removedTrksPtSq += (Xrmu_2->track()->pt())*(Xrmu_2->track()->pt()) ;
                        continue; }
                        if (i == trackPos->track().key()) { removedTrksPtSq += (trackPos->track()->pt())*(trackPos->track()->pt()) ;
                          continue; }
                          if (i == trackNeg->track().key()) { removedTrksPtSq += (trackNeg->track()->pt())*(trackNeg->track()->pt()) ;
                            continue; }

                            reco::TrackRef trk_now(pvtracks, i) ;
                            TransientTrack transientTrack = theTTBuilder->build( trk_now );
                            transientTrack.setBeamSpot( beamSpot );
                            XLess.push_back( transientTrack );
                          }
                          if ( removedTrksPtSq > 0. ) {
                            X_pvs = revertex.makeVertices(XLess, *pvbeamspot, iSetup) ; // list of PV
                          } else
                          if (Debug_) std::cout <<"\n\\\\\\\\\\\\\\\\\\\\ excluded tracks pT^2 = 0 \\\\\\\\\\\\\\\\\\\\\n" <<std::endl ;
                          if ( !X_pvs.empty() ) {
                            XLessPV = Vertex(X_pvs.front());
                            XLessPV_tracksPtSq->push_back( vertexHigherPtSquared.sumPtSquared(XLessPV) ) ;
                            XLessPV_4tracksPtSq->push_back( removedTrksPtSq ) ;
                            if (Debug_) {
                              std::cout <<"\nXLessPV_z = " <<XLessPV.position().z() <<std::endl ;
                              std::cout <<"XLessPV_tracks = " <<XLessPV.tracksSize() <<std::endl ;
                              std::cout <<"XLessPV_tracksPtSq = " <<vertexHigherPtSquared.sumPtSquared(XLessPV) <<std::endl ;
                              std::cout <<"XLessPV_removedTracksPtSq = " <<removedTrksPtSq <<std::endl ;
                              std::cout <<"X_pvs->size() = " <<X_pvs.size() <<std::endl ;
                              std::cout <<"priVtx_tracks = " <<priVtx_tracks <<std::endl ;
                              std::cout <<"tracksPtSq_pV = " <<tracksPtSq_pV <<std::endl ;
                              std::cout <<"recVtxs->size() = " <<recVtxs->size() <<std::endl ;
                            }
                          }
                        }
                      }


                      PriVtxXLess_n->push_back( X_pvs.size() ) ;
                      xLessPvs.push_back( XLessPV);
                      // PriVtxXLess_X->push_back( XLessPV.position().x() ) ;
                      // PriVtxXLess_Y->push_back( XLessPV.position().y() ) ;
                      // PriVtxXLess_Z->push_back( XLessPV.position().z() ) ;
                      // PriVtxXLess_EX->push_back( XLessPV.xError() ) ;
                      // PriVtxXLess_EY->push_back( XLessPV.yError() ) ;
                      // PriVtxXLess_EZ->push_back( XLessPV.zError() ) ;
                      // PriVtxXLess_CL->push_back( ChiSquaredProbability( (double)(XLessPV.chi2()), (double)(XLessPV.ndof())) );
                      // PriVtxXLess_Chi2->push_back( XLessPV.chi2() ) ;
                      // PriVtxXLess_tracks->push_back( XLessPV.tracksSize() ) ;

                      /// dxy, dz, dxyE, dzE for kaons from BS
                      math::XYZPoint BSVtx;
                      BSVtx = theBeamSpotVtx.position();
                      kaonPos_dxy_BS->push_back( trackPos->track()->dxy(BSVtx) );
                      kaonPos_dz_BS->push_back( trackPos->track()->dz(BSVtx) );
                      kaonNeg_dxy_BS->push_back( trackNeg->track()->dxy(BSVtx) );
                      kaonNeg_dz_BS->push_back( trackNeg->track()->dz(BSVtx) );


                      ////////////////// Lifetime wrt B0LessPV for B0 //////////////////
                      X_v2e = XLessPV.error();
                      X_vXYe = X_v1e.matrix() + X_v2e.matrix();
                      /// 2D
                      X_pvtx.SetXYZ(XLessPV.position().x(), XLessPV.position().y(), 0) ;
                      X_vdiff = X_vtx - X_pvtx ;
                      X_cosAlpha = X_vdiff.Dot(X_pperp)/(X_vdiff.Perp()*X_pperp.Perp());
                      X_lxy = X_vdiff.Perp();
                      X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
                      X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
                      X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(XLessPV) ) ;
                      X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      /// 3D
                      X_pvtx3D.SetXYZ(XLessPV.position().x(), XLessPV.position().y(), XLessPV.position().z());
                      X_vdiff3D = X_vtx3D - X_pvtx3D;
                      X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/( X_vdiff3D.Mag()*X_pperp3D.Mag() );
                      X_lxyz = X_vdiff3D.Mag();
                      X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                      X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();

                      xCosAlphaXLessPV->push_back( X_cosAlpha ) ; xCosAlpha3DXLessPV->push_back( X_cosAlpha3D ) ;
                      xCTauXLessPV->push_back( X_ctau ) ; xCTauXLessPVE->push_back( X_ctauErr ) ;
                      xLxyXLessPV->push_back( X_lxy ) ; xLxyXLessPVE->push_back( X_lxyErr ) ;
                      xLxyzXLessPV->push_back( X_lxyz ) ; xLxyzXLessPVE->push_back( X_lxyzErr ) ;

                      /// dxy, dz, dxyE, dzE for kaons from B0LessPV
                      math::XYZPoint XLessPVvtx;
                      XLessPVvtx = XLessPV.position();
                      kaonPos_dxy_XLessPV->push_back( trackPos->track()->dxy(XLessPVvtx) );
                      kaonPos_dz_XLessPV->push_back( trackPos->track()->dz(XLessPVvtx) );
                      kaonNeg_dxy_XLessPV->push_back( trackNeg->track()->dxy(XLessPVvtx) );
                      kaonNeg_dz_XLessPV->push_back( trackNeg->track()->dz(XLessPVvtx) );

                      kaonPos_dxyE->push_back( trackPos->track()->dxyError() );
                      kaonPos_dzE->push_back( trackPos->track()->dzError() );
                      kaonNeg_dxyE->push_back( trackNeg->track()->dxyError() );
                      kaonNeg_dzE->push_back( trackNeg->track()->dzError() );


                      /// Find the PV among the original offlinePV with the largest B0_cos(alpha)
                      Vertex theCosAlphaV = thePrimaryVtx ;
                      float maxCosAlpha = -1. ;

                      for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
                        X_pvtx.SetXYZ(itv->position().x(), itv->position().y(), 0) ;
                        X_vdiff = X_vtx - X_pvtx ;
                        float cosAlpha_temp = X_vdiff.Dot(X_pperp) / (X_vdiff.Perp()*X_pperp.Perp()) ; // Perp() == Mag() when z = 0

                        if ( cosAlpha_temp > maxCosAlpha ) {
                          maxCosAlpha = cosAlpha_temp ;
                          theCosAlphaV = Vertex(*itv) ;
                        }
                      }

                      PriVtx_XCosAlpha_n->push_back( recVtxs->size() ) ;
                      xCosAlphaPVs.push_back(theCosAlphaV);
                      // PriVtx_XCosAlpha_X->push_back( theCosAlphaV.position().x() ) ;
                      // PriVtx_XCosAlpha_Y->push_back( theCosAlphaV.position().y() ) ;
                      // PriVtx_XCosAlpha_Z->push_back( theCosAlphaV.position().z() ) ;
                      // PriVtx_XCosAlpha_EX->push_back( theCosAlphaV.xError() ) ;
                      // PriVtx_XCosAlpha_EY->push_back( theCosAlphaV.yError() ) ;
                      // PriVtx_XCosAlpha_EZ->push_back( theCosAlphaV.zError() ) ;
                      // PriVtx_XCosAlpha_CL->push_back( ChiSquaredProbability((double)(theCosAlphaV.chi2()), (double)(theCosAlphaV.ndof())) ) ;
                      // PriVtx_XCosAlpha_Chi2->push_back( theCosAlphaV.chi2() ) ;
                      // PriVtx_XCosAlpha_tracks->push_back( theCosAlphaV.tracksSize() ) ;


                      /// Find the PV among the original offlinePV with the largest B0_cos(alpha) 3D
                      Vertex theCosAlpha3DV = thePrimaryVtx ;
                      float maxCosAlpha3D = -1. ;

                      for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
                        X_pvtx3D.SetXYZ(itv->position().x(), itv->position().y(), itv->position().z()) ;
                        X_vdiff3D = X_vtx3D - X_pvtx3D ;
                        float cosAlpha_temp3D = X_vdiff3D.Dot(X_pperp3D) / (X_vdiff3D.Mag()*X_pperp3D.Mag()) ;

                        if ( cosAlpha_temp3D > maxCosAlpha3D ) {
                          maxCosAlpha3D = cosAlpha_temp3D ;
                          theCosAlpha3DV = Vertex(*itv) ;
                        }
                      }


                      ////////////////// Lifetime wrt PV with largest B0_cos(alpha) candidate //////////////////
                      X_v2e = theCosAlpha3DV.error();
                      X_vXYe = X_v1e.matrix() + X_v2e.matrix();
                      /// 2D
                      X_pvtx.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), 0) ;
                      X_vdiff = X_vtx - X_pvtx ;
                      X_cosAlpha = X_vdiff.Dot(X_pperp)/(X_vdiff.Perp()*X_pperp.Perp()); ;
                      X_lxy = X_vdiff.Perp();
                      X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
                      X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
                      X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
                      X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      X_lxy = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
                      /// 3D
                      X_pvtx3D.SetXYZ(theCosAlpha3DV.position().x(), theCosAlpha3DV.position().y(), theCosAlpha3DV.position().z()) ;
                      X_vdiff3D = X_vtx3D - X_pvtx3D ;
                      X_cosAlpha3D =  maxCosAlpha3D ;
                      X_lxyz = X_vdiff3D.Mag();
                      X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                      X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();

                      xCosAlphaPVCosAlpha3D->push_back( X_cosAlpha ) ; xCosAlpha3DPVCosAlpha3D->push_back( X_cosAlpha3D ) ;
                      xCTauPVCosAlpha3D->push_back( X_ctau ) ; xCTauPVCosAlpha3DE->push_back( X_ctauErr ) ;
                      xLxyPVCosAlpha3D->push_back( X_lxy ) ; xLxyPVCosAlpha3DE->push_back( X_lxyErr ) ;
                      xLxyzPVCosAlpha3D->push_back( X_lxyz ) ; xLxyzPVCosAlpha3DE->push_back( X_lxyzErr ) ;


                      /// Find the PV among the B0lessPV with the largest B0_cos(alpha)
                      Vertex theXLessCosAlphaV = thePrimaryVtx ;
                      maxCosAlpha = -1. ;

                      for (std::vector<TransientVertex>::iterator itv = X_pvs.begin(), itvend = X_pvs.end(); itv != itvend; ++itv) {
                        X_pvtx.SetXYZ(itv->position().x(), itv->position().y(), 0) ;
                        X_vdiff = X_vtx - X_pvtx ;
                        float cosAlpha_temp = X_vdiff.Dot(X_pperp) / (X_vdiff.Perp()*X_pperp.Perp()) ; // Perp() == Mag() when z = 0

                        if ( cosAlpha_temp > maxCosAlpha ) {
                          maxCosAlpha = cosAlpha_temp ;
                          theXLessCosAlphaV = Vertex(*itv) ;
                        }
                      }

                      PriVtxXLess_XCosAlpha_n->push_back( X_pvs.size() ) ;
                      xCosAlphaXLessPVs.push_back( theXLessCosAlphaV) ;
                      // PriVtxXLess_XCosAlpha_Y->push_back( theXLessCosAlphaV.position().y() ) ;
                      // PriVtxXLess_XCosAlpha_Z->push_back( theXLessCosAlphaV.position().z() ) ;
                      // PriVtxXLess_XCosAlpha_EX->push_back( theXLessCosAlphaV.xError() ) ;
                      // PriVtxXLess_XCosAlpha_EY->push_back( theXLessCosAlphaV.yError() ) ;
                      // PriVtxXLess_XCosAlpha_EZ->push_back( theXLessCosAlphaV.zError() ) ;
                      // PriVtxXLess_XCosAlpha_CL->push_back( ChiSquaredProbability((double)(theXLessCosAlphaV.chi2()), (double)(theXLessCosAlphaV.ndof())) ) ;
                      // PriVtxXLess_XCosAlpha_Chi2->push_back( theXLessCosAlphaV.chi2() ) ;
                      // PriVtxXLess_XCosAlpha_tracks->push_back( theXLessCosAlphaV.tracksSize() ) ;


                      ////////////////// Lifetime wrt PV with largest B0_cos(alpha) candidate //////////////////
                      X_v2e = theCosAlphaV.error();
                      X_vXYe = X_v1e.matrix() + X_v2e.matrix();
                      /// 2D
                      X_pvtx.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), 0) ;
                      X_vdiff = X_vtx - X_pvtx ;
                      X_cosAlpha =  maxCosAlpha ;
                      X_lxy = X_vdiff.Perp();
                      X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
                      X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
                      X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
                      X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      X_lxy = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
                      /// 3D
                      X_pvtx3D.SetXYZ(theCosAlphaV.position().x(), theCosAlphaV.position().y(), theCosAlphaV.position().z());
                      X_vdiff3D = X_vtx3D - X_pvtx3D;
                      X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/( X_vdiff3D.Mag()*X_pperp3D.Mag() );
                      X_lxyz = X_vdiff3D.Mag();
                      X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                      X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();

                      xCosAlphaPVCosAlpha->push_back( X_cosAlpha ) ; xCosAlpha3DPVCosAlpha->push_back( X_cosAlpha3D ) ;
                      xCTauPVCosAlpha->push_back( X_ctau ) ; xCTauPVCosAlphaE->push_back( X_ctauErr ) ;
                      xLxyPVCosAlpha->push_back( X_lxy ) ; xLxyPVCosAlphaE->push_back( X_lxyErr ) ;
                      xLxyzPVCosAlpha->push_back( X_lxyz ) ; xLxyzPVCosAlphaE->push_back( X_lxyzErr ) ;

                      PriVtx_XCosAlpha3D_n->push_back( recVtxs->size() ) ;
                      xCosAlpha3DPVs.push_back( theCosAlpha3DV) ;
                      // PriVtx_XCosAlpha3D_X->push_back( theCosAlpha3DV.position().x() ) ;
                      // PriVtx_XCosAlpha3D_Y->push_back( theCosAlpha3DV.position().y() ) ;
                      // PriVtx_XCosAlpha3D_Z->push_back( theCosAlpha3DV.position().z() ) ;
                      // PriVtx_XCosAlpha3D_EX->push_back( theCosAlpha3DV.xError() ) ;
                      // PriVtx_XCosAlpha3D_EY->push_back( theCosAlpha3DV.yError() ) ;
                      // PriVtx_XCosAlpha3D_EZ->push_back( theCosAlpha3DV.zError() ) ;
                      // PriVtx_XCosAlpha3D_CL->push_back( ChiSquaredProbability((double)(theCosAlpha3DV.chi2()), (double)(theCosAlpha3DV.ndof())) ) ;
                      // PriVtx_XCosAlpha3D_Chi2->push_back( theCosAlpha3DV.chi2() ) ;
                      // PriVtx_XCosAlpha3D_tracks->push_back( theCosAlpha3DV.tracksSize() ) ;


                      ////////////////// Lifetime wrt B0LessPV with largest B0_cos(alpha) 3D candidate
                      X_v2e = theXLessCosAlphaV.error();
                      X_vXYe = X_v1e.matrix() + X_v2e.matrix();
                      /// 2D
                      X_pvtx.SetXYZ(theXLessCosAlphaV.position().x(), theXLessCosAlphaV.position().y(), 0) ;
                      X_vdiff = X_vtx - X_pvtx ;
                      X_cosAlpha =  maxCosAlpha ;
                      X_lxy = X_vdiff.Perp();
                      X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
                      X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp() ;
                      X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(theXLessCosAlphaV) ) ;
                      X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      X_lxy = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
                      /// 3D
                      X_pvtx3D.SetXYZ(theXLessCosAlphaV.position().x(), theXLessCosAlphaV.position().y(), theXLessCosAlphaV.position().z());
                      X_vdiff3D = X_vtx3D - X_pvtx3D;
                      X_cosAlpha3D = X_vdiff3D.Dot(X_pperp3D)/( X_vdiff3D.Mag()*X_pperp3D.Mag() );
                      X_lxyz = X_vdiff3D.Mag();
                      X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                      X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();

                      xCosAlphaXLessPVCosAlpha->push_back( X_cosAlpha ) ; xCosAlpha3DXLessPVCosAlpha->push_back( X_cosAlpha3D ) ;
                      xCTauXLessPVCosAlpha->push_back( X_ctau ) ; xCTauXLessPVCosAlphaE->push_back( X_ctauErr ) ;
                      xLxyXLessPVCosAlpha->push_back( X_lxy ) ; xLxyXLessPVCosAlphaE->push_back( X_lxyErr ) ;
                      xLxyzXLessPVCosAlpha->push_back( X_lxyz ) ; xLxyzXLessPVCosAlphaE->push_back( X_lxyzErr ) ;


                      /// Find the PV among the B0lessPV with the largest B0_cos(alpha) 3D
                      Vertex theXLessCosAlpha3DV = thePrimaryVtx ;
                      maxCosAlpha3D = -1. ;

                      for (std::vector<TransientVertex>::iterator itv = X_pvs.begin(), itvend = X_pvs.end(); itv != itvend; ++itv) {
                        X_pvtx3D.SetXYZ(itv->position().x(), itv->position().y(), itv->position().z()) ;
                        X_vdiff3D = X_vtx3D - X_pvtx3D ;
                        float cosAlpha_temp3D = X_vdiff3D.Dot(X_pperp3D) / (X_vdiff3D.Mag()*X_pperp3D.Mag()) ;

                        if ( cosAlpha_temp3D > maxCosAlpha3D ) {
                          maxCosAlpha3D = cosAlpha_temp3D ;
                          theXLessCosAlpha3DV = Vertex(*itv) ;
                        }
                      }

                      PriVtxXLess_XCosAlpha3D_n->push_back( X_pvs.size() ) ;
                      xCosAlpha3DXLessPVs.push_back(theXLessCosAlpha3DV);
                      //
                      // PriVtxXLess_XCosAlpha3D_X->push_back( theXLessCosAlpha3DV.position().x() ) ;
                      // PriVtxXLess_XCosAlpha3D_Y->push_back( theXLessCosAlpha3DV.position().y() ) ;
                      // PriVtxXLess_XCosAlpha3D_Z->push_back( theXLessCosAlpha3DV.position().z() ) ;
                      // PriVtxXLess_XCosAlpha3D_EX->push_back( theXLessCosAlpha3DV.xError() ) ;
                      // PriVtxXLess_XCosAlpha3D_EY->push_back( theXLessCosAlpha3DV.yError() ) ;
                      // PriVtxXLess_XCosAlpha3D_EZ->push_back( theXLessCosAlpha3DV.zError() ) ;
                      // PriVtxXLess_XCosAlpha3D_CL->push_back( ChiSquaredProbability((double)(theXLessCosAlpha3DV.chi2()), (double)(theXLessCosAlpha3DV.ndof())) ) ;
                      // PriVtxXLess_XCosAlpha3D_Chi2->push_back( theXLessCosAlpha3DV.chi2() ) ;
                      // PriVtxXLess_XCosAlpha3D_tracks->push_back( theXLessCosAlpha3DV.tracksSize() ) ;


                      ////////////////// Lifetime wrt B0LessPV with largest B0_cos(alpha) 3D candidate
                      X_v2e = theXLessCosAlpha3DV.error();
                      X_vXYe = X_v1e.matrix() + X_v2e.matrix();
                      /// 2D
                      X_pvtx.SetXYZ(theXLessCosAlpha3DV.position().x(), theXLessCosAlpha3DV.position().y(), 0) ;
                      X_vdiff = X_vtx - X_pvtx ;
                      X_cosAlpha = X_vdiff.Dot(X_pperp)/(X_vdiff.Perp()*X_pperp.Perp());
                      X_lxy = X_vdiff.Perp();
                      X_vDiff[0] = X_vdiff.x(); X_vDiff[1] = X_vdiff.y(); X_vDiff[2] = 0 ; // needed by Similarity method
                      X_lxyErr = sqrt(ROOT::Math::Similarity(X_vDiff,X_vXYe)) / X_vdiff.Perp();
                      X_distXY = X_vdistXY.distance( Vertex(*XCand_vertex_fromMCFit), Vertex(theCosAlphaV) ) ;
                      X_ctau = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      X_ctauErr = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYe)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      X_lxy = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
                      /// 3D
                      X_pvtx3D.SetXYZ(theXLessCosAlpha3DV.position().x(), theXLessCosAlpha3DV.position().y(), theXLessCosAlpha3DV.position().z()) ;
                      X_vdiff3D = X_vtx3D - X_pvtx3D ;
                      X_cosAlpha3D =  maxCosAlpha3D ;
                      X_lxyz = X_vdiff3D.Mag();
                      X_vDiff3D[0] = X_vdiff3D.x(); X_vDiff3D[1] = X_vdiff3D.y(); X_vDiff3D[2] = X_vdiff3D.z() ;
                      X_lxyzErr = sqrt(ROOT::Math::Similarity(X_vDiff3D,X_vXYe)) / X_vdiff3D.Mag();
                      X_lxy = X_vdiff3D.Dot(X_pperp) / X_pperp.Mag() ;

                      xCosAlphaXLessPVCosAlpha3D->push_back( X_cosAlpha ) ; xCosAlpha3DXLessPVCosAlpha3D->push_back( X_cosAlpha3D ) ;
                      xCTauXLessPVCosAlpha3D->push_back( X_ctau ) ; xCTauXLessPVCosAlpha3DE->push_back( X_ctauErr ) ;
                      xLxyXLessPVCosAlpha3D->push_back( X_lxy ) ; xLxyXLessPVCosAlpha3DE->push_back( X_lxyErr ) ;
                      xLxyzXLessPVCosAlpha3D->push_back( X_lxyz ) ; xLxyzXLessPVCosAlpha3DE->push_back( X_lxyzErr ) ;


                      Vertex theOtherV = thePrimaryVtx;
                      if (resolveAmbiguity_) {
                        float minDz = 999999. ;
                        if (!addXlessPrimaryVertex_) {
                          for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv)
                          {
                            float deltaZ = fabs((*XCand_vertex_fromMCFit).position().z() - itv->position().z()) ;
                            if ( deltaZ < minDz ) {
                              minDz = deltaZ;
                              thePrimaryVtx = Vertex(*itv);
                              theOtherV = thePrimaryVtx;
                            }
                          }
                        } else {
                          for (std::vector<TransientVertex>::iterator itv2 = X_pvs.begin(), itvend2 = X_pvs.end(); itv2 != itvend2; ++itv2)
                          {
                            float deltaZ = fabs((*XCand_vertex_fromMCFit).position().z() - itv2->position().z()) ;
                            if ( deltaZ < minDz ) {
                              minDz = deltaZ;
                              Vertex XLessPV = Vertex(*itv2);
                              thePrimaryVtx = XLessPV;
                              theOtherV = XLessPV;
                            }
                          }
                        }
                      }

                      Vertex TheOtherVertex3D = thePrimaryVtx;
                      if (Debug_) std::cout<<" choose PV ="<< std::endl;
                      Int_t theXCorrPV_multiplicity = -1 ;
                      if (resolveAmbiguity_) {
                        float minDz = 999999.;
                        if (!addXlessPrimaryVertex_) {
                          theXCorrPV_multiplicity = recVtxs->size() ;
                          for (VertexCollection::const_iterator itv = recVtxs->begin(), itvend = recVtxs->end(); itv != itvend; ++itv) {
                            float deltaZ = fabs((*XCand_vertex_fromMCFit).position().z() - itv->position().z()) ;
                            if ( deltaZ < minDz ) {
                              minDz = deltaZ;
                              TheOtherVertex3D = Vertex(*itv);
                            }
                          }
                        } else {
                          theXCorrPV_multiplicity = X_pvs.size() ;
                          for (std::vector<TransientVertex>::iterator itv2 = X_pvs.begin(), itvend2 = X_pvs.end(); itv2 != itvend2; ++itv2) {
                            VertexDistance3D a3d;
                            float deltaZ   = a3d.distance(Vertex(*itv2), Vertex(*XCand_vertex_fromMCFit)).value();
                            if ( deltaZ < minDz ) {
                              minDz = deltaZ;
                              Vertex XLessPV = Vertex(*itv2);
                              TheOtherVertex3D = XLessPV;
                              //std::cout<<" z(X) - z(vtx) min="<<minDz<<std::endl;
                            }

                          }
                        }
                      }

                      PriVtxXCorr_n->push_back( theXCorrPV_multiplicity ) ;
                      corrPVs.push_back(thePrimaryVtx);
                      //
                      // PriVtxXCorr_X->push_back( thePrimaryVtx.position().x() ) ;
                      // PriVtxXCorr_Y->push_back( thePrimaryVtx.position().y() ) ;
                      // PriVtxXCorr_Z->push_back( thePrimaryVtx.position().z() ) ;
                      // PriVtxXCorr_EX->push_back( thePrimaryVtx.xError() ) ;
                      // PriVtxXCorr_EY->push_back( thePrimaryVtx.yError() ) ;
                      // PriVtxXCorr_EZ->push_back( thePrimaryVtx.zError() ) ;
                      // PriVtxXCorr_CL->push_back( ChiSquaredProbability( (double)(thePrimaryVtx.chi2()), (double)(thePrimaryVtx.ndof())) );
                      // PriVtxXCorr_Chi2->push_back( thePrimaryVtx.chi2() ) ;
                      // PriVtxXCorr_tracks->push_back( thePrimaryVtx.tracksSize() ) ;


                      ////////////////// Lifetime wrt PV with smaller longitudinal X impact parameter for B0  //////////////////
                      X_pvtx.SetXYZ(theOtherV.position().x(), theOtherV.position().y(), 0);
                      X_vdiff = X_vtx - X_pvtx;
                      X_cosAlpha = X_vdiff.Dot(X_pperp) / (X_vdiff.Perp()*X_pperp.Perp());
                      X_distXY = X_vdistXY.distance(Vertex(*XCand_vertex_fromMCFit), Vertex(theOtherV));
                      double X_ctauPVX = X_distXY.value() * X_cosAlpha * XCand_fromMCFit->currentState().mass() / X_pperp.Perp();
                      GlobalError X_v1eX = (Vertex(*XCand_vertex_fromMCFit)).error();
                      GlobalError X_v2eX = theOtherV.error();
                      AlgebraicSymMatrix33 X_vXYeX = X_v1eX.matrix() + X_v2eX.matrix();
                      double ctauErrPVX = sqrt(ROOT::Math::Similarity(X_v3pperp,X_vXYeX)) * XCand_fromMCFit->currentState().mass() / (X_pperp.Perp2());
                      float lxyPVX = X_vdiff.Dot(X_pperp) / X_pperp.Mag() ;
                      float lxyzPVX = X_vdiff3D.Dot(X_pperp3D) / X_pperp3D.Mag() ;
                      xCosAlphaPVX->push_back(X_cosAlpha);
                      xCTauPVX->push_back(X_ctauPVX); xCTauPVXE->push_back(ctauErrPVX);
                      xLxyPVX->push_back(lxyPVX);
                      xLxyzPVX->push_back(lxyzPVX);
                      VertexDistance3D a3d;
                      float Dist3DPV     = a3d.distance(TheOtherVertex3D, Vertex(*XCand_vertex_fromMCFit)).value() ;
                      float Dist3DPV_err = a3d.distance(TheOtherVertex3D, Vertex(*XCand_vertex_fromMCFit)).error() ;
                      xCTauPVX_3D->push_back(Dist3DPV);
                      xCTauPVX_3D_err->push_back(Dist3DPV_err);
                      //std::cout << Dist3DPV << " " << Dist3DPV_err << std::endl;
                      X_MuMuIdx->push_back(nMuMu-1);
                      // X_ka1Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), trackPos));
                      // X_ka2Idx->push_back(std::distance(theKaonRefittedPATTrackHandle->begin(), trackNeg));
                      nX++;
                      xDaughters.clear();
                      xDaughters_unref.clear();

                      ////////////////// flag for checking the Kaons from PV or not PV //////////////////
                      /// flag for kaonPos
                      std::vector<TransientTrack> vertexTrackskaonPos;
                      //std::cout << "\nthePrimaryVtx.tracksSize() = " << thePrimaryVtx.tracksSize() << std::endl;
                      //std::cout << "thePrimaryVtx.nTracks() = " << thePrimaryVtx.nTracks() << std::endl;
                      for ( std::vector<TrackBaseRef >::const_iterator iTrack = thePrimaryVtx.tracks_begin(); iTrack != thePrimaryVtx.tracks_end(); ++iTrack) {

                        TrackRef trackRefkaonPos = iTrack->castTo<TrackRef>();
                        //std::cout << "\ntrackRefkaonPos = " << trackRefkaonPos << std::endl;
                        //std::cout <<"before match" ;
                        if ( (trackPos->track().key() == trackRefkaonPos.key()) ) {
                          std::cout << "\ninside match" << std::endl;
                          TransientTrack kaonPosTT(trackRefkaonPos, &(*bFieldHandle) );
                          vertexTrackskaonPos.push_back(kaonPosTT);
                        }
                      }
                      //std::cout << "\nvertexTrackskaonPos.size() = " << vertexTrackskaonPos.size() << std::endl;
                      if (vertexTrackskaonPos.size()==0)
                      kaonPosFromPV->push_back(false);
                      else
                      kaonPosFromPV->push_back(true);

                      /// flag for kaonNeg
                      std::vector<TransientTrack> vertexTrackskaonNeg;
                      //std::cout << "\nthePrimaryVtx.tracksSize() = " << thePrimaryVtx.tracksSize() << std::endl;
                      //std::cout << "thePrimaryVtx.nTracks() = " << thePrimaryVtx.nTracks() << std::endl;
                      for ( std::vector<TrackBaseRef >::const_iterator iTrack = thePrimaryVtx.tracks_begin(); iTrack != thePrimaryVtx.tracks_end(); ++iTrack) {

                        TrackRef trackRefkaonNeg = iTrack->castTo<TrackRef>();
                        //std::cout << "\ntrackRefkaonNeg = " << trackRefkaonNeg << std::endl;
                        //std::cout <<"before match" ;
                        if (  (trackNeg->track().key() == trackRefkaonNeg.key()) ) {
                          TransientTrack kaonNegTT(trackRefkaonNeg, &(*bFieldHandle) );
                          vertexTrackskaonNeg.push_back(kaonNegTT);
                        }
                      }
                      //std::cout << "\nvertexTrackskaonPos.size() = " << vertexTrackskaonPos.size() << std::endl;
                      if (vertexTrackskaonNeg.size()==0)
                      kaonNegFromPV->push_back(false);
                      else
                      kaonNegFromPV->push_back(true);

                    } // 2nd loop over track (look for k2)
                  } // 1st loop over track (look for k1)
                } // 2nd loop over muons (look for mu-)
              } //first loop over muons (look for mu+)
            } // if (thePATMuonHandle->size() >= 2  && hasRequestedTrigger) {
            } // if (doMC || doData)
            // AT THE END OF THE EVENT fill the tree and clear the vectors
            // ===========================================================

            if (nX > 0)
            mumukktree->Fill() ;

            /// trigger stuff
            trigRes->clear(); trigNames->clear(); L1TT->clear(); MatchTriggerNames->clear();
            /// event numbers
            runNum = 0; evtNum = 0; lumiNum = 0;
            /// counters for x(4140)
            nMu = 0; nMuMu = 0; nX = 0; nKK = 0;
            nX_pre0 = 0; nX_pre1 = 0; nX_pre2 = 0; nX_pre3 = 0; nX_pre4 = 0; nX_pre5 = 0; nX_pre6 = 0; nX_pre7 = 0; nX_pre8 = 0; nX_pre9 = 0; nX_pre10 = 0; nX_pre11 = 0; nX_pre12 = 0; nX_pre13 = 0; nX_pre14 = 0; nX_pre15 = 0;
            //nX = 0;
            /// indices
            mu1Idx->clear(); mu2Idx->clear();
            ka1Idx->clear(); ka2Idx->clear();
            X_MuMuIdx->clear(); X_ka1Idx->clear(); X_ka2Idx->clear();

            /// MC Analysis
            if (doMC) {
              // Gen Primary Vertex
              n_genEvtVtx = 0;
              genEvtVtx_X->clear(); genEvtVtx_Y->clear(); genEvtVtx_Z->clear();
              genEvtVtx_particles->clear();
              n_XAncestors->clear();
              nMCAll = 0, nMCX = 0; //nMCXVtx = 0;
              // Gen Primary Vertex
              PriVtxGen_X->clear(); PriVtxGen_Y->clear(); PriVtxGen_Z->clear();
              PriVtxGen_EX->clear(); PriVtxGen_EY->clear(); PriVtxGen_EZ->clear();
              PriVtxGen_Chi2->clear(); PriVtxGen_CL->clear(); PriVtxGen_Ndof->clear();
              PriVtxGen_tracks->clear();

              MCPdgIdAll->clear(); MCDanNumAll->clear();
              MCJPsiPx->clear(); MCJPsiPy->clear(); MCJPsiPz->clear();
              MCmupPx->clear(); MCmupPy->clear(); MCmupPz->clear();
              MCmumPx->clear(); MCmumPy->clear(); MCmumPz->clear();
              MCPhiPx->clear(); MCPhiPy->clear(); MCPhiPz->clear();
              MCkpPx->clear(); MCkpPy->clear(); MCkpPz->clear();
              MCkmPx->clear(); MCkmPy->clear(); MCkmPz->clear();
              //MCpionPx->clear(); MCpionPy->clear(); MCpionPz->clear();
              //MCkaonPx->clear(); MCkaonPy->clear(); MCkaonPz->clear();
              //MCpionCh->clear(); MCkaonCh->clear();
              MCPx->clear(); MCPy->clear(); MCPz->clear();
            }
            if (Debug_) std::cout <<"after MC stuff clear" <<std::endl ;
            /// Primary Vertex
            n_pV = 0;
            tracksPtSq_pV = 0 ;
            // priVtx_X = 0; priVtx_Y = 0; priVtx_Z = 0 ;
            // priVtx_XE = 0; priVtx_YE = 0; priVtx_ZE = 0 ;
            // priVtx_NormChi2 = 0; priVtx_Chi2 = 0; priVtx_CL = 0; priVtx_tracks = 0;
            /// MuMu cand & KK cand
            MuMuMass->clear(); MuMuVtx_CL->clear(); MuMuVtx_Chi2->clear();
            MuMuPx->clear(); MuMuPy->clear(); MuMuPz->clear();
            MuMuDecayVtx_X->clear(); MuMuDecayVtx_Y->clear(); MuMuDecayVtx_Z->clear();
            MuMuDecayVtx_XE->clear(); MuMuDecayVtx_YE->clear(); MuMuDecayVtx_ZE->clear();
            MuMuMuonTrigMatch->clear();

            if (Debug_) std::cout <<"after mumus stuff clear" <<std::endl ;

            KKMass->clear(); KKPx->clear(); KKPy->clear(); KKPz->clear();
            KKVtx_CL->clear(); KKVtx_Chi2->clear();
            KKDecayVtx_X->clear(); KKDecayVtx_Y->clear(); KKDecayVtx_Z->clear();
            KKDecayVtx_XE->clear(); KKDecayVtx_YE->clear(); KKDecayVtx_ZE->clear();

            if (Debug_) std::cout <<"after kks stuff clear" <<std::endl ;

            /// muons from JPsi (MuMu) fit & kaons from Phi (KK) fit
            muPos_MuMu_Px->clear(); muPos_MuMu_Py->clear(); muPos_MuMu_Pz->clear(); muPos_MuMu_Chi2->clear(); muPos_MuMu_NDF->clear();
            muNeg_MuMu_Px->clear(); muNeg_MuMu_Py->clear(); muNeg_MuMu_Pz->clear(); muNeg_MuMu_Chi2->clear(); muNeg_MuMu_NDF->clear();
            MuMuType->clear();
            kaonPos_KK_Px->clear(); kaonPos_KK_Py->clear(); kaonPos_KK_Pz->clear(); kaonPos_KK_Chi2->clear(); kaonPos_KK_NDF->clear();
            kaonNeg_KK_Px->clear(); kaonNeg_KK_Py->clear();  kaonNeg_KK_Pz->clear(); kaonNeg_KK_Chi2->clear(); kaonNeg_KK_NDF->clear();
            // DR_MuMu_K1->clear(); DR_MuMu_K2->clear(); DR_MuMuKK_K1->clear(); DR_MuMuKK_K2->clear();
            if (Debug_) std::cout <<"after kaons stuff clear" <<std::endl ;
            /// Primary Vertex with "MuMu correction"
            mumuLessPvs_n.clear();
            PriVtxMuMuCorr_X->clear(); PriVtxMuMuCorr_Y->clear(); PriVtxMuMuCorr_Z->clear();
            PriVtxMuMuCorr_EX->clear(); PriVtxMuMuCorr_EY->clear(); PriVtxMuMuCorr_EZ->clear();
            PriVtxMuMuCorr_Chi2->clear(); PriVtxMuMuCorr_CL->clear(); PriVtxMuMuCorr_tracks->clear();
            nTrk->clear();
            /// X candidates
            xMass->clear(); xVtx_CL->clear(); xVtx_Chi2->clear();
            xPx->clear(); xPy->clear(); xPz->clear();
            xPxE->clear(); xPyE->clear(); xPzE->clear();
            xDecayVtx_X->clear(); xDecayVtx_Y->clear(); xDecayVtx_Z->clear();
            xDecayVtx_XE->clear(); xDecayVtx_YE->clear(); xDecayVtx_ZE->clear();
            if (Debug_) std::cout <<"after x cands stuff clear" <<std::endl ;
            /// Muons and tracks after X candidates fit
            mu1Px_MuMuKK->clear(); mu1Py_MuMuKK->clear(); mu1Pz_MuMuKK->clear(); mu1E_MuMuKK->clear();
            mu2Px_MuMuKK->clear(); mu2Py_MuMuKK->clear(); mu2Pz_MuMuKK->clear(); mu2E_MuMuKK->clear();
            k1Px_MuMuKK->clear(); k1Py_MuMuKK->clear(); k1Pz_MuMuKK->clear(); k1E_MuMuKK->clear();
            kaonPos_nsigdedx->clear(); kaonPos_dedx->clear(); kaonPos_dedxMass->clear(); kaonPos_theo->clear(); kaonPos_sigma->clear();
            kaonPos_dedx_byHits->clear(); kaonPos_dedxErr_byHits->clear(); kaonPos_saturMeas_byHits->clear(); kaonPos_Meas_byHits->clear();
            k2Px_MuMuKK->clear(); k2Py_MuMuKK->clear(); k2Pz_MuMuKK->clear(); k2E_MuMuKK->clear();
            kaonNeg_nsigdedx->clear(); kaonNeg_dedx->clear(); kaonNeg_dedxMass->clear(); kaonNeg_theo->clear(); kaonNeg_sigma->clear();
            kaonNeg_dedx_byHits->clear(); kaonNeg_dedxErr_byHits->clear(); kaonNeg_saturMeas_byHits->clear(); kaonNeg_Meas_byHits->clear();
            if (Debug_) std::cout <<"after mu tracks stuff clear" <<std::endl ;
            /// Primary Vertex with largest B0_cos(alpha)
            PriVtxXLess_n->clear();
            PriVtxXLess_X->clear(); PriVtxXLess_Y->clear(); PriVtxXLess_Z->clear();
            PriVtxXLess_EX->clear(); PriVtxXLess_EY->clear(); PriVtxXLess_EZ->clear();
            PriVtxXLess_Chi2->clear(); PriVtxXLess_CL->clear(); PriVtxXLess_tracks->clear();
            if (Debug_) std::cout <<"after pvxless stuff clear" <<std::endl ;
            XLessPV_tracksPtSq->clear(); XLessPV_4tracksPtSq->clear();
            PriVtx_XCosAlpha_n->clear();
            PriVtx_XCosAlpha_X->clear(); PriVtx_XCosAlpha_Y->clear(); PriVtx_XCosAlpha_Z->clear();
            PriVtx_XCosAlpha_EX->clear(); PriVtx_XCosAlpha_EY->clear(); PriVtx_XCosAlpha_EZ->clear();
            PriVtx_XCosAlpha_Chi2->clear(); PriVtx_XCosAlpha_CL->clear(); PriVtx_XCosAlpha_tracks->clear();
            PriVtxXLess_XCosAlpha_n->clear();
            PriVtxXLess_XCosAlpha_X->clear(); PriVtxXLess_XCosAlpha_Y->clear(); PriVtxXLess_XCosAlpha_Z->clear();
            PriVtxXLess_XCosAlpha_EX->clear(); PriVtxXLess_XCosAlpha_EY->clear(); PriVtxXLess_XCosAlpha_EZ->clear();
            PriVtxXLess_XCosAlpha_Chi2->clear(); PriVtxXLess_XCosAlpha_CL->clear(); PriVtxXLess_XCosAlpha_tracks->clear();
            if (Debug_) std::cout <<"after cosA  stuff clear" <<std::endl ;
            /// Primary Vertex with "B0 correction"
            PriVtxXCorr_n->clear();
            PriVtxXCorr_X->clear(); PriVtxXCorr_Y->clear(); PriVtxXCorr_Z->clear();
            PriVtxXCorr_EX->clear(); PriVtxXCorr_EY->clear(); PriVtxXCorr_EZ->clear();
            PriVtxXCorr_Chi2->clear(); PriVtxXCorr_CL->clear(); PriVtxXCorr_tracks->clear();
            if (Debug_) std::cout <<"after pvs stuff clear" <<std::endl ;
            /// Lifetime variables for B0
            xCosAlphaBS->clear(); xCosAlpha3DBS->clear(); xCTauBS->clear(); xCTauBSE->clear(); xLxyBS->clear(); xLxyBSE->clear(); xLxyzBS->clear(); xLxyzBSE->clear();
            xCosAlphaPV->clear(); xCosAlpha3DPV->clear(); xCTauPV->clear(); xCTauPVE->clear(); xLxyPV->clear(); xLxyPVE->clear(); xLxyzPV->clear(); xLxyzPVE->clear();
            xCosAlphaPVCosAlpha->clear(); xCosAlpha3DPVCosAlpha->clear(); xCTauPVCosAlpha->clear(); xCTauPVCosAlphaE->clear(); xLxyPVCosAlpha->clear(); xLxyPVCosAlphaE->clear(); xLxyzPVCosAlpha->clear(); xLxyzPVCosAlphaE->clear();
            xCosAlphaPVCosAlpha3D->clear(); xCosAlpha3DPVCosAlpha3D->clear(); xCTauPVCosAlpha3D->clear(); xCTauPVCosAlpha3DE->clear(); xLxyPVCosAlpha3D->clear(); xLxyPVCosAlpha3DE->clear(); xLxyzPVCosAlpha3D->clear(); xLxyzPVCosAlpha3DE->clear();
            xCosAlphaXLessPV->clear(); xCosAlpha3DXLessPV->clear(); xCTauXLessPV->clear() ; xCTauXLessPVE->clear() ; xLxyXLessPV->clear() ; xLxyXLessPVE->clear() ; xLxyzXLessPV->clear() ; xLxyzXLessPVE->clear() ;
            xCosAlphaXLessPVCosAlpha->clear(); xCosAlpha3DXLessPVCosAlpha->clear(); xCTauXLessPVCosAlpha->clear() ; xCTauXLessPVCosAlphaE->clear() ; xLxyXLessPVCosAlpha->clear() ; xLxyXLessPVCosAlphaE->clear() ; xLxyzXLessPVCosAlpha->clear() ; xLxyzXLessPVCosAlphaE->clear() ;
            xCosAlphaXLessPVCosAlpha3D->clear(); xCosAlpha3DXLessPVCosAlpha3D->clear(); xCTauXLessPVCosAlpha3D->clear() ; xCTauXLessPVCosAlpha3DE->clear() ; xLxyXLessPVCosAlpha3D->clear() ; xLxyXLessPVCosAlpha3DE->clear() ; xLxyzXLessPVCosAlpha3D->clear() ; xLxyzXLessPVCosAlpha3DE->clear() ;
            xCosAlphaPVX->clear(); xCTauPVX->clear(); xCTauPVXE->clear(); xLxyPVX->clear(); xLxyzPVX->clear();
            xCTauPVX_3D->clear(); xCTauPVX_3D_err->clear();
            if (Debug_) std::cout <<"after other ctau stuff clear" <<std::endl ;

            /// dxy, dz, dxyE, dzE for kaons from PV, BS, B0LessPV
            kaonPos_dxy_PV->clear(); kaonPos_dz_PV->clear(); kaonNeg_dxy_PV->clear(); kaonNeg_dz_PV->clear();
            kaonPos_dxy_BS->clear(); kaonPos_dz_BS->clear(); kaonNeg_dxy_BS->clear(); kaonNeg_dz_BS->clear();
            kaonPos_dxy_XLessPV->clear(); kaonPos_dz_XLessPV->clear(); kaonNeg_dxy_XLessPV->clear(); kaonNeg_dz_XLessPV->clear();
            kaonPos_dxyE->clear(); kaonPos_dzE->clear(); kaonNeg_dxyE->clear(); kaonNeg_dzE->clear();

            kaonPosFromPV->clear(); kaonNegFromPV->clear();

            if (Debug_) std::cout <<"before muon stuff clear" <<std::endl ;
            /// muons
            muPx->clear(); muPy->clear(); muPz->clear(); muCharge->clear();
            muD0->clear(); muDz->clear(); muChi2->clear(); muGlChi2->clear();
            mufHits->clear(); muFirstBarrel->clear(); muFirstEndCap->clear(); muD0E->clear() ;  muDzVtxErr->clear() ; muKey->clear() ;
            muIsGlobal->clear(); muIsPF->clear();
            muDzVtx->clear(); muDxyVtx->clear(); muGlMatchedStation->clear(); muGlDzVtx->clear(); muGlDxyVtx->clear();
            nMatchedStations->clear();
            muNDF->clear(); muGlNDF->clear(); muPhits->clear(); muShits->clear(); muGlMuHits->clear(); muType->clear();
            muQual->clear(); muTrack->clear(); muNOverlap->clear(); muNSharingSegWith->clear();

            if (Debug_) std::cout <<"after muon stuff clear" <<std::endl ;
            /// tracks
            tracks.clear();
            // trNotRef->clear(); trRef->clear();
            // trPx->clear(); trPy->clear(); trPz->clear(); trE->clear();
            // trNDF->clear(); trPhits->clear(); trShits->clear(); trChi2->clear();
            // trD0->clear(); trD0E->clear(); trCharge->clear();
            // trQualityHighPurity->clear(); trQualityTight->clear();
            // trfHits->clear(); trFirstBarrel->clear(); trFirstEndCap->clear();
            trDzVtx->clear(); trDxyVtx->clear();
            tr_nsigdedx->clear(); tr_dedx->clear(); tr_dedxMass->clear(); tr_theo->clear(); tr_sigma->clear();
            tr_dedx_byHits->clear(); tr_dedxErr_byHits->clear(); tr_saturMeas_byHits->clear(); tr_Meas_byHits->clear();

            if (Debug_) std::cout <<"end of branches clear" <<std::endl ;
          }
          //}/// analyze
          /// ------------ method called once each job just before starting event loop  ------------
          void MuMuKKPAT::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
          {
          }
          void MuMuKKPAT::beginJob()
          {
            edm::Service<TFileService> fs;

            mumukktree = fs->make<TTree>("JPsi/Phi Tree", "MuMuKK Data");

            mumukktree->Branch("TrigRes", &trigRes);
            mumukktree->Branch("TrigNames", &trigNames);
            mumukktree->Branch("MatchTriggerNames", &MatchTriggerNames);
            mumukktree->Branch("L1TrigRes", &L1TT);
            mumukktree->Branch("evtNum", &evtNum,"evtNum/i");
            mumukktree->Branch("runNum", &runNum,"runNum/i");
            mumukktree->Branch("lumiNum", &lumiNum, "lumiNum/i");

            mumukktree->Branch("pV", "reco::Vertex", &pV);
            mumukktree->Branch("n_pV", &n_pV, "n_pV/f");
            mumukktree->Branch("tracksPtSq_pV", &tracksPtSq_pV, "tracksPtSq_pV/f");

            // mumukktree->Branch("priVtx_X", &priVtx_X, "priVtx_X/f");
            // mumukktree->Branch("priVtx_Y", &priVtx_Y, "priVtx_Y/f");
            // mumukktree->Branch("priVtx_Z", &priVtx_Z, "priVtx_Z/f");
            // mumukktree->Branch("priVtx_XE", &priVtx_XE, "priVtx_XE/f");
            // mumukktree->Branch("priVtx_YE", &priVtx_YE, "priVtx_YE/f");
            // mumukktree->Branch("priVtx_ZE", &priVtx_ZE, "priVtx_ZE/f");
            // mumukktree->Branch("priVtx_NormChi2",&priVtx_NormChi2, "priVtx_NormChi2/f");
            // mumukktree->Branch("priVtx_Chi2",&priVtx_Chi2, "priVtx_Chi2/f");
            // mumukktree->Branch("priVtx_CL",&priVtx_CL, "priVtx_CL/f");
            // mumukktree->Branch("priVtx_tracks", &priVtx_tracks, "priVtx_tracks/i");

            /// MC Analysis
            if (doMC) {
              // Gen Primary Vertex
              mumukktree->Branch("genEvtVtx_X", &genEvtVtx_X);
              mumukktree->Branch("genEvtVtx_Y", &genEvtVtx_Y);
              mumukktree->Branch("genEvtVtx_Z", &genEvtVtx_Z);
              mumukktree->Branch("genEvtVtx_particles", &genEvtVtx_particles);
              mumukktree->Branch("n_XAncestors", &n_XAncestors);
              mumukktree->Branch("nMCAll", &nMCAll, "nMCAll/i");
              mumukktree->Branch("MCPdgIdAll", &MCPdgIdAll);
              mumukktree->Branch("MCDanNumAll", &MCDanNumAll);
              mumukktree->Branch("nMCX",&nMCX,"nMCX/i");
              // Gen Primary Vertex
              mumukktree->Branch("PriVtxGen_X",&PriVtxGen_X);
              mumukktree->Branch("PriVtxGen_Y",&PriVtxGen_Y);
              mumukktree->Branch("PriVtxGen_Z",&PriVtxGen_Z);
              mumukktree->Branch("PriVtxGen_EX",&PriVtxGen_EX);
              mumukktree->Branch("PriVtxGen_EY",&PriVtxGen_EY);
              mumukktree->Branch("PriVtxGen_EZ",&PriVtxGen_EZ);
              mumukktree->Branch("PriVtxGen_Chi2",&PriVtxGen_Chi2);
              mumukktree->Branch("PriVtxGen_CL",&PriVtxGen_CL);
              mumukktree->Branch("PriVtxGen_Ndof",&PriVtxGen_Ndof);
              mumukktree->Branch("PriVtxGen_tracks",&PriVtxGen_tracks);
              mumukktree->Branch("MCJPsiPx",&MCJPsiPx);
              mumukktree->Branch("MCJPsiPy",&MCJPsiPy);
              mumukktree->Branch("MCJPsiPz",&MCJPsiPz);
              mumukktree->Branch("MCmupPx",&MCmupPx);
              mumukktree->Branch("MCmupPy",&MCmupPy);
              mumukktree->Branch("MCmupPz",&MCmupPz);
              mumukktree->Branch("MCmumPx",&MCmumPx);
              mumukktree->Branch("MCmumPy",&MCmumPy);
              mumukktree->Branch("MCmumPz",&MCmumPz);
              mumukktree->Branch("MCPhiPx",&MCPhiPx);
              mumukktree->Branch("MCPhiPy",&MCPhiPy);
              mumukktree->Branch("MCPhiPz",&MCPhiPz);
              mumukktree->Branch("MCkpPx",&MCkpPx);
              mumukktree->Branch("MCkpPy",&MCkpPy);
              mumukktree->Branch("MCkpPz",&MCkpPz);
              mumukktree->Branch("MCkmPx",&MCkmPx);
              mumukktree->Branch("MCkmPy",&MCkmPy);
              mumukktree->Branch("MCkmPz",&MCkmPz);
              //mumukktree->Branch("MCpionPx",&MCpionPx);
              //mumukktree->Branch("MCpionPy",&MCpionPy);
              //mumukktree->Branch("MCpionPz",&MCpionPz);
              //mumukktree->Branch("MCpionCh",&MCpionCh);
              //mumukktree->Branch("MCkaonPx",&MCkaonPx);
              //mumukktree->Branch("MCkaonPy",&MCkaonPy);
              //mumukktree->Branch("MCkaonPz",&MCkaonPz);
              //mumukktree->Branch("MCkaonCh",&MCkaonCh);
              mumukktree->Branch("MCPx", &MCPx);
              mumukktree->Branch("MCPy", &MCPy);
              mumukktree->Branch("MCPz", &MCPz);
            }
            /// generic tracks
            mumukktree->Branch("tracks", "reco::TrackCollection", &tracks);
            // mumukktree->Branch("trNotRef", &trNotRef);
            // mumukktree->Branch("trRef", &trRef);
            // mumukktree->Branch("trackPx", &trPx);
            // mumukktree->Branch("trackPy", &trPy);
            // mumukktree->Branch("trackPz", &trPz);
            // mumukktree->Branch("trackEnergy", &trE);
            // mumukktree->Branch("trackNDF", &trNDF);
            // mumukktree->Branch("trackPhits", &trPhits);
            // mumukktree->Branch("trackShits", &trShits);
            // mumukktree->Branch("trackChi2", &trChi2);
            // mumukktree->Branch("trackD0", &trD0);
            // mumukktree->Branch("trackD0Err", &trD0E);
            // mumukktree->Branch("trackCharge", &trCharge);
            // mumukktree->Branch("TrackHighPurity", &trQualityHighPurity);
            // mumukktree->Branch("TrackTight", &trQualityTight);
            // mumukktree->Branch("trackfHits", &trfHits);
            // mumukktree->Branch("trackFirstBarrel", &trFirstBarrel);
            // mumukktree->Branch("trackFirstEndCap", &trFirstEndCap);
            mumukktree->Branch("trackDzVtx", &trDzVtx);
            mumukktree->Branch("trackDxyVtx", &trDxyVtx);
            mumukktree->Branch("tr_nsigdedx", &tr_nsigdedx);
            mumukktree->Branch("tr_dedx", &tr_dedx);
            mumukktree->Branch("tr_dedxMass", &tr_dedxMass);
            mumukktree->Branch("tr_theo", &tr_theo);
            mumukktree->Branch("tr_sigma", &tr_sigma);
            mumukktree->Branch("tr_dedx_byHits", &tr_dedx_byHits );
            mumukktree->Branch("tr_dedxErr_byHits", &tr_dedxErr_byHits );
            mumukktree->Branch("tr_saturMeas_byHits", &tr_saturMeas_byHits );
            mumukktree->Branch("tr_Meas_byHits", &tr_Meas_byHits );
            /// Generic muons
            mumukktree->Branch("nMu", &nMu, "nMu/i");
            mumukktree->Branch("muPx",&muPx);
            mumukktree->Branch("muPy",&muPy);
            mumukktree->Branch("muPz",&muPz);
            mumukktree->Branch("muCharge", &muCharge);
            mumukktree->Branch("muD0",&muD0);
            mumukktree->Branch("muDz",&muDz);
            mumukktree->Branch("muChi2",&muChi2);
            mumukktree->Branch("muNDF",&muNDF);
            mumukktree->Branch("muPhits",&muPhits);
            mumukktree->Branch("muShits",&muShits);
            mumukktree->Branch("muLayersTr",&muLayersTr);
            mumukktree->Branch("muLayersPix",&muLayersPix);
            mumukktree->Branch("muD0E",&muD0E);
            mumukktree->Branch("muDzVtxErr",&muDzVtxErr);
            mumukktree->Branch("muKey",&muKey);
            mumukktree->Branch("muIsGlobal",&muIsGlobal);
            mumukktree->Branch("muIsPF",&muIsPF);
            mumukktree->Branch("muGlMuHits",&muGlMuHits);
            mumukktree->Branch("muGlChi2",&muGlChi2);
            mumukktree->Branch("muGlNDF",&muGlNDF);
            mumukktree->Branch("muGlMatchedStation",&muGlMatchedStation);
            mumukktree->Branch("muGlDzVtx", &muGlDzVtx);
            mumukktree->Branch("muGlDxyVtx", &muGlDxyVtx);
            mumukktree->Branch("nMatchedStations", &nMatchedStations);
            mumukktree->Branch("muType",&muType);
            mumukktree->Branch("muQual",&muQual);
            mumukktree->Branch("muTrack",&muTrack);
            mumukktree->Branch("muNOverlap",&muNOverlap);
            mumukktree->Branch("muNSharingSegWith",&muNSharingSegWith);
            mumukktree->Branch("mufHits", &mufHits);
            mumukktree->Branch("muFirstBarrel", &muFirstBarrel);
            mumukktree->Branch("muFirstEndCap", &muFirstEndCap);
            mumukktree->Branch("muDzVtx", &muDzVtx);
            mumukktree->Branch("muDxyVtx", &muDxyVtx);
            /// MuMu cand
            mumukktree->Branch("nMuMu",&nMuMu,"nMuMu/i");
            mumukktree->Branch("ref_Jpsi", "reco::CompositeCandidateCollection", &ref_Jpsi);
            mumukktree->Branch("MuMuMass",&MuMuMass);
            mumukktree->Branch("MuMuPx",&MuMuPx);
            mumukktree->Branch("MuMuPy",&MuMuPy);
            mumukktree->Branch("MuMuPz",&MuMuPz);
            mumukktree->Branch("MuMuVtx_CL",&MuMuVtx_CL);
            mumukktree->Branch("MuMuVtx_Chi2",&MuMuVtx_Chi2);
            mumukktree->Branch("MuMuDecayVtx_X",&MuMuDecayVtx_X);
            mumukktree->Branch("MuMuDecayVtx_Y",&MuMuDecayVtx_Y);
            mumukktree->Branch("MuMuDecayVtx_Z",&MuMuDecayVtx_Z);
            mumukktree->Branch("MuMuDecayVtx_XE",&MuMuDecayVtx_XE);
            mumukktree->Branch("MuMuDecayVtx_YE",&MuMuDecayVtx_YE);
            mumukktree->Branch("MuMuDecayVtx_ZE",&MuMuDecayVtx_ZE);
            /// muons from JPsi (MuMu) fit
            mumukktree->Branch("mu1Idx",&mu1Idx);
            mumukktree->Branch("mu2Idx",&mu2Idx);
            mumukktree->Branch("mu1Px_MuMu",&muPos_MuMu_Px);
            mumukktree->Branch("mu1Py_MuMu",&muPos_MuMu_Py);
            mumukktree->Branch("mu1Pz_MuMu",&muPos_MuMu_Pz);
            mumukktree->Branch("mu1Chi2_MuMu",&muPos_MuMu_Chi2);
            mumukktree->Branch("mu1NDF_MuMu",&muPos_MuMu_NDF);
            mumukktree->Branch("mu2Px_MuMu",&muNeg_MuMu_Px);
            mumukktree->Branch("mu2Py_MuMu",&muNeg_MuMu_Py);
            mumukktree->Branch("mu2Pz_MuMu",&muNeg_MuMu_Pz);
            mumukktree->Branch("mu2Chi2_MuMu",&muNeg_MuMu_Chi2);
            mumukktree->Branch("mu2NDF_MuMu",&muNeg_MuMu_NDF);
            mumukktree->Branch("MuMuType",&MuMuType);
            mumukktree->Branch("MuMuMuonTrigMatch",&MuMuMuonTrigMatch);
            /// Primary Vertex with "MuMu correction"
            mumukktree->Branch("mumuLessPvs_n", &mumuLessPvs_n);
            mumukktree->Branch("PriVtxMuMuCorr_X", &PriVtxMuMuCorr_X);
            mumukktree->Branch("PriVtxMuMuCorr_Y", &PriVtxMuMuCorr_Y);
            mumukktree->Branch("PriVtxMuMuCorr_Z", &PriVtxMuMuCorr_Z);
            mumukktree->Branch("PriVtxMuMuCorr_EX", &PriVtxMuMuCorr_EX);
            mumukktree->Branch("PriVtxMuMuCorr_EY", &PriVtxMuMuCorr_EY);
            mumukktree->Branch("PriVtxMuMuCorr_EZ", &PriVtxMuMuCorr_EZ);
            mumukktree->Branch("PriVtxMuMuCorr_Chi2", &PriVtxMuMuCorr_Chi2);
            mumukktree->Branch("PriVtxMuMuCorr_CL", &PriVtxMuMuCorr_CL);
            mumukktree->Branch("PriVtxMuMuCorr_tracks", &PriVtxMuMuCorr_tracks);
            mumukktree->Branch("nTrk_afterMuMu", &nTrk);
            /// KK cand
            mumukktree->Branch("nKK",&nKK,"nKK/i");
            mumukktree->Branch("KKMass",&KKMass);
            mumukktree->Branch("KKPx",&KKPx);
            mumukktree->Branch("KKPy",&KKPy);
            mumukktree->Branch("KKPz",&KKPz);
            mumukktree->Branch("KKVtx_CL",&KKVtx_CL);
            mumukktree->Branch("KKVtx_Chi2",&KKVtx_Chi2);
            mumukktree->Branch("KKDecayVtx_X",&KKDecayVtx_X);
            mumukktree->Branch("KKDecayVtx_Y",&KKDecayVtx_Y);
            mumukktree->Branch("KKDecayVtx_Z",&KKDecayVtx_Z);
            mumukktree->Branch("KKDecayVtx_XE",&KKDecayVtx_XE);
            mumukktree->Branch("KKDecayVtx_YE",&KKDecayVtx_YE);
            mumukktree->Branch("KKDecayVtx_ZE",&KKDecayVtx_ZE);
            /// kaons from Phi (KK) fit
            mumukktree->Branch("ka1Idx",&ka1Idx);
            mumukktree->Branch("ka2Idx",&ka2Idx);
            mumukktree->Branch("ka1Px_KK",&kaonPos_KK_Px);
            mumukktree->Branch("ka1Py_KK",&kaonPos_KK_Py);
            mumukktree->Branch("ka1Pz_KK",&kaonPos_KK_Pz);
            mumukktree->Branch("ka1Chi2_KK",&kaonPos_KK_Chi2);
            mumukktree->Branch("ka1NDF_KK",&kaonPos_KK_NDF);
            mumukktree->Branch("ka2Px_KK",&kaonNeg_KK_Px);
            mumukktree->Branch("ka2Py_KK",&kaonNeg_KK_Py);
            mumukktree->Branch("ka2Pz_KK",&kaonNeg_KK_Pz);
            mumukktree->Branch("ka2Chi2_KK",&kaonNeg_KK_Chi2);
            mumukktree->Branch("ka2NDF_KK",&kaonNeg_KK_NDF);
            // mumukktree->Branch("DR_MuMu_K1",&DR_MuMu_K1);
            // mumukktree->Branch("DR_MuMu_K2",&DR_MuMu_K2);
            // mumukktree->Branch("DR_MuMuKK_K1",&DR_MuMuKK_K1);
            // mumukktree->Branch("DR_MuMuKK_K2",&DR_MuMuKK_K2);
            /// counters for X
            mumukktree->Branch("nX",&nX,"nX/i");
            mumukktree->Branch("nX_pre0",&nX_pre0,"nX_pre0/i");
            mumukktree->Branch("nX_pre1",&nX_pre1,"nX_pre1/i");
            mumukktree->Branch("nX_pre2",&nX_pre2,"nX_pre2/i");
            mumukktree->Branch("nX_pre3",&nX_pre3,"nX_pre3/i");
            mumukktree->Branch("nX_pre4",&nX_pre4,"nX_pre4/i");
            mumukktree->Branch("nX_pre5",&nX_pre5,"nX_pre5/i");
            mumukktree->Branch("nX_pre6",&nX_pre6,"nX_pre6/i");
            mumukktree->Branch("nX_pre7",&nX_pre7,"nX_pre7/i");
            mumukktree->Branch("nX_pre8",&nX_pre8,"nX_pre8/i");
            mumukktree->Branch("nX_pre9",&nX_pre9,"nX_pre9/i");
            mumukktree->Branch("nX_pre10",&nX_pre10,"nX_pre10/i");
            mumukktree->Branch("nX_pre11",&nX_pre11,"nX_pre11/i");
            mumukktree->Branch("nX_pre12",&nX_pre12,"nX_pre12/i");
            mumukktree->Branch("nX_pre13",&nX_pre13,"nX_pre13/i");
            mumukktree->Branch("nX_pre14",&nX_pre14,"nX_pre14/i");
            mumukktree->Branch("nX_pre15",&nX_pre15,"nX_pre15/i");
            /// B0 cand
            mumukktree->Branch("XMass",&xMass);
            mumukktree->Branch("XPx",&xPx);
            mumukktree->Branch("XPy",&xPy);
            mumukktree->Branch("XPz",&xPz);
            mumukktree->Branch("XPxE",&xPxE);
            mumukktree->Branch("XPyE",&xPyE);
            mumukktree->Branch("XPzE",&xPzE);
            mumukktree->Branch("XVtx_CL",&xVtx_CL);
            mumukktree->Branch("XVtx_Chi2",&xVtx_Chi2);
            mumukktree->Branch("XDecayVtx_X",&xDecayVtx_X);
            mumukktree->Branch("XDecayVtx_Y",&xDecayVtx_Y);
            mumukktree->Branch("XDecayVtx_Z",&xDecayVtx_Z);
            mumukktree->Branch("XDecayVtx_XE",&xDecayVtx_XE);
            mumukktree->Branch("XDecayVtx_YE",&xDecayVtx_YE);
            mumukktree->Branch("XDecayVtx_ZE",&xDecayVtx_ZE);
            mumukktree->Branch("XCosAlphaBS", &xCosAlphaBS);
            mumukktree->Branch("XCosAlpha3DBS", &xCosAlpha3DBS);
            mumukktree->Branch("XCTauBS", &xCTauBS);
            mumukktree->Branch("XCTauBSE", &xCTauBSE);
            mumukktree->Branch("XLxyBS", &xLxyBS);
            mumukktree->Branch("XLxyBSE", &xLxyBSE);
            mumukktree->Branch("XLxyzBS", &xLxyzBS);
            mumukktree->Branch("XLxyzBSE", &xLxyzBSE);
            mumukktree->Branch("XCosAlphaPV", &xCosAlphaPV);
            mumukktree->Branch("XCosAlpha3DPV", &xCosAlpha3DPV);
            mumukktree->Branch("XCTauPV", &xCTauPV);
            mumukktree->Branch("XCTauPVE", &xCTauPVE);
            mumukktree->Branch("XLxyPV", &xLxyPV);
            mumukktree->Branch("XLxyPVE", &xLxyPVE);
            mumukktree->Branch("XLxyzPV", &xLxyzPV);
            mumukktree->Branch("XLxyzPVE", &xLxyzPVE);
            /// Primary Vertex with largest B0_cos(alpha)
            mumukktree->Branch("PriVtx_XCosAlpha_n",&PriVtx_XCosAlpha_n);
            mumukktree->Branch("PriVtx_XCosAlpha_X",&PriVtx_XCosAlpha_X);
            mumukktree->Branch("PriVtx_XCosAlpha_Y",&PriVtx_XCosAlpha_Y);
            mumukktree->Branch("PriVtx_XCosAlpha_Z",&PriVtx_XCosAlpha_Z);
            mumukktree->Branch("PriVtx_XCosAlpha_EX",&PriVtx_XCosAlpha_EX);
            mumukktree->Branch("PriVtx_XCosAlpha_EY",&PriVtx_XCosAlpha_EY);
            mumukktree->Branch("PriVtx_XCosAlpha_EZ",&PriVtx_XCosAlpha_EZ);
            mumukktree->Branch("PriVtx_XCosAlpha_Chi2",&PriVtx_XCosAlpha_Chi2);
            mumukktree->Branch("PriVtx_XCosAlpha_CL",&PriVtx_XCosAlpha_CL);
            mumukktree->Branch("PriVtx_XCosAlpha_tracks",&PriVtx_XCosAlpha_tracks);
            mumukktree->Branch("XCosAlphaPVCosAlpha", &xCosAlphaPVCosAlpha);
            mumukktree->Branch("XCosAlpha3DPVCosAlpha", &xCosAlpha3DPVCosAlpha);
            mumukktree->Branch("XCTauPVCosAlpha", &xCTauPVCosAlpha);
            mumukktree->Branch("XCTauPVCosAlphaE", &xCTauPVCosAlphaE);
            mumukktree->Branch("XLxyPVCosAlpha", &xLxyPVCosAlpha);
            mumukktree->Branch("XLxyPVCosAlphaE", &xLxyPVCosAlphaE);
            mumukktree->Branch("XLxyzPVCosAlpha", &xLxyzPVCosAlpha);
            mumukktree->Branch("XLxyzPVCosAlphaE", &xLxyzPVCosAlphaE);
            mumukktree->Branch("PriVtx_XCosAlpha3D_n",&PriVtx_XCosAlpha3D_n);
            mumukktree->Branch("PriVtx_XCosAlpha3D_X",&PriVtx_XCosAlpha3D_X);
            mumukktree->Branch("PriVtx_XCosAlpha3D_Y",&PriVtx_XCosAlpha3D_Y);
            mumukktree->Branch("PriVtx_XCosAlpha3D_Z",&PriVtx_XCosAlpha3D_Z);
            mumukktree->Branch("PriVtx_XCosAlpha3D_EX",&PriVtx_XCosAlpha3D_EX);
            mumukktree->Branch("PriVtx_XCosAlpha3D_EY",&PriVtx_XCosAlpha3D_EY);
            mumukktree->Branch("PriVtx_XCosAlpha3D_EZ",&PriVtx_XCosAlpha3D_EZ);
            mumukktree->Branch("PriVtx_XCosAlpha3D_Chi2",&PriVtx_XCosAlpha3D_Chi2);
            mumukktree->Branch("PriVtx_XCosAlpha3D_CL",&PriVtx_XCosAlpha3D_CL);
            mumukktree->Branch("PriVtx_XCosAlpha3D_tracks",&PriVtx_XCosAlpha3D_tracks);
            mumukktree->Branch("XCosAlphaPVCosAlpha3D", &xCosAlphaPVCosAlpha3D);
            mumukktree->Branch("XCosAlpha3DPVCosAlpha3D", &xCosAlpha3DPVCosAlpha3D);
            mumukktree->Branch("XCTauPVCosAlpha3D", &xCTauPVCosAlpha3D);
            mumukktree->Branch("XCTauPVCosAlpha3DE", &xCTauPVCosAlpha3DE);
            mumukktree->Branch("XLxyPVCosAlpha3D", &xLxyPVCosAlpha3D);
            mumukktree->Branch("XLxyPVCosAlpha3DE", &xLxyPVCosAlpha3DE);
            mumukktree->Branch("XLxyzPVCosAlpha3D", &xLxyzPVCosAlpha3D);
            mumukktree->Branch("XLxyzPVCosAlpha3DE", &xLxyzPVCosAlpha3DE);

            mumukktree->Branch("XLessPV_tracksPtSq",&XLessPV_tracksPtSq);
            mumukktree->Branch("XLessPV_4tracksPtSq",&XLessPV_4tracksPtSq);
            mumukktree->Branch("PriVtxXLess_n",&PriVtxXLess_n);
            mumukktree->Branch("PriVtxXLess_X",&PriVtxXLess_X);
            mumukktree->Branch("PriVtxXLess_Y",&PriVtxXLess_Y);
            mumukktree->Branch("PriVtxXLess_Z",&PriVtxXLess_Z);
            mumukktree->Branch("PriVtxXLess_EX",&PriVtxXLess_EX);
            mumukktree->Branch("PriVtxXLess_EY",&PriVtxXLess_EY);
            mumukktree->Branch("PriVtxXLess_EZ",&PriVtxXLess_EZ);
            mumukktree->Branch("PriVtxXLess_Chi2",&PriVtxXLess_Chi2);
            mumukktree->Branch("PriVtxXLess_CL",&PriVtxXLess_CL);
            mumukktree->Branch("PriVtxXLess_tracks",&PriVtxXLess_tracks);
            mumukktree->Branch("XCosAlphaXLessPV", &xCosAlphaXLessPV);
            mumukktree->Branch("XCosAlpha3DXLessPV", &xCosAlpha3DXLessPV);
            mumukktree->Branch("XCTauXLessPV", &xCTauXLessPV);
            mumukktree->Branch("XCTauXLessPVE", &xCTauXLessPVE);
            mumukktree->Branch("XLxyXLessPV", &xLxyXLessPV);
            mumukktree->Branch("XLxyXLessPVE", &xLxyXLessPVE);
            mumukktree->Branch("XLxyzXLessPV", &xLxyzXLessPV);
            mumukktree->Branch("XLxyzXLessPVE", &xLxyzXLessPVE);
            mumukktree->Branch("PriVtxXLess_XCosAlpha_n",&PriVtxXLess_XCosAlpha_n);
            mumukktree->Branch("PriVtxXLess_XCosAlpha_X",&PriVtxXLess_XCosAlpha_X);
            mumukktree->Branch("PriVtxXLess_XCosAlpha_Y",&PriVtxXLess_XCosAlpha_Y);
            mumukktree->Branch("PriVtxXLess_XCosAlpha_Z",&PriVtxXLess_XCosAlpha_Z);
            mumukktree->Branch("PriVtxXLess_XCosAlpha_EX",&PriVtxXLess_XCosAlpha_EX);
            mumukktree->Branch("PriVtxXLess_XCosAlpha_EY",&PriVtxXLess_XCosAlpha_EY);
            mumukktree->Branch("PriVtxXLess_XCosAlpha_EZ",&PriVtxXLess_XCosAlpha_EZ);
            mumukktree->Branch("PriVtxXLess_XCosAlpha_Chi2",&PriVtxXLess_XCosAlpha_Chi2);
            mumukktree->Branch("PriVtxXLess_XCosAlpha_CL",&PriVtxXLess_XCosAlpha_CL);
            mumukktree->Branch("PriVtxXLess_XCosAlpha_tracks",&PriVtxXLess_XCosAlpha_tracks);
            mumukktree->Branch("XCosAlphaXLessPVCosAlpha", &xCosAlphaXLessPVCosAlpha);
            mumukktree->Branch("XCosAlpha3DXLessPVCosAlpha", &xCosAlpha3DXLessPVCosAlpha);
            mumukktree->Branch("XCTauXLessPVCosAlpha", &xCTauXLessPVCosAlpha);
            mumukktree->Branch("XCTauXLessPVCosAlphaE", &xCTauXLessPVCosAlphaE);
            mumukktree->Branch("XLxyXLessPVCosAlpha", &xLxyXLessPVCosAlpha);
            mumukktree->Branch("XLxyXLessPVCosAlphaE", &xLxyXLessPVCosAlphaE);
            mumukktree->Branch("XLxyzXLessPVCosAlpha", &xLxyzXLessPVCosAlpha);
            mumukktree->Branch("XLxyzXLessPVCosAlphaE", &xLxyzXLessPVCosAlphaE);
            mumukktree->Branch("PriVtxXLess_XCosAlpha3D_n",&PriVtxXLess_XCosAlpha3D_n);
            mumukktree->Branch("PriVtxXLess_XCosAlpha3D_X",&PriVtxXLess_XCosAlpha3D_X);
            mumukktree->Branch("PriVtxXLess_XCosAlpha3D_Y",&PriVtxXLess_XCosAlpha3D_Y);
            mumukktree->Branch("PriVtxXLess_XCosAlpha3D_Z",&PriVtxXLess_XCosAlpha3D_Z);
            mumukktree->Branch("PriVtxXLess_XCosAlpha3D_EX",&PriVtxXLess_XCosAlpha3D_EX);
            mumukktree->Branch("PriVtxXLess_XCosAlpha3D_EY",&PriVtxXLess_XCosAlpha3D_EY);
            mumukktree->Branch("PriVtxXLess_XCosAlpha3D_EZ",&PriVtxXLess_XCosAlpha3D_EZ);
            mumukktree->Branch("PriVtxXLess_XCosAlpha3D_Chi2",&PriVtxXLess_XCosAlpha3D_Chi2);
            mumukktree->Branch("PriVtxXLess_XCosAlpha3D_CL",&PriVtxXLess_XCosAlpha3D_CL);
            mumukktree->Branch("PriVtxXLess_XCosAlpha3D_tracks",&PriVtxXLess_XCosAlpha3D_tracks);
            mumukktree->Branch("XCosAlphaXLessPVCosAlpha3D", &xCosAlphaXLessPVCosAlpha3D);
            mumukktree->Branch("XCosAlpha3DXLessPVCosAlpha3D", &xCosAlpha3DXLessPVCosAlpha3D);
            mumukktree->Branch("XCTauXLessPVCosAlpha3D", &xCTauXLessPVCosAlpha3D);
            mumukktree->Branch("XCTauXLessPVCosAlpha3DE", &xCTauXLessPVCosAlpha3DE);
            mumukktree->Branch("XLxyXLessPVCosAlpha3D", &xLxyXLessPVCosAlpha3D);
            mumukktree->Branch("XLxyXLessPVCosAlpha3DE", &xLxyXLessPVCosAlpha3DE);
            mumukktree->Branch("XLxyzXLessPVCosAlpha3D", &xLxyzXLessPVCosAlpha3D);
            mumukktree->Branch("XLxyzXLessPVCosAlpha3DE", &xLxyzXLessPVCosAlpha3DE);
            /// Primary Vertex with "B0 correction"
            mumukktree->Branch("PriVtxXCorr_n",&PriVtxXCorr_n);
            mumukktree->Branch("PriVtxXCorr_X",&PriVtxXCorr_X);
            mumukktree->Branch("PriVtxXCorr_Y",&PriVtxXCorr_Y);
            mumukktree->Branch("PriVtxXCorr_Z",&PriVtxXCorr_Z);
            mumukktree->Branch("PriVtxXCorr_EX",&PriVtxXCorr_EX);
            mumukktree->Branch("PriVtxXCorr_EY",&PriVtxXCorr_EY);
            mumukktree->Branch("PriVtxXCorr_EZ",&PriVtxXCorr_EZ);
            mumukktree->Branch("PriVtxXCorr_Chi2",&PriVtxXCorr_Chi2);
            mumukktree->Branch("PriVtxXCorr_CL",&PriVtxXCorr_CL);
            mumukktree->Branch("PriVtxXCorr_tracks",&PriVtxXCorr_tracks);
            /// Lifetime variables for B0
            mumukktree->Branch("XCosAlphaPVX", &xCosAlphaPVX);
            mumukktree->Branch("XCTauPVX", &xCTauPVX);
            mumukktree->Branch("XCTauPVXE", &xCTauPVXE);
            mumukktree->Branch("XLxyPVX", &xLxyPVX);
            mumukktree->Branch("XLxyzPVX", &xLxyzPVX);
            mumukktree->Branch("XCTauPVX_3D", &xCTauPVX_3D);
            mumukktree->Branch("XCTauPVX_3D_err", &xCTauPVX_3D_err);
            /// dxy, dz, dxyE, dzE for kaons from PV, BS, B0LessPV
            mumukktree->Branch("kaonPos_dxy_PV", &kaonPos_dxy_PV);
            mumukktree->Branch("kaonPos_dz_PV", &kaonPos_dz_PV);
            mumukktree->Branch("kaonNeg_dxy_PV", &kaonNeg_dxy_PV);
            mumukktree->Branch("kaonNeg_dz_PV", &kaonNeg_dz_PV);
            mumukktree->Branch("kaonPos_dxy_BS", &kaonPos_dxy_BS);
            mumukktree->Branch("kaonPos_dz_BS", &kaonPos_dz_BS);
            mumukktree->Branch("kaonNeg_dxy_BS", &kaonNeg_dxy_BS);
            mumukktree->Branch("kaonNeg_dz_BS", &kaonNeg_dz_BS);
            mumukktree->Branch("kaonPos_dxy_XLessPV", &kaonPos_dxy_XLessPV);
            mumukktree->Branch("kaonPos_dz_XLessPV", &kaonPos_dz_XLessPV);
            mumukktree->Branch("kaonNeg_dxy_XLessPV", &kaonNeg_dxy_XLessPV);
            mumukktree->Branch("kaonNeg_dz_XLessPV", &kaonNeg_dz_XLessPV);
            mumukktree->Branch("kaonPos_dxyE", &kaonPos_dxyE);
            mumukktree->Branch("kaonPos_dzE", &kaonPos_dzE);
            mumukktree->Branch("kaonNeg_dxyE", &kaonNeg_dxyE);
            mumukktree->Branch("kaonNeg_dzE", &kaonNeg_dzE);

            mumukktree->Branch("XMuMuIdx", &X_MuMuIdx);
            mumukktree->Branch("XkaonPosIdx", &X_ka1Idx);
            mumukktree->Branch("XkaonNegIdx", &X_ka2Idx);

            mumukktree->Branch("kaonPosFromPV",&kaonPosFromPV);
            mumukktree->Branch("kaonNegFromPV",&kaonNegFromPV );

            /// Muons and tracks after X candidates fit
            mumukktree->Branch("Muon1Px_MuMuKK", &mu1Px_MuMuKK);
            mumukktree->Branch("Muon1Py_MuMuKK", &mu1Py_MuMuKK);
            mumukktree->Branch("Muon1Pz_MuMuKK", &mu1Pz_MuMuKK);
            mumukktree->Branch("Muon1E_MuMuKK", &mu1E_MuMuKK);
            mumukktree->Branch("Muon2Px_MuMuKK", &mu2Px_MuMuKK);
            mumukktree->Branch("Muon2Py_MuMuKK", &mu2Py_MuMuKK);
            mumukktree->Branch("Muon2Pz_MuMuKK", &mu2Pz_MuMuKK);
            mumukktree->Branch("Muon2E_MuMuKK", &mu2E_MuMuKK);
            mumukktree->Branch("kaonPosPx_MuMuKK", &k1Px_MuMuKK);
            mumukktree->Branch("kaonPosPy_MuMuKK", &k1Py_MuMuKK);
            mumukktree->Branch("kaonPosPz_MuMuKK", &k1Pz_MuMuKK);
            mumukktree->Branch("kaonPosE_MuMuKK", &k1E_MuMuKK);
            mumukktree->Branch("kaonPos_nsigdedx", &kaonPos_nsigdedx);
            mumukktree->Branch("kaonPos_dedx", &kaonPos_dedx);
            mumukktree->Branch("kaonPos_dedxMass", &kaonPos_dedxMass);
            mumukktree->Branch("kaonPos_theo", &kaonPos_theo);
            mumukktree->Branch("kaonPos_sigma", &kaonPos_sigma);
            mumukktree->Branch("kaonPos_dedx_byHits", &kaonPos_dedx_byHits);
            mumukktree->Branch("kaonPos_dedxErr_byHits", &kaonPos_dedxErr_byHits);
            mumukktree->Branch("kaonPos_saturMeas_byHits", &kaonPos_saturMeas_byHits);
            mumukktree->Branch("kaonPos_Meas_byHits", &kaonPos_Meas_byHits);
            mumukktree->Branch("kaonNegPx_MuMuKK", &k2Px_MuMuKK);
            mumukktree->Branch("kaonNegPy_MuMuKK", &k2Py_MuMuKK);
            mumukktree->Branch("kaonNegPz_MuMuKK", &k2Pz_MuMuKK);
            mumukktree->Branch("kaonNegE_MuMuKK", &k2E_MuMuKK);
            mumukktree->Branch("kaonNeg_nsigdedx", &kaonNeg_nsigdedx);
            mumukktree->Branch("kaonNeg_dedx", &kaonNeg_dedx);
            mumukktree->Branch("kaonNeg_dedxMass", &kaonNeg_dedxMass);
            mumukktree->Branch("kaonNeg_theo", &kaonNeg_theo);
            mumukktree->Branch("kaonNeg_sigma", &kaonNeg_sigma);
            mumukktree->Branch("kaonNeg_dedx_byHits", &kaonNeg_dedx_byHits);
            mumukktree->Branch("kaonNeg_dedxErr_byHits", &kaonNeg_dedxErr_byHits);
            mumukktree->Branch("kaonNeg_saturMeas_byHits", &kaonNeg_saturMeas_byHits);
            mumukktree->Branch("kaonNeg_Meas_byHits", &kaonNeg_Meas_byHits);

          }/// begin Job

          /// ------------ method called once each job just after ending the event loop  ------------
          void MuMuKKPAT::endJob() {
            mumukktree->GetDirectory()->cd();
            mumukktree->Write();
          }/// endjob


          bool MuMuKKPAT::isAbHadron(int pdgID) {

            if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
            return false;

          }

          bool MuMuKKPAT::isAMixedbHadron(int pdgID, int momPdgID) {

            if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) ||
            (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0))
            return true;
            return false;

          }

          std::pair<int, float> MuMuKKPAT::findCandMCInfo(reco::GenParticleRef genCand) {

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

          double MuMuKKPAT::getSigmaOfLogdEdx(double logde)
          {
            return 0.3;
          }

          float MuMuKKPAT::getEnergyLoss(const reco::TrackRef & track)
          {
            if (iexception_dedx==1) return 9999.;
            const reco::DeDxDataValueMap & eloss = *energyLoss;
            return eloss[track].dEdx();
          }

          double MuMuKKPAT::nsigmaofdedx(const reco::TrackRef & track, double & theo, double & sigma)
          {

            // no usable dE/dx if p > 2
            double nsigma = 99 ;
            if (iexception_dedx==1) return nsigma ;

            double m  = 0.13957;
            double bg = track->p() / m;

            theo = getLogdEdx(bg);


            int nhitr = track->numberOfValidHits();
            double meas = log(getEnergyLoss(track));
            sigma = getSigmaOfLogdEdx(theo) * pow(nhitr,-0.65);
            if (sigma>0)
            nsigma = (meas-theo) / sigma ;
            return nsigma;
          }


          double MuMuKKPAT::getLogdEdx(double bg)
          {
            const double a =  3.25 ;
            const double b =  0.288;
            const double c = -0.852;

            double beta = bg/sqrt(bg*bg + 1);
            double dedx = log( a/(beta*beta) + b * log(bg) + c );

            return dedx;

          }


          double MuMuKKPAT::GetMass(const reco::TrackRef & track){
            double P = track->p();
            double C = 2.625;
            double K = 2.495;
            double I = getEnergyLoss(track);
            return sqrt((I-C)/K)*P;
          }


          template<typename T>
          bool MuMuKKPAT::isBetterMuon(const T &mu1, const T &mu2) const {
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

          bool MuMuKKPAT::isSameMuon(const reco::Muon &mu1, const reco::Muon &mu2) const {
            return (& mu1 == & mu2) ||
            //(mu1.originalObjectRef() == mu2.originalObjectRef()) ||
            (mu1.reco::Muon::innerTrack().isNonnull() ?
            mu1.reco::Muon::innerTrack() == mu2.reco::Muon::innerTrack() :
            mu1.reco::Muon::outerTrack() == mu2.reco::Muon::outerTrack());
          }


          /// define this as a plug-in
          DEFINE_FWK_MODULE(MuMuKKPAT);

          // rsync -vut --existing src/MuMuKKPAT.cc semrat@lxplus.cern.ch:/afs/cern.ch/user/s/semrat/scratch0/CMSSW_5_3_22/src/X4140/MuMuKKPAT/src/MuMuKKPAT.cc
