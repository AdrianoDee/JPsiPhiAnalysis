/**
   \file
   Declaration of DiMuonDiTrakProducer

   \author Alberto Sanchez-Hernandez
   \date 2 Mar 2014
*/

#ifndef __DiMuonDiTrakProducer_h_
#define __DiMuonDiTrakProducer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include <TLorentzVector.h>
#include <vector>
#include <tuple>

#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "../interface/DiMuonVtxReProducer.h"

#include "CommonTools/UtilAlgos/interface/PhysObjectMatcher.h"
#include "CommonTools/UtilAlgos/interface/MCMatchSelector.h"
#include "CommonTools/UtilAlgos/interface/MatchByDRDPt.h"
#include "CommonTools/UtilAlgos/interface/MatchLessByDPt.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

template<typename T>
struct GreaterByVProb {
  bool operator()( const T & t1, const T & t2 ) const {
    return t1.userFloat("vProb") > t2.userFloat("vProb");
  }
};

/**
   Create a HF candidate by mathing DiMuon(chi,psi,etc.) and a track (K, pi, etc.)
 */

class DiMuonDiTrakProducer : public edm::EDProducer {

 public:
  explicit DiMuonDiTrakProducer(const edm::ParameterSet& ps);

 private:

  void produce(edm::Event& event, const edm::EventSetup& esetup) override;

  void endJob() override;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> DiMuonCollection_;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> TrakCollection_;
  double TrackPtCut_,MaxDeltaR_, MaxDeltaR_;
  edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> TrackGenMap_;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> TriggerCollection_;
  edm::EDGetTokenT<edm::TriggerResults> TriggerResults_;
  std::vector<double> DiMuonMassCuts_;
  std::vector<double> TrakTrakMassCuts_;
  std::vector<double> DiMuonDiTrakMassCuts_;
  std::vector<double> MassTraks_;
  double JPsiMass_,PhiMass_;
  bool OnlyBest_;
  std::string product_name_;
  std::vector<std::string>  HLTFilters_;
  bool IsMC_;
  bool AddSameSig_;
  bool doPionRefit_;

  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
  edm::EDGetTokenT<pat::MuonCollection> allMuons_;
  edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> genMap_;
  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  bool IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu);
  pat::CompositeCandidate makeDiMuonTTCandidate(const pat::CompositeCandidate& DiMuon,
						    const pat::CompositeCandidate& tt);
  pat::CompositeCandidate makeTTCandidate(const pat::PackedCandidate& track1,
                                                const pat::PackedCandidate& track2);
  const pat::CompositeCandidate makeTTTriggerCandidate(
                                            const pat::TriggerObjectStandAlone& trackP,
                                            const pat::TriggerObjectStandAlone& trackN
                                          );
  const pat::CompositeCandidate makeTTTriggerMixedCandidate(
                                            const pat::PackedCandidate& trackP,
                                            const pat::TriggerObjectStandAlone& trackN
                                          );
  std::tuple<int, float, float> findJpsiMCInfo(reco::GenParticleRef genParticle);
  bool isSameTrack(reco::Track t1, reco::Track t2);
  bool isAbHadron(int pdgID);
  bool isAMixedbHadron(int pdgID, int momPdgID);
  bool isTheCandidate(reco::GenParticleRef genY);

  bool MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  float DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  float DeltaPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);

  float maxDeltaR;
  float maxDPtRel;

  int candidates;
  int nevents;
  int ndimuon;
  int nreco;
};

#endif // __DiMuonDiTrakProducer_h_
