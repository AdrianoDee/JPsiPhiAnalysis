/**
   \file
   Declaration of FiveTracksProducer

   \author Alberto Sanchez-Hernandez
   \date 2 Mar 2014
*/

#ifndef __DiMuonDiTrakProducerFit_h_
#define __DiMuonDiTrakProducerFit_h_

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

#include <TLorentzVector.h>
#include <vector>
#include <tuple>

#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "../interface/DiMuonVtxReProducer.h"

/**
   Create a HF candidate by mathing DiMuon(chi,psi,etc.) and a track (K, pi, etc.)
 */

class FiveTracksProducer : public edm::EDProducer {

 public:
  explicit FiveTracksProducer(const edm::ParameterSet& ps);

 private:

  void produce(edm::Event& event, const edm::EventSetup& esetup) override;

  void endJob() override;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> DiMuonDiTrakCollection_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> TrakCollection_;
  double trakPtCut_;
  edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> TriggerCollection_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  std::vector<double> FiveTrakMassCuts_;
  bool OnlyBest_;
  std::vector<std::string>  HLTFilters_;
  bool isMC_;
  bool addMCTruth_;
  bool doDoubleConstant_;
  bool addSameSig_;

  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  bool IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu);
  pat::CompositeCandidate makeDiMuonTTCandidate(const pat::CompositeCandidate& DiMuon,
						    const pat::CompositeCandidate& tt);
  pat::CompositeCandidate makeFiveCandidate(const pat::PackedCandidate& trak1,
                                                const pat::PackedCandidate& trak2);
  pat::CompositeCandidate makePsi2SCandidate(const pat::CompositeCandidate& dimuon,
                                             const pat::CompositeCandidate& t1,
                                             const pat::CompositeCandidate& t2
                                           );
pat::CompositeCandidate makeFiveCandidateMixed(
                                            const pat::CompositeCandidate& dimuon,
                                            const pat::PackedCandidate& trakP,
                                            const pat::PackedCandidate& trakN,
                                            const pat::PackedCandidate& trak3,
                                            double massOne,
                                            double massTwo,
                                            double massThree
                                          );
pat::CompositeCandidate makeFiveCandidateMixed(
                                              const pat::CompositeCandidate& dimuon,
                                              const pat::CompositeCandidate& trakP,
                                              const pat::CompositeCandidate& trakN,
                                              const pat::CompositeCandidate& trak3
                                            );
  pat::CompositeCandidate makeFiveCandidate(
                                            const pat::CompositeCandidate& dimuonditrak,
                                            const pat::PackedCandidate& trakN
                                          );

  std::tuple<int, float, float> findJpsiMCInfo(reco::GenParticleRef genParticle);
  bool isSameTrack(reco::Track t1, reco::Track t2);
  bool IsTheSame(const pat::PackedCandidate& t1, const pat::PackedCandidate& t2);
  bool isAbHadron(int pdgID);
  bool isAMixedbHadron(int pdgID, int momPdgID);
  bool isTheCandidate(reco::GenParticleRef genY);

  bool MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  float DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);

  float maxDeltaR;
  float maxDPtRel;
  double kaonmass, pionmass, trackmass, psi2smass;

  int nevents;
  uint ncombo,ncomboneg,ncomboneu;
};

#endif // __DiMuonDiTrakProducerFit_h_
