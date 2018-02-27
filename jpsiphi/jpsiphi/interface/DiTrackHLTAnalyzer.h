/**
   \file
   Declaration of DiTrackHLTAnalyzer

   \author Alberto Sanchez-Hernandez
   \date 2 Mar 2014
*/

#ifndef __DiTrackHLTAnalyzer_h_
#define __DiTrackHLTAnalyzer_h_

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "CommonTools/UtilAlgos/interface/DeltaR.h"
#include "CommonTools/UtilAlgos/interface/MatchByDRDPt.h"
#include <TLorentzVector.h>
#include <vector>
#include <cmath>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

/**
   Create a HF candidate by mathing DiMuon(chi,psi,etc.) and a track (K, pi, etc.)
 */

class DiTrackHLTAnalyzer : public edm::EDAnalyzer {

 public:
  explicit DiTrackHLTAnalyzer(const edm::ParameterSet& ps);
  ~DiTrackHLTAnalyzer() override;
  static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

 private:

  void beginJob() override;
  void analyze(edm::Event& event, const edm::EventSetup& esetup) override;
  void endJob() override;
  void beginRun(edm::Run const &, edm::EventSetup const &) override;
  void endRun(edm::Run const &, edm::EventSetup const &) override;
  void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;
  void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;
  UInt_t getTriggerBits(const edm::Event &);
  UInt_t isTriggerMatched(const pat::CompositeCandidate *diTrig_Candidate);

  std::string file_name;

  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> TrakCollection_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> TriggerCollection_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  std::vector<double> TrakTrakMassCuts_;
  std::vector<double> MassTraks_;
  bool OnlyBest_;
  std::string TTTrigger_name_; //TTCandidate_name_
  std::vector<std::string>  HLTs_;
  std::vector <std::string> HLTFilters_;


  reco::Candidate::LorentzVector convertVector(const math::XYZTLorentzVectorF& v);
  bool IsTheSame(const pat::PackedCandidate& t1, const pat::PackedCandidate& t2);
  const pat::CompositeCandidate makeTTTriggerCandidate(const pat::TriggerObjectStandAlone& t1,
						    const pat::TriggerObjectStandAlone& t2);
  const pat::CompositeCandidate makeTTCandidate(const pat::PackedCandidate& t1,
                                                const pat::PackedCandidate& t2);

  bool MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);
  float DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2);

  int candidates;
  int nevents;
  int ndimuon;
  int nreco;
  float maxDeltaR;
  float maxDPtRel;

  UInt_t    run;
  ULong64_t event;
  UInt_t    lumiblock;

  UInt_t    trigger;
  UInt_t    tMatch;

  UInt_t charge;

  TLorentzVector ditrak_p4;
  TLorentzVector trakP_p4;
  TLorentzVector trakN_p4;

  TLorentzVector ditrig_p4;
  TLorentzVector trigP_p4;
  TLorentzVector trigN_p4;

  UInt_t numPrimaryVertices;

  TTree *ditrak_tree;

  UInt_t nditrak, ntraks;

};

#endif // __DiTrackHLTAnalyzer_h_
