#ifndef JpsiPhiAnalysis_DiTrak_h
#define JpsiPhiAnalysis_DiTrak_h


// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/Utils/interface/PtComparator.h"

// DataFormat includes
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"



//
// class decleration
//

class DiTrakPAT : public edm::EDProducer {
 public:
  explicit DiTrakPAT(const edm::ParameterSet&);
  ~DiTrakPAT() override;

 private:
  void beginJob() override ;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endJob() override ;


  // ----------member data ---------------------------
 private:

  edm::EDGetTokenT<std::vector<pat::PackedCandidate>> traks_;
  std::vector<double> ditrakMassCuts_;
  std::vector<double> massTraks_;

};

#endif
