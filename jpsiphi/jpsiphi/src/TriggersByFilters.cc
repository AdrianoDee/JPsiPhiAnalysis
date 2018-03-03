// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

//
// class declaration
//

class TriggersByFilters:public edm::EDAnalyzer {
      public:
	explicit TriggersByFilters(const edm::ParameterSet &);
	~TriggersByFilters() override;

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

  private:

  	void beginJob() override;
  	void analyze(const edm::Event &, const edm::EventSetup &) override;
  	void endJob() override;

  	void beginRun(edm::Run const &, edm::EventSetup const &) override;
  	void endRun(edm::Run const &, edm::EventSetup const &) override;
  	void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;
  	void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;

	// ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggers_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults;
  std::vector<std::string>  HLTFilters_;

	UInt_t    run;
	ULong64_t event;
  UInt_t    lumiblock;

  UInt_t    ntrigs;
  UInt_t    trigger;
  UInt_t    filter;

  UInt_t charge;
  TLorentzVector trig_p4;
	UInt_t numPrimaryVertices;

	TTree *trig_tree;




};

//
// constructors and destructor
//

TriggersByFilters::TriggersByFilters(const edm::ParameterSet & iConfig):
triggers_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
triggerResults(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters"))
{
  edm::Service < TFileService > fs;
  trig_tree = fs->make < TTree > ("TriggerTree", "Tree of Trigger");

  trig_tree->Branch("run",      &run,      "run/i");
  trig_tree->Branch("event",    &event,    "event/l");
  trig_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  trig_tree->Branch("ntrigs",   &ntrigs,   "ntrigs/I");
  trig_tree->Branch("trigger",  &trigger,  "trigger/I");
  trig_tree->Branch("filter",   &filter,   "filter/I");

  trig_tree->Branch("charge",   &charge,   "charge/I");
  trig_tree->Branch("trig_p4", "TLorentzVector", &trig_p4);
  trig_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");

}

TriggersByFilters::~TriggersByFilters() {}

UInt_t TriggersByFilters::getTriggerBits(const edm::Event& iEvent, const edm::Handle< edm::TriggerResults >& triggerResults_handle) {

  UInt_t trigger = 0;
  const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

     unsigned int NTRIGGERS = HLTs_.size();

     for (unsigned int i = 0; i < NTRIGGERS; i++) {
        for (int version = 1; version < 20; version++) {
           std::stringstream ss;
           ss << HLTs_[i] << "_v" << version;
           unsigned int bit = names.triggerIndex(edm::InputTag(ss.str()).label());
           if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
              trigger += (1<<i);
              break;
           }
        }
     }

   return trigger;
}

// ------------ method called for each event  ------------
void TriggersByFilters::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trigs;
  iEvent.getByToken(triggers_,trigs);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults , triggerResults_handle);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  trig_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  trigger = 0;

  if (triggerResults_handle.isValid())
    trigger = getTriggerBits(iEvent,triggerResults_handle);
  else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;

  ntrigs = trigs->size();

  for ( size_t iTrigObj = 0; iTrigObj < trigs->size(); ++iTrigObj ) {

    pat::TriggerObjectStandAlone unPackedTrigger( trigs->at( iTrigObj ) );

    if(unPackedTrigger.charge()==0) continue

    unPackedTrigger.unpackPathNames( names );
    unPackedTrigger.unpackFilterLabels(iEvent,*triggerResults_handle);

    bool filtered = false;
    UInt_t thisFilter = 0;

    for (size_t i = 0; i < HLTFilters_.size(); i++)
    {
      if(unPackedTrigger.hasFilterLabel(HLTFilters_[i]))
        {
          thisFilter += (1<<i);
          filtered = true;
        }
    }

    if(!filtered) continue;

    filter = thisFilter;
    trig_p4.SetPtEtaPhiM(unPackedTrigger.pt(),unPackedTrigger.eta(),unPackedTrigger.phi(),unPackedTrigger.mass());

    trig_tree->Fill();

  }

}

// ------------ method called once each job just before starting event loop  ------------
void TriggersByFilters::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void TriggersByFilters::endJob() {}

// ------------ method called when starting to processes a run  ------------
void TriggersByFilters::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void TriggersByFilters::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void TriggersByFilters::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void TriggersByFilters::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TriggersByFilters::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggersByFilters);
