// -*- C++ -*-
//
// Package:    DiMuonRootupler
// Class:      DiMuonRootupler
//
// Description: Dimuon(mu+ mu-)  producer
//
// Author:  Adriano Di Florio
//    based on : Alberto Sanchez Hernandez Onia2MuMu code
//

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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

//
// class declaration
//

class DiMuonRootupler:public edm::EDAnalyzer {
      public:
	explicit DiMuonRootupler(const edm::ParameterSet &);
	~DiMuonRootupler() override;

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

      private:
        UInt_t getTriggerBits(const edm::Event &);
        bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
        const  reco::Candidate* GetAncestor(const reco::Candidate *);

	void beginJob() override;
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	void endJob() override;

	void beginRun(edm::Run const &, edm::EventSetup const &) override;
	void endRun(edm::Run const &, edm::EventSetup const &) override;
	void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;
	void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;

	// ----------member data ---------------------------
	std::string file_name;
	edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  int  pdgid_;
  std::vector<double> DimuonMassCuts_;
	bool isMC_;
  bool OnlyBest_;
  bool OnlyGen_;
  std::vector<std::string>  HLTs_;

	UInt_t    run;
	ULong64_t event;
  UInt_t    lumiblock;
  UInt_t    ndimuon;
  UInt_t    nmuons;
  UInt_t    trigger;
  UInt_t    tMatch;
  Int_t     charge;

	TLorentzVector dimuon_p4;
	TLorentzVector highMuon_p4;
	TLorentzVector lowMuon_p4;

  Float_t MassErr;
  Float_t vProb;
  Float_t DCA;
  Float_t ppdlPV;
  Float_t ppdlErrPV;
  Float_t ppdlBS;
  Float_t ppdlErrBS;
  Float_t cosAlpha;
  Float_t lxyPV;
  Float_t lxyBS;

  Float_t dimuon_m, dimuon_pt, dimuon_eta, dimuon_phi, dimuon_p;

	UInt_t numPrimaryVertices;

	TTree *dimuon_tree;

  Int_t mother_pdgId;
  Int_t dimuon_pdgId;
  TLorentzVector gen_mother_p4;
	TLorentzVector gen_dimuon_p4;
	TLorentzVector gen_highMuon_p4;
	TLorentzVector gen_muonM_p4;

  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;

};

//
// constructors and destructor
//

DiMuonRootupler::DiMuonRootupler(const edm::ParameterSet & iConfig):
dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuons"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
pdgid_(iConfig.getParameter<uint32_t>("dimuon_pdgid")),
DimuonMassCuts_(iConfig.getParameter<std::vector<double>>("dimuon_mass_cuts")),
isMC_(iConfig.getParameter<bool>("isMC")),
OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs"))
{
  edm::Service < TFileService > fs;
  dimuon_tree = fs->make < TTree > ("dimuonTree", "Tree of DiMuon");

  dimuon_tree->Branch("run",      &run,      "run/i");
  dimuon_tree->Branch("event",    &event,    "event/l");
  dimuon_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  if (!OnlyGen_) {
    dimuon_tree->Branch("ndimuon",    &ndimuon,    "ndimuon/i");
    dimuon_tree->Branch("nmuons",   &nmuons,   "nmuons/i");
    dimuon_tree->Branch("trigger",  &trigger,  "trigger/i");
    dimuon_tree->Branch("charge",   &charge,   "charge/I");

    dimuon_tree->Branch("tMatch",    &tMatch,      "tMatch/I");

    dimuon_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
    dimuon_tree->Branch("highMuon_p4",  "TLorentzVector", &highMuon_p4);
    dimuon_tree->Branch("lowMuon_p4",  "TLorentzVector", &lowMuon_p4);

    dimuon_tree->Branch("dimuon_m",      &dimuon_m,    "dimuon_m/F");
    dimuon_tree->Branch("dimuon_pt",     &dimuon_pt,   "dimuon_pt/F");
    dimuon_tree->Branch("dimuon_eta",    &dimuon_eta,  "dimuon_eta/F");
    dimuon_tree->Branch("dimuon_phi",    &dimuon_phi,  "dimuon_phi/F");
    dimuon_tree->Branch("dimuon_p",     &dimuon_p,     "dimuon_p/F");

    dimuon_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
    dimuon_tree->Branch("vProb",     &vProb,      "vProb/F");
    dimuon_tree->Branch("DCA",       &DCA,        "DCA/F");
    dimuon_tree->Branch("ppdlPV",    &ppdlPV,     "ppdlPV/F");
    dimuon_tree->Branch("ppdlErrPV", &ppdlErrPV,  "ppdlErrPV/F");
    dimuon_tree->Branch("ppdlBS",    &ppdlBS,     "ppdlBS/F");
    dimuon_tree->Branch("ppdlErrBS", &ppdlErrBS,  "ppdlErrBS/F");
    dimuon_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");
    dimuon_tree->Branch("lxyPV",     &lxyPV,      "lxyPV/F");
    dimuon_tree->Branch("lxyBS",     &lxyBS,      "lxyBS/F");

    dimuon_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");
  }

  if (isMC_ || OnlyGen_) {
     // std::cout << "DiMuonRootupler::DiMuonRootupler: Dimuon id " << pdgid_ << std::endl;
     dimuon_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
     dimuon_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
     dimuon_tree->Branch("gen_mother_p4", "TLorentzVector",  &gen_mother_p4);
     dimuon_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
     dimuon_tree->Branch("gen_highMuon_p4",  "TLorentzVector",  &gen_highMuon_p4);
     dimuon_tree->Branch("gen_lowMuon_p4",  "TLorentzVector",  &gen_muonM_p4);
  }

  genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
  packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
}

DiMuonRootupler::~DiMuonRootupler() {}

//
// member functions
//

const reco::Candidate* DiMuonRootupler::GetAncestor(const reco::Candidate* p) {
   if (p->numberOfMothers()) {
      if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
      else return p->mother(0);
   }
   return p;
}

//Check recursively if any ancestor of particle is the given one
bool DiMuonRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}

/* Grab Trigger information. Save it in variable trigger, trigger is an uint between 0 and 256, in binary it is:
   (pass 2)(pass 1)(pass 0)
   ex. 7 = pass 0, 1 and 2
   ex. 6 = pass 1, 2
   ex. 1 = pass 0
*/

UInt_t DiMuonRootupler::getTriggerBits(const edm::Event& iEvent ) {

  UInt_t trigger = 0;

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  if (triggerResults_handle.isValid()) {
     const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
     unsigned int NTRIGGERS = HLTs_.size();

     for (unsigned int i = 0; i < NTRIGGERS; i++) {
        for (int version = 1; version < 20; version++) {
           std::stringstream ss;
           ss << HLTs_[i] << "_v" << version;
           unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
           if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
              trigger += (1<<i);
              break;
           }
        }
     }
   } else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;

   return trigger;
}

// ------------ method called for each event  ------------
void DiMuonRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  iEvent.getByToken(dimuon_Label,dimuons);

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();
  trigger = getTriggerBits(iEvent);

  dimuon_pdgId = 0;
  mother_pdgId = 0;
  ndimuon  = 0;
  nmuons = 0;

  dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  lowMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_mother_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muonM_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  // Pruned particles are the one containing "important" stuff
  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genCands_, pruned);

  // Packed particles are all the status 1. The navigation to pruned is possible (the other direction should be made by hand)
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packCands_,  packed);

  if ( (isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid() ) {
    for (size_t i=0; i<pruned->size(); i++) {
      const reco::Candidate *adimuon = &(*pruned)[i];
      if ( (abs(adimuon->pdgId()) == pdgid_) && (adimuon->status() == 2) ) {
        int foundit = 1;
        dimuon_pdgId = adimuon->pdgId();
        for ( size_t j=0; j<packed->size(); j++ ) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
          const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
          const reco::Candidate * d = &(*packed)[j];
          if ( motherInPrunedCollection != nullptr && (d->pdgId() ==  13 ) && isAncestor(adimuon , motherInPrunedCollection) ) {
            gen_muonM_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
            foundit++;
          }
          if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(adimuon , motherInPrunedCollection) ) {
            gen_highMuon_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
            foundit++;
          }
          if ( foundit == 3 ) break;
        }
        if ( foundit == 3 ) {
          gen_dimuon_p4 = gen_muonM_p4 + gen_highMuon_p4;   // this should take into account FSR
          mother_pdgId  = GetAncestor(adimuon)->pdgId();
          auto mother = GetAncestor(adimuon);
          gen_mother_p4.SetPtEtaPhiM(mother->pt(),mother->eta(),mother->phi(),mother->mass());
          break;
        } else dimuon_pdgId = 0;
      }  // if ( p_id
    } // for (size
    // if ( dimuon_pdgId ) std::cout << "DiMuonRootupler: found the given decay " << run << "," << event << std::endl; // sanity check
  }  // end if isMC

  float DimuonMassMax_ = DimuonMassCuts_[1];
  float DimuonMassMin_ = DimuonMassCuts_[0];

  bool already_stored = false;
  if ( ! OnlyGen_ ) { // we will look for dimuons, then for muons
    if ( dimuons.isValid() && !dimuons->empty()) {
      for ( pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuons->begin(); dimuonCand != dimuons->end(); ++dimuonCand ) {
        if (dimuonCand->mass() > DimuonMassMin_ && dimuonCand->mass() < DimuonMassMax_ && dimuonCand->charge() == 0) {
          dimuon_p4.SetPtEtaPhiM(dimuonCand->pt(),dimuonCand->eta(),dimuonCand->phi(),dimuonCand->mass());

          dimuon_m    = dimuonCand->mass();
          dimuon_pt   = dimuonCand->pt();
          dimuon_eta  = dimuonCand->eta();
          dimuon_phi  = dimuonCand->phi();
          dimuon_p    = dimuonCand->p();

          reco::Candidate::LorentzVector vP = dimuonCand->daughter("highMuon")->p4();
          reco::Candidate::LorentzVector vM = dimuonCand->daughter("lowMuon")->p4();
          if ( dimuonCand->daughter("highMuon")->pt() < dimuonCand->daughter("lowMuon")->pt() ) {
              vP = dimuonCand->daughter("lowMuon")->p4();
              vM = dimuonCand->daughter("highMuon")->p4();
          }
          highMuon_p4.SetPtEtaPhiM(vP.pt(),vP.eta(),vP.phi(),vP.mass());
          lowMuon_p4.SetPtEtaPhiM(vM.pt(),vM.eta(),vM.phi(),vM.mass());
          MassErr = -1.0;
          if (dimuonCand->hasUserFloat("MassErr"))    MassErr = dimuonCand->userFloat("MassErr");
          vProb = dimuonCand->userFloat("vProb");
          DCA = -1.;
          if (dimuonCand->hasUserFloat("DCA"))  DCA = dimuonCand->userFloat("DCA");
          ppdlPV = dimuonCand->userFloat("ppdlPV");
          ppdlErrPV = dimuonCand->userFloat("ppdlErrPV");
          ppdlBS = dimuonCand->userFloat("ppdlBS");
          ppdlErrBS = dimuonCand->userFloat("ppdlErrBS");
          cosAlpha = dimuonCand->userFloat("cosAlpha");
          tMatch = dimuonCand->userInt("isTriggerMatched");
          charge = dimuonCand->charge();
          TVector3 pperp(dimuonCand->px(),dimuonCand->py(),0);
          lxyPV = ppdlPV * pperp.Perp() / dimuonCand->mass();
          lxyBS = ppdlBS * pperp.Perp() / dimuonCand->mass();
          ndimuon++;
          if (OnlyBest_) break;
          else {
            dimuon_tree->Fill();   // be aware, we are storing all combinations
            already_stored = true;
          }
        }
      }
    } //..else {
      //std::cout << "DiMuonRootupler: (" << run << "," << event << ") -> ";

  }  // !OnlyGen_

  if ( !already_stored ) {  // we have to make sure, we are not double storing an combination
    if ( !isMC_ ) {
      if ( ndimuon > 0 ) dimuon_tree->Fill();   // if not MC filter out
    } else dimuon_tree->Fill();
  }
}

// ------------ method called once each job just before starting event loop  ------------
void DiMuonRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DiMuonRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DiMuonRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void DiMuonRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DiMuonRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DiMuonRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DiMuonRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonRootupler);
