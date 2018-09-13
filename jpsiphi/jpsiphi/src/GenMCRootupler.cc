// // -*- C++ -*-
// //
// // Package:    GenMCRootupler
// // Class:      GenMCRootupler
// //
// // Description: Dimuon(mu+ mu-)  producer
// //
// // Author:  Adriano Di Florio
// //    based on : Alberto Sanchez Hernandez Onia2MuMu code
// //
//
// // system include files
// #include <memory>
//
// // user include files
// #include "FWCore/Framework/Dnterface/Frameworkfwd.h"
// #include "FWCore/Framework/Dnterface/EDAnalyzer.h"
// #include "FWCore/Framework/Dnterface/Event.h"
// #include "FWCore/Framework/Dnterface/MakerMacros.h"
// #include "FWCore/ParameterSet/Dnterface/ParameterSet.h"
// #include "FWCore/ServiceRegistry/Dnterface/Service.h"
//
// #include "DataFormats/PatCandidates/Dnterface/CompositeCandidate.h"
// #include "DataFormats/PatCandidates/Dnterface/PackedGenParticle.h"
// #include "DataFormats/PatCandidates/Dnterface/Muon.h"
// #include "DataFormats/Candidate/Dnterface/Candidate.h"
// #include "DataFormats/HepMCCandidate/Dnterface/GenParticleFwd.h"
// #include "DataFormats/VertexReco/Dnterface/VertexFwd.h"
//
// #include "DataFormats/Common/Dnterface/TriggerResults.h"
// #include "FWCore/Common/Dnterface/TriggerNames.h"
//
// #include "CommonTools/UtilAlgos/Dnterface/TFileService.h"
// #include "TLorentzVector.h"
// #include "TTree.h"
//
// //
// // class declaration
// //
//
// class GenMCRootupler:public edm::EDAnalyzer {
//       public:
// 	explicit GenMCRootupler(const edm::ParameterSet &);
// 	~GenMCRootupler() override;
//
// 	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
//
//       private:
//         UInt_t getTriggerBits(const edm::Event &);
//         bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
//         const  reco::Candidate* GetAncestor(const reco::Candidate *);
//
// 	void beginJob() override;
// 	void analyze(const edm::Event &, const edm::EventSetup &) override;
// 	void endJob() override;
//
// 	void beginRun(edm::Run const &, edm::EventSetup const &) override;
// 	void endRun(edm::Run const &, edm::EventSetup const &) override;
// 	void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;
// 	void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;
//
// 	// ----------member data ---------------------------
// 	std::string file_name;
// 	edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
//   edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
//   edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
//   std::vector<int>  PdgIds_;
//   std::vector<double> DimuonMassCuts_;
// 	bool isMC_;
//   bool OnlyBest_;
//   bool OnlyGen_;
//   std::vector<std::string>  HLTs_;
//
// 	Double_t    run;
// 	ULong64_t event;
//   UInt_t    lumiblock;
//   UInt_t    ndimuon;
//   UInt_t    nmuons;
//   UInt_t    trigger;
//   UInt_t    tMatch;
//   Int_t     charge;
//
//   std::vector < TLorentzVector > gen_dau_p4, gen_gda_p4;
//
// 	TLorentzVector dimuon_p4;
// 	TLorentzVector highMuon_p4;
// 	TLorentzVector lowMuon_p4;
//
//   Float_t MassErr;
//   Float_t vProb;
//   Float_t DCA;
//   Float_t ppdlPV;
//   Float_t ppdlErrPV;
//   Float_t ppdlBS;
//   Float_t ppdlErrBS;
//   Float_t cosAlpha;
//   Float_t lxyPV;
//   Float_t lxyBS;
//
// 	UInt_t numPrimaryVertices;
//
// 	TTree *gen_tree;
//
//   Int_t mother_pdgId;
//   Int_t dimuon_pdgId;
//   TLorentzVector gen_mother_p4;
// 	TLorentzVector gen_dimuon_p4;
// 	TLorentzVector gen_highMuon_p4;
// 	TLorentzVector gen_muonM_p4;
//
//   edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
//   edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
//
// };
//
// //
// // constructors and destructor
// //
//
// GenMCRootupler::GenMCRootupler(const edm::ParameterSet & iConfig):
// PdgIds_(iConfig.getParameter<std::vector<int>>("PdgIds"))
// {
//   edm::Service < TFileService > fs;
//   gen_tree = fs->make < TTree > ("dimuonTree", "Tree of DiMuon");
//
//   gen_tree->Branch("run",      &run,      "run/D");
//   gen_tree->Branch("event",    &event,    "event/l");
//   gen_tree->Branch("lumiblock",&lumiblock,"lumiblock/D");
//
//   gen_tree->Branch("trigger",  &trigger,  "trigger/D");
//   gen_tree->Branch("gen_p4", "TLorentzVector", &dimuon_p4);
//
//   gen_tree->Branch("status",   &status,   "status/D");
//   gen_tree->Branch("charge",   &charge,   "charge/D");
//   gen_tree->Branch("pdgId",   &pdgId,   "pdgId/D");
//
//   gen_tree->Branch("pt",   &pt,   "pt/D");
//   gen_tree->Branch("eta",   &eta,   "eta/D");
//   gen_tree->Branch("phi",   &phi,   "phi/D");
//   gen_tree->Branch("mass",   &mass,   "mass/D");
//
//   gen_tree->Branch("ndaughter",   &ndaughter,   "ndaughter/D");
//   gen_tree->Branch("ngdaughter",   &ngdaughter,   "ngdaughter/D");
//
//   //Up to 4 daughtes
//   gen_tree->Branch("gen_dau1_p4", "TLorentzVector", &dimuon_p4);
//   gen_tree->Branch("gen_dau2_p4", "TLorentzVector", &dimuon_p4);
//   gen_tree->Branch("gen_dau3_p4", "TLorentzVector", &dimuon_p4);
//   gen_tree->Branch("gen_dau4_p4", "TLorentzVector", &dimuon_p4);
//
//   //Upt to 8 gdaughters
//
//   if (!OnlyGen_) {
//     gen_tree->Branch("ndimuon",    &ndimuon,    "ndimuon/D");
//     gen_tree->Branch("nmuons",   &nmuons,   "nmuons/D");
//     gen_tree->Branch("trigger",  &trigger,  "trigger/D");
//     gen_tree->Branch("charge",   &charge,   "charge/D");
//
//     gen_tree->Branch("tMatch",    &tMatch,      "tMatch/D");
//
//     gen_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
//     gen_tree->Branch("highMuon_p4",  "TLorentzVector", &highMuon_p4);
//     gen_tree->Branch("lowMuon_p4",  "TLorentzVector", &lowMuon_p4);
//
//     gen_tree->Branch("MassErr",   &MassErr,    "MassErr/F");
//     gen_tree->Branch("vProb",     &vProb,      "vProb/F");
//     gen_tree->Branch("DCA",       &DCA,        "DCA/F");
//     gen_tree->Branch("ppdlPV",    &ppdlPV,     "ppdlPV/F");
//     gen_tree->Branch("ppdlErrPV", &ppdlErrPV,  "ppdlErrPV/F");
//     gen_tree->Branch("ppdlBS",    &ppdlBS,     "ppdlBS/F");
//     gen_tree->Branch("ppdlErrBS", &ppdlErrBS,  "ppdlErrBS/F");
//     gen_tree->Branch("cosAlpha",  &cosAlpha,   "cosAlpha/F");
//     gen_tree->Branch("lxyPV",     &lxyPV,      "lxyPV/F");
//     gen_tree->Branch("lxyBS",     &lxyBS,      "lxyBS/F");
//
//     gen_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/D");
//   }
//
//   if (isMC_ || OnlyGen_) {
//      // std::cout << "GenMCRootupler::GenMCRootupler: Dimuon id " << PdgIds_ << std::endl;
//      gen_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/D");
//      gen_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/D");
//      gen_tree->Branch("gen_mother_p4", "TLorentzVector",  &gen_mother_p4);
//      gen_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
//      gen_tree->Branch("gen_highMuon_p4",  "TLorentzVector",  &gen_highMuon_p4);
//      gen_tree->Branch("gen_lowMuon_p4",  "TLorentzVector",  &gen_muonM_p4);
//   }
//
//   genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
//   packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
// }
//
// GenMCRootupler::~GenMCRootupler() {}
//
// //
// // member functions
// //
//
// const reco::Candidate* GenMCRootupler::GetAncestor(const reco::Candidate* p) {
//    if (p->numberOfMothers()) {
//       if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
//       else return p->mother(0);
//    }
//    return p;
// }
//
// //Check recursively if any ancestor of particle is the given one
// bool GenMCRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
//    if (ancestor == particle ) return true;
//    for (size_t i=0; i< particle->numberOfMothers(); i++) {
//       if (isAncestor(ancestor, particle->mother(i))) return true;
//    }
//    return false;
// }
//
// /* Grab Trigger information. Save it in variable trigger, trigger is an uint between 0 and 256, in binary it is:
//    (pass 2)(pass 1)(pass 0)
//    ex. 7 = pass 0, 1 and 2
//    ex. 6 = pass 1, 2
//    ex. 1 = pass 0
// */
//
// UInt_t GenMCRootupler::getTriggerBits(const edm::Event& iEvent ) {
//
//   UInt_t trigger = 0;
//
//   edm::Handle< edm::TriggerResults > triggerResults_handle;
//   iEvent.getByToken( triggerResults_Label , triggerResults_handle);
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
//
// // ------------ method called for each event  ------------
// void GenMCRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {
//
//   run       = iEvent.id().run();
//   event     = iEvent.id().event();
//   lumiblock = iEvent.id().luminosityBlock();
//
//   trigger = 0;
//
//   dimuon_pdgId = 0;
//   mother_pdgId = 0;
//   ndimuon  = 0;
//   nmuons = 0;
//
//   dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
//   highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
//   lowMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
//   gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
//   gen_mother_p4.SetPtEtaPhiM(0.,0.,0.,0.);
//   gen_highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
//   gen_muonM_p4.SetPtEtaPhiM(0.,0.,0.,0.);
//
//   // Pruned particles are the one containing "important" stuff
//   edm::Handle<reco::GenParticleCollection> pruned;
//   iEvent.getByToken(genCands_, pruned);
//
//
//   if (pruned.isValid() )
//   {
//     for (size_t i=0; i<pruned->size(); i++)
//     {
//       const reco::Candidate *xcand = &(*pruned)[i];
//       int thePdg = xcand->pdgId();
//
//       if ( (std::find(PdgIds_.begin(),PdgIds_.end(),thePdg) != PdgIds_.end()) && ( (xcand->status() == 2) || (xcand->status() >=11)) )
//       {
//         int nDau = xcand->numberOfDaughters();
//         if( nDau > 0)
//         {
//           std::cout << " >>>>>> Particle : " << thePdg << std::endl;
//           std::cout << " >>> no. daug: " << nDau << std::endl;
//
//           for(int jj = 0; jj<nDau;jj++)
//           {
//             std::cout << " > "<< jj << " - " << xcand->daughter(jj)->pdgId() << std::endl;
//           }
//         }
//
//       }
//     }
//   }
//
// }
//
// // ------------ method called once each job just before starting event loop  ------------
// void GenMCRootupler::beginJob() {}
//
// // ------------ method called once each job just after ending the event loop  ------------
// void GenMCRootupler::endJob() {}
//
// // ------------ method called when starting to processes a run  ------------
// void GenMCRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}
//
// // ------------ method called when ending the processing of a run  ------------
// void GenMCRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}
//
// // ------------ method called when starting to processes a luminosity block  ------------
// void GenMCRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}
//
// // ------------ method called when ending the processing of a luminosity block  ------------
// void GenMCRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}
//
// // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
// void GenMCRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
// 	//The following says we do not know what parameters are allowed so do no validation
// 	// Please change this to state exactly what you do use, even if it is no parameters
// 	edm::ParameterSetDescription desc;
// 	desc.setUnknown();
// 	descriptions.addDefault(desc);
// }
//
// //define this as a plug-in
// DEFINE_FWK_MODULE(GenMCRootupler);
