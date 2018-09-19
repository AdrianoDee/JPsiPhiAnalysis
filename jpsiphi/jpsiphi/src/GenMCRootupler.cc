// -*- C++ -*-
//
// Package:    GenMCRootupler
// Class:      GenMCRootupler
//
// Author:  Adriano Di Florio
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

class GenMCRootupler:public edm::EDAnalyzer {
public:
  explicit GenMCRootupler(const edm::ParameterSet &);
  ~GenMCRootupler() override;

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


  std::vector<int>  PdgIds_,GoodDaughters_,GoodGDaughters_;
  int MaxNumOfDaughters_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;

  std::vector<double> DimuonMassCuts_;
  bool isMC_;
  bool OnlyBest_;
  bool OnlyGen_;
  std::vector<std::string>  HLTs_;


  UInt_t   run, event, lumiblock, trigger;
  Double_t charge, pdgId, pt, eta, phi, p;
  Double_t mass, status, isPrompt;
  Double_t isGood;

  std::vector < TLorentzVector > gen_dau_p4, dummyP4;
  std::vector < std::vector < TLorentzVector > > gen_gdau_p4;

  std::vector < Double_t > gen_dau_pt, gen_dau_eta, gen_dau_phi, gen_dau_p;
  std::vector < Double_t > gen_dau_mass, gen_dau_pdg, gen_dau_status, dummy;
  std::vector < Double_t > gen_dau_ndaughter,gen_dau_isPrompt, gen_dau_charge;
  std::vector < std::vector < Double_t > > gen_gdau_pt, gen_gdau_eta, gen_gdau_phi, gen_gdau_charge, gen_gdau_isPrompt;
  std::vector < std::vector < Double_t > > gen_gdau_p, gen_gdau_mass, gen_gdau_pdg, gen_gdau_status;

  TLorentzVector gen_p4;

  UInt_t numPrimaryVertices, ndaughter, ngdaughter;

  TTree *gen_tree;

  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;

};

//
// constructors and destructor
//

GenMCRootupler::GenMCRootupler(const edm::ParameterSet & iConfig):
PdgIds_(iConfig.getParameter<std::vector<int>>("PdgIds")),
GoodDaughters_(iConfig.getParameter<std::vector<int>>("GoodDaughters")),
GoodGDaughters_(iConfig.getParameter<std::vector<int>>("GoodGDaughters")),
MaxNumOfDaughters_(iConfig.existsAs<int>("MaxDaughters") ? iConfig.getParameter<double>("MaxDaughters") : 4),
triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
HLTs_(iConfig.getParameter<std::vector<std::string>>("HLTs"))
{
  edm::Service < TFileService > fs;
  gen_tree = fs->make < TTree > ("genTree", "Tree of Gen Particles");

  gen_tree->Branch("run",      &run,      "run/i");
  gen_tree->Branch("event",    &event,    "event/i");
  gen_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");

  gen_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");

  gen_tree->Branch("trigger",  &trigger,  "trigger/i");
  gen_tree->Branch("gen_p4", "TLorentzVector", &gen_p4);

  gen_tree->Branch("status",    &status,  "status/D");
  gen_tree->Branch("charge",    &charge,  "charge/D");
  gen_tree->Branch("pdgId",     &pdgId,   "pdgId/D");

  gen_tree->Branch("pt",        &pt,      "pt/D");
  gen_tree->Branch("eta",       &eta,     "eta/D");
  gen_tree->Branch("phi",       &phi,     "phi/D");
  gen_tree->Branch("mass",      &mass,    "mass/D");

  gen_tree->Branch("p",         &p,       "p/D");

  gen_tree->Branch("isPrompt",  &isPrompt,"isPrompt/D");

  gen_tree->Branch("ndaughter",   &ndaughter,   "ndaughter/i");
  gen_tree->Branch("ngdaughter",  &ngdaughter,  "ngdaughter/i");

  gen_tree->Branch("isGood",  &isGood,  "isGood/D");

  TLorentzVector zero;
  zero.SetPtEtaPhiM(-10.0,-10.0,-10.0,-10.0);

  //Up to n daughtes (n**2 gran daughters)

  for (int i = 0; i < MaxNumOfDaughters_; i++)
  {
    dummyP4.clear();

    gen_dau_p4.push_back(zero);

    for (int j = 0; j < MaxNumOfDaughters_; j++)
    dummyP4.push_back(zero);

    gen_gdau_p4.push_back(dummyP4);
  }

  std::string name, var;


  for (int i = 0; i < MaxNumOfDaughters_; i++)
  {
    dummy.clear();

    gen_dau_pt.push_back(-1.0);
    gen_dau_eta.push_back(-5.0);
    gen_dau_phi.push_back(-10.0);
    gen_dau_p.push_back(-1.0);
    gen_dau_mass.push_back(-1.0);
    gen_dau_pdg.push_back(0.0);
    gen_dau_status.push_back(-10.0);
    gen_dau_charge.push_back(-10.0);
    gen_dau_ndaughter.push_back(-10.0);
    gen_dau_isPrompt.push_back(-10.0);

    for (int j = 0; j < MaxNumOfDaughters_; j++)
      dummy.push_back(-100.0);

    gen_gdau_pt.push_back(dummy);
    gen_gdau_eta.push_back(dummy);
    gen_gdau_phi.push_back(dummy);
    gen_gdau_p.push_back(dummy);
    gen_gdau_mass.push_back(dummy);
    gen_gdau_pdg.push_back(dummy);
    gen_gdau_status.push_back(dummy);
    gen_gdau_charge.push_back(dummy);

  }

  for (int i = 0; i < MaxNumOfDaughters_; i++)
  {
    name = "gen_dau_" + std::to_string(i) + "_p4";
    gen_tree->Branch(name.c_str(), "TLorentzVector", &gen_dau_p4[i]);


    name = "gen_dau_" + std::to_string(i) + "_pdg"; var = name + "/D";
    gen_tree->Branch(name.c_str(),&gen_dau_pdg[i],var.c_str());

    name = "gen_dau_" + std::to_string(i) + "_charge"; var = name + "/D";
    gen_tree->Branch(name.c_str(),&gen_dau_charge[i],var.c_str());

    name = "gen_dau_" + std::to_string(i) + "_status"; var = name + "/D";
    gen_tree->Branch(name.c_str(),&gen_dau_status[i],var.c_str());


    name = "gen_dau_" + std::to_string(i) + "_pt"; var = name + "/D";
    gen_tree->Branch(name.c_str(),&gen_dau_pt[i],var.c_str());

    name = "gen_dau_" + std::to_string(i) + "_eta"; var = name + "/D";
    gen_tree->Branch(name.c_str(),&gen_dau_eta[i],var.c_str());

    name = "gen_dau_" + std::to_string(i) + "_phi"; var = name + "/D";
    gen_tree->Branch(name.c_str(),&gen_dau_phi[i],var.c_str());

    name = "gen_dau_" + std::to_string(i) + "_mass"; var = name + "/D";
    gen_tree->Branch(name.c_str(),&gen_dau_mass[i],var.c_str());

    name = "gen_dau_" + std::to_string(i) + "_p"; var = name + "/D";
    gen_tree->Branch(name.c_str(),&gen_dau_p[i],var.c_str());


    name = "gen_dau_" + std::to_string(i) + "_isPrompt"; var = name + "/D";
    gen_tree->Branch(name.c_str(),&gen_dau_isPrompt[i],var.c_str());

    name = "gen_dau_" + std::to_string(i) + "_ndaughter"; var = name + "/D";
    gen_tree->Branch(name.c_str(),&gen_dau_ndaughter[i],var.c_str());

    for (int j = 0; j < MaxNumOfDaughters_; j++)
    {
      name = "gen_gdau_" + std::to_string(i) + "_" + std::to_string(j) + "_p4";
      gen_tree->Branch(name.c_str(), "TLorentzVector", &gen_gdau_p4[i][j]);


      name = "gen_gdau_" + std::to_string(i) + "_" + std::to_string(j) + "_pdg"; var = name + "/D";
      gen_tree->Branch(name.c_str(),&gen_gdau_pdg[i][j],var.c_str());

      name = "gen_gdau_" + std::to_string(i) + "_" + std::to_string(j) + "_charge"; var = name + "/D";
      gen_tree->Branch(name.c_str(),&gen_gdau_charge[i][j],var.c_str());

      name = "gen_gdau_" + std::to_string(i) + "_" + std::to_string(j) + "_status"; var = name + "/D";
      gen_tree->Branch(name.c_str(),&gen_gdau_status[i][j],var.c_str());


      name = "gen_gdau_" + std::to_string(i) + "_" + std::to_string(j) + "_pt"; var = name + "/D";
      gen_tree->Branch(name.c_str(),&gen_gdau_pt[i][j],var.c_str());

      name = "gen_gdau_" + std::to_string(i) + "_" + std::to_string(j) + "_eta"; var = name + "/D";
      gen_tree->Branch(name.c_str(),&gen_gdau_eta[i][j],var.c_str());

      name = "gen_gdau_" + std::to_string(i) + "_" + std::to_string(j) + "_phi"; var = name + "/D";
      gen_tree->Branch(name.c_str(),&gen_gdau_phi[i][j],var.c_str());

      name = "gen_gdau_" + std::to_string(i) + "_" + std::to_string(j) + "_m"; var = name + "/D";
      gen_tree->Branch(name.c_str(),&gen_gdau_mass[i][j],var.c_str());

      name = "gen_gdau_" + std::to_string(i) + "_" + std::to_string(j) + "_p"; var = name + "/D";
      gen_tree->Branch(name.c_str(),&gen_gdau_p[i][j],var.c_str());


      // name = "gen_gdau_" + std::to_string(i) + "_" + std::to_string(j) + "_p"; var = name + "/D";
      // gen_tree->Branch(name.c_str(),&gen_gdau_isPrompt[i][j],var.c_str());

    }

  }

  genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
  packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
}

GenMCRootupler::~GenMCRootupler() {}

//
// member functions
//

const reco::Candidate* GenMCRootupler::GetAncestor(const reco::Candidate* p) {
  if (p->numberOfMothers()) {
    if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
    else return p->mother(0);
  }
  return p;
}

//Check recursively if any ancestor of particle is the given one
bool GenMCRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
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

UInt_t GenMCRootupler::getTriggerBits(const edm::Event& iEvent ) {

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
void GenMCRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  trigger = 0;

  // dimuon_pdgId = 0;
  // mother_pdgId = 0;
  // ndimuon  = 0;
  // nmuons = 0;

  // dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  // highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  // lowMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  // gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  // gen_mother_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  // gen_highMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  // gen_muonM_p4.SetPtEtaPhiM(0.,0.,0.,0.);

  // Pruned particles are the one containing "important" stuff
  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genCands_, pruned);

  edm::Handle<std::vector<pat::PackedGenParticle>> packed;
  iEvent.getByToken(packCands_,  packed);

  // std::map < int, int > prunedMap, packedMap;
  // std::map < int, int > parCountPruned,parCoutPacked;
  // parCountPruned[13] = 0;
  // parCountPruned[211] = 0;
  // parCountPruned[321] = 0;
  // parCoutPacked[13] = 0;
  // parCoutPacked[211] = 0;
  // parCoutPacked[321] = 0

  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  numPrimaryVertices = 0;
  if (primaryVertices_handle.isValid()) numPrimaryVertices = (int) primaryVertices_handle->size();
  trigger = getTriggerBits(iEvent);

  std::sort(GoodDaughters_.begin(),GoodDaughters_.end());
  std::sort(GoodGDaughters_.begin(),GoodGDaughters_.end());

  if (pruned.isValid() )
  {
    std::vector<int> theDaughters,theGDaughters;

    for (size_t i=0; i<pruned->size(); i++)
    {
      bool thisIsGood;
      const reco::GenParticle *xcand = &(*pruned)[i];
      int thePdg = xcand->pdgId();

      theDaughters.clear();
      theGDaughters.clear();

      if ( (std::find(PdgIds_.begin(),PdgIds_.end(),thePdg) != PdgIds_.end()))
      {

        status        = (Double_t)xcand->status();
        charge        = (Double_t)xcand->charge();
        pdgId         = (Double_t)thePdg;

        pt            = (Double_t)xcand->pt();
        eta           = (Double_t)xcand->eta();
        phi           = (Double_t)xcand->phi();
        mass          = (Double_t)xcand->mass();
        mass          = (Double_t)xcand->p();

        gen_p4.SetPtEtaPhiM(pt,eta,phi,mass);

        isPrompt      = (Double_t)xcand->isPromptDecayed();
        ndaughter     = (UInt_t)xcand->numberOfDaughters();

        ngdaughter    = 0;
        thisIsGood        = (ndaughter==(UInt_t)(GoodDaughters_.size()));

        int maxD = -std::max(-MaxNumOfDaughters_,-(int)ndaughter);

        for(int jj = 0; jj<maxD;jj++)
        {
          auto thisDau = xcand->daughter(jj);

          gen_dau_pt[jj]     = (Double_t)thisDau->pt();
          gen_dau_eta[jj]    = (Double_t)thisDau->eta();
          gen_dau_phi[jj]    = (Double_t)thisDau->phi();
          gen_dau_mass[jj]   = (Double_t)thisDau->mass();

          gen_dau_p4[jj].SetPtEtaPhiM(gen_dau_pt[jj],gen_dau_eta[jj],gen_dau_phi[jj],gen_dau_mass[jj]);

          gen_dau_p[jj]      = thisDau->p();

          gen_dau_charge[jj] = (Double_t)thisDau->charge();
          gen_dau_pdg[jj]    = (Double_t)thisDau->pdgId();
          gen_dau_status[jj] = (Double_t)thisDau->status();

          theDaughters.push_back(gen_dau_pdg[jj]);

          gen_dau_isPrompt[jj]      = 0.0;//(Double_t)thisDau->isPromptDecayed();
          gen_dau_ndaughter[jj]     = (Double_t)thisDau->numberOfDaughters();

          ngdaughter += thisDau->numberOfDaughters();


          int maxGD = -std::max(-MaxNumOfDaughters_,-(int)gen_dau_ndaughter[jj]);

          for(int kk = 0; kk<maxGD;kk++)
          {
            auto thisGDau = thisDau->daughter(kk);

            gen_gdau_pt[jj][kk]     = (Double_t)thisGDau->pt();
            gen_gdau_eta[jj][kk]    = (Double_t)thisGDau->eta();
            gen_gdau_phi[jj][kk]    = (Double_t)thisGDau->phi();
            gen_gdau_mass[jj][kk]   = (Double_t)thisGDau->mass();

            gen_gdau_p4[jj][kk].SetPtEtaPhiM(gen_gdau_pt[jj][kk],gen_gdau_eta[jj][kk],gen_gdau_phi[jj][kk],gen_gdau_mass[jj][kk]);

            gen_gdau_p[jj][kk]      = thisGDau->p();

            gen_gdau_charge[jj][kk] = (Double_t)thisGDau->charge();
            gen_gdau_pdg[jj][kk]    = (Double_t)thisGDau->pdgId();
            gen_gdau_status[jj][kk] = (Double_t)thisGDau->status();

            theGDaughters.push_back(gen_gdau_pdg[jj][kk]);

            // gen_gdau_isPrompt[jj][kk]      = 0.0;//(Double_t)thisGDau->isPromptDecayed();

          }//gdaugher cycle

        }//daugher cycle

        std::sort(theDaughters.begin(),theDaughters.end());
        std::sort(theGDaughters.begin(),theGDaughters.end());

        std::vector < int > dComparison, gdComparison;

        std::set_intersection(theDaughters.begin(), theDaughters.end(), GoodDaughters_.begin(), GoodDaughters_.end(), std::back_inserter(dComparison));
        std::set_intersection(theGDaughters.begin(), theGDaughters.end(), GoodGDaughters_.begin(), GoodGDaughters_.end(), std::back_inserter(gdComparison));

        thisIsGood = (theDaughters.size()  == GoodDaughters_.size())  && thisIsGood;
        thisIsGood = (theGDaughters.size() == GoodGDaughters_.size()) && thisIsGood;
        thisIsGood = (theDaughters.size()  == dComparison.size())     && thisIsGood;
        thisIsGood = (theGDaughters.size() == gdComparison.size())    && thisIsGood;

        // std::cout << " >>>>>> " << i << " Particle : " << thePdg << std::endl;
        // std::cout << " >>> no. daug  : " << ndaughter << std::endl;
        // std::cout << " >>> no. gdaug : " << ngdaughter << std::endl;
        // std::cout << " >>> status    : " << xcand->status() << std::endl;
        // std::cout << " >>> ancestor  : " << GetAncestor(xcand)->pdgId() << std::endl;
        // std::cout << " >>> Daughters : ";
        // for (size_t i = 0; i < theDaughters.size(); i++)  std::cout << theDaughters[i] << " - ";
        // std::cout << std::endl;
        // std::cout << " >>> GDaughters : ";
        // for (size_t i = 0; i < theGDaughters.size(); i++)  std::cout << theGDaughters[i] << " - ";
        // std::cout << std::endl;
        //
        // std::cout << " >>> Good Daughters : ";
        // for (size_t i = 0; i < GoodDaughters_.size(); i++)  std::cout << GoodDaughters_[i] << " - ";
        // std::cout << std::endl;
        // std::cout << " >>> Good GDaughters : ";
        // for (size_t i = 0; i < GoodGDaughters_.size(); i++)  std::cout << GoodGDaughters_[i] << " - ";
        // std::cout << std::endl;

        isGood = (Double_t)thisIsGood;

        gen_tree->Fill();

        // std::cout << " >>> Is Good : " << thisIsGood << std::endl;
        // std::cout << " >>> dau     : " << (theDaughters.size()  == GoodDaughters_.size()) << std::endl;
        // std::cout << " >>> gdau    : " << (theGDaughters.size() == GoodGDaughters_.size()) << std::endl;
        // std::cout << " >>> dau in  : " << (theDaughters.size()  == gdComparison.size()) << std::endl;
        // std::cout << " >>> gdau in : " << (theGDaughters.size() == gdComparison.size()) << std::endl;
      }

    }
}
// if (packed.isValid() )
// {
//   std::<< "Packed is valid "<<std::endl;
//   for (size_t i=0; i<packed->size(); i++)
//   {
//     const reco::Candidate *xcand = &(*packed)[jj];
//     int thePdg = xcand->pdgId();
//
//     if ( (std::find(PdgIds_.begin(),PdgIds_.end(),thePdg) != PdgIds_.end()) && ( (xcand->status() == 2) || (xcand->status() >=11)) )
//     {
//       int nDau = xcand->numberOfDaughters();
//       if( nDau > 0)
//       {
//         std::cout << " >>>>>> " << i << " Packed Particle : " << thePdg << std::endl;
//         std::cout << " >>> no. daug : " << nDau << std::endl;
//         std::cout << " >>> status   : " << xcand->status() << std::endl;
//         std::cout << " >>> ancestor : " << GetAncestor(xcand)->pdgId() << std::endl;
//
//         packedMap[xcand->status()] = 1.0;
//
//         for(int jj = 0; jj<nDau;jj++)
//         {
//           std::cout << " > "<< jj << " - " << xcand->daughter(jj)->pdgId() << std::endl;
//         }
//       }
//
//     }
//   }
// }

// std::cout << "Pruned ";
// for ( const auto &myPair : prunedMap ) {
//   std::cout << myPair.first << " - ";
// }
// std::cout << "\n";
//
// std::cout << "Packed ";
// for ( const auto &myPair : packedMap ) {
//   std::cout << myPair.first << " - ";
// }
// std::cout << "\n";

}

// ------------ method called once each job just before starting event loop  ------------
void GenMCRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void GenMCRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
void GenMCRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void GenMCRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void GenMCRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void GenMCRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenMCRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenMCRootupler);
