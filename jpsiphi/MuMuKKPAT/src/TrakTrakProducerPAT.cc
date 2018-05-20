// -*- C++ -*-
//
// Package:    TrakTrakProducerPAT
// Class:      TrakTrakProducerPAT
//
/**\class TrakTrakProducerPAT TrakTrakProducerPAT.cc myProducers/TrakTrakProducerPAT/src/TrakTrakProducerPAT.cc

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
#include "TLorentzVector.h"

#include "../interface/TrakTrakProducerPAT.h"
#include "../interface/VertexReProducer.h"
//#include "DataFormats/Candidate/interface/OverlapChecker.h"

///
/// constants, enums and typedefs
///

typedef math::Error<3>::type CovarianceMatrix;

//float small_sigma = kaon_mass*1.e-6; /// SEMRA

///
/// static data member definitions
///

///
/// constructors and destructor
///

TrakTrakProducerPAT::TrakTrakProducerPAT(const edm::ParameterSet& iConfig):

hlTriggerResults_(iConfig.getUntrackedParameter<edm::InputTag>("HLTriggerResults",edm::InputTag("TriggerResults::HLT")) ),
inputGEN_(iConfig.getUntrackedParameter<edm::InputTag>("inputGEN",edm::InputTag("genParticles"))),
vtxSample_(iConfig.getUntrackedParameter<std::string>("vtxSample_",std::string("offlinePrimaryVertices"))),

phiMassCuts_(iConfig.getParameter<std::vector<double>>("PhiMassCuts")),

doData_( iConfig.getUntrackedParameter<bool>("DoDataAnalysis", true) ),
doMC_( iConfig.getUntrackedParameter<bool>("DoMonteCarloTree", true) ),

addCommonVertex_(iConfig.getParameter<bool>("addCommonVertex")),
resolveAmbiguity_(iConfig.getParameter<bool>("resolvePileUpAmbiguity")),

addMCTruth_(iConfig.getParameter<bool>("addMCTruth")),
MCParticle_( iConfig.getUntrackedParameter<int>("MonteCarloParticleId", 20443) ), /// 20443 X, 100443 Psi(2S), 9120443 X from B / decide later for X(4140)
MCExclusiveDecay_( iConfig.getUntrackedParameter<bool>("MonteCarloExclusiveDecay", true) ),
MCMother_( iConfig.getUntrackedParameter<int>("MonteCarloMotherId", 511) ), /// 511 B0 (=anti-B0), 531 B0 / decide later MCMotherId for X(4140)
MCDaughtersN_( iConfig.getUntrackedParameter<int>(" MonteCarloDaughtersN", 3) ), /// will be same

TrMinPixHits_(iConfig.getUntrackedParameter<int>("MinNumTrPixHits", 0)),
TrMinSiHits_(iConfig.getUntrackedParameter<int>("MinNumTrSiHits", 0)),
TrMaxNormChi_(iConfig.getUntrackedParameter<double>("MaxTrNormChi2", 1000)),
TrMaxD0_(iConfig.getUntrackedParameter<double>("MaxTrD0", 1000)),

TriggerCut_(iConfig.getUntrackedParameter<bool>("TriggerCut",true)),
HLTPaths_(iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching")),
FiltersForMatching_(iConfig.getUntrackedParameter<std::vector<std::string> >("FiltersForMatching")),
Debug_(iConfig.getUntrackedParameter<bool>("Debug_Output",true)),
SameSign_(iConfig.getUntrackedParameter<bool>("SameSign",false))

{
  // revtxtrks_ = "generalTracks"; //if that is not true, we will raise an exception
  // revtxbs_ = "offlineBeamSpot";
  // genCands_ = "genParticles";

  produces<pat::CompositeCandidateCollection>( "DiTrakCandidates" ).setBranchAlias( "DiTrakCandidates");

  /// now do what ever initialization is needed


}

TrakTrakProducerPAT::~TrakTrakProducerPAT()
{
  /// do anything here that needs to be done at desctruction time
  /// (e.g. close files, deallocate resources etc.)

}

UInt_t TrakTrakProducerPAT::getTriggerBits(const edm::Event& iEvent ) {

  UInt_t trigger = 0;

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByLabel( hlTriggerResults_ , triggerResults_handle);

  if (triggerResults_handle.isValid()) {
     const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
     unsigned int NTRIGGERS = HLTPaths_.size();

     for (unsigned int i = 0; i < NTRIGGERS; i++) {
        for (int version = 1; version < 20; version++) {
           std::stringstream ss;
           ss << HLTPaths_[i] << "_v" << version;
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


///
/// member functions
///

/// ------------ method called to for each event  ------------
void TrakTrakProducerPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  const ParticleMass kaon_mass = 0.493667; //pdg mass
  float small_sigma = kaon_mass*1.e-6; /// SEMRA

  using namespace edm;
  using namespace std;
  using namespace reco;

  int evtNum = iEvent.id().event();

  std::auto_ptr<pat::CompositeCandidateCollection> oniaOutput(new pat::CompositeCandidateCollection);

  bool decayChainOK = false;

  ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  /// first get HLT results
  unsigned int triggerBit = -1;

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

    if (ntrigs==0)
      std::cout << "No trigger name given in TriggerResults of the input " << std::endl;

    triggerBit = getTriggerBits(iEvent);

  } /// end valid trigger


  Vertex thePrimaryVtx, theBeamSpotVtx;
  math::XYZPoint RefVtx;

  // Int_t thePrimaryVtx_multiplicity = -1 ;

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
    // thePrimaryVtx_multiplicity = recVtxs->size() ;

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
    // thePrimaryVtx_multiplicity = 1 ;
  }

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  RefVtx = thePrimaryVtx.position(); /// reference primary vertex choosen

  /// /// /// /// /// /// /// /// /// /// /// /// /// ///
  /// Paricles
  /// /// /// /// /// /// /// /// /// /// /// /// /// ///

  Handle< std::vector<pat::GenericParticle> > thePATTrackHandle;
  iEvent.getByLabel("cleanPatTrackCands", thePATTrackHandle); /// container of tracks with pion mass hypothesis
  Handle< std::vector<pat::GenericParticle> > theKaonRefittedPATTrackHandle;
  iEvent.getByLabel("cleanPatTrackKaonCands", theKaonRefittedPATTrackHandle);
  Handle< std::vector<pat::Muon> > thePATMuonHandle;
  iEvent.getByLabel("patMuonsWithTrigger", thePATMuonHandle);

  // //Trigger Event for matching NOT USEFULL IN 2012, no track matched
  //
  // edm::Handle<edm::hlTriggerEvent_> hltEvent;
  // iEvent.getByLabel(hlTriggerEvent_, hltEvent);
  //
  // std::vector < TriggerObjectRefVector > theFilteredTriggers;
  //
  // //Filtering trigger objects by filter name
  //
  // for (size_t i = 0; i < FiltersForMatching_.size(); ++i)
  // {
  //   TriggerObjectRefVector thisTriggers = hltEvent.filterObjects(FiltersForMatching_[i]);
  //
  //   if(thisTriggers.isNull())
  //     continue;
  //   else
  //   {
  //     theFilteredTriggers.push_back(thisTriggers);
  //   }
  // }
  //
  // for (size_t i = 0; i < theFilteredTriggers.size(); i++) {
  //   /* code */
  // }
  //
  // for ( std::vector<pat::GenericParticle>::const_iterator trackPos = theKaonRefittedPATTrackHandle->begin(); trackPos != theKaonRefittedPATTrackHandle->end(); ++trackPos )
  // {
  //
  // }
  //
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
  //           bool trktrkOK = false;
  //           bool kkOK = false;
  //           //bool pionOK = false, kaonOK = false;
  //
  //           for (int j=0; j<dauNum; ++j) {
  //             const Candidate *dau = p.daughter(j);
  //             if (Debug_) std::cout << "dauPdgId = " << dau->pdgId() << std::endl;
  //
  //             /// check if one of B0 daughters is a psi(nS) whitch has 2 muons as daughters /// SEMRA ask again !!!
  //             int trktrkId = 0 ;
  //             if (skipJPsi) /// SEMRA cleaned skipPsi2S
  //             if (Debug_) std::cout <<"Skipping J/psi!" <<std::endl ; /// SEMRA cleaned skipPsi2S
  //             //else if (skipPsi2S) /// SEMRA
  //             //  trktrkId = 443 ; /// SEMRA (JPsi ID)
  //
  //             if ( ((skipJPsi) && (dau->pdgId() == trktrkId)) ||
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
  //               if ( muNum == 2 ) trktrkOK = true ;
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
  //       if (Debug_) std::cout << "trktrkOK = " << trktrkOK << ", kkOK = " << kkOK << std::endl;
  //       if ( trktrkOK && kkOK ) {
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

    bool isEventWithInvalidTrack = false;

    if (Debug_) std::cout << "Starting event with " << thePATMuonHandle->size() << " muons" << std::endl;

    if ((thePATMuonHandle->size() > 10000) || (thePATTrackHandle->size() > 10000))
      std::cout << "Too many Muons: " << thePATMuonHandle->size() << std::endl;
    else //if (thePATMuonHandle->size() >= 2) { // check
      if (thePATTrackHandle->size() >= 2  && (triggerBit > 0 || !TriggerCut_)) {

        if (Debug_) std::cout <<"============================  Evt: " <<evtNum <<" accept event with 2 mu and trigger ==============================================" <<std::endl;

        float posTrackDzVtx, posTrackDxyVtx, negTrackDzVtx, negTrackDxyVtx;

        /// get TrTr cands
        for ( std::vector<pat::GenericParticle>::const_iterator trackPos = theKaonRefittedPATTrackHandle->begin(); trackPos != theKaonRefittedPATTrackHandle->end(); ++trackPos ) {

          if (trackPos->track().isNull()) continue;

          if (trackPos->charge()<=0.0 && !SameSign_) continue;

          if (trackPos->pt()<=0.7) continue;

          if ( !(trackPos->track()->hitPattern().trackerLayersWithMeasurement()) ) {
            isEventWithInvalidTrack = true;
            if (Debug_) std::cout <<"evt:" <<evtNum <<" problem with trackerLayersWithMeasurement" <<std::endl;
            continue ;
          }
          if ( !(trackPos->track()->hitPattern().pixelLayersWithMeasurement()) ) {
            if (Debug_) std::cout <<"evt:" <<evtNum <<" problem with pixelLayersWithMeasurement" <<std::endl;
            isEventWithInvalidTrack = true;
            continue ;
          }

          /// cuts on trackpos
          if (trackPos->track()->hitPattern().numberOfValidPixelHits() < TrMinPixHits_
              || trackPos->track()->hitPattern().numberOfValidStripHits() < TrMinSiHits_
              || trackPos->track()->chi2()/trackPos->track()->ndof() > TrMaxNormChi_
              || fabs(trackPos->track()->dxy(RefVtx)) > TrMaxD0_)

            continue;

          posTrackDzVtx = trackPos->track()->dz(RefVtx);
          posTrackDxyVtx = trackPos->track()->dxy(RefVtx);

          ////////////////// check for muon2 //////////////////
          for ( std::vector<pat::GenericParticle>::const_iterator trackNeg = trackPos+1; trackNeg != theKaonRefittedPATTrackHandle->end(); ++trackNeg ){

            if (trackNeg->track().isNull()) continue;

            if (trackNeg->track().key() == trackPos->track().key())
            continue ;

            if(trackNeg->charge() >= 0 && !SameSign_) continue ;

            if(trackNeg->pt() <= 0.7) continue ;

            if ( !(trackNeg->track()->hitPattern().trackerLayersWithMeasurement()) ) {
              isEventWithInvalidTrack = true;
              if (Debug_) std::cout <<"evt:" <<evtNum <<" problem with trackerLayersWithMeasurement" <<std::endl;
              continue ;
            }
            if ( !(trackNeg->track()->hitPattern().pixelLayersWithMeasurement()) ) {
              if (Debug_) std::cout <<"evt:" <<evtNum <<" problem with pixelLayersWithMeasurement" <<std::endl;
              isEventWithInvalidTrack = true;
              continue ;
            }


            /// cuts on track 2
            if (trackNeg->track()->hitPattern().numberOfValidPixelHits() < TrMinPixHits_
                || trackNeg->track()->hitPattern().numberOfValidStripHits() < TrMinSiHits_
                || trackNeg->track()->chi2()/trackNeg->track()->ndof() > TrMaxNormChi_
                || fabs(trackNeg->track()->dxy(RefVtx)) > TrMaxD0_)

                continue;

            negTrackDzVtx = trackNeg->bestTrackRef()->dz(RefVtx);
            negTrackDxyVtx = trackNeg->bestTrackRef()->dxy(RefVtx);

            //Vertex fit
            TransientTrack kaonPosTT( trackPos->track(), &(*bFieldHandle) );
            TransientTrack kaonNegTT( trackNeg->track(), &(*bFieldHandle) );
            KinematicParticleFactoryFromTransientTrack pFactory;

            /// initial chi2 and ndf before kinematic fits
            float chi = 0., ndf = 0.;

            std::vector<RefCountedKinematicParticle> kaons;
            kaons.push_back( pFactory.particle( kaonPosTT, kaon_mass, chi, ndf, small_sigma));
            kaons.push_back( pFactory.particle( kaonNegTT, kaon_mass, chi, ndf, small_sigma));
            KinematicParticleVertexFitter KKFitter;
            RefCountedKinematicTree KKVertexFitTree;
            KKVertexFitTree = KKFitter.fit(kaons);

            if (!KKVertexFitTree->isValid())
            continue ;

            KKVertexFitTree->movePointerToTheTop();
            RefCountedKinematicParticle KKCand_fromFit = KKVertexFitTree->currentParticle();
            RefCountedKinematicVertex KKCand_vertex_fromFit = KKVertexFitTree->currentDecayVertex();

            float trktrkVProb = ChiSquaredProbability((float)( KKCand_vertex_fromFit->chiSquared()),(float)( KKCand_vertex_fromFit->degreesOfFreedom()));

            if (trktrkVProb < 0.001)
            continue;

            float trktrkChi2 = KKCand_vertex_fromFit->chiSquared();
            float trktrkNDof = KKCand_vertex_fromFit->degreesOfFreedom();

            KKVertexFitTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle kaonPosCand_fromFit = KKVertexFitTree->currentParticle();
            KKVertexFitTree->movePointerToTheNextChild();
            RefCountedKinematicParticle kaonNegCand_fromFit = KKVertexFitTree->currentParticle();

            double ditrack_vx_fit = KKCand_vertex_fromFit->position().x();
            double ditrack_vy_fit = KKCand_vertex_fromFit->position().y();
            double ditrack_vz_fit = KKCand_vertex_fromFit->position().z();


            ////////////////// fill the KK vectors //////////////////
            if (KKCand_fromFit->currentState().mass() < phiMassCuts_[0]  ||  KKCand_fromFit->currentState().mass() > phiMassCuts_[1])
            continue ;

            float ditrack_ma_fit = KKCand_fromFit->currentState().mass();
            int   ditrack_ch_fit = KKCand_fromFit->currentState().particleCharge();
            float ditrack_px_fit = KKCand_fromFit->currentState().kinematicParameters().momentum().x();
            float ditrack_py_fit = KKCand_fromFit->currentState().kinematicParameters().momentum().y();
            float ditrack_pz_fit = KKCand_fromFit->currentState().kinematicParameters().momentum().z();
            float ditrack_en_fit = sqrt(ditrack_ma_fit*ditrack_ma_fit+ditrack_px_fit*ditrack_px_fit+ditrack_py_fit*ditrack_py_fit+ditrack_pz_fit*ditrack_pz_fit);

            reco::CompositeCandidate reco_ref_Phi(ditrack_ch_fit,math::XYZTLorentzVector(ditrack_px_fit,ditrack_py_fit,ditrack_pz_fit,ditrack_en_fit),
                                                     math::XYZPoint(ditrack_vx_fit,ditrack_vy_fit,ditrack_vz_fit),443);
            pat::CompositeCandidate pat_ref_Phi(reco_ref_Phi);

            //////////////////// For Lifetimes Calculations ////////////////////
            TVector3 TrTr_vtx((*KKCand_vertex_fromFit).position().x(), (*KKCand_vertex_fromFit).position().y(), 0) ;
            TVector3 TrTr_pperp(KKCand_fromFit->currentState().globalMomentum().x(), KKCand_fromFit->currentState().globalMomentum().y(), 0);
            TVector3 TrTr_vtx3D((*KKCand_vertex_fromFit).position().x(), (*KKCand_vertex_fromFit).position().y(), (*KKCand_vertex_fromFit).position().z()) ;
            TVector3 TrTr_pperp3D(KKCand_fromFit->currentState().globalMomentum().x(),KKCand_fromFit->currentState().globalMomentum().y(), KKCand_fromFit->currentState().globalMomentum().z());


            float kaonPos_ma_fit = kaonPosCand_fromFit->currentState().mass();
            int   kaonPos_ch_fit = kaonPosCand_fromFit->currentState().particleCharge();
            float kaonPos_px_fit = kaonPosCand_fromFit->currentState().kinematicParameters().momentum().x();
            float kaonPos_py_fit = kaonPosCand_fromFit->currentState().kinematicParameters().momentum().y();
            float kaonPos_pz_fit = kaonPosCand_fromFit->currentState().kinematicParameters().momentum().z();
            float kaonPos_en_fit = sqrt(kaonPos_ma_fit*kaonPos_ma_fit+kaonPos_px_fit*kaonPos_px_fit+kaonPos_py_fit*kaonPos_py_fit+kaonPos_pz_fit*kaonPos_pz_fit);

            reco::CompositeCandidate reco_ref_PK(kaonPos_ch_fit,math::XYZTLorentzVector(kaonPos_px_fit,kaonPos_py_fit,kaonPos_pz_fit,kaonPos_en_fit),
                                                     math::XYZPoint(ditrack_vx_fit,ditrack_vy_fit,ditrack_vz_fit),-13);
            pat::CompositeCandidate pat_ref_PK(reco_ref_PK);

            float kaonNeg_ma_fit = kaonNegCand_fromFit->currentState().mass();
            int   kaonNeg_ch_fit = kaonNegCand_fromFit->currentState().particleCharge();
            float kaonNeg_px_fit = kaonNegCand_fromFit->currentState().kinematicParameters().momentum().x();
            float kaonNeg_py_fit = kaonNegCand_fromFit->currentState().kinematicParameters().momentum().y();
            float kaonNeg_pz_fit = kaonNegCand_fromFit->currentState().kinematicParameters().momentum().z();
            float kaonNeg_en_fit = sqrt(kaonNeg_ma_fit*kaonNeg_ma_fit+kaonNeg_px_fit*kaonNeg_px_fit+kaonNeg_py_fit*kaonNeg_py_fit+kaonNeg_pz_fit*kaonNeg_pz_fit);

            reco::CompositeCandidate reco_ref_NK(kaonNeg_ch_fit,math::XYZTLorentzVector(kaonNeg_px_fit,kaonNeg_py_fit,kaonNeg_pz_fit,kaonNeg_en_fit),
                                                     math::XYZPoint(ditrack_vx_fit,ditrack_vy_fit,ditrack_vz_fit),13);

            pat::CompositeCandidate pat_ref_NK(reco_ref_NK);


            pat::CompositeCandidate trktrkCandidate;

            // ---- define and set candidate's 4momentum  ----
            ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > trktrk = trackPos->p4() + trackNeg->p4();
            TLorentzVector trP, trN,trktrkP4;

            trP.SetXYZM(trackPos->bestTrackRef()->px(),trackPos->bestTrackRef()->py(),trackPos->bestTrackRef()->pz(),kaon_mass);
            trN.SetXYZM(trackNeg->bestTrackRef()->px(),trackNeg->bestTrackRef()->py(),trackNeg->bestTrackRef()->pz(),kaon_mass);
            // LorentzVector trktrk;

            trktrkP4 = trP + trN;

            trktrkCandidate.setP4(trktrk);
            trktrkCandidate.setCharge(trackPos->charge() + trackNeg->charge());

            float deltaRTrTr = reco::deltaR2(trackPos->eta(),trackPos->phi(),trackNeg->eta(),trackNeg->phi());

            trktrkCandidate.addDaughter(*trackPos, "trackPos");
            trktrkCandidate.addDaughter(*trackNeg,"trackNeg");

            pat_ref_Phi.addUserFloat("deltaR",deltaRTrTr);
            pat_ref_Phi.addUserFloat("trktrkP4",trktrkP4.M());

            pat_ref_Phi.addDaughter(pat_ref_PK, "ref_kaonPos");
            pat_ref_Phi.addDaughter(pat_ref_NK, "ref_kaonNeg");
            pat_ref_Phi.addDaughter(trktrkCandidate, "trktrkCandidate");
            // pat_ref_Phi.addDaughter(*pat_ref_PM, "trackPos");
            // pat_ref_Phi.addDaughter(*pat_ref_NM, "trackNeg");

            pat_ref_Phi.addUserFloat("trackPosDzVtx",posTrackDzVtx);
            pat_ref_Phi.addUserFloat("trackPosDxyVtx",posTrackDxyVtx);

            pat_ref_Phi.addUserFloat("trackNegDzVtx",negTrackDzVtx);
            pat_ref_Phi.addUserFloat("trackNegDxyVtx",negTrackDxyVtx);


            pat_ref_Phi.addUserFloat("trPos_Chi2", kaonPosCand_fromFit->chiSquared());
            pat_ref_Phi.addUserFloat("trPos_NDF",  kaonPosCand_fromFit->degreesOfFreedom());
            pat_ref_Phi.addUserFloat("trNeg_Chi2", kaonNegCand_fromFit->chiSquared());
            pat_ref_Phi.addUserFloat("trNeg_NDF",  kaonNegCand_fromFit->degreesOfFreedom());

            pat_ref_Phi.addUserFloat("VProb", trktrkVProb);
            pat_ref_Phi.addUserFloat("Chi2",  trktrkChi2);
            pat_ref_Phi.addUserFloat("NDof",  trktrkNDof);
            pat_ref_Phi.addUserFloat("SS",  float(SameSign_));

            AlgebraicVector3 TrTr_v3pperp ;
            TrTr_v3pperp[0] = TrTr_pperp.x(); TrTr_v3pperp[1] = TrTr_pperp.y(); TrTr_v3pperp[2] = 0.;

            TVector3 TrTr_pvtx, TrTr_pvtx3D, TrTr_vdiff, TrTr_vdiff3D ;
            float TrTr_cosAlpha, TrTr_cosAlpha3D, TrTr_ctau ;
            VertexDistanceXY TrTr_vdistXY ;
            Measurement1D TrTr_distXY ;
            GlobalError TrTr_v1e = (Vertex(*KKCand_vertex_fromFit)).error(), TrTr_v2e;
            AlgebraicSymMatrix33 TrTr_vXYe ;
            float TrTr_ctauErr ;
            float TrTr_lxy, TrTr_lxyErr, TrTr_lxyz, TrTr_lxyzErr ;
            ROOT::Math::SVector<double, 3> TrTr_vDiff, TrTr_vDiff3D ; // needed by Similarity method

            ////////////////// Lifetime wrt PV for TrTr //////////////////
            TrTr_v2e = thePrimaryVtx.error();
            TrTr_vXYe = TrTr_v1e.matrix() + TrTr_v2e.matrix() ;

            /// 2D
            TrTr_pvtx.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), 0) ;
            TrTr_vdiff = TrTr_vtx - TrTr_pvtx ;
            TrTr_cosAlpha = TrTr_vdiff.Dot(TrTr_pperp) / (TrTr_vdiff.Perp()*TrTr_pperp.Perp()) ;
            TrTr_lxy = TrTr_vdiff.Perp();
            TrTr_vDiff[0] = TrTr_vdiff.x(); TrTr_vDiff[1] = TrTr_vdiff.y(); TrTr_vDiff[2] = 0 ; // needed by Similarity method
            TrTr_lxyErr = sqrt(ROOT::Math::Similarity(TrTr_vDiff,TrTr_vXYe)) / TrTr_vdiff.Perp();
            TrTr_distXY = TrTr_vdistXY.distance(Vertex(*KKCand_vertex_fromFit), Vertex(thePrimaryVtx));
            TrTr_ctau = TrTr_distXY.value() * TrTr_cosAlpha * KKCand_fromFit->currentState().mass() / TrTr_pperp.Perp();
            TrTr_ctauErr = sqrt(ROOT::Math::Similarity(TrTr_v3pperp,TrTr_vXYe)) * KKCand_fromFit->currentState().mass() / (TrTr_pperp.Perp2()) ;

            /// 3D
            TrTr_pvtx3D.SetXYZ(thePrimaryVtx.position().x(), thePrimaryVtx.position().y(), thePrimaryVtx.position().z());
            TrTr_vdiff3D = TrTr_vtx3D - TrTr_pvtx3D;
            TrTr_cosAlpha3D = TrTr_vdiff3D.Dot(TrTr_pperp3D)/(TrTr_vdiff3D.Mag()*TrTr_pperp3D.Mag());
            TrTr_lxyz = TrTr_vdiff3D.Mag();
            TrTr_vDiff3D[0] = TrTr_vdiff3D.x(); TrTr_vDiff3D[1] = TrTr_vdiff3D.y(); TrTr_vDiff3D[2] = TrTr_vdiff3D.z() ;
            TrTr_lxyzErr = sqrt(ROOT::Math::Similarity(TrTr_vDiff3D,TrTr_vXYe)) / TrTr_vdiff3D.Mag();

            pat_ref_Phi.addUserFloat("cosAlpha",    TrTr_cosAlpha);
            pat_ref_Phi.addUserFloat("cosAlpha3D",  TrTr_cosAlpha3D);
            pat_ref_Phi.addUserFloat("ctau",        TrTr_ctau);
            pat_ref_Phi.addUserFloat("ctauErr",     TrTr_ctauErr);
            pat_ref_Phi.addUserFloat("lxy",         TrTr_lxy);
            pat_ref_Phi.addUserFloat("lxyErr",      TrTr_lxyErr);
            pat_ref_Phi.addUserFloat("lxyz",        TrTr_lxyz);
            pat_ref_Phi.addUserFloat("lxyzErr",     TrTr_lxyzErr);

            ////////////////// Lifetime wrt BS for TrTr //////////////////
            TrTr_v2e = theBeamSpotVtx.error();
            TrTr_vXYe = TrTr_v1e.matrix() + TrTr_v2e.matrix();

            /// 2D
            TrTr_pvtx.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), 0);
            TrTr_vdiff = TrTr_vtx - TrTr_pvtx;
            TrTr_cosAlpha = TrTr_vdiff.Dot(TrTr_pperp)/(TrTr_vdiff.Perp()*TrTr_pperp.Perp());
            TrTr_lxy = TrTr_vdiff.Perp();
            TrTr_vDiff[0] = TrTr_vdiff.x(); TrTr_vDiff[1] = TrTr_vdiff.y(); TrTr_vDiff[2] = 0 ; // needed by Similarity method
            TrTr_lxyErr = sqrt(ROOT::Math::Similarity(TrTr_vDiff,TrTr_vXYe)) / TrTr_vdiff.Perp();
            TrTr_distXY = TrTr_vdistXY.distance(Vertex(*KKCand_vertex_fromFit), Vertex(theBeamSpotVtx));
            TrTr_ctau = TrTr_distXY.value() * TrTr_cosAlpha * (KKCand_fromFit->currentState().mass() / TrTr_pperp.Perp()) ;
            TrTr_ctauErr = sqrt(ROOT::Math::Similarity(TrTr_v3pperp,TrTr_vXYe)) * KKCand_fromFit->currentState().mass()/TrTr_pperp.Perp2();

            /// 3D
            TrTr_pvtx3D.SetXYZ(theBeamSpotVtx.position().x(), theBeamSpotVtx.position().y(), theBeamSpotVtx.position().z());
            TrTr_vdiff3D = TrTr_vtx3D - TrTr_pvtx3D;
            TrTr_cosAlpha3D = TrTr_vdiff3D.Dot(TrTr_pperp3D)/(TrTr_vdiff3D.Mag()*TrTr_pperp3D.Mag());
            TrTr_lxyz = TrTr_vdiff3D.Mag();
            TrTr_vDiff3D[0] = TrTr_vdiff3D.x(); TrTr_vDiff3D[1] = TrTr_vdiff3D.y(); TrTr_vDiff3D[2] = TrTr_vdiff3D.z() ;
            TrTr_lxyzErr = sqrt(ROOT::Math::Similarity(TrTr_vDiff3D,TrTr_vXYe)) / TrTr_vdiff3D.Mag();

            pat_ref_Phi.addUserFloat("BS_cosAlpha",    TrTr_cosAlpha);
            pat_ref_Phi.addUserFloat("BS_cosAlpha3D",  TrTr_cosAlpha3D);
            pat_ref_Phi.addUserFloat("BS_ctau",        TrTr_ctau);
            pat_ref_Phi.addUserFloat("BS_ctauErr",     TrTr_ctauErr);
            pat_ref_Phi.addUserFloat("BS_lxy",         TrTr_lxy);
            pat_ref_Phi.addUserFloat("BS_lxyErr",      TrTr_lxyErr);
            pat_ref_Phi.addUserFloat("BS_lxyz",        TrTr_lxyz);
            pat_ref_Phi.addUserFloat("BS_lxyzErr",     TrTr_lxyzErr);

            pat_ref_Phi.addUserFloat("isEventWithInvalidTrack",    float(isEventWithInvalidTrack));

            // TrTrDecayVtx_XE->push_back( sqrt( KKCand_vertex_fromFit->error().cxx()) );
            // TrTrDecayVtx_YE->push_back( sqrt( KKCand_vertex_fromFit->error().cyy()) );
            // TrTrDecayVtx_ZE->push_back( sqrt( KKCand_vertex_fromFit->error().czz()) );

            // pat_ref_Phi.addUserInt("isTriggerMatchedPos",isTriggerMatched(&(*trackPos),&(*trackNeg)));
            // pat_ref_Phi.addUserInt("isTriggerMatchedNeg",isTriggerMatched(&(*trackPos),&(*trackNeg)));

            kaons.clear();

            oniaOutput->push_back(pat_ref_Phi);

          } // 2nd loop over tracks (look for tr-)
        } //first loop over tracks (look for tr+)
      } // if (thePATTrackHandle->size() >= 2  && hasRequestedTrigger) {
      } // if (doMC_ || doData_)
      // AT THE END OF THE EVENT fill the tree and clear the vectors
      // ===========================================================

      std::sort(oniaOutput->begin(),oniaOutput->end(),vPComparator_);
      //std::cout << "MuMu candidates count : " << oniaOutput->size() << std::endl;
      //iEvent.put(std::move(oniaOutput));
      if(Debug_) std::cout << "No. ditracks: " << oniaOutput->size() << std::endl;
      iEvent.put( oniaOutput, "DiTrakCandidates" );


    }
    //}/// produce
    /// ------------ method called once each job just before starting event loop  ------------
    void TrakTrakProducerPAT::beginRun(edm::Run& iRun, edm::EventSetup const& iSetup)
    {
    }
    void TrakTrakProducerPAT::beginJob()
    {


    }/// begin Job

        // ------------ method called when ending the processing of a run  ------------
    void
    TrakTrakProducerPAT::endRun(edm::Run&, edm::EventSetup const&)
    {
    }

    // ------------ method called when starting to processes a luminosity block  ------------
    void
    TrakTrakProducerPAT::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
    {
    }

    // ------------ method called when ending the processing of a luminosity block  ------------
    void
    TrakTrakProducerPAT::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
    {
    }


    /// ------------ method called once each job just after ending the event loop  ------------
    void TrakTrakProducerPAT::endJob() {

    }/// endjob

    void
    TrakTrakProducerPAT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


    bool TrakTrakProducerPAT::isAbHadron(int pdgID) {

      if (abs(pdgID) == 511 || abs(pdgID) == 521 || abs(pdgID) == 531 || abs(pdgID) == 5122) return true;
      return false;

    }

    bool TrakTrakProducerPAT::isAMixedbHadron(int pdgID, int momPdgID) {

      if ((abs(pdgID) == 511 && abs(momPdgID) == 511 && pdgID*momPdgID < 0) ||
      (abs(pdgID) == 531 && abs(momPdgID) == 531 && pdgID*momPdgID < 0))
      return true;
      return false;

    }

    std::pair<int, float> TrakTrakProducerPAT::findCandMCInfo(reco::GenParticleRef genCand) {

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

    // float TrakTrakProducerPAT::getSigmaOfLogdEdx(float logde)
    // {
    //   return 0.3;
    // }
    //
    // float TrakTrakProducerPAT::getEnergyLoss(const reco::TrackRef & track)
    // {
    //   if (iexception_dedx==1) return 9999.;
    //   const reco::DeDxDataValueMap & eloss = *energyLoss;
    //   return eloss[track].dEdx();
    // }
    //
    // float TrakTrakProducerPAT::nsigmaofdedx(const reco::TrackRef & track, float & theo, float & sigma)
    // {
    //
    //   // no usable dE/dx if p > 2
    //   float nsigma = 99 ;
    //   if (iexception_dedx==1) return nsigma ;
    //
    //   float m  = 0.13957;
    //   float bg = track->p() / m;
    //
    //   theo = getLogdEdx(bg);
    //
    //
    //   int nhitr = track->numberOfValidHits();
    //   float meas = log(getEnergyLoss(track));
    //   sigma = getSigmaOfLogdEdx(theo) * pow(nhitr,-0.65);
    //   if (sigma>0)
    //   nsigma = (meas-theo) / sigma ;
    //   return nsigma;
    // }
    //
    //
    // float TrakTrakProducerPAT::getLogdEdx(float bg)
    // {
    //   const float a =  3.25 ;
    //   const float b =  0.288;
    //   const float c = -0.852;
    //
    //   float beta = bg/sqrt(bg*bg + 1);
    //   float dedx = log( a/(beta*beta) + b * log(bg) + c );
    //
    //   return dedx;
    //
    // }
    //
    //
    // float TrakTrakProducerPAT::GetMass(const reco::TrackRef & track){
    //   float P = track->p();
    //   float C = 2.625;
    //   float K = 2.495;
    //   float I = getEnergyLoss(track);
    //   return sqrt((I-C)/K)*P;
    // }


    /// define this as a plug-in
    DEFINE_FWK_MODULE(TrakTrakProducerPAT);

    // rsync -vut --existing src/TrakTrakProducerPAT.cc semrat@lxplus.cern.ch:/afs/cern.ch/user/s/semrat/scratch0/CMSSW_5_3_22/src/X4140/TrakTrakProducerPAT/src/TrakTrakProducerPAT.cc
