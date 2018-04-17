// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "../interface/X4140FilterPAT.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"


//
// class declaration
//

class X4140FilterPAT : public edm::EDProducer {
   public:
      explicit X4140FilterPAT(const edm::ParameterSet&);
      ~X4140FilterPAT();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
X4140FilterPAT::X4140FilterPAT(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


X4140FilterPAT::~X4140FilterPAT()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
X4140FilterPAT::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  Handle< vector<pat::Muon> > thePATMuonHandle;
  iEvent.getByLabel("cleanPatMuons", thePATMuonHandle);

  Handle< vector<pat::GenericParticle> > thePATTrackHandle;
  iEvent.getByLabel("cleanPatTrackCands", thePATTrackHandle);

  if (thePATTrackHandle->size() > 10000) return false;
  if ((thePATMuonHandle->size()) * (thePATTrackHandle->size()) > 20000) return false;

  if (thePATMuonHandle->size()>=1) {
    for ( std::vector<pat::Muon>::const_iterator iMuonP = thePATMuonHandle->begin();
          iMuonP != thePATMuonHandle->end(); ++iMuonP ) {
      //check for mu+ first
      if (iMuonP->charge() != 1) continue;
      TrackRef muTrackP = iMuonP->track();
      if ( muTrackP.isNull() ) {
        cout << "continue due to no track ref" << endl;
        continue;
      }
      const reco::Muon* mup = dynamic_cast<const reco::Muon * >(iMuonP->originalObject());

      //next check for mu-
      for ( std::vector<pat::Muon>::const_iterator iMuonM = thePATMuonHandle->begin();
            iMuonM != thePATMuonHandle->end(); ++iMuonM ) {
        if (iMuonM->charge() != -1) continue;
        TrackRef muTrackM = iMuonM->track();
        if ( muTrackM.isNull() ) {
          cout << "continue from no track ref" << endl;
          continue;
        }
        const reco::Muon* mum = dynamic_cast<const reco::Muon * >(iMuonM->originalObject());
        if(!muon::overlap(*mup, *mum)) {
          cout << "\n================================= NEW EVENT PASSING THE FILTER ========================================\n" << endl;
          return true;
        }
      }// 2nd loop over muons
    }//1st loop over muons

  }//if two muons
   return false;
}

// ------------ method called once each job just before starting event loop  ------------
void
X4140FilterPAT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
X4140FilterPAT::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool
X4140FilterPAT::beginRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool
X4140FilterPAT::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool
X4140FilterPAT::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool
X4140FilterPAT::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
X4140FilterPAT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(X4140FilterPAT);
