// system include files
#include <memory>

// FW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// DataFormat includes
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include <CommonTools/UtilAlgos/interface/StringCutObjectSelector.h>

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TMath.h"

template<typename T>
struct EqualByPosition {
  bool operator() ( T  v1, T  v2 ) {
    return (v1.x() == v2.x() && v1.y() == v2.y() && v1.z() == v2.z());
  }
};

class TrackFilter : public edm::EDProducer {
 public:
  explicit TrackFilter(const edm::ParameterSet&);
  ~TrackFilter() override {};
 private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  int GetPVFromOnia(edm::Event&);
 private:
  edm::EDGetTokenT<std::vector<pat::GenericParticle>> theTracks_;
  edm::EDGetTokenT<reco::VertexCollection> thePVs_;
  edm::EDGetTokenT<pat::CompositeCandidateCollection> theOnias_;
  StringCutObjectSelector<pat::GenericParticle> TrackSelection_;
  EqualByPosition<reco::Vertex> vertexComparator_;
};

TrackFilter::TrackFilter(const edm::ParameterSet& iConfig):
  theTracks_(consumes<pat::GenericParticleCollection>(iConfig.getParameter<edm::InputTag>("TrackTag"))),
  thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertexTag"))),
  theOnias_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("OniaTag"))),
  TrackSelection_(iConfig.existsAs<std::string>("TrackSelection") ? iConfig.getParameter<std::string>("TrackSelection") : "")
{  
    produces<pat::GenericParticleCollection>();  
}

int TrackFilter::GetPVFromOnia(edm::Event& iEvent) {
  int ivertex = -1;
  edm::Handle<pat::CompositeCandidateCollection> onias_;
  iEvent.getByToken(theOnias_, onias_);
  edm::Handle<reco::VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);
  //std::cout << "TrackFilter::GetPVFromOnia: onia size " << onias_->size() << " pv size " << priVtxs->size() << std::endl;
  if (onias_.isValid() && !onias_->empty() && !priVtxs->empty()) {
    const pat::CompositeCandidate *ionia = &(onias_->at(0));
    if (ionia) {
      const reco::Vertex *ipv = ionia->userData<reco::Vertex>("PVwithmuons");
      ivertex = 0;               // if not vertex available in onia the defakt one is used
      if (ipv) {
        //std::cout << "TrackFilter::GetPVFromOnia: finding pv, from PVwithmuons " << ipv->x() << " " << ipv->y() << " " << ipv->z() << std::endl;
        int index_v = -1;
        for (std::vector<reco::Vertex>::const_iterator iv = priVtxs->begin(), ivend = priVtxs->end(); iv != ivend; ++iv) {
          index_v++;
          if (vertexComparator_(*iv,*ipv)) {
            //std::cout << index_v << ": " << iv->x() << " " << iv->y() << " " << iv->z() << std::endl;
	    ivertex = index_v;
	    break;
          }  	
        } 
	if (index_v < 0) std::cout << "TrackFilter::GetPVFromOnia: *** non matching PV to Onia PV found in PVCollection, using iPV=0" << std::endl;
      } else std::cout << "TrackFilter::GetPVFromOnia: *** no PV associate to Onia object found in composite candidate, using iPV=0" << std::endl;
    } //else std::cout << "TrackFilter::GetPVFromOnia: *** no Onia combination pass selection" << std::endl; 
  } else { if (!priVtxs->empty()) ivertex = 0; }  // in case no onia is present we will use the PV==0, so a filter has to be applied before
  return ivertex;
}

// ------------ method called to produce the data  ------------
void TrackFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {  

  edm::Handle<pat::GenericParticleCollection> kaons;
  iEvent.getByToken(theTracks_,kaons);

  edm::Handle<reco::VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  
  int npv = priVtxs->size(); 
  int which_vertex = GetPVFromOnia(iEvent);
  if (which_vertex < 0)   std::cout << "TrackFilter::GetPVFromOnia: *** Event does not have enough data to continue == -1" << std::endl; 
  if (which_vertex > npv) std::cout << "TrackFilter::GetPVFromOnia: *** Onia vertex index above bondaries " << which_vertex << " " << npv << std::endl;

  int ntracks = kaons->size();
  std::unique_ptr<pat::GenericParticleCollection> filetered(new pat::GenericParticleCollection);

  if (ntracks>0 && npv > 0 && which_vertex >= 0) {
     for (size_t i=0; i<(size_t)ntracks; i++) {
	double dzmin = 99999.;
	int    j_pv   = -1;
	pat::GenericParticle tk = kaons->at(i);
        if (!TrackSelection_(tk)) continue;
	reco::TransientTrack t_tk = theTTBuilder->build(tk.track());
	for (size_t j=0;j<(size_t)npv; j++) {
	    reco::Vertex t_pv = priVtxs->at(j);
	    GlobalPoint vert(t_pv.x(), t_pv.y(), t_pv.z());
	    TrajectoryStateClosestToPoint traj = t_tk.trajectoryStateClosestToPoint(vert);
	    double dz = TMath::Abs(traj.perigeeParameters().longitudinalImpactParameter());
	    if (dz < dzmin) {
	       dzmin = dz;
	       j_pv  = j; 
	    } 	
        }
	if (j_pv>-1 && j_pv == which_vertex) filetered->push_back(tk);  
     }
  }
  iEvent.put(std::move(filetered));
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackFilter);
