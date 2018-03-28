//#include "../interface/DiMuonDiTrakML.h"
#include "../interface/DiMuonVtxReProducer.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackerSingleRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include <TH2F.h>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"


DiMuonDiTrakML::DiMuonDiTrakML(const edm::ParameterSet& iConfig):
muons_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("Muons"))),
traks_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("Tracks"))),
thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
DiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonCuts")),
DiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiTrakCuts"))
// thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
// thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
// DiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonCuts")),
// DiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiTrakCuts")),
// DiMuonDiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonDiTrakCuts")),
// massCands_(iConfig.getParameter<std::vector<double>>("CandsMasses"))
{
  // produces<pat::CompositeCandidateCollection>();
  muon_mass = 0.1056583715;
}


DiMuonDiTrakML::~DiMuonDiTrakML()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


void
DiMuonDiTrakML::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::vector<int> pixelDets{0,1,2,3,14,15,16,29,30,31};

  std::vector<double> mmMasses;
  mmMasses.push_back( 0.1056583715 );
  mmMasses.push_back( 0.1056583715 );

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(traks_,tracks);

  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muons_,muons);

  reco::Vertex thePrimaryV;

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  const MagneticField* field = magneticField.product();

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;

  edm::Handle<reco::VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);
  if ( priVtxs->begin() != priVtxs->end() ) {
    thePrimaryV = reco::Vertex(*(priVtxs->begin()));
  }
  else {
    thePrimaryV = reco::Vertex(bs.position(), bs.covariance3D());
  }


  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  KalmanVertexFitter vtxFitter(true);

  int max = 0;

  float vProb, vNDF, vChi2, minDz = 999999.;
  float vProb_mumu, vNDF_mumu, vChi2_mumu;
  float cosAlpha, ctauPV, ctauErrPV, dca;
  float l_xy, lErr_xy;

  for(reco::MuonCollection::const_iterator mPos = muons->begin();mPos != muons->end(); ++mPos )
  {
    if(mPos->charge()<=0.0) continue;
    if (!(mPos->bestTrackRef().isNonnull())) continue;
    if (!(mPos->innerTrack().isNonnull())) continue;

    for(reco::MuonCollection::const_iterator mNeg = muons->begin();mNeg != muons->end(); ++mNeg )
    {
      if(mNeg->charge()>=0.0) continue;
      if (!(mNeg->bestTrackRef().isNonnull())) continue;
      if (!(mNeg->innerTrack().isNonnull())) continue;

      std::vector<TransientVertex> vDiMuon;

      // Candidate::LorentzVector mumu = mNeg->p4() + mPos->p4();

      TLorentzVector mu1, mu2,mumuP4;

      mu1.SetXYZM(mNeg->track()->px(),mNeg->track()->py(),mNeg->track()->pz(),muon_mass);
      mu2.SetXYZM(mPos->track()->px(),mPos->track()->py(),mPos->track()->pz(),muon_mass);

      mumuP4=mu1+mu2;
      // mumucand.setP4(mumu);
      // mumucand.setCharge(mNeg->charge()+mPos->charge());

      if(mumuP4.M() < DiMuonMassCuts_[1]) continue;
      if(mumuP4.M() > DiMuonMassCuts_[0]) continue;

      std::vector<reco::TransientTrack> mm_ttks;

      mm_ttks.push_back(theTTBuilder->build(mNeg->track()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
      mm_ttks.push_back(theTTBuilder->build(mPos->track()));

      TransientVertex mumuVertex = vtxFitter.vertex(mm_ttks);
      CachingVertex<5> VtxForInvMass = vtxFitter.vertex( mm_ttks );

      Measurement1D MassWErr(mPos->mass(),-9999.);
      if ( field->nominalValue() > 0 )
          MassWErr = massCalculator.invariantMass( VtxForInvMass, mmMasses );
      else
          mumuVertex = TransientVertex();                      // with no arguments it is invalid

      if (!(mumuVertex.isValid()))
          continue;

      vChi2_mumu = mumuVertex.totalChiSquared();
      vNDF_mumu  = mumuVertex.degreesOfFreedom();
      vProb_mumu = TMath::Prob(vChi2_mumu,(int)vNDF_mumu);

      if(vProb_mumu < 0.0) continue;

      for(reco::TrackCollection::const_iterator posTrack = tracks->begin();posTrack != tracks->end(); ++posTrack )
      {
        if(!(posTrack->extra())) continue;
        if(posTrack->charge()<=0.0) continue;
        // if(!(posTrack->isNonnull())) continue;

        for(reco::TrackCollection::const_iterator negTrack = tracks->begin();negTrack != tracks->end(); ++negTrack )
        {
          if(!(negTrack->extra())) continue;
          if(negTrack->charge()>=0.0) continue;
          // if(!(negTrack->isNonnull())) continue;


        }


      }


    }

  }

// 	int padHalfSize = 8;
// 	int padSize = padHalfSize*2;
//   int maxpix = 0;
//
//   for(reco::TrackCollection::const_iterator itTrack = tracks->begin();itTrack != tracks->end(); ++itTrack )
//     {
// 	int noPixels= 0,noStripOne = 0, noStripTwo = 0;
// 	int counter = 0;
// 	float clusterSize = 0.0;
// //	std::cout<<"On "<< itTrack->found() ;
// 	for ( trackingRecHit_iterator recHit = (*itTrack).recHitsBegin();recHit != (*itTrack).recHitsEnd(); ++recHit )
// 	{
// 		counter++;
// /*		 if(!(*recHit))
//                  continue;
//
//                  if (!((*recHit)->isValid()))
//                  continue;
//
//                  if(!((*recHit)->hasPositionAndError()))
//                  continue;
// */
// 		TrackerSingleRecHit const * hit= dynamic_cast<TrackerSingleRecHit const *>(*recHit);
// //    BaseTrackerRecHit const * bhit = dynamic_cast<BaseTrackerRecHit const *>(recHit);
//
// 		DetId detid = (*recHit)->geographicalId();
// 		unsigned int subdetid = detid.subdetId();
//
//
// 	        //if(!(siPix))
// 	        //continue;
// 	        //
//
//     	 if(detid.det() != DetId::Tracker) continue;
// //	 if (!((subdetid==1) || (subdetid==2))) continue;
// // 	 if()
// 		if (dynamic_cast<SiPixelRecHit const *>(hit))
// 		{		noPixels++;
// 				clusterSize += float(dynamic_cast<SiPixelRecHit const *>(hit)->cluster()->size());
// 				clusterSize /= float(counter);
//
// 		auto clust = dynamic_cast<SiPixelRecHit const *>(hit)->cluster();
// 		TH2F hClust("hClust","hClust",
//               padSize,
//               clust->x()-padHalfSize,
//               clust->x()+padHalfSize,
//               padSize,
//               clust->y()-padHalfSize,
//               clust->y()+padHalfSize);
//
// 		for (int nx = 0; nx < padSize; ++nx)
//               for (int ny = 0; ny < padSize; ++ny)
//               hClust.SetBinContent(nx,ny,0.0);
//
//               for (int k = 0; k < clust->size(); ++k)
//               hClust.SetBinContent(hClust.FindBin((float)clust->pixel(k).x, (float)clust->pixel(k).y),(float)clust->pixel(k).adc);
//
// 		for (int ny = padSize; ny>0; --ny)
//               {
//                 for(int nx = 0; nx<padSize; nx++)
//                 {
//                   int n = (ny+2)*(padSize + 2) - 2 -2 - nx - padSize; //see TH2 reference for clarification
//
//        //          std::cout << hClust.GetBinContent(n) << " ";
// 		}
//               }
// 	//	std::cout << std::endl;
//
// 		}
// 		if (dynamic_cast<SiStripRecHit1D const *>(hit))
// 		noStripOne++;
//
// 	if (dynamic_cast<SiStripRecHit2D const *>(hit))
// 			noStripTwo++;
//
// 	}
// //	std::cout << " n. pixels = " << noPixels<< " 1DStrips = " << noStripOne << " 2DStrips = " << noStripTwo<< " clustsize : "<< clusterSize <<std::endl;
//
// 	     max = std::max(max,int(itTrack->found()));
// 		maxpix = std::max(maxpix,noPixels);
// }
//   std::cout<<"Max = " << max<< " Max pixels " << maxpix << std::endl;
  //loop on
  // std::sort(mmttCollection->begin(),mmttCollection->end(),vPComparator_);
  // iEvent.put(std::move(mmttCollection));

}


// ------------ method called once each job just before starting event loop  ------------
void
DiMuonDiTrakML::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DiMuonDiTrakML::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakML);
