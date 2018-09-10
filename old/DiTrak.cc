// #include "../interface/DiTrak.h"
//
// DiTrakPAT::DiTrakPAT(const edm::ParameterSet& iConfig):
// traks_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Traks"))),
// TriggerCollection_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
// thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpot"))),
// thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"))),
// ditrakSelection_(iConfig.existsAs<std::string>("DiTrakCuts") ? iConfig.getParameter<std::string>("DiTrakCuts") : ""),
// massTraks_(iConfig.getParameter<std::vector<double>>("TraksMasses"))
// {
//   produces<pat::CompositeCandidateCollection>();
// }
//
//
// DiTrakPAT::~DiTrakPAT()
// {
//
//   // do anything here that needs to be done at desctruction time
//   // (e.g. close files, deallocate resources etc.)
//
// }
//
// const pat::CompositeCandidate DiTrakPAT::makeTTCandidate(
//                                           const pat::PackedCandidate& trakP,
//                                           const pat::PackedCandidate& trakN
//                                          ){
//
//   pat::CompositeCandidate trktrkcand;
//   trktrkcand.addDaughter(trakP,"trakP");
//   trktrkcand.addDaughter(trakN,"trakN");
//   trktrkcand.setCharge(trakP.charge()+trakN.charge());
//
//   double m_trakP = massTraks_[0];
//   math::XYZVector mom_trakP = trakP.momentum();
//   double e_trakP = sqrt(m_trakP*m_trakP + mom_trakP.Mag2());
//   math::XYZTLorentzVector p4_trakP = math::XYZTLorentzVector(mom_trakP.X(),mom_trakP.Y(),mom_trakP.Z(),e_trakP);
//   double m_trakN = massTraks_[1];
//   math::XYZVector mom_trakN = trakN.momentum();
//   double e_trakN = sqrt(m_trakN*m_trakN + mom_trakN.Mag2());
//   math::XYZTLorentzVector p4_trakN = math::XYZTLorentzVector(mom_trakN.X(),mom_trakN.Y(),mom_trakN.Z(),e_trakN);
//   reco::Candidate::LorentzVector vTT = p4_trakP + p4_trakN;
//   trktrkcand.setP4(vTT);
//
//   return trktrkcand;
// }
//
//
// void
// DiTrakPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
// {
//   using namespace edm;
//   using namespace std;
//   using namespace reco;
//   typedef Candidate::LorentzVector LorentzVector;
//
//   std::unique_ptr<pat::CompositeCandidateCollection> trakCollection(new pat::CompositeCandidateCollection);
//
//   edm::Handle<std::vector<pat::PackedCandidate> > traks;
//   iEvent.getByToken(traks_,traks);
//
//   edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerColl;
//   iEvent.getByToken(TriggerCollection_,triggerColl);
//
//   pat::TriggerObjectStandAloneCollection filteredColl;
//
//   Vertex thePrimaryV;
//
//   ESHandle<MagneticField> magneticField;
//   iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
//   const MagneticField* field = magneticField.product();
//
//   edm::Handle<BeamSpot> theBeamSpot;
//   iEvent.getByToken(thebeamspot_,theBeamSpot);
//   BeamSpot bs = *theBeamSpot;
//
//   edm::Handle<VertexCollection> priVtxs;
//   iEvent.getByToken(thePVs_, priVtxs);
//   if ( priVtxs->begin() != priVtxs->end() ) {
//     thePrimaryV = Vertex(*(priVtxs->begin()));
//   }
//   else {
//     thePrimaryV = Vertex(bs.position(), bs.covariance3D());
//   }
//
//   edm::Handle< edm::TriggerResults > triggerResults_handle;
//   iEvent.getByToken( triggerResults_Label , triggerResults_handle);
//
//   const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );
//
//   for ( size_t iTrigObj = 0; iTrigObj < triggerColl->size(); ++iTrigObj ) {
//
//     pat::TriggerObjectStandAlone unPackedTrigger( triggerColl->at( iTrigObj ) );
//
//     unPackedTrigger.unpackPathNames( names );
//     unPackedTrigger.unpackFilterLabels(iEvent,*triggerResults_handle);
//
//     bool filtered = false;
//     UInt_t thisFilter = 0;
//
//     for (size_t i = 0; i < HLTFilters_.size(); i++)
//       if(unPackedTrigger.hasFilterLabel(HLTFilters_[i]))
//         {
//           thisFilter += (1<<i);
//           filtered = true;
//         }
//
//     if(filtered)
//     {
//       filteredColl.push_back(unPackedTrigger);
//       filterResults.push_back(thisFilter);
//     }
//   }
//
//   for (std::vector<pat::PackedCandidate>::const_iterator trak = trakColl->begin(), trakend=trakColl->end(); trak!= trakend; ++trak)
//   {
//     bool matched = false;
//     for (std::vector<pat::TriggerObjectStandAlone>::const_iterator trigger = filteredColl.begin(), triggerEnd=filteredColl.end(); trigger!= triggerEnd; ++trigger)
//     {
//       if(MatchByDRDPt(*trak,*trigger))
//       {
//         if(matched)
//         {
//           if(DeltaR(*trak,matchedColl.back()) > DeltaR(*trak,*trigger))
//           {
//             matchedColl.pop_back();
//             matchedColl.push_back(*trigger);
//
//           }
//         }
//
//         if(!matched)
//           {
//             filteredTracks.push_back(*trak);
//             matchedColl.push_back(*trigger);
//           }
//
//         matched = true;
//       }
//     }
//   }
//
//
//
//   edm::ESHandle<TransientTrackBuilder> theTTBuilder;
//   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
//   KalmanVertexFitter vtxFitter(true);
//
//   // ParticleMass trakP_mass = massTraks_[0];
//   // ParticleMass trakN_mass = massTraks_[1];
//
//   vector<double> ttMasses;
//   ttMasses.push_back(massTraks_[0]);
//   ttMasses.push_back(massTraks_[1]);
//
//   // float trakP_sigma = trakP_mass*1.e-6;
//   // float trakN_sigma = trakN_mass*1.e-6;
//
//   float vProb, vNDF, vChi2, minDz = 999999.;
//   float cosAlpha, ctauPV, ctauErrPV, dca;
//   float l_xy, lErr_xy;
//   for (size_t i = 0; i < traks->size(); i++)
//   {
//     auto posTrack = traks->at(i);
//
//     if(posTrack.charge() <= 0 ) continue;
//     if(posTrack.pt()<0.5) continue;
//     if(fabs(posTrack.pdgId())!=211) continue;
//     if(!posTrack.hasTrackDetails()) continue;
//
//     for (size_t j = 0; j < traks->size(); j++){
//
//       vProb = -1.0; vNDF = -1.0; vChi2 = -1.0;
//       cosAlpha = -1.0; ctauPV = -1.0; ctauErrPV = -1.0;
//       dca = -1.0; minDz = 999999.; dca = 1E20;
//
//       if (i == j) continue;
//
//       auto negTrack = traks->at(j);
//
//       if(negTrack.charge() >= 0 ) continue;
//       if(negTrack.pt()<0.5) continue;
//       if(fabs(negTrack.pdgId())!=211) continue;
//       if(!negTrack.hasTrackDetails()) continue;
//
//       pat::CompositeCandidate trktrkcand = makeTTCandidate(posTrack,negTrack);
//       vector<TransientVertex> pvs;
//
//       if(!ditrakSelection_(trktrkcand)) continue;
//
//       vector<TransientTrack> tt_ttks;
//       tt_ttks.push_back(theTTBuilder->build(negTrack.bestTrack()));  // pass the reco::Track, not  the reco::TrackRef (which can be transient)
//       tt_ttks.push_back(theTTBuilder->build(posTrack.bestTrack()));
//
//       TransientVertex ttVertex = vtxFitter.vertex(tt_ttks);
//       CachingVertex<5> VtxForInvMass = vtxFitter.vertex( tt_ttks );
//
//       LorentzVector trktrk = posTrack.p4() + negTrack.p4();
//
//       Measurement1D MassWErr(posTrack.mass(),-9999.);
//       if ( field->nominalValue() > 0 )
//           MassWErr = massCalculator.invariantMass( VtxForInvMass, ttMasses );
//       else
//           ttVertex = TransientVertex();                      // with no arguments it is invalid
//
//       if (!(ttVertex.isValid()))
//           continue;
//
//       vChi2 = ttVertex.totalChiSquared();
//       vNDF  = ttVertex.degreesOfFreedom();
//       vProb = TMath::Prob(vChi2,(int)vNDF);
//
//       //Vertex parameters
//       TVector3 vtx,vtx3D;
//       TVector3 pvtx,pvtx3D;
//       VertexDistanceXY vdistXY;
//
//       vtx.SetXYZ(ttVertex.position().x(),ttVertex.position().y(),0);
//       vtx3D.SetXYZ(ttVertex.position().x(),ttVertex.position().y(),ttVertex.position().z());
//       TVector3 pperp(trktrk.px(), trktrk.py(), 0);
//       TVector3 pperp3D(trktrk.px(), trktrk.py(), trktrk.pz());
//       AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
//       AlgebraicVector3 vpperp3D(pperp.x(),pperp.y(),pperp.z());
//
//       //Resolving pileup ambiguity with two trak min distance
//       TwoTrackMinimumDistance ttmd;
//       bool status = ttmd.calculate( GlobalTrajectoryParameters(
//         GlobalPoint(ttVertex.position().x(), ttVertex.position().y(), ttVertex.position().z()),
//         GlobalVector(trktrkcand.px(),trktrkcand.py(),trktrkcand.pz()),TrackCharge(0),&(*magneticField)),
//         GlobalTrajectoryParameters(
//           GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
//           GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(0),&(*magneticField)));
//
//       float extrapZ=-9E20;
//
//       if (status) extrapZ=ttmd.points().first.z();
//
//       for(VertexCollection::const_iterator itv = priVtxs->begin(), itvend = priVtxs->end(); itv != itvend; ++itv)
//       {
//           float deltaZ = fabs(extrapZ - itv->position().z()) ;
//           if ( deltaZ < minDz ) {
//               minDz = deltaZ;
//               thePrimaryV = Vertex(*itv);
//             }
//         }
//
//       //Distance of Closest Approach
//       TrajectoryStateClosestToPoint mu1TS = tt_ttks[0].impactPointTSCP();
//       TrajectoryStateClosestToPoint mu2TS = tt_ttks[1].impactPointTSCP();
//
//       if (mu1TS.isValid() && mu2TS.isValid()) {
//         ClosestApproachInRPhi cApp;
//         cApp.calculate(mu1TS.theState(), mu2TS.theState());
//         if (cApp.status() ) dca = cApp.distance();
//       }
//
//       //Lifetime calculations
//       pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
//       TVector3 vdiff = vtx - pvtx;
//       cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
//
//       Measurement1D distXY = vdistXY.distance(Vertex(ttVertex), thePrimaryV);
//       ctauPV = distXY.value()*cosAlpha * trktrkcand.mass()/pperp.Perp();
//
//       GlobalError v1e = (Vertex(ttVertex)).error();
//       GlobalError v2e = thePrimaryV.error();
//       AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
//       ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*trktrkcand.mass()/(pperp.Perp2());
//
//       AlgebraicVector3 vDiff;
//       vDiff[0] = vdiff.x(); vDiff[1] = vdiff.y(); vDiff[2] = 0 ;
//       l_xy = vdiff.Perp();
//       lErr_xy = sqrt(ROOT::Math::Similarity(vDiff,vXYe)) / vdiff.Perp();
//
//
//       trktrkcand.addUserFloat("vNChi2",vChi2/vNDF);
//       trktrkcand.addUserFloat("vProb",vProb);
//       trktrkcand.addUserFloat("DCA", dca );
//       trktrkcand.addUserFloat("MassErr",MassWErr.error());
//       trktrkcand.addUserFloat("ctauPV",ctauPV);
//       trktrkcand.addUserFloat("ctauErrPV",ctauErrPV);
//       trktrkcand.addUserFloat("lxy",l_xy);
//       trktrkcand.addUserFloat("lErrxy",lErr_xy);
//       trktrkcand.addUserFloat("cosAlpha",cosAlpha);
//       trktrkcand.addUserData("thePV",Vertex(thePrimaryV));
//       trktrkcand.addUserData("theVertex",Vertex(ttVertex));
//
//       trakCollection->push_back(trktrkcand);
//
//
//     } // loop over second track
//   }
//
//   std::sort(trakCollection->begin(),trakCollection->end(),vPComparator_);
//   iEvent.put(std::move(trakCollection));
//
// }
//
//
// // ------------ method called once each job just before starting event loop  ------------
// void
// DiTrakPAT::beginJob()
// {
// }
//
// // ------------ method called once each job just after ending the event loop  ------------
// void
// DiTrakPAT::endJob() {
// }
//
// //define this as a plug-in
// DEFINE_FWK_MODULE(DiTrakPAT);
