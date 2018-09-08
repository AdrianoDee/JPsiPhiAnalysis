// #include "../interface/DiMuonDiTrakProducerAOD.h"
//
// float DiMuonDiTrakProducerAOD::DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
// {
//    float p1 = t1.phi();
//    float p2 = t2.phi();
//    float e1 = t1.eta();
//    float e2 = t2.eta();
//    auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);
//
//    return sqrt((e1-e2)*(e1-e2) + dp*dp);
// }
//
// bool DiMuonDiTrakProducerAOD::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
// {
//   return (fabs(t1.pt()-t2.pt())/t2.pt()<maxDPtRel &&
// 	DeltaR(t1,t2) < maxDeltaR);
// }
//
//
// DiMuonDiTrakProducerAOD::DiMuonDiTrakProducerAOD(const edm::ParameterSet& iConfig):
//   DiMuonCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiMuon"))),
//   TrakCollection_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
//   TriggerEvent_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
//   triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
//   DiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonMassCuts")),
//   TrakTrakMassCuts_(iConfig.getParameter<std::vector<double>>("TrakTrakMassCuts")),
//   DiMuonDiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonDiTrakMassCuts")),
//   MassTraks_(iConfig.getParameter<std::vector<double>>("MassTraks")),
//   OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
//   product_name_(iConfig.getParameter<std::string>("Product")),
//   HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters")),
//   isMC_(iConfig.getParameter<bool>("IsMC"))
// {
//   produces<pat::CompositeCandidateCollection>(product_name_);
//   candidates = 0;
//   nevents = 0;
//   ndimuon = 0;
//   nreco = 0;
//
//   maxDeltaR = 0.01;
//   maxDPtRel = 2.0;
//
// }
//
// void DiMuonDiTrakProducerAOD::produce(edm::Event& iEvent, const edm::EventSetup& esetup){
//
//   std::unique_ptr<pat::CompositeCandidateCollection> DiMuonTTCandColl(new pat::CompositeCandidateCollection);
//
//   edm::Handle<pat::CompositeCandidateCollection> dimuon;
//   iEvent.getByToken(DiMuonCollection_,dimuon);
//
//   edm::Handle<std::vector<pat::PackedCandidate> > trak;
//   iEvent.getByToken(TrakCollection_,trak);
//
//   edm::Handle<trigger::TriggerEvent> triggerEvent;
//   iEvent.getByToken(TriggerEvent_,triggerEvent);
//
//   edm::Handle< edm::TriggerResults > triggerResults_handle;
//   iEvent.getByToken( triggerResults_Label , triggerResults_handle);
//
//   const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );
//
//   uint ncombo = 0;
//   float DiMuonMassMax_ = DiMuonMassCuts_[1];
//   float DiMuonMassMin_ = DiMuonMassCuts_[0];
//   float TrakTrakMassMax_ = TrakTrakMassCuts_[1];
//   float TrakTrakMassMin_ = TrakTrakMassCuts_[0];
//   float DiMuonDiTrakMassMax_ = DiMuonDiTrakMassCuts_[1];
//   float DiMuonDiTrakMassMin_ = DiMuonDiTrakMassCuts_[0];
//
//   const trigger::size_type nFilters(triggerEvent->sizeFilters());
//   const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
//
//   pat::TriggerObjectCollection filteredColl,matchedColl;
//   std::vector < UInt_t > iFilterKeys, filterResults,filters;
//   std::vector< std::vector <size_t> > filterIndices;
//
//   for (size_t i = 0; i < HLTFilters_.size()
//   {
//     for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter)
//     {
//     std::string filterTag = triggerEvent->filterTag(iFilter).encode();
//   }
//
//   }
//
//     for (size_t i = 0; i < HLTFilters_.size(); i++)
//       //if(unPackedTrigger.hasFilterLabel(HLTFilters_[i]))
//       if ( ( std::string(filterTag) == HLTFilters_[i]) )
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
//
//
//     {
//       trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
//       const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
//
//       for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey)
//       {
//         trigger::size_type objKey = objectKeys.at(iKey);
//         const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);
//
//         HLTObjCand hltObj;
//
//         hltObj.filterTag = filterTag;
//
//         hltObj.pt  = triggerObj.pt();
//         hltObj.eta = triggerObj.eta();
//         hltObj.phi = triggerObj.phi();
//         hltObj.id  = triggerObj.id();
//
//
//
//         if (isTag)       event_.hltTag.objects.push_back(hltObj);
//         else             event_.hlt   .objects.push_back(hltObj);
//       }
//     }
//   }
//
//   for ( size_t iTrigObj = 0; iTrigObj < triggerEvent->size(); ++iTrigObj ) {
//
//     pat::TriggerObjectStandAlone unPackedTrigger( triggerEvent->at( iTrigObj ) );
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
//   for (std::vector<pat::PackedCandidate>::const_iterator t = trak->begin(), trakend=trak->end(); t!= trakend; ++t)
//   {
//     bool matched = false;
//     for (std::vector<pat::TriggerObjectStandAlone>::const_iterator trigger = filteredColl.begin(), triggerEnd=filteredColl.end(); trigger!= triggerEnd; ++trigger)
//   for ( size_t iTrigObj = 0; iTrigObj < filteredColl.size(); ++iTrigObj )
//     {
//       auto thisTrig = filteredColl.at(iTrigObj);
//       if(MatchByDRDPt(*t,filteredColl[iTrigObj]))
//       {
//         if(matched)
//         {
//           if(DeltaR(*t,matchedColl.back()) > DeltaR(*t,thisTrig))
//           {
//             filters.pop_back();
//             filters.push_back(filterResults[iTrigObj]);
//             matchedColl.pop_back();
//             matchedColl.push_back(thisTrig);
//
//           }
//         }
//
//         if(!matched)
//           {
//             filters.push_back(filterResults[iTrigObj]);
//             matchedColl.push_back(thisTrig);
//           }
//
//         matched = true;
//       }
//     }
//
//     if (!matched)
//       filters.push_back(0);
//   }
//
// // Note: Dimuon cand are sorted by decreasing vertex probability then first is associated with "best" dimuon
//   for (pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuon->begin(); dimuonCand != dimuon->end(); ++dimuonCand){
//      if ( dimuonCand->mass() < DiMuonMassMax_  && dimuonCand->mass() > DiMuonMassMin_ ) {
//        const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon1"));
//        const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon2"));
//
// // loop on track candidates, make DiMuonT candidate, positive charge
//        // for (std::vector<pat::PackedCandidate>::const_iterator posTrack = trak->begin(), trakend=trak->end(); posTrack!= trakend; ++posTrack){
//        for (size_t i = 0; i < trak->size(); i++) {
//          auto posTrack = trak->at(i);
//
//          if(posTrack.charge()==0) continue;
//          if(posTrack.pt()<0.5) continue;
// 	       if(!isMC_ and fabs(posTrack.pdgId())!=211) continue;
// 	       if(!(posTrack.trackHighPurity())) continue;
//
//          if ( IsTheSame(posTrack,*pmu1) || IsTheSame(posTrack,*pmu2) || posTrack.charge() < 0 ) continue;
//
// // loop over second track candidate, negative charge
//          // for (std::vector<pat::PackedCandidate>::const_iterator negTrack = trak->begin(); negTrack!= trakend; ++negTrack){
//          for (size_t j = 0; j < trak->size(); j++) {
//            auto negTrack = trak->at(j);
//
//            if(negTrack.charge()==0) continue;
//            if(negTrack.pt()<0.5) continue;
//
//   	       if(!isMC_ and fabs(negTrack.pdgId())!=211) continue;
//   	       if(!(negTrack.trackHighPurity())) continue;
//
//            if (i == j) continue;
//            if ( IsTheSame(negTrack,*pmu1) || IsTheSame(negTrack,*pmu2) || negTrack.charge() > 0 ) continue;
//
//            pat::CompositeCandidate TTCand = makeTTCandidate(posTrack, negTrack);
//
//            if ( TTCand.mass() < TrakTrakMassMax_ && TTCand.mass() > TrakTrakMassMin_ ) {
//
//            pat::CompositeCandidate DiMuonTTCand = makeDiMuonTTCandidate(*dimuonCand, *&TTCand);
//
//            DiMuonTTCand.addUserInt("tPMatch",filters[i]);
//            DiMuonTTCand.addUserInt("tNMatch",filters[j]);
//
//            if ( DiMuonTTCand.mass() < DiMuonDiTrakMassMax_ && DiMuonTTCand.mass() > DiMuonDiTrakMassMin_) {
//
//              DiMuonTTCandColl->push_back(DiMuonTTCand);
//              candidates++;
//              ncombo++;
//            }
//         }
//
//          }
//          } // loop over second track
//        }   // loop on track candidates
//        if (OnlyBest_) break;
//      }
//
//   if ( ncombo != DiMuonTTCandColl->size() ) std::cout <<"ncombo ("<<ncombo<< ") != DiMuonTT ("<<DiMuonTTCandColl->size()<<")"<< std::endl;
//   if ( !dimuon->empty() )  ndimuon++;
//   if ( ncombo > 0 ) nreco++;
//   iEvent.put(std::move(DiMuonTTCandColl),product_name_);
//   nevents++;
// }
//
// void DiMuonDiTrakProducerAOD::endJob(){
//   std::cout << "###########################" << std::endl;
//   std::cout << "DiMuonDiTrak Candidate producer report:" << std::endl;
//   std::cout << "###########################" << std::endl;
//   std::cout << "Found " << nevents << " Events" << std::endl;
//   std::cout << "Events with DiMuon candidates " << ndimuon << std::endl;
//   std::cout << "Events with DiMuonDiTrak candidates " << nreco << std::endl;
//   std::cout << "###########################" << std::endl;
//   std::cout << "Found " << candidates << " DiMuonDiTrak candidates." << std::endl;
//   std::cout << "###########################" << std::endl;
// }
//
// bool DiMuonDiTrakProducerAOD::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
//   double DeltaEta = fabs(mu.eta()-tk.eta());
//   double DeltaP   = fabs(mu.p()-tk.p());
//   if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
//   return false;
// }
//
// pat::CompositeCandidate DiMuonDiTrakProducerAOD::makeDiMuonTTCandidate(
//                                           const pat::CompositeCandidate& dimuon,
// 				          const pat::CompositeCandidate& tt
//                                          ){
//
//   pat::CompositeCandidate DiMuonTCand;
//   DiMuonTCand.addDaughter(dimuon,"dimuon");
//   DiMuonTCand.addDaughter(tt,"ditrak");
//   DiMuonTCand.setVertex(dimuon.vertex());
//   DiMuonTCand.setCharge(tt.charge());
//
//   reco::Candidate::LorentzVector vDiMuonT = dimuon.p4() + tt.p4();
//   DiMuonTCand.setP4(vDiMuonT);
//
//   return DiMuonTCand;
//
// }
//
// pat::CompositeCandidate DiMuonDiTrakProducerAOD::makeTTCandidate(
//                                           const pat::PackedCandidate& trakP,
//                                           const pat::PackedCandidate& trakN
//                                          ){
//
//   pat::CompositeCandidate TTCand;
//   TTCand.addDaughter(trakP,"trakP");
//   TTCand.addDaughter(trakN,"trakN");
//   TTCand.setCharge(trakP.charge()+trakN.charge());
//
//   double m_kaon1 = MassTraks_[0];
//   math::XYZVector mom_kaon1 = trakP.momentum();
//   double e_kaon1 = sqrt(m_kaon1*m_kaon1 + mom_kaon1.Mag2());
//   math::XYZTLorentzVector p4_kaon1 = math::XYZTLorentzVector(mom_kaon1.X(),mom_kaon1.Y(),mom_kaon1.Z(),e_kaon1);
//   double m_kaon2 = MassTraks_[1];
//   math::XYZVector mom_kaon2 = trakN.momentum();
//   double e_kaon2 = sqrt(m_kaon2*m_kaon2 + mom_kaon2.Mag2());
//   math::XYZTLorentzVector p4_kaon2 = math::XYZTLorentzVector(mom_kaon2.X(),mom_kaon2.Y(),mom_kaon2.Z(),e_kaon2);
//   reco::Candidate::LorentzVector vTT = p4_kaon1 + p4_kaon2;
//   TTCand.setP4(vTT);
//
//   return TTCand;
// }
//
//
// reco::Candidate::LorentzVector DiMuonDiTrakProducerAOD::convertVector(const math::XYZTLorentzVectorF& v){
//
//   return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
// }
// //define this as a plug-in
// DEFINE_FWK_MODULE(DiMuonDiTrakProducerAOD);
