#include "../interface/DiMuonDiTrakProducer.h"

float DiMuonDiTrakProducer::DeltaR(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
   float p1 = t1.phi();
   float p2 = t2.phi();
   float e1 = t1.eta();
   float e2 = t2.eta();
   auto dp=std::abs(p1-p2); if (dp>float(M_PI)) dp-=float(2*M_PI);

   return sqrt((e1-e2)*(e1-e2) + dp*dp);
}

bool DiMuonDiTrakProducer::MatchByDRDPt(const pat::PackedCandidate t1, const pat::TriggerObjectStandAlone t2)
{
  return (fabs(t1.pt()-t2.pt())/t2.pt()<maxDPtRel &&
	DeltaR(t1,t2) < maxDeltaR);
}


std::pair<int, float>
DiMuonDiTrakProducer::findJpsiMCInfo(reco::GenParticleRef genJpsi) {

  // std::cout << "findJpsiMCInfo 1 " << std::endl;
  int momJpsiID = 0;
  float trueLife = -99.;

  if (genJpsi->numberOfMothers()>0) {

    // std::cout << "findJpsiMCInfo 1 " << std::endl;

    TVector3 trueVtx(0.0,0.0,0.0);
    TVector3 trueP(0.0,0.0,0.0);
    TVector3 trueVtxMom(0.0,0.0,0.0);

    trueVtx.SetXYZ(genJpsi->vertex().x(),genJpsi->vertex().y(),genJpsi->vertex().z());
    trueP.SetXYZ(genJpsi->momentum().x(),genJpsi->momentum().y(),genJpsi->momentum().z());

    bool aBhadron = false;
    reco::GenParticleRef Jpsimom = genJpsi->motherRef();       // find mothers
    // std::cout << "findJpsiMCInfo 1 " << std::endl;
    if (Jpsimom.isNull()) {
      std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
      return result;
    } else {
      reco::GenParticleRef Jpsigrandmom = Jpsimom->motherRef();
      if (isAbHadron(Jpsimom->pdgId())) {
        if (Jpsigrandmom.isNonnull() && isAMixedbHadron(Jpsimom->pdgId(),Jpsigrandmom->pdgId())) {
          momJpsiID = Jpsigrandmom->pdgId();
          trueVtxMom.SetXYZ(Jpsigrandmom->vertex().x(),Jpsigrandmom->vertex().y(),Jpsigrandmom->vertex().z());
        } else {
          momJpsiID = Jpsimom->pdgId();
          trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
        }
        aBhadron = true;
      } else {
        if (Jpsigrandmom.isNonnull() && isAbHadron(Jpsigrandmom->pdgId())) {
          reco::GenParticleRef JpsiGrandgrandmom = Jpsigrandmom->motherRef();
          if (JpsiGrandgrandmom.isNonnull() && isAMixedbHadron(Jpsigrandmom->pdgId(),JpsiGrandgrandmom->pdgId())) {
            momJpsiID = JpsiGrandgrandmom->pdgId();
            trueVtxMom.SetXYZ(JpsiGrandgrandmom->vertex().x(),JpsiGrandgrandmom->vertex().y(),JpsiGrandgrandmom->vertex().z());
          } else {
            momJpsiID = Jpsigrandmom->pdgId();
            trueVtxMom.SetXYZ(Jpsigrandmom->vertex().x(),Jpsigrandmom->vertex().y(),Jpsigrandmom->vertex().z());
          }
          aBhadron = true;
        }
      }
      if (!aBhadron) {
        momJpsiID = Jpsimom->pdgId();
        trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
      }
    }

    TVector3 vdiff = trueVtx - trueVtxMom;
    //trueLife = vdiff.Perp()*3.09688/trueP.Perp();
    trueLife = vdiff.Perp()*genJpsi->mass()/trueP.Perp();
  }
  std::pair<int, float> result = std::make_pair(momJpsiID, trueLife);
  return result;

}

DiMuonDiTrakProducer::DiMuonDiTrakProducer(const edm::ParameterSet& iConfig):
  DiMuonCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DiMuon"))),
  TrakCollection_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
  TriggerCollection_(consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  DiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonMassCuts")),
  TrakTrakMassCuts_(iConfig.getParameter<std::vector<double>>("TrakTrakMassCuts")),
  DiMuonDiTrakMassCuts_(iConfig.getParameter<std::vector<double>>("DiMuonDiTrakMassCuts")),
  MassTraks_(iConfig.getParameter<std::vector<double>>("MassTraks")),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  product_name_(iConfig.getParameter<std::string>("Product")),
  HLTFilters_(iConfig.getParameter<std::vector<std::string>>("Filters")),
  isMC_(iConfig.getParameter<bool>("IsMC")),
  AddMCTruth_(iConfig.getParameter<bool>("AddMCTruth"))
{
  produces<pat::CompositeCandidateCollection>(product_name_);
  candidates = 0;
  nevents = 0;
  ndimuon = 0;
  nreco = 0;

  maxDeltaR = 0.01;
  maxDPtRel = 2.0;

}

void DiMuonDiTrakProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::unique_ptr<pat::CompositeCandidateCollection> DiMuonTTCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> dimuon;
  iEvent.getByToken(DiMuonCollection_,dimuon);

  edm::Handle<std::vector<pat::PackedCandidate> > trak;
  iEvent.getByToken(TrakCollection_,trak);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trig;
  iEvent.getByToken(TriggerCollection_,trig);

  edm::Handle< edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken( triggerResults_Label , triggerResults_handle);

  const edm::TriggerNames & names = iEvent.triggerNames( *triggerResults_handle );

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  uint ncombo = 0;
  float DiMuonMassMax_ = DiMuonMassCuts_[1];
  float DiMuonMassMin_ = DiMuonMassCuts_[0];
  float TrakTrakMassMax_ = TrakTrakMassCuts_[1];
  float TrakTrakMassMin_ = TrakTrakMassCuts_[0];
  float DiMuonDiTrakMassMax_ = DiMuonDiTrakMassCuts_[1];
  float DiMuonDiTrakMassMin_ = DiMuonDiTrakMassCuts_[0];

  pat::TriggerObjectStandAloneCollection filteredColl,matchedColl;
  std::vector < UInt_t > filterResults,filters;

  for ( size_t iTrigObj = 0; iTrigObj < trig->size(); ++iTrigObj ) {

    pat::TriggerObjectStandAlone unPackedTrigger( trig->at( iTrigObj ) );

    unPackedTrigger.unpackPathNames( names );
    unPackedTrigger.unpackFilterLabels(iEvent,*triggerResults_handle);

    bool filtered = false;
    UInt_t thisFilter = 0;

    for (size_t i = 0; i < HLTFilters_.size(); i++)
      if(unPackedTrigger.hasFilterLabel(HLTFilters_[i]))
        {
          thisFilter += (1<<i);
          filtered = true;
        }

    if(filtered)
    {
      filteredColl.push_back(unPackedTrigger);
      filterResults.push_back(thisFilter);
    }
  }

  for (std::vector<pat::PackedCandidate>::const_iterator t = trak->begin(), trakend=trak->end(); t!= trakend; ++t)
  {
    bool matched = false;
    for (std::vector<pat::TriggerObjectStandAlone>::const_iterator trigger = filteredColl.begin(), triggerEnd=filteredColl.end(); trigger!= triggerEnd; ++trigger)
  for ( size_t iTrigObj = 0; iTrigObj < filteredColl.size(); ++iTrigObj )
    {
      auto thisTrig = filteredColl.at(iTrigObj);
      if(MatchByDRDPt(*t,filteredColl[iTrigObj]))
      {
        if(matched)
        {
          if(DeltaR(*t,matchedColl.back()) > DeltaR(*t,thisTrig))
          {
            filters.pop_back();
            filters.push_back(filterResults[iTrigObj]);
            matchedColl.pop_back();
            matchedColl.push_back(thisTrig);

          }
        }

        if(!matched)
          {
            filters.push_back(filterResults[iTrigObj]);
            matchedColl.push_back(thisTrig);
          }

        matched = true;
      }
    }

    if (!matched)
      filters.push_back(0);
  }

// Note: Dimuon cand are sorted by decreasing vertex probability then first is associated with "best" dimuon
  for (pat::CompositeCandidateCollection::const_iterator dimuonCand = dimuon->begin(); dimuonCand != dimuon->end(); ++dimuonCand){
     if ( dimuonCand->mass() < DiMuonMassMax_  && dimuonCand->mass() > DiMuonMassMin_ ) {

       if(dimuonCand->userFloat("vProb")<0.0)
         continue;

       const reco::Vertex thePrimaryV = *dimuonCand->userData<reco::Vertex>("PVwithmuons");

       const pat::Muon *pmu1 = dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon1"));
       const pat::Muon *pmu2 = dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon2"));

// loop on track candidates, make DiMuonT candidate, positive charge
       // for (std::vector<pat::PackedCandidate>::const_iterator posTrack = trak->begin(), trakend=trak->end(); posTrack!= trakend; ++posTrack){
       for (size_t i = 0; i < trak->size(); i++) {
         auto posTrack = trak->at(i);

         if(posTrack.charge()==0) continue;
         if(posTrack.pt()<0.5) continue;
	       if(!isMC_ and fabs(posTrack.pdgId())!=211) continue;
	       if(!(posTrack.trackHighPurity())) continue;
         if(!(posTrack.hasTrackDetails())) continue;

         if ( IsTheSame(posTrack,*pmu1) || IsTheSame(posTrack,*pmu2) || posTrack.charge() < 0 ) continue;

// loop over second track candidate, negative charge
         // for (std::vector<pat::PackedCandidate>::const_iterator negTrack = trak->begin(); negTrack!= trakend; ++negTrack){
         for (size_t j = 0; j < trak->size(); j++) {
           auto negTrack = trak->at(j);

           if(negTrack.charge()==0) continue;
           if(negTrack.pt()<0.5) continue;

  	       if(!isMC_ and fabs(negTrack.pdgId())!=211) continue;
  	       if(!(negTrack.trackHighPurity())) continue;
           if(!(negTrack.hasTrackDetails())) continue;

           if (i == j) continue;
           if ( IsTheSame(negTrack,*pmu1) || IsTheSame(negTrack,*pmu2) || negTrack.charge() > 0 ) continue;

           pat::CompositeCandidate TTCand = makeTTCandidate(posTrack, negTrack);


           if ( !(TTCand.mass() < TrakTrakMassMax_ && TTCand.mass() > TrakTrakMassMin_) ) continue;

           pat::CompositeCandidate DiMuonTTCand = makeDiMuonTTCandidate(*dimuonCand, *&TTCand);

           if ( !(DiMuonTTCand.mass() < DiMuonDiTrakMassMax_ && DiMuonTTCand.mass() > DiMuonDiTrakMassMin_)) continue;

           float refittedMass = -1.0, mumuVtxCL = -1.0;

           const ParticleMass muonMass(0.1056583);
           float muonSigma = muonMass*1E-6;
           const ParticleMass trakMass1(MassTraks_[0]);
           float trakSigma1 = trakMass1*1E-6;
           const ParticleMass trakMass2(MassTraks_[1]);
           float trakSigma2 = trakMass2*1E-6;

           std::vector<reco::TransientTrack> xTracks;
           KinematicParticleFactoryFromTransientTrack pFactory;
           std::vector<RefCountedKinematicParticle> xParticles;

           float kinChi = 0.;
           float kinNdf = 0.;

           xTracks.push_back((*theB).build(*(pmu1->innerTrack()))); // K+
           xTracks.push_back((*theB).build(*(pmu2->innerTrack()))); // K+
           xTracks.push_back((*theB).build(*(posTrack.bestTrack()))); // K+
           xTracks.push_back((*theB).build(*(negTrack.bestTrack()))); // K+

           xParticles.push_back(pFactory.particle(xTracks[0],muonMass,kinChi,kinNdf,muonSigma));
           xParticles.push_back(pFactory.particle(xTracks[1],muonMass,kinChi,kinNdf,muonSigma));
           xParticles.push_back(pFactory.particle(xTracks[2],trakMass1,kinChi,kinNdf,trakSigma1));
           xParticles.push_back(pFactory.particle(xTracks[3],trakMass2,kinChi,kinNdf,trakSigma2));

           KinematicParticleVertexFitter kFitter;
           RefCountedKinematicTree xVertexFitTree;
           xVertexFitTree = kFitter.fit(xParticles);

           if (xVertexFitTree->isEmpty()) continue;

           xVertexFitTree->movePointerToTheTop();
           RefCountedKinematicParticle fitX = xVertexFitTree->currentParticle();
           RefCountedKinematicVertex fitXVertex = xVertexFitTree->currentDecayVertex();

           double x_ma_fit = 14000.;
           double x_vp_fit = -9999.;
           double x_x2_fit = 10000.;
           double x_ndof_fit = 10000.;

           if (!fitX->currentState().isValid()) continue;

           x_ma_fit = fitX->currentState().mass();
           x_x2_fit = fitXVertex->chiSquared();
           x_vp_fit = ChiSquaredProbability(x_x2_fit,
                                                (double)(fitXVertex->degreesOfFreedom()));
           x_ndof_fit = (double)(fitXVertex->degreesOfFreedom());

           TVector3 vtx;
           TVector3 pvtx;
           VertexDistanceXY vdistXY;
           int   x_ch_fit = x->charge();
           double x_px_fit = fitX->currentState().kinematicParameters().momentum().x();
           double x_py_fit = fitX->currentState().kinematicParameters().momentum().y();
           double x_pz_fit = fitX->currentState().kinematicParameters().momentum().z();
           double x_en_fit = sqrt(x_ma_fit*x_ma_fit+x_px_fit*x_px_fit+
                                     x_py_fit*x_py_fit+x_pz_fit*x_pz_fit);
           double x_vx_fit = fitXVertex->position().x();
	         double x_vy_fit = fitXVertex->position().y();
           double x_vz_fit = fitXVertex->position().z();

           vtx.SetXYZ(x_vx_fit,x_vy_fit,0);
           TVector3 pperp(x_px_fit, x_py_fit, 0);
           AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
           pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
           TVector3 vdiff = vtx - pvtx;
           double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
           Measurement1D distXY = vdistXY.distance(reco::Vertex(*fitXVertex), thePrimaryV);
           double ctauPV = distXY.value()*cosAlpha * x_ma_fit/pperp.Perp();
           GlobalError v1e = (reco::Vertex(*fitXVertex)).error();
           GlobalError v2e = thePrimaryV.error();
           AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
           double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*x_ma_fit/(pperp.Perp2());


           DiMuonTTCand.addUserInt("tPMatch",filters[i]);
           DiMuonTTCand.addUserInt("tNMatch",filters[j]);

           DiMuonTTCand.addUserInt("mass_rf",x_ma_fit);
           DiMuonTTCand.addUserFloat("vProb",x_vp_fit);
           DiMuonTTCand.addUserFloat("vChi2",x_x2_fit);
           DiMuonTTCand.addUserFloat("nDof",x_ndof_fit);
           DiMuonTTCand.addUserFloat("cosAlpha",cosAlpha);
           DiMuonTTCand.addUserFloat("ctauPV",ctauPV);
           DiMuonTTCand.addUserFloat("ctauErrPV",ctauErrPV);



             if (addMCTruth_) {
               reco::GenParticleRef genMu1 = pmu1->genParticleRef();
               reco::GenParticleRef genMu2 = pmu2->genParticleRef();
               reco::GenParticleRef genKaon1 = posTrack->genParticleRef();
               reco::GenParticleRef genKaon2 = negTrack->genParticleRef();

               if (genMu1.isNonnull() && genMu2.isNonnull() && genKaon1.isNonnull() && genKaon2.isNonnull()) {
                 if (genMu1->numberOfMothers()>0 && genMu2->numberOfMothers()>0 && genKaon1->numberOfMothers()>0 && gengenKaon1Mu2->numberOfMothers()>0){
                   reco::GenParticleRef mumu_mom1 = genMu1->motherRef();
                   reco::GenParticleRef mumu_mom2 = genMu2->motherRef();
                   if (mumu_mom1.isNonnull() && (mumu_mom1 == mumu_mom2)) {

                     std::pair<int, float> MCinfo = findJpsiMCInfo(mumu_mom1);
                     DiMuonTTCand.addUserInt("jpsiPDGId",MCinfo.first);
                     DiMuonTTCand.addUserFloat("jpsiPpdlTrue",MCinfo.second);
                   } else {
                     DiMuonTTCand.addUserInt("jpsiPDGId",0);
                     DiMuonTTCand.addUserFloat("jpsiPpdlTrue",-99.);
                   }

                   reco::GenParticleRef tktk_mom1 = genKaon1->motherRef();
                   reco::GenParticleRef tktk_mom2 = genKaon2->motherRef();
                   if (tktk_mom1.isNonnull() && (tktk_mom1 == tktk_mom2)) {

                     std::pair<int, float> MCinfo = findJpsiMCInfo(tktk_mom1);
                     DiMuonTTCand.addUserInt("phiPDGId",MCinfo.first);
                     DiMuonTTCand.addUserFloat("phiPpdlTrue",MCinfo.second);
                   } else {
                     DiMuonTTCand.addUserInt("phiPDGId",0);
                     DiMuonTTCand.addUserFloat("phiPpdlTrue",-99.);
                   }

                   reco::GenParticleRef x_mom1 = tktk_mom1->motherRef();
                   reco::GenParticleRef x_mom2 = mumu_mom1->motherRef();

                   if (x_mom1.isNonnull() && (x_mom1 == x_mom2)) {
                     DiMuonTTCand.setGenParticleRef(x_mom1); // set
                     DiMuonTTCand.embedGenParticle();      // and embed
                     DiMuonTTCand.addUserInt("xPDGId",x_mom1.pdgId());
                   } else {
                     DiMuonTTCand.addUserInt("xPDGId",0);
                   }


                 }
              } else {
                 DiMuonTTCand.addUserInt("xPDGId",0);
                 DiMuonTTCand.addUserInt("phiPDGId",0);
                 DiMuonTTCand.addUserFloat("phiPpdlTrue",-99.);
                 DiMuonTTCand.addUserInt("jpsiPDGId",0);
                 DiMuonTTCand.addUserFloat("jpsiPpdlTrue",-99.);
               }
             }


             DiMuonTTCandColl->push_back(DiMuonTTCand);
             candidates++;
             ncombo++;


         }
         } // loop over second track
       }   // loop on track candidates
       if (OnlyBest_) break;
     }

  if ( ncombo != DiMuonTTCandColl->size() ) std::cout <<"ncombo ("<<ncombo<< ") != DiMuonTT ("<<DiMuonTTCandColl->size()<<")"<< std::endl;
  if ( !dimuon->empty() )  ndimuon++;
  if ( ncombo > 0 ) nreco++;
  iEvent.put(std::move(DiMuonTTCandColl),product_name_);
  nevents++;
}

void DiMuonDiTrakProducer::endJob(){
  std::cout << "###########################" << std::endl;
  std::cout << "DiMuonDiTrak Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "Events with DiMuon candidates " << ndimuon << std::endl;
  std::cout << "Events with DiMuonDiTrak candidates " << nreco << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " DiMuonDiTrak candidates." << std::endl;
  std::cout << "###########################" << std::endl;
}

bool DiMuonDiTrakProducer::IsTheSame(const pat::PackedCandidate& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

pat::CompositeCandidate DiMuonDiTrakProducer::makeDiMuonTTCandidate(
                                          const pat::CompositeCandidate& dimuon,
				          const pat::CompositeCandidate& tt
                                         ){

  pat::CompositeCandidate DiMuonTCand;
  DiMuonTCand.addDaughter(dimuon,"dimuon");
  DiMuonTCand.addDaughter(tt,"ditrak");
  DiMuonTCand.setVertex(dimuon.vertex());
  DiMuonTCand.setCharge(tt.charge());

  reco::Candidate::LorentzVector vDiMuonT = dimuon.p4() + tt.p4();
  DiMuonTCand.setP4(vDiMuonT);

  return DiMuonTCand;

}

pat::CompositeCandidate DiMuonDiTrakProducer::makeTTCandidate(
                                          const pat::PackedCandidate& trakP,
                                          const pat::PackedCandidate& trakN
                                         ){

  pat::CompositeCandidate TTCand;
  TTCand.addDaughter(trakP,"trakP");
  TTCand.addDaughter(trakN,"trakN");
  TTCand.setCharge(trakP.charge()+trakN.charge());

  double m_kaon1 = MassTraks_[0];
  math::XYZVector mom_kaon1 = trakP.momentum();
  double e_kaon1 = sqrt(m_kaon1*m_kaon1 + mom_kaon1.Mag2());
  math::XYZTLorentzVector p4_kaon1 = math::XYZTLorentzVector(mom_kaon1.X(),mom_kaon1.Y(),mom_kaon1.Z(),e_kaon1);
  double m_kaon2 = MassTraks_[1];
  math::XYZVector mom_kaon2 = trakN.momentum();
  double e_kaon2 = sqrt(m_kaon2*m_kaon2 + mom_kaon2.Mag2());
  math::XYZTLorentzVector p4_kaon2 = math::XYZTLorentzVector(mom_kaon2.X(),mom_kaon2.Y(),mom_kaon2.Z(),e_kaon2);
  reco::Candidate::LorentzVector vTT = p4_kaon1 + p4_kaon2;
  TTCand.setP4(vTT);

  return TTCand;
}


reco::Candidate::LorentzVector DiMuonDiTrakProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonDiTrakProducer);
