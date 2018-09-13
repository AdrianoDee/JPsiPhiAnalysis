
#include "../interface/DoubleDiMuonProducer.h"

DoubleDiMuonProducer::DoubleDiMuonProducer(const edm::ParameterSet& iConfig):
  HighDiMuonCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("HighDiMuonCollection"))),
  LowDiMuonCollection_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("LowDiMuonCollection"))),
  thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
  thePVs_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexTag"))),
  HighDiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("HighDiMuonMassCuts")),
  LowDiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("LowDiMuonMassCuts")),
  DoubleDiMuonMassCuts_(iConfig.getParameter<std::vector<double>>("DoubleDiMuonMassCuts")),
  JPsiMass_(iConfig.getParameter<double>("JPsiMass")),
  PhiMass_(iConfig.getParameter<double>("PhiMass")),
  addMCTruth_(iConfig.getParameter<bool>("AddMCTruth")),
  addSameSig_(iConfig.getParameter<bool>("AddSS")),
  doDoubleConstant_(iConfig.getParameter<bool>("DoDouble"))
{
  produces<pat::CompositeCandidateCollection>("FourMuonCandidates");
  candidates = 0;
  nevents = 0;
  nLdM = 0;
  nHdM = 0;
}

std::tuple<int, float, float>
DoubleDiMuonProducer::findJpsiMCInfo(reco::GenParticleRef genJpsi) {

  // std::cout << "findJpsiMCInfo 1 " << std::endl;
  int momJpsiID = 0;
  float trueLife = -99.;
  float isPrompt = -99.;
  if (genJpsi->numberOfMothers()>0) {

    // std::cout << "findJpsiMCInfo 1 " << std::endl;

    TVector3 trueVtx(0.0,0.0,0.0);
    TVector3 trueP(0.0,0.0,0.0);
    TVector3 trueVtxMom(0.0,0.0,0.0);

    trueVtx.SetXYZ(genJpsi->vertex().x(),genJpsi->vertex().y(),genJpsi->vertex().z());
    trueP.SetXYZ(genJpsi->momentum().x(),genJpsi->momentum().y(),genJpsi->momentum().z());

    reco::GenParticleRef Jpsimom = genJpsi->motherRef();       // find mothers
    // std::cout << "findJpsiMCInfo 1 " << std::endl;
    if (Jpsimom.isNull()) {
      std::tuple<int, float, float> result = std::make_tuple(momJpsiID, trueLife,isPrompt);
      return result;
    } else
    {
    momJpsiID = Jpsimom->pdgId();
    isPrompt = Jpsimom->isPromptDecayed();
    trueVtxMom.SetXYZ(Jpsimom->vertex().x(),Jpsimom->vertex().y(),Jpsimom->vertex().z());
    TVector3 vdiff = trueVtx - trueVtxMom;
    trueLife = vdiff.Perp()*genJpsi->mass()/trueP.Perp();
  }
}
  std::tuple<int,float,float> result = std::make_tuple(momJpsiID, trueLife,isPrompt);
  return result;

}

void DoubleDiMuonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::unique_ptr<pat::CompositeCandidateCollection> DoubleDiMuonCandColl(new pat::CompositeCandidateCollection);

  edm::Handle<pat::CompositeCandidateCollection> highDiMuon;
  iEvent.getByToken(HighDiMuonCollection_,highDiMuon);

  edm::Handle<pat::CompositeCandidateCollection> lowDiMuon;
  iEvent.getByToken(LowDiMuonCollection_,lowDiMuon);

  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_,theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;
  reco::Vertex theBeamSpotV = reco::Vertex(bs.position(), bs.covariance3D());

  edm::Handle<reco::VertexCollection> priVtxs;
  iEvent.getByToken(thePVs_, priVtxs);

  edm::ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  float HighDiMuonMassMax_ = HighDiMuonMassCuts_[1];
  float HighDiMuonMassMin_ = HighDiMuonMassCuts_[0];
  float LowDiMuonMassMax_ = LowDiMuonMassCuts_[1];
  float LowDiMuonMassMin_ = LowDiMuonMassCuts_[0];

  float DoubleDiMuonMassMax_ = DoubleDiMuonMassCuts_[1];
  float DoubleDiMuonMassMin_ = DoubleDiMuonMassCuts_[0];

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

// Note: Dimuon cand are sorted by decreasing vertex probability then first is associated with "best" dimuon
  //Looking for J/Psi

  for (pat::CompositeCandidateCollection::const_iterator highCand = highDiMuon->begin(); highCand != highDiMuon->end(); ++highCand){

     if ( highCand->mass() < HighDiMuonMassMax_  && highCand->mass() > HighDiMuonMassMin_ ) {
       const pat::Muon *jPsiMuHigh = dynamic_cast<const pat::Muon*>(highCand->daughter("highMuon"));
       const pat::Muon *jPsiMuLow = dynamic_cast<const pat::Muon*>(highCand->daughter("lowMuon"));


       for (pat::CompositeCandidateCollection::const_iterator lowCand = lowDiMuon->begin(); lowCand != lowDiMuon->end(); ++lowCand){

          if ( lowCand->mass() < LowDiMuonMassMax_  && lowCand->mass() > LowDiMuonMassMin_ ) {

            const pat::Muon *phiMuHigh = dynamic_cast<const pat::Muon*>(lowCand->daughter("highMuon"));
            const pat::Muon *phiMuLow = dynamic_cast<const pat::Muon*>(lowCand->daughter("lowMuon"));

            if( phiMuHigh == phiMuLow || phiMuHigh == jPsiMuHigh || phiMuHigh == jPsiMuLow ) continue;
            if( phiMuLow == jPsiMuHigh || phiMuLow == jPsiMuLow ) continue;
            if( jPsiMuHigh == jPsiMuLow ) continue;

            if(lowCand->userFloat("vProb")<0.0)
              continue;

            if(highCand->userFloat("vProb")<0.0)
              continue;

            pat::CompositeCandidate FourMuonCandidate = makeCandidate(*lowCand, *highCand);

            if(!addSameSig_ && FourMuonCandidate.charge() != 0.0) continue;

            if ( FourMuonCandidate.mass() < DoubleDiMuonMassMax_ && FourMuonCandidate.mass() > DoubleDiMuonMassMin_)
              {
                candidates++;

                float kinChi = 0.;
                float kinNdf = 0.;

                const ParticleMass muonMass(0.1056583);
                float muonSigma = muonMass*1E-6;

                //Normal kinematic vertex fit
                std::vector<reco::TransientTrack> xTracks;
                KinematicParticleFactoryFromTransientTrack pFactory;
                std::vector<RefCountedKinematicParticle> xParticles;

                xTracks.push_back((*theB).build(*(jPsiMuHigh->innerTrack()))); // µ
                xTracks.push_back((*theB).build(*(jPsiMuLow->innerTrack()))); // µ
                xTracks.push_back((*theB).build(*(phiMuHigh->innerTrack()))); // µ
                xTracks.push_back((*theB).build(*(phiMuLow->innerTrack()))); // µ

                xParticles.push_back(pFactory.particle(xTracks[0],muonMass,kinChi,kinNdf,muonSigma));
                xParticles.push_back(pFactory.particle(xTracks[1],muonMass,kinChi,kinNdf,muonSigma));
                xParticles.push_back(pFactory.particle(xTracks[2],muonMass,kinChi,kinNdf,muonSigma));
                xParticles.push_back(pFactory.particle(xTracks[3],muonMass,kinChi,kinNdf,muonSigma));

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

                if (!(fitX->currentState().isValid())) continue;


                x_ma_fit = fitX->currentState().mass();
                x_x2_fit = fitXVertex->chiSquared();
                x_vp_fit = ChiSquaredProbability(x_x2_fit,
                                                     (double)(fitXVertex->degreesOfFreedom()));
                x_ndof_fit = (double)(fitXVertex->degreesOfFreedom());

                FourMuonCandidate.addUserFloat("mass_rf",x_ma_fit);
                FourMuonCandidate.addUserFloat("vProb",x_vp_fit);
                FourMuonCandidate.addUserFloat("vChi2",x_x2_fit);
                FourMuonCandidate.addUserFloat("nDof",x_ndof_fit);
                // FourMuonCandidate.addUserInt("pId",i);
                // FourMuonCandidate.addUserInt("mId",j);

                //////////////////////////////////////////////////////////////////////////////
                //PV Selection(s)

                std::vector <double> vertexWeight,sumPTPV,cosAlpha,ctauPV,ctauErrPV;
                std::vector <int> countTksOfPV;

                TVector3 vtx, vdiff, pvtx;
                VertexDistanceXY vdistXY;
                reco::Vertex thePrimaryV,thePrimaryVDZ;
                TwoTrackMinimumDistance ttmd;

                double x_px_fit = fitX->currentState().kinematicParameters().momentum().x();
                double x_py_fit = fitX->currentState().kinematicParameters().momentum().y();
                double x_pz_fit = fitX->currentState().kinematicParameters().momentum().z();
                // double x_en_fit = sqrt(x_ma_fit*x_ma_fit+x_px_fit*x_px_fit+x_py_fit*x_py_fit+x_pz_fit*x_pz_fit);
                double x_vx_fit = fitXVertex->position().x();
     	          double x_vy_fit = fitXVertex->position().y();
                double x_vz_fit = fitXVertex->position().z();
                vtx.SetXYZ(x_vx_fit,x_vy_fit,0);

                bool status = ttmd.calculate( GlobalTrajectoryParameters(
                  GlobalPoint(x_vx_fit,x_vy_fit,x_vz_fit),
                  GlobalVector(x_px_fit,x_py_fit,x_pz_fit),TrackCharge(0),&(*magneticField)),
                  GlobalTrajectoryParameters(
                    GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()),
                    GlobalVector(bs.dxdz(), bs.dydz(), 1.),TrackCharge(0),&(*magneticField)));
                float extrapZ=-9E20;

                if (status) extrapZ=ttmd.points().first.z();

                TVector3 pperp(x_px_fit, x_py_fit, 0);
                AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);

                reco::VertexCollection verteces;
                std::vector<int> vKeys;
                verteces.push_back(theBeamSpotV);
                vKeys.push_back(0);
                float minDz = 999999.;
                double maxCosAlpha = -1.0;
                if ( !(priVtxs->begin() != priVtxs->end()) )
                {

                  thePrimaryV = reco::Vertex(*(priVtxs->begin()));
                  thePrimaryVDZ = reco::Vertex(*(priVtxs->begin()));
                  verteces.push_back(thePrimaryV);
                  verteces.push_back(thePrimaryVDZ);
                  vKeys.push_back(0);
                  vKeys.push_back(0);
                }else
                {

                  reco::Vertex p,pz;
                  int thisp,thispz;
                  for(size_t pV = 0; pV<priVtxs->size();++pV)
                  {
                    auto thisPV = priVtxs->at(pV);

                    pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
                    vdiff = vtx - pvtx;
                    double thisCosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
                    if(thisCosAlpha>maxCosAlpha)
                    {
                      thePrimaryV = reco::Vertex(thisPV);
                      maxCosAlpha = thisCosAlpha;
                      p = reco::Vertex(thisPV);
                      thisp = pV;
                    }

                    float deltaZ = fabs(extrapZ - thisPV.position().z()) ;
                    if ( deltaZ < minDz ) {
                      minDz = deltaZ;
                      pz = reco::Vertex(thisPV);
                      thispz =pV;
                    }
                  }
                  verteces.push_back(p);
                  vKeys.push_back(thisp);
                  verteces.push_back(pz);
                  vKeys.push_back(thispz);
                }

                //CosAlpha for the 3 pvx
                for(size_t i = 0; i < verteces.size(); i++)
                {
                  pvtx.SetXYZ(verteces[i].position().x(),verteces[i].position().y(),0);
                  vdiff = vtx - pvtx;
                  cosAlpha.push_back(vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp()));
                  Measurement1D distXY = vdistXY.distance(reco::Vertex(*fitXVertex), verteces[i]);
                  ctauPV.push_back(distXY.value()*cosAlpha[i] * x_ma_fit/pperp.Perp());
                  GlobalError v1e = (reco::Vertex(*fitXVertex)).error();
                  GlobalError v2e = verteces[i].error();
                  AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
                  ctauErrPV.push_back(sqrt(ROOT::Math::Similarity(vpperp,vXYe))*x_ma_fit/(pperp.Perp2()));
                }

                float candRef = -1.0, cand_const_ref = -1.0;

                FourMuonCandidate.addUserData("bestPV",reco::Vertex(thePrimaryV));

                FourMuonCandidate.addUserFloat("vtxX",x_vx_fit);
                FourMuonCandidate.addUserFloat("vtxY",x_vy_fit);

                FourMuonCandidate.addUserFloat("vtxX",x_vx_fit);
                FourMuonCandidate.addUserFloat("vtxY",x_vy_fit);
                FourMuonCandidate.addUserFloat("vtxZ",x_vz_fit);

                FourMuonCandidate.addUserFloat("cosAlphaBS",cosAlpha[0]);
                FourMuonCandidate.addUserFloat("ctauPVBS",ctauPV[0]);
                FourMuonCandidate.addUserFloat("ctauErrPVBS",ctauErrPV[0]);

                FourMuonCandidate.addUserFloat("cosAlpha",cosAlpha[1]);
                FourMuonCandidate.addUserFloat("ctauPV",ctauPV[1]);
                FourMuonCandidate.addUserFloat("ctauErrPV",ctauErrPV[1]);

                FourMuonCandidate.addUserFloat("cosAlphaDZ",cosAlpha[2]);
                FourMuonCandidate.addUserFloat("ctauPVDZ",ctauPV[2]);
                FourMuonCandidate.addUserFloat("ctauErrPVDZ",ctauErrPV[2]);

                ///DCAs
                std::vector<float> DCAs;
                for(size_t i = 0; i < xTracks.size();++i)
                {
                  for(size_t j = i+1; j < xTracks.size();++j)
                  {
                    TrajectoryStateClosestToPoint TS1 = xTracks[i].impactPointTSCP();
                    TrajectoryStateClosestToPoint TS2 = xTracks[j].impactPointTSCP();
                    float dca = 1E20;
                    if (TS1.isValid() && TS2.isValid()) {
                      ClosestApproachInRPhi cApp;
                      cApp.calculate(TS1.theState(), TS2.theState());
                      if (cApp.status() ) dca = cApp.distance();
                    }
                    DCAs.push_back(dca);
                  }
                }

                FourMuonCandidate.addUserFloat("dca_mp1mp2",DCAs[0]);
                FourMuonCandidate.addUserFloat("dca_mp1mj1",DCAs[1]);
                FourMuonCandidate.addUserFloat("dca_mp1mj2",DCAs[2]);
                FourMuonCandidate.addUserFloat("dca_mp2mj1",DCAs[3]);
                FourMuonCandidate.addUserFloat("dca_mp2mj2",DCAs[4]);
                FourMuonCandidate.addUserFloat("dca_mj1mj2",DCAs[5]);

                //JPsi Mass Constrained Fit
                KinematicConstrainedVertexFitter vertexFitter;
                MultiTrackKinematicConstraint *jpsi_mtc = new  TwoTrackMassKinematicConstraint(JPsiMass_);
                RefCountedKinematicTree FourMuTree = vertexFitter.fit(xParticles,jpsi_mtc);

                if (!FourMuTree->isEmpty()) {

                   FourMuTree->movePointerToTheTop();
                   RefCountedKinematicParticle fitFourMuT = FourMuTree->currentParticle();
                   RefCountedKinematicVertex FourMuDecayVertex = FourMuTree->currentDecayVertex();
            // Get FourMu reffited
                   double dimuontt_ma_fit = 14000.;
                   double dimuontt_vp_fit = -9999.;
                   double dimuontt_x2_fit = 10000.;
                   double dimuontt_ndof_fit = 10000.;

                   if (fitFourMuT->currentState().isValid()) {
                     dimuontt_ma_fit = fitFourMuT->currentState().mass();
                     dimuontt_x2_fit = FourMuDecayVertex->chiSquared();
                     dimuontt_vp_fit = ChiSquaredProbability(dimuontt_x2_fit,
                                                          (double)(FourMuDecayVertex->degreesOfFreedom()));
                     dimuontt_ndof_fit = (double)(FourMuDecayVertex->degreesOfFreedom());
                   }

                   if ( dimuontt_vp_fit > 0.0 ) {

                        TVector3 vtx;
                        TVector3 pvtx;
                        VertexDistanceXY vdistXY;
                        int   dimuontt_ch_fit = FourMuonCandidate.charge();
                        double dimuontt_px_fit = fitFourMuT->currentState().kinematicParameters().momentum().x();
                        double dimuontt_py_fit = fitFourMuT->currentState().kinematicParameters().momentum().y();
                        double dimuontt_pz_fit = fitFourMuT->currentState().kinematicParameters().momentum().z();
                        double dimuontt_en_fit = sqrt(dimuontt_ma_fit*dimuontt_ma_fit+dimuontt_px_fit*dimuontt_px_fit+
                                                  dimuontt_py_fit*dimuontt_py_fit+dimuontt_pz_fit*dimuontt_pz_fit);
                        double dimuontt_vx_fit = FourMuDecayVertex->position().x();
                        double dimuontt_vy_fit = FourMuDecayVertex->position().y();
                        double dimuontt_vz_fit = FourMuDecayVertex->position().z();

                        vtx.SetXYZ(dimuontt_vx_fit,dimuontt_vy_fit,0);
                        TVector3 pperp(dimuontt_px_fit, dimuontt_py_fit, 0);
                        AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
                        pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
                        TVector3 vdiff = vtx - pvtx;
                        double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
                        Measurement1D distXY = vdistXY.distance(reco::Vertex(*FourMuDecayVertex), thePrimaryV);
                        double ctauPV = distXY.value()*cosAlpha * dimuontt_ma_fit/pperp.Perp();
                        GlobalError v1e = (reco::Vertex(*FourMuDecayVertex)).error();
                        GlobalError v2e = thePrimaryV.error();
                        AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
                        double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*dimuontt_ma_fit/(pperp.Perp2());

                        reco::CompositeCandidate recoFourMu_rf(dimuontt_ch_fit,math::XYZTLorentzVector(dimuontt_px_fit,dimuontt_py_fit,dimuontt_pz_fit,dimuontt_en_fit),
                                                           math::XYZPoint(dimuontt_vx_fit,dimuontt_vy_fit,dimuontt_vz_fit),531);

                        pat::CompositeCandidate FourMuonCandidate_rf(recoFourMu_rf);

                        FourMuonCandidate.addUserFloat("vProb_ref",dimuontt_vp_fit);
                        FourMuonCandidate.addUserFloat("vChi2_ref",dimuontt_x2_fit);
                        FourMuonCandidate.addUserFloat("nDof_ref",dimuontt_ndof_fit);
                        FourMuonCandidate.addUserFloat("cosAlpha_ref",cosAlpha);
                        FourMuonCandidate.addUserFloat("ctauPV_ref",ctauPV);
                        FourMuonCandidate.addUserFloat("ctauErrPV_ref",ctauErrPV);

            // get first muon
                        bool child = FourMuTree->movePointerToTheFirstChild();
                        RefCountedKinematicParticle fitMu1 = FourMuTree->currentParticle();
                        if (!child) break;
                        float m1_ma_fit = fitMu1->currentState().mass();
                        int   m1_ch_fit = fitMu1->currentState().particleCharge();
                        float m1_px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
                        float m1_py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
                        float m1_pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
                        float m1_en_fit = sqrt(m1_ma_fit*m1_ma_fit+m1_px_fit*m1_px_fit+m1_py_fit*m1_py_fit+m1_pz_fit*m1_pz_fit);
                        reco::CompositeCandidate recoMu1(m1_ch_fit,math::XYZTLorentzVector(m1_px_fit,m1_py_fit,m1_pz_fit,m1_en_fit),
                                                         math::XYZPoint(dimuontt_vx_fit,dimuontt_vy_fit,dimuontt_vz_fit),13);
                        pat::CompositeCandidate patMu1(recoMu1);
            // get second muon
                        child = FourMuTree->movePointerToTheNextChild();
                        RefCountedKinematicParticle fitMu2 = FourMuTree->currentParticle();
                        if (!child) break;
                        float m2_ma_fit = fitMu2->currentState().mass();
                        int   m2_ch_fit = fitMu2->currentState().particleCharge();
                        float m2_px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
                        float m2_py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
                        float m2_pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
                        float m2_en_fit = sqrt(m2_ma_fit*m2_ma_fit+m2_px_fit*m2_px_fit+m2_py_fit*m2_py_fit+m2_pz_fit*m2_pz_fit);
                        reco::CompositeCandidate recoMu2(m2_ch_fit,math::XYZTLorentzVector(m2_px_fit,m2_py_fit,m2_pz_fit,m2_en_fit),
                                                         math::XYZPoint(dimuontt_vx_fit,dimuontt_vy_fit,dimuontt_vz_fit),13);
                        pat::CompositeCandidate patMu2(recoMu2);

            // get kaon
                        child = FourMuTree->movePointerToTheNextChild();
                        RefCountedKinematicParticle fitMu3 = FourMuTree->currentParticle();
                        if (!child) break;
                        float m3_ma_fit = fitMu3->currentState().mass();
                        int   m3_ch_fit = fitMu3->currentState().particleCharge();
                        float m3_px_fit = fitMu3->currentState().kinematicParameters().momentum().x();
                        float m3_py_fit = fitMu3->currentState().kinematicParameters().momentum().y();
                        float m3_pz_fit = fitMu3->currentState().kinematicParameters().momentum().z();
                        float m3_en_fit = sqrt(m3_ma_fit*m3_ma_fit+m3_px_fit*m3_px_fit+m3_py_fit*m3_py_fit+m3_pz_fit*m3_pz_fit);
                        reco::CompositeCandidate recoMu3(m3_ch_fit,math::XYZTLorentzVector(m3_px_fit,m3_py_fit,m3_pz_fit,m3_en_fit),
                                                         math::XYZPoint(dimuontt_vx_fit,dimuontt_vy_fit,dimuontt_vz_fit),321);
                        pat::CompositeCandidate patMu3(recoMu3);

            // get kaon2
                        child = FourMuTree->movePointerToTheNextChild();
                        RefCountedKinematicParticle fitMu4 = FourMuTree->currentParticle();
                        if (!child) break;
                        float m4_ma_fit = fitMu4->currentState().mass();
                        int   m4_ch_fit = fitMu4->currentState().particleCharge();
                        float m4_px_fit = fitMu4->currentState().kinematicParameters().momentum().x();
                        float m4_py_fit = fitMu4->currentState().kinematicParameters().momentum().y();
                        float m4_pz_fit = fitMu4->currentState().kinematicParameters().momentum().z();
                        float m4_en_fit = sqrt(m4_ma_fit*m4_ma_fit+m4_px_fit*m4_px_fit+m4_py_fit*m4_py_fit+m4_pz_fit*m4_pz_fit);
                        reco::CompositeCandidate recoMu4(m4_ch_fit,math::XYZTLorentzVector(m4_px_fit,m4_py_fit,m4_pz_fit,m4_en_fit),
                                                         math::XYZPoint(dimuontt_vx_fit,dimuontt_vy_fit,dimuontt_vz_fit),321);
                        pat::CompositeCandidate patMu4(recoMu4);

            // Define phi from two kaons
                        pat::CompositeCandidate phi;
                        phi.addDaughter(patMu1,"highMuon");
                        phi.addDaughter(patMu2,"lowMuon");
                        phi.setP4(patMu2.p4()+patMu1.p4());

                        // Define psi from two muons
                        pat::CompositeCandidate psi;
                        psi.addDaughter(patMu3,"highMuon");
                        psi.addDaughter(patMu4,"lowMuon");
                        psi.setP4(patMu3.p4()+patMu4.p4());


                        candRef = 1.0;
                        FourMuonCandidate_rf.addDaughter(phi,"phi");
                        FourMuonCandidate_rf.addDaughter(psi,"jpsi");
                        FourMuonCandidate.addDaughter(FourMuonCandidate_rf,"ref_cand");
                      }
                   }

                xParticles.clear();
                xParticles.push_back(pFactory.particle(xTracks[2],muonMass,kinChi,kinNdf,muonSigma));
                xParticles.push_back(pFactory.particle(xTracks[3],muonMass,kinChi,kinNdf,muonSigma));
                xParticles.push_back(pFactory.particle(xTracks[0],muonMass,kinChi,kinNdf,muonSigma));
                xParticles.push_back(pFactory.particle(xTracks[1],muonMass,kinChi,kinNdf,muonSigma));

                if(doDoubleConstant_)
                {
                  std::vector<RefCountedKinematicParticle> JPsiParticles;
                  std::vector<reco::TransientTrack> JFourMurTk;
                  JFourMurTk.push_back(xTracks[2]);
                  JFourMurTk.push_back(xTracks[3]);

                  JPsiParticles.push_back(pFactory.particle(JFourMurTk[0],muonMass,float(0),float(0),muonSigma));
                  JPsiParticles.push_back(pFactory.particle(JFourMurTk[1],muonMass,float(0),float(0),muonSigma));

                  KinematicParticleVertexFitter fitter;
                  KinematicParticleFitter csFitterJPsi;
                  RefCountedKinematicTree jpsiVertexFitTree;
                  jpsiVertexFitTree = fitter.fit(JPsiParticles);


                  if (jpsiVertexFitTree->isValid())
                  {

                    const ParticleMass jpsi_mass(JPsiMass_);
                    float jpsi_sigma = 1E-6;

                    KinematicConstraint * jpsi_c = new MassKinematicConstraint(jpsi_mass,jpsi_sigma);

                    jpsiVertexFitTree->movePointerToTheTop();
                    jpsiVertexFitTree = csFitterJPsi.fit(jpsi_c,jpsiVertexFitTree);

                    if (jpsiVertexFitTree->isValid())
                    {

                      jpsiVertexFitTree->movePointerToTheTop();
                    	 RefCountedKinematicParticle fitJPsi = jpsiVertexFitTree->currentParticle();

                      //Phi
                      std::vector<RefCountedKinematicParticle> allFourMuDaughters;

                      allFourMuDaughters.push_back(pFactory.particle(xTracks[0],muonMass,kinChi,kinNdf,muonSigma));
                      allFourMuDaughters.push_back(pFactory.particle(xTracks[1],muonMass,kinChi,kinNdf,muonSigma));
                      allFourMuDaughters.push_back(fitJPsi);

                      KinematicConstrainedVertexFitter vertexFitter;
                      MultiTrackKinematicConstraint *phi_mtc = new  TwoTrackMassKinematicConstraint(PhiMass_);
                      RefCountedKinematicTree PsiPhiTree = vertexFitter.fit(allFourMuDaughters,phi_mtc);

                      if (!PsiPhiTree->isEmpty()) {
                         PsiPhiTree->movePointerToTheTop();
                         RefCountedKinematicParticle fitFourMuT = PsiPhiTree->currentParticle();
                         RefCountedKinematicVertex FourMuDecayVertex = PsiPhiTree->currentDecayVertex();
                  // Get FourMu reffited
                         double dimuontt_ma_fit = 14000.;
                         double dimuontt_vp_fit = -9999.;
                         double dimuontt_x2_fit = 10000.;
                         double dimuontt_ndof_fit = 10000.;

                         if (fitFourMuT->currentState().isValid()) {
                           dimuontt_ma_fit = fitFourMuT->currentState().mass();
                           dimuontt_x2_fit = FourMuDecayVertex->chiSquared();
                           dimuontt_vp_fit = ChiSquaredProbability(dimuontt_x2_fit,
                                                                (double)(FourMuDecayVertex->degreesOfFreedom()));
                           dimuontt_ndof_fit = (double)(FourMuDecayVertex->degreesOfFreedom());
                         }

                         if ( dimuontt_vp_fit > 0.0 ) {
                              TVector3 vtx;
                              TVector3 pvtx;
                              VertexDistanceXY vdistXY;
                              int   dimuontt_ch_fit = FourMuonCandidate.charge();
                              double dimuontt_px_fit = fitFourMuT->currentState().kinematicParameters().momentum().x();
                              double dimuontt_py_fit = fitFourMuT->currentState().kinematicParameters().momentum().y();
                              double dimuontt_pz_fit = fitFourMuT->currentState().kinematicParameters().momentum().z();
                              double dimuontt_en_fit = sqrt(dimuontt_ma_fit*dimuontt_ma_fit+dimuontt_px_fit*dimuontt_px_fit+
                                                        dimuontt_py_fit*dimuontt_py_fit+dimuontt_pz_fit*dimuontt_pz_fit);
                              double dimuontt_vx_fit = FourMuDecayVertex->position().x();
                  	           double dimuontt_vy_fit = FourMuDecayVertex->position().y();
                              double dimuontt_vz_fit = FourMuDecayVertex->position().z();

                              vtx.SetXYZ(dimuontt_vx_fit,dimuontt_vy_fit,0);
                              TVector3 pperp(dimuontt_px_fit, dimuontt_py_fit, 0);
                              AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
                              pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
                              TVector3 vdiff = vtx - pvtx;
                              double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
                              Measurement1D distXY = vdistXY.distance(reco::Vertex(*FourMuDecayVertex), thePrimaryV);
                              double ctauPV = distXY.value()*cosAlpha * dimuontt_ma_fit/pperp.Perp();
                              GlobalError v1e = (reco::Vertex(*FourMuDecayVertex)).error();
                              GlobalError v2e = thePrimaryV.error();
                              AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
                              double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*dimuontt_ma_fit/(pperp.Perp2());

                  	           reco::CompositeCandidate recoFourMu_rf(dimuontt_ch_fit,math::XYZTLorentzVector(dimuontt_px_fit,dimuontt_py_fit,dimuontt_pz_fit,dimuontt_en_fit),
                                                                 math::XYZPoint(dimuontt_vx_fit,dimuontt_vy_fit,dimuontt_vz_fit),531);

                              pat::CompositeCandidate FourMuonCandidate_rf(recoFourMu_rf);

                              FourMuonCandidate.addUserFloat("vProb_const_ref",dimuontt_vp_fit);
                              FourMuonCandidate.addUserFloat("vChi2_const_ref",dimuontt_x2_fit);
                              FourMuonCandidate.addUserFloat("nDof_const_ref",dimuontt_ndof_fit);
                              FourMuonCandidate.addUserFloat("cosAlpha_const_ref",cosAlpha);
                              FourMuonCandidate.addUserFloat("ctauPV_const_ref",ctauPV);
                              FourMuonCandidate.addUserFloat("ctauErrPV_const_ref",ctauErrPV);

                              cand_const_ref = 1.0;
                              // FourMuonCandidate_rf.addDaughter(patJPsi_rf,"dimuon");
                  	           // FourMuonCandidate_rf.addDaughter(phi,"ditrak");
                              FourMuonCandidate.addDaughter(FourMuonCandidate_rf,"ref_const_cand");
                            }
                  	      }

                    }
                  }
                }

                FourMuonCandidate.addUserFloat("has_ref",candRef);
                FourMuonCandidate.addUserFloat("has_const_ref",cand_const_ref);


                DoubleDiMuonCandColl->push_back(FourMuonCandidate);
              }
            }
          }
        }
      }
     // if (OnlyBest_) break;

     if ( !(highDiMuon->empty()) )  nHdM++;
     if ( !(lowDiMuon->empty()) )  nLdM++;

     iEvent.put(std::move(DoubleDiMuonCandColl),"FourMuonCandidates");
     nevents++;
  }


void DoubleDiMuonProducer::endJob(){
  std::cout << "#########################################" << std::endl;
  std::cout << "DoubleDiMuon Candidate producer report:" << std::endl;
  std::cout << "#########################################" << std::endl;
  std::cout << "Found " << nevents << " Events" << std::endl;
  std::cout << "Events with LowDiMuon  candidates " << nLdM << std::endl;
  std::cout << "Events with HighDiMuon candidates " << nHdM << std::endl;
  std::cout << "#########################################" << std::endl;
  std::cout << "Found " << candidates << " DoubleDiMuon candidates." << std::endl;
  std::cout << "#########################################" << std::endl;
}


const pat::CompositeCandidate DoubleDiMuonProducer::makeCandidate(const pat::CompositeCandidate& lowDiMuon,
  const pat::CompositeCandidate& higDiMuon){
    pat::CompositeCandidate xCand;
    xCand.addDaughter(lowDiMuon,"phi");
    xCand.addDaughter(higDiMuon,"jpsi");
    reco::Candidate::LorentzVector vX = lowDiMuon.p4() + higDiMuon.p4();
    xCand.setP4(vX);
    return xCand;
  }

reco::Candidate::LorentzVector DoubleDiMuonProducer::convertVector(const math::XYZTLorentzVectorF& v){

  return reco::Candidate::LorentzVector(v.x(),v.y(), v.z(), v.t());
}
//define this as a plug-in
DEFINE_FWK_MODULE(DoubleDiMuonProducer);
