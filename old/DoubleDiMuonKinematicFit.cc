// -*- C++ -*-
//
// Package:    DoubleDiMuonKinematicFit
// Class:      DoubleDiMuonKinematicFit
//
/**\class DoubleDiMuonKinematicFit Ponia/OniaTrak/src/DoubleDiMuonKinematicFit.cc

 Description: performs vertex kinematical fit for DiMuon + DiTrak candidates

**/
//
// Original Author: Adriano Di Florio
//         Created: Based on Alberto Sanchez-Hernandez PsiTkTk code

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

///For kinematic fit:
#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>

#include <boost/foreach.hpp>
#include <string>

//
// class declaration
//



class DoubleDiMuonKinematicFit : public edm::EDProducer {
   public:
      explicit DoubleDiMuonKinematicFit(const edm::ParameterSet&);
      ~DoubleDiMuonKinematicFit() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginJob() override ;
      void produce(edm::Event&, const edm::EventSetup&) override;
      void endJob() override ;

      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      std::tuple<int, float, float> findJpsiMCInfo(reco::GenParticleRef genJpsi);

// ----------member data ---------------------------
  edm::EDGetTokenT<pat::CompositeCandidateCollection> doubledimuon_cand_;
  double mass_low_dimuon,mass_hig_dimuon;
  std::vector<double> DoubleDiMuonMassCuts_;
  std::string product_name_;

  template<typename T>
  struct GreaterByVProb {
     typedef T first_argument_type;
     typedef T second_argument_type;
     bool operator()( const T & t1, const T & t2 ) const {
        return t1.userFloat("vProb") > t2.userFloat("vProb");
     }
  };
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
DoubleDiMuonKinematicFit::DoubleDiMuonKinematicFit(const edm::ParameterSet& iConfig) {
  doubledimuon_cand_   = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("DoubleDiMuonCollection"));
  mass_low_dimuon         = iConfig.getParameter<double>("LowMassConstraint");
  mass_hig_dimuon         = iConfig.getParameter<double>("HighMassConstraint");
  DoubleDiMuonMassCuts_ = iConfig.getParameter<std::vector<double>>("DoubleDiMuonMassCuts");
  product_name_ = iConfig.getParameter<std::string>("product_name");

// kinematic refit collections
  produces<pat::CompositeCandidateCollection>(product_name_);

// now do what ever other initialization is needed
}

DoubleDiMuonKinematicFit::~DoubleDiMuonKinematicFit() {
// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

std::tuple<int, float, float>
DoubleDiMuonKinematicFit::findJpsiMCInfo(reco::GenParticleRef genJpsi) {

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

// ------------ method called to produce the data  ------------
void DoubleDiMuonKinematicFit::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // Grab paramenters
  edm::Handle<pat::CompositeCandidateCollection> DoubleDiMuHandle;
  iEvent.getByToken(doubledimuon_cand_, DoubleDiMuHandle);

// Kinematic refit collection
  std::unique_ptr< pat::CompositeCandidateCollection > DoubleDiMuRefitColl(new pat::CompositeCandidateCollection);

// Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  int indexDoubleDiMu=-1;
  for (pat::CompositeCandidateCollection::const_iterator oniat=DoubleDiMuHandle->begin(); oniat!=DoubleDiMuHandle->end(); ++oniat) {

    const pat::Muon* jpsiMu1 = dynamic_cast<const pat::Muon*>(oniat->daughter("higdimuon")->daughter("muon1"));
    const pat::Muon* jpsiMu2 = dynamic_cast<const pat::Muon*>(oniat->daughter("higdimuon")->daughter("muon2"));
    const pat::Muon* phiMu1  = dynamic_cast<const pat::Muon*>(oniat->daughter("lowdimuon")->daughter("muon1"));
    const pat::Muon* phiMu2  = dynamic_cast<const pat::Muon*>(oniat->daughter("lowdimuon")->daughter("muon2"));

    int phiGenPdgId, xGenPdgId;
    float phiPpdlTrue, xGenIsPrompt;

    if (true) {

      phiGenPdgId  = -1;
      xGenPdgId    = -1;
      phiPpdlTrue  = -1.0;
      xGenIsPrompt = -1.0;

      reco::GenParticleRef genMu1 = phiMu1->genParticleRef();
      reco::GenParticleRef genMu2 = phiMu2->genParticleRef();
      reco::GenParticleRef genMu3 = jpsiMu1->genParticleRef();
      reco::GenParticleRef genMu4 = jpsiMu2->genParticleRef();
      // reco::GenParticleRef genKaon1 = posTrack.genParticleRef();
      // reco::GenParticleRef genKaon2 = negTrack.genParticleRef();

      if (genMu1.isNonnull() && genMu2.isNonnull() && genMu3.isNonnull() && genMu4.isNonnull() ) {
        if (genMu1->numberOfMothers()>0 && genMu2->numberOfMothers()>0){
          reco::GenParticleRef mumu_mom1 = genMu1->motherRef();
          reco::GenParticleRef mumu_mom2 = genMu2->motherRef();

          if (mumu_mom1.isNonnull() && (mumu_mom1 == mumu_mom2)) {

            std::tuple<int,float,float> MCinfo = findJpsiMCInfo(mumu_mom1);
            phiGenPdgId  = mumu_mom1->pdgId();
            phiPpdlTrue  = std::get<1>(MCinfo);
            xGenPdgId    = std::get<0>(MCinfo);
            xGenIsPrompt = std::get<2>(MCinfo);

          }
        }
     }
    }

    indexDoubleDiMu++;
    reco::TrackRef HighPhiTk[4]={
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("higdimuon")->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("higdimuon")->daughter("muon2") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("lowdimuon")->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("lowdimuon")->daughter("muon2") ) )->innerTrack()
    };

    reco::TrackRef HighTk[2]={
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("higdimuon")->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("higdimuon")->daughter("muon2") ) )->innerTrack()
    };

    reco::TrackRef PhiTk[2]={
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("lowdimuon")->daughter("muon1") ) )->innerTrack(),
      ( dynamic_cast<const pat::Muon*>(oniat->daughter("lowdimuon")->daughter("muon2") ) )->innerTrack()
    };

    std::vector<reco::TransientTrack> HighMuTT;
    HighMuTT.push_back((*theB).build(&HighTk[0]));
    HighMuTT.push_back((*theB).build(&HighTk[1]));

    std::vector<reco::TransientTrack> LowMuTT;
    LowMuTT.push_back((*theB).build(&PhiTk[0]));
    LowMuTT.push_back((*theB).build(&PhiTk[1]));

    const pat::CompositeCandidate *highDiMCand = dynamic_cast<const pat::CompositeCandidate *>(oniat->daughter("higdimuon"));
    const reco::Vertex thePrimaryV = *highDiMCand->userData<reco::Vertex>("PVwithmuons");

    std::vector<reco::TransientTrack> doubeDiMuTT;
    doubeDiMuTT.push_back((*theB).build(&HighPhiTk[0]));
    doubeDiMuTT.push_back((*theB).build(&HighPhiTk[1]));
    doubeDiMuTT.push_back((*theB).build(&HighPhiTk[2]));
    doubeDiMuTT.push_back((*theB).build(&HighPhiTk[3]));

    KinematicParticleFactoryFromTransientTrack pFactory;

    const ParticleMass muonMass(0.1056583);
    float muonSigma = muonMass*1E-6;

    std::vector<RefCountedKinematicParticle> allDoubleDiMuDaughters;

    std::vector<RefCountedKinematicParticle> highDiMParticles;
    highDiMParticles.push_back(pFactory.particle(HighMuTT[0],muonMass,float(0),float(0),muonSigma));
    highDiMParticles.push_back(pFactory.particle(HighMuTT[1],muonMass,float(0),float(0),muonSigma));

    std::vector<RefCountedKinematicParticle> lowDiMParticles;
    lowDiMParticles.push_back(pFactory.particle(LowMuTT[0],muonMass,float(0),float(0),muonSigma));
    lowDiMParticles.push_back(pFactory.particle(LowMuTT[1],muonMass,float(0),float(0),muonSigma));

    KinematicParticleVertexFitter fitter;

    RefCountedKinematicTree highDiMVertexFitTree;
    highDiMVertexFitTree = fitter.fit(highDiMParticles);

    if(highDiMVertexFitTree->isValid())
    {

      const ParticleMass highDiM_mass(mass_hig_dimuon);
      float highDiM_sigma = 1E-6 * highDiM_mass;

      KinematicParticleFitter csFitterHigh;
      KinematicConstraint * highDiM_c = new MassKinematicConstraint(highDiM_mass,highDiM_sigma);

      highDiMVertexFitTree->movePointerToTheTop();
      highDiMVertexFitTree = csFitterHigh.fit(highDiM_c,highDiMVertexFitTree);

      if (highDiMVertexFitTree->isValid()) {

        RefCountedKinematicTree lowDiMVertexFitTree;
        lowDiMVertexFitTree = fitter.fit(lowDiMParticles);

        highDiMVertexFitTree->movePointerToTheTop();
      	RefCountedKinematicParticle fitHigh = highDiMVertexFitTree->currentParticle();

        allDoubleDiMuDaughters.push_back(pFactory.particle (LowMuTT[0], muonMass, float(0), float(0), muonSigma));
        allDoubleDiMuDaughters.push_back(pFactory.particle (LowMuTT[1], muonMass, float(0), float(0), muonSigma));
      	allDoubleDiMuDaughters.push_back(fitHigh);

        KinematicConstrainedVertexFitter constVertexFitter;

      	MultiTrackKinematicConstraint *lowDiM_mtc = new  TwoTrackMassKinematicConstraint(mass_low_dimuon);
      	RefCountedKinematicTree DoubleDiMuTree = constVertexFitter.fit(allDoubleDiMuDaughters,lowDiM_mtc);


        if (!DoubleDiMuTree->isEmpty())
        {
          DoubleDiMuTree->movePointerToTheTop();
          RefCountedKinematicParticle fitDoubleDiMu = DoubleDiMuTree->currentParticle();
          RefCountedKinematicVertex DoubleDiMuDecayVertex = DoubleDiMuTree->currentDecayVertex();

          if (fitDoubleDiMu->currentState().isValid())
          {


            float DoubleDiMuM_fit  = fitDoubleDiMu->currentState().mass();

            if ( DoubleDiMuM_fit < DoubleDiMuonMassCuts_[0] || DoubleDiMuM_fit > DoubleDiMuonMassCuts_[1])
            continue;

	          float DoubleDiMuPx_fit = fitDoubleDiMu->currentState().kinematicParameters().momentum().x();
            float DoubleDiMuPy_fit = fitDoubleDiMu->currentState().kinematicParameters().momentum().y();
	          float DoubleDiMuPz_fit = fitDoubleDiMu->currentState().kinematicParameters().momentum().z();
            float DoubleDiMuVtxX_fit = DoubleDiMuDecayVertex->position().x();
            float DoubleDiMuVtxY_fit = DoubleDiMuDecayVertex->position().y();
            float DoubleDiMuVtxZ_fit = DoubleDiMuDecayVertex->position().z();
            float DoubleDiMu_x2_fit = DoubleDiMuDecayVertex->chiSquared();
            float DoubleDiMu_nDof_fit = DoubleDiMuDecayVertex->degreesOfFreedom();
            float DoubleDiMuVtxP_fit = ChiSquaredProbability((double)(DoubleDiMu_x2_fit),
                                                       (double)(DoubleDiMu_nDof_fit));
            float DoubleDiMu_en_fit = sqrt(DoubleDiMuM_fit*DoubleDiMuM_fit+DoubleDiMuPx_fit*DoubleDiMuPx_fit+
                                                       DoubleDiMuPy_fit*DoubleDiMuPy_fit+DoubleDiMuPz_fit*DoubleDiMuPz_fit);

            if ( DoubleDiMuVtxP_fit < 0.0 )
            continue;

            TVector3 vtx;
            TVector3 pvtx;
            VertexDistanceXY vdistXY;

            vtx.SetXYZ(DoubleDiMuVtxX_fit,DoubleDiMuVtxY_fit,0);

            TVector3 pperp(DoubleDiMuPx_fit, DoubleDiMuPy_fit, 0);
            AlgebraicVector3 vpperp(pperp.x(),pperp.y(),0);
            pvtx.SetXYZ(thePrimaryV.position().x(),thePrimaryV.position().y(),0);
            TVector3 vdiff = vtx - pvtx;
            double cosAlpha = vdiff.Dot(pperp)/(vdiff.Perp()*pperp.Perp());
            Measurement1D distXY = vdistXY.distance(reco::Vertex(*DoubleDiMuDecayVertex), thePrimaryV);
            double ctauPV = distXY.value()*cosAlpha * DoubleDiMuM_fit/pperp.Perp();
            GlobalError v1e = (reco::Vertex(*DoubleDiMuDecayVertex)).error();
            GlobalError v2e = thePrimaryV.error();
            AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
            double ctauErrPV = sqrt(ROOT::Math::Similarity(vpperp,vXYe))*DoubleDiMuM_fit/(pperp.Perp2());

            reco::CompositeCandidate recoDoubleDiMu(0, math::XYZTLorentzVector(DoubleDiMuPx_fit, DoubleDiMuPy_fit, DoubleDiMuPz_fit,DoubleDiMu_en_fit), math::XYZPoint(DoubleDiMuVtxX_fit,
                                              DoubleDiMuVtxY_fit, DoubleDiMuVtxZ_fit), 531);


	          pat::CompositeCandidate patDoubleDiMu(recoDoubleDiMu);
            patDoubleDiMu.addUserFloat("vProb",DoubleDiMuVtxP_fit);
            patDoubleDiMu.addUserFloat("vChi2",DoubleDiMu_x2_fit);
            patDoubleDiMu.addUserFloat("nDof",DoubleDiMu_nDof_fit);
            patDoubleDiMu.addUserFloat("cosAlpha",cosAlpha);
            patDoubleDiMu.addUserFloat("ctauPV",ctauPV);
            patDoubleDiMu.addUserFloat("ctauErrPV",ctauErrPV);

            patDoubleDiMu.addUserInt("bIndex",indexDoubleDiMu);

            patDoubleDiMu.addUserInt("phiGenPdgId",phiGenPdgId);
            patDoubleDiMu.addUserFloat("phiPpdlTrue",phiPpdlTrue);
            patDoubleDiMu.addUserInt("xGenPdgId",xGenPdgId);
            patDoubleDiMu.addUserFloat("xGenIsPrompt",xGenIsPrompt);

            //get first muon
            bool child = DoubleDiMuTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle fitPhiMu1 = DoubleDiMuTree->currentParticle();
            if (!child) break;

            float lowDiMMu1M_fit  = fitPhiMu1->currentState().mass();
            float lowDiMMu1Q_fit  = fitPhiMu1->currentState().particleCharge();
            float lowDiMMu1Px_fit = fitPhiMu1->currentState().kinematicParameters().momentum().x();
	          float lowDiMMu1Py_fit = fitPhiMu1->currentState().kinematicParameters().momentum().y();
            float lowDiMMu1Pz_fit = fitPhiMu1->currentState().kinematicParameters().momentum().z();
	          reco::CompositeCandidate recoPhiMu1(lowDiMMu1Q_fit, math::XYZTLorentzVector(lowDiMMu1Px_fit, lowDiMMu1Py_fit, lowDiMMu1Pz_fit,
                                             sqrt(lowDiMMu1M_fit*lowDiMMu1M_fit + lowDiMMu1Px_fit*lowDiMMu1Px_fit + lowDiMMu1Py_fit*lowDiMMu1Py_fit +
                                             lowDiMMu1Pz_fit*lowDiMMu1Pz_fit)), math::XYZPoint(DoubleDiMuVtxX_fit, DoubleDiMuVtxY_fit, DoubleDiMuVtxZ_fit), 13);
	          pat::CompositeCandidate patPhiMu1(recoPhiMu1);

            //get second muon
            child = DoubleDiMuTree->movePointerToTheNextChild();
	          RefCountedKinematicParticle fitPhiMu2 = DoubleDiMuTree->currentParticle();
            if (!child) break;

	          float lowDiMMu2M_fit  = fitPhiMu2->currentState().mass();
            float lowDiMMu2Q_fit  = fitPhiMu2->currentState().particleCharge();
            float lowDiMMu2Px_fit = fitPhiMu2->currentState().kinematicParameters().momentum().x();
	          float lowDiMMu2Py_fit = fitPhiMu2->currentState().kinematicParameters().momentum().y();
	          float lowDiMMu2Pz_fit = fitPhiMu2->currentState().kinematicParameters().momentum().z();
            reco::CompositeCandidate recoPhiMu2(lowDiMMu2Q_fit, math::XYZTLorentzVector(lowDiMMu2Px_fit, lowDiMMu2Py_fit, lowDiMMu2Pz_fit,
                                             sqrt(lowDiMMu2M_fit*lowDiMMu2M_fit + lowDiMMu2Px_fit*lowDiMMu2Px_fit + lowDiMMu2Py_fit*lowDiMMu2Py_fit +
                                             lowDiMMu2Pz_fit*lowDiMMu2Pz_fit)), math::XYZPoint(DoubleDiMuVtxX_fit, DoubleDiMuVtxY_fit, DoubleDiMuVtxZ_fit), 13);
            pat::CompositeCandidate patPhiMu2(recoPhiMu2);

            pat::CompositeCandidate lowDiMRefit;
            lowDiMRefit.addDaughter(patPhiMu1,"muon1");
            lowDiMRefit.addDaughter(patPhiMu2,"muon2");
            lowDiMRefit.setP4(patPhiMu1.p4()+patPhiMu2.p4());

            child = DoubleDiMuTree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitHigh = DoubleDiMuTree->currentParticle();
	          if (!child) break;

            float highDiMM_fit  = fitHigh->currentState().mass();
            float highDiMPx_fit = fitHigh->currentState().kinematicParameters().momentum().x();
	          float highDiMPy_fit = fitHigh->currentState().kinematicParameters().momentum().y();
            float highDiMPz_fit = fitHigh->currentState().kinematicParameters().momentum().z();

            patDoubleDiMu.addUserFloat("highDiMM_fit",highDiMM_fit);
            patDoubleDiMu.addUserFloat("highDiMPx_fit",highDiMPx_fit);
            patDoubleDiMu.addUserFloat("highDiMPy_fit",highDiMPy_fit);
            patDoubleDiMu.addUserFloat("highDiMPz_fit",highDiMPz_fit);

            child = highDiMVertexFitTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle fitjPsiMu1 = highDiMVertexFitTree->currentParticle();
            if (!child) break;

            float highDiMMu1M_fit  = fitjPsiMu1->currentState().mass();
            float highDiMMu1Q_fit  = fitjPsiMu1->currentState().particleCharge();
            float highDiMMu1Px_fit = fitjPsiMu1->currentState().kinematicParameters().momentum().x();
            float highDiMMu1Py_fit = fitjPsiMu1->currentState().kinematicParameters().momentum().y();
            float highDiMMu1Pz_fit = fitjPsiMu1->currentState().kinematicParameters().momentum().z();
            reco::CompositeCandidate recojPsiMu1(highDiMMu1Q_fit, math::XYZTLorentzVector(highDiMMu1Px_fit, highDiMMu1Py_fit, highDiMMu1Pz_fit,
                                             sqrt(highDiMMu1M_fit*highDiMMu1M_fit + highDiMMu1Px_fit*highDiMMu1Px_fit + highDiMMu1Py_fit*highDiMMu1Py_fit +
                                             highDiMMu1Pz_fit*highDiMMu1Pz_fit)), math::XYZPoint(DoubleDiMuVtxX_fit, DoubleDiMuVtxY_fit, DoubleDiMuVtxZ_fit), 13);
            pat::CompositeCandidate patjPsiMu1(recojPsiMu1);

            //get second muon
            child = highDiMVertexFitTree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitjPsiMu2 = highDiMVertexFitTree->currentParticle();
            if (!child) break;

            float highDiMMu2M_fit  = fitjPsiMu2->currentState().mass();
            float highDiMMu2Q_fit  = fitjPsiMu2->currentState().particleCharge();
            float highDiMMu2Px_fit = fitjPsiMu2->currentState().kinematicParameters().momentum().x();
            float highDiMMu2Py_fit = fitjPsiMu2->currentState().kinematicParameters().momentum().y();
            float highDiMMu2Pz_fit = fitjPsiMu2->currentState().kinematicParameters().momentum().z();
            reco::CompositeCandidate recojPsiMu2(highDiMMu2Q_fit, math::XYZTLorentzVector(highDiMMu2Px_fit, highDiMMu2Py_fit, highDiMMu2Pz_fit,
                                             sqrt(highDiMMu2M_fit*highDiMMu2M_fit + highDiMMu2Px_fit*highDiMMu2Px_fit + highDiMMu2Py_fit*highDiMMu2Py_fit +
                                             highDiMMu2Pz_fit*highDiMMu2Pz_fit)), math::XYZPoint(DoubleDiMuVtxX_fit, DoubleDiMuVtxY_fit, DoubleDiMuVtxZ_fit), 13);
            pat::CompositeCandidate patjPsiMu2(recojPsiMu2);

            pat::CompositeCandidate highDiMRefit;
            highDiMRefit.addDaughter(patjPsiMu1,"muon1");
            highDiMRefit.addDaughter(patjPsiMu2,"muon2");
            highDiMRefit.setP4(patjPsiMu1.p4()+patjPsiMu2.p4());



            patDoubleDiMu.addDaughter(highDiMRefit,"lowdimuon");
            patDoubleDiMu.addDaughter(lowDiMRefit,"higdimuon");

            DoubleDiMuRefitColl->push_back(patDoubleDiMu);
          }

        }

      }

    }
  }


// now sort by vProb
  DoubleDiMuonKinematicFit::GreaterByVProb<pat::CompositeCandidate> vPComparator;
  std::sort(DoubleDiMuRefitColl->begin(),DoubleDiMuRefitColl->end(),vPComparator);
  iEvent.put(std::move(DoubleDiMuRefitColl),product_name_);
}

// ------------ method called once each job just before starting event loop  ------------
void DoubleDiMuonKinematicFit::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void DoubleDiMuonKinematicFit::endJob() {}

// ------------ method called when starting to processes a run  ------------
void DoubleDiMuonKinematicFit::beginRun(edm::Run&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void DoubleDiMuonKinematicFit::endRun(edm::Run&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void DoubleDiMuonKinematicFit::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void DoubleDiMuonKinematicFit::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DoubleDiMuonKinematicFit::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DoubleDiMuonKinematicFit);
