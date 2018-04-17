// -*- C++ -*-
//
// Package:    DiMuonPAT
// Class:      DiMuonPAT
// 
/**\class DiMuonPAT DiMuonPAT.cc myFilters/DiMuonPAT/src/DiMuonPAT.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Warren Clarida
//         Created:  Mon Apr 19 10:52:47 CDT 2010
// $Id: X4140FilterPAT.h,v 1.2 2010/07/19 13:15:47 wclarida Exp $
//
//

#ifndef _XFilterPAT_h
#define _XFilterPAT_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

//
// class declaration
//
/*
class X4140FilterPAT : public edm::EDFilter {
   public:
      explicit X4140FilterPAT(const edm::ParameterSet&);
      ~X4140FilterPAT();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------

};*/
#endif
