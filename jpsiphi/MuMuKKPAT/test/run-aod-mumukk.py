import FWCore.ParameterSet.Config as cms

process = cms.Process('NTUPLE')

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
	#,SkipEvent = cms.untracked.vstring('ProductNotFound')
)
# import of standard configurations
process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.suppressInfo = cms.untracked.vstring( "mkcands" )
process.MessageLogger.suppressWarning = cms.untracked.vstring( "mkcands" )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Input source
process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32( 0 ), #with 11976 Processing run: 201707 lumi: 281 event: 383901681
                            fileNames = cms.untracked.vstring()
)


sourceFiles = cms.untracked.vstring('file:00054BD0-5668-E211-8091-00215E21DC7E.root')
process.PoolSource.fileNames = sourceFiles


process.source.inputCommands = cms.untracked.vstring(
        "keep *",
        "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__RECO",
        "drop *_MEtoEDMConverter_*_*"
	)

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32( 200 )
	)

#Output size of CRAB jobs ~200MB usually works well. (max 300-500 Mb according to Cesare)

process.load('Configuration.Geometry.GeometryIdeal_cff') # 53x
process.load("Configuration.StandardSequences.GeometryExtended_cff") # from Lucia
process.load("Configuration.StandardSequences.Reconstruction_cff") # from Lucia
process.load("Configuration.StandardSequences.MagneticField_cff") # for using TransientTrackBuilder
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff") # for using TransientTrackBuilder
#process.GlobalTag.globaltag = 'FT_53_V6_AN3::All' # for using TransientTrackBuilder
process.GlobalTag.globaltag = 'START53_V19F::All'
#process.GlobalTag.globaltag = 'START53_V7C::All'
#process.GlobalTag.globaltag = 'FT_R_53_V18::All' #Global tag for 2012B data
#process.GlobalTag.globaltag = 'FT53_V21A_AN6::All' #Global tag for 2012C data

process.load('Configuration/EventContent/EventContent_cff')
#
#  Load common sequences
#
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

###############################################################################################
################################## good collisions ############################################

#### 44x
#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                                      vertexCollection = cms.InputTag('offlinePrimaryVertices'),
#                                                      minimumNDOF = cms.uint32(4) ,
#                                                      maxAbsZ = cms.double(24),
#                                                      maxd0 = cms.double(2)
#                                           )

## 53x

pvSelection = cms.PSet(
        minNdof = cms.double( 4. )
        , maxZ    = cms.double( 24. )
        , maxRho  = cms.double( 2. )
)

process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter", # checks for fake PVs automatically
                                                  filterParams = pvSelection,
                                                  filter       = cms.bool( False ), # use only as producer
                                                  src          = cms.InputTag( 'offlinePrimaryVertices' )
                                          )

process.primaryVertexFilter = process.goodOfflinePrimaryVertices.clone( filter = True )


process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  #debugOn = cms.untracked.bool(True),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                          )


# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
patMuons.embedTrack = cms.bool(True)
patMuons.embedPickyMuon = cms.bool(False)
patMuons.embedTpfmsMuon = cms.bool(False)

# Prune generated particles to muons and their parents
process.genMuons = cms.EDProducer("GenParticlePruner",
                                  src = cms.InputTag("genParticles"),
                                  select = cms.vstring(
                                          "drop  *  ",                     # this is the default
                                          "++keep abs(pdgId) = 13",        # keep muons and their parents
                                          "drop pdgId == 21 && status = 2" # remove intermediate qcd spam carrying no flavour info
                                  )
 )

########################################################################
##Trigger Matching for Muons
########################################################################

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import  addMCinfo, useExistingPATMuons, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution, addDiMuonTriggers
    # with some customization
if MC:
        addMCinfo(process)
        # since we match inner tracks, keep the matching tight and make it one-to-one
        process.muonMatch.maxDeltaR = 0.05
        process.muonMatch.resolveByMatchQuality = True

addDiMuonTriggers(process)
useExistingPATMuons(process,'cleanPatMuons',addL1Info=False)
changeTriggerProcessName(process, 'HLT')
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
useL1MatchingWindowForSinglets(process)

process.muonL1Info.maxDeltaR     = 0.3
process.muonL1Info.fallbackToME1 = True
process.muonMatchHLTL1.maxDeltaR     = 0.3
process.muonMatchHLTL1.fallbackToME1 = True
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL2.maxDPtRel = 10.0
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0

from PhysicsTools.PatAlgos.tools.trackTools import *

########################################################################
##Tracks - adding tracks refitted with different mass
from RecoTracker.TrackProducer.TrackRefitters_cff import *
from TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi import *
#process.RungeKuttaTrackerPropagatorForMuons = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone( Mass = cms.double(0.10565837), ComponentName = cms.string('RungeKuttaTrackerPropagatorForMuons') )
#process.refittedGeneralTracksMuon = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForMuons" )
process.RungeKuttaTrackerPropagatorForPions = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone( Mass = cms.double(0.13957), ComponentName = cms.string('RungeKuttaTrackerPropagatorForPions') )
process.refittedGeneralTracksPion = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForPions" )

#######################################################################
##Pions

makeTrackCandidates( process,                                # patAODTrackCands
                     label = 'TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
                     tracks = cms.InputTag('generalTracks'), # input track collection
                     #tracks = cms.InputTag('refittedGeneralTracksMuon'), # input track collection               // AP changed from generalTracks
                     #tracks = cms.InputTag('refittedGeneralTracksPion'), # input track collection               // AP changed from generalTracks
                     #particleType = 'mu+',                   # particle type (for assigning a mass) # not working, everything is a pion
                     particleType = 'pi+',                   # particle type (for assigning a mass) # not working, everything is a pion
                     preselection = 'pt > 0.5',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
                     #selection = 'pt > 0.35',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     #selection = 'p > 0.5',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     selection = 'pt > 0.5 && p > 0.5',     # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
		             isolation = {},                         # Isolations to use ('source':deltaR; set to {} for None)
                     isoDeposits = [],
                     mcAs = None                           # Replicate MC match as the one used for Muons
             );                                    # you can specify more than one collection for this

l1cands = getattr(process, 'patTrackCands')
l1cands.addGenMatch = False

###################################################
##Kaons

#process.RungeKuttaTrackerPropagator.Mass = cms.double(0.493677)
process.RungeKuttaTrackerPropagatorForKaons = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone(
        Mass = cms.double(0.493677), ComponentName = cms.string('RungeKuttaTrackerPropagatorForKaons') )
process.refittedGeneralTracksKaon = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForKaons" )

makeTrackCandidates( process,                                        # patAODTrackCands
                     label = 'TrackKaonCands',                       # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
                     #tracks = cms.InputTag('refittedGeneralTracksKaon'), # input track collection               // AP changed from generalTracks
                     tracks = cms.InputTag('generalTracks'), # input track collection               // AP changed from generalTracks
                     particleType = 'K+',                            # particle type (for assigning a mass)  // AP changed from pi to K # not working, everything is a pion
                     #particleType = 'pi+',                            # particle type (for assigning a mass)  // AP changed from pi to K # not working, everything is a pion
                     #particleType = 'mu+',                            # particle type (for assigning a mass)  // AP changed from pi to K # not working, everything is a pion
                     preselection = 'pt > 0.5',                      # preselection cut on candidates. Only methods of 'reco::Candidate' are available
                     #selection = 'pt > 0.35',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     #selection = 'p > 0.5',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     selection = 'pt > 0.5 && p > 0.5',     # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
		             isolation = {},                                 # Isolations to use ('source':deltaR; set to {} for None)
                     isoDeposits = [],
                     #mcAs = 'muon'                                   # Replicate MC match as the one used for Muons # AP "=None"  ??
                     mcAs = None                                    # Replicate MC match as the one used for Muons
             );                                                      # you can specify more than one collection for this

l1cands = getattr(process, 'patTrackKaonCands')
l1cands.addGenMatch = False

###################################################
##Filtering the event

#process.PATfilter = cms.EDFilter("X3872FilterPAT")
process.PATfilter = cms.EDFilter("X4140FilterPAT")


###################################################
##The Producers

process.psitomumu = cms.EDProducer("MuMuProducerPAT",
                                 HLTriggerResults = cms.untracked.InputTag("DiMuonCandidates"),
                                 inputGEN  = cms.untracked.InputTag("genParticles"),
                                 VtxSample   = cms.untracked.string('offlinePrimaryVertices'),

                                 JPsiMassCuts = cms.untracked.vdouble((2.95,3.25)),
                                 PsiMassCuts = cms.untracked.vdouble((3.45,3.85)),

                                 DoDataAnalysis = cms.untracked.bool( True ),
                                 DoMonteCarloTree = cms.untracked.bool( False ),

                                 addCommonVertex = cms.untracked.bool( False ),
                                 resolvePileUpAmbiguity = cms.untracked.bool( False ),
                                 addMuMulessPrimaryVertex = cms.untracked.bool( True ),

                                 addMCTruth = cms.untracked.bool( False ),
                                 MonteCarloParticleId = cms.untracked.int32(20443),
                                 MonteCarloExclusiveDecay = cms.untracked.bool( False ),
                                 MonteCarloMotherId = cms.untracked.int32( 511 ),
                                 MonteCarloDaughtersN = cms.untracked.int32( 3 ), # 3 for exclusive B0->psi'KPi

                                 MinNumMuPixHits = cms.untracked.int32(1),
                                 MinNumMuSiHits = cms.untracked.int32(8),
                                 MaxMuNormChi2 = cms.untracked.double(7),
                                 MaxMuD0 = cms.untracked.double(10.0),

                                 Debug_Output = cms.untracked.bool(True), # true
                                 ##
                                 ##  use the correct trigger path
                                 ##
                                 TriggersForMatching = cms.untracked.vstring(
                                         #2012 displaced J/psi = Alessandra
                                         "HLT_DoubleMu4_Jpsi_Displaced",
                    					 "HLT_Dimuon8_Jpsi",
                                 ),

                                 FiltersForMatching = cms.untracked.vstring(
                                         "hltDisplacedmumuFilterDoubleMu4Jpsi",
                                         "hltVertexmumuFilterDimuon8Jpsi")

                         )

process.psitomumu = cms.EDProducer("DiMuonRootupler",
                                 dimuons = cms.untracked.InputTag("genParticles"),
                                 primaryVertices = cms.untracked.InputTag("genParticles"),
                                 TriggerResults = cms.untracked.InputTag("genParticles"),

                                 dimuon_pdgid = cms.untracked.int32( 511 ),
                                 dimuon_mass_cuts = cms.untracked.vdouble((2.95,3.25)),
                                 isMC = cms.untracked.bool( True ),
                                 OnlyBest = cms.untracked.bool( True ),
                                 OnlyGen = cms.untracked.bool( True ),
                                 TriggersForMatching = cms.untracked.vstring(
                                         #2012 displaced J/psi = Alessandra
                                         "HLT_DoubleMu4_Jpsi_Displaced",
                    					 "HLT_Dimuon8_Jpsi",
                                 ),

                                 FiltersForMatching = cms.untracked.vstring(
                                         "hltDisplacedmumuFilterDoubleMu4Jpsi",
                                         "hltVertexmumuFilterDimuon8Jpsi")

                         )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('set_below.root')
)
if (not MC) :
    process.TFileService.fileName = cms.string('MuOniaParked_Run2012B_MuMuKKPAT_ntpl.root')
elif MC :
    process.TFileService.fileName = cms.string('BsToJPsiPhi_Summer12DR53X_MuMuPiKPAT_ntpl.root')



# turn off MC matching for the process
from PhysicsTools.PatAlgos.tools.coreTools import *
# old: removeMCMatching(process, ['All'], outputInProcess = False)
removeMCMatching(process,['All'],"",None,[])

process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetPartons)

process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)

process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('onia2MuMuPAT.root'),
                               outputCommands = cms.untracked.vstring('drop *',
                                             #'keep *_genMuons_*_Onia2MuMuPAT', # generated muons and parents
                                             'keep patMuons_patMuonsWithTrigger_*_NTUPLE', # All PAT muons including general tracks and matches to triggers
                                                              )
                       )

process.filter = cms.Sequence(
        process.goodOfflinePrimaryVertices
        + process.primaryVertexFilter
        + process.noscraping
)
#44x process.filter = cms.Sequence(process.primaryVertexFilter+process.noscraping)

process.ntup = cms.Path(
        #process.refittedGeneralTracksPion *
        #process.refittedGeneralTracksMuon *
        #process.refittedGeneralTracksKaon *
        #process.offlineBeamSpot * process.TrackRefitter * process.dedxHitInfo
        #process.dedxHarmonic2Kaon *
        process.offlineBeamSpot #* process.dedxHitInfo
        * process.filter
        * process.patDefaultSequence
        * process.patMuonsWithTriggerSequence
        * process.PATfilter
        * process.mkcands
)

process.schedule = cms.Schedule(process.ntup)

# rsync -vut --existing test/runMuMuKKPAT_dataOrMC.py semrat@lxplus.cern.ch:/afs/cern.ch/user/s/semrat/scratch0/CMSSW_5_3_22/src/X4140/MuMuKKPAT/test/runMuMuKKPAT_dataOrMC.py
