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

#MC = False
MC = True
if MC :
        official = True
MCMotherId = 531

if MCMotherId == 511 :
    MCExclusiveDecay = True
elif MCMotherId == 531 :
    MCExclusiveDecay = False

# Input source
process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32( 0 ), #with 11976 Processing run: 201707 lumi: 281 event: 383901681
                            fileNames = cms.untracked.vstring()
)

if (not MC) :
    sourceFiles = cms.untracked.vstring('/store/data/Run2012B/MuOniaParked/AOD/22Jan2013-v1/20001/EE35A843-3E69-E211-9292-00215E2226AC.root')
elif MC :
    readFiles = cms.untracked.vstring()
    sourceFiles.extend( [
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/00AC758F-4558-E611-AF39-02163E014BCA.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/188EC304-1F58-E611-B1E5-FA163E7918CE.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/3C06A515-2258-E611-98CA-FA163E9AE9D0.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/4018D9C2-2058-E611-A698-02163E00F483.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/502525F1-1A58-E611-9A06-FA163EBE010E.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/5264DCBB-2058-E611-B18C-FA163E99F0C3.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/56C0264B-1E58-E611-ADF7-FA163E881C27.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/729EFC9E-1D58-E611-9D8F-FA163E28C01D.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/74D6D987-4558-E611-84FC-FA163E7846B2.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/842E1AA9-1C58-E611-B214-02163E00E9D1.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/982A5A01-1F58-E611-81E1-FA163E918327.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/A8201DEB-1B58-E611-8D87-FA163E2ADA19.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/AEF65ABB-1858-E611-B346-0CC47A07FA1E.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/BC358A48-1958-E611-B192-FA163ECBD6A3.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/CA04869D-1D58-E611-B451-FA163EA543D3.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/DABCF18D-1C58-E611-9FA6-001E67E95B68.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/E2FBF0EE-1A58-E611-B579-FA163E784361.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/50000/E6960119-2358-E611-9108-FA163E6265F2.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/008D511F-4458-E611-B9C6-FA163E28C01D.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/0278351E-1058-E611-B38C-FA163E437AF6.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/04AAC8B0-1558-E611-9BB6-FA163EC5AB84.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/0C3ADC4E-1258-E611-94B1-FA163EB99C8C.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/0C5F248A-1958-E611-8518-001E67579498.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/167F9524-1658-E611-A520-FA163E954418.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/18D52F20-FB57-E611-BF0E-008CFAFBE8F2.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/1C355DE5-1758-E611-9EF6-FA163E97FA4D.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/1C7A3DDD-1D58-E611-8168-FA163E4CCB09.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/1E0979EF-4358-E611-9B2B-FA163EE20B5D.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/208FFFDB-1658-E611-9BD0-FA163E5CC901.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/224454F7-1258-E611-8179-FA163E824933.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/3083217B-1A58-E611-959F-FA163E665D2C.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/4A9CFFB3-1158-E611-B2BE-02163E014D37.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/62AB1A60-FB57-E611-BC23-F04DA275BF11.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/86364225-FB57-E611-80AA-008CFAFC043A.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/865EF4AD-1558-E611-8E0D-FA163EC0B67E.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/86A26E4D-1C58-E611-9D00-FA163EE6D545.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/8CD8CCF4-4358-E611-BD05-7845C4FC3AC1.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/94FC71A1-1B58-E611-8B86-FA163E7CDD49.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/9CB277BD-1358-E611-982D-FA163E6CDEDC.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/A0846A7B-1458-E611-B381-0025903D0B32.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/A23AD6DC-1658-E611-860B-FA163E83B5F4.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/B09E5782-1958-E611-8689-FA163E023B23.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/B6670C17-1058-E611-92C3-FA163EA1151A.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/CAC24AF8-1758-E611-9696-02163E014AFC.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/CEC6E26C-0D58-E611-90BC-02163E016486.root',
    '/store/mc/Summer12DR53X/Y4140ToJpsiPhi_BFilter_TuneCUEP8M1_8TeV-Pythia8/AODSIM/PU_RD2_START53_V19F-v2/60000/EEF0A85B-0E58-E611-B261-FA163E8D3CEB.root' ] );


process.PoolSource.fileNames = sourceFiles;


process.source.inputCommands = cms.untracked.vstring(
        "keep *",
        "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__RECO",
        "drop *_MEtoEDMConverter_*_*"
	)

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32( -1 ) # 256Kb in 2' for 100 events, 1Mb in 7' for 1k events, 6Mb in 50' for 8650 events, 11Mb in 66' for 10k events, 100Mb in 14h for 150k events, 1.4Gb in 4 days for 1.2M events of official MC
        #input = cms.untracked.int32( 1000 ) # 310Kb in 3' for 1k events of private MC
        #input = cms.untracked.int32( 100 ) # = 20Mb in 2h for 15k events, 2Mb in 10' for 1k events of Run2012C/MuOniaParked/AOD/22Jan2013-v1
	#input = cms.untracked.int32( 1000 ) # = 3Mb for 6546 events, 85Kb for 100, 800kb for 1k events of BsToPsiMuMu
	#input = cms.untracked.int32( 24000 ) # = 870Kb # timeout after 24500 for Run2012A/MuOnia
	#input = cms.untracked.int32( -1 ) # = 5718Kb # timeout after 3700 for Run2012A/MuOnia
	)

#Output size of CRAB jobs ~200MB usually works well. (max 300-500 Mb according to Cesare)

process.load('Configuration.Geometry.GeometryIdeal_cff') # 53x

process.load("Configuration.StandardSequences.GeometryExtended_cff") # from Lucia
process.load("Configuration.StandardSequences.Reconstruction_cff") # from Lucia

process.load("Configuration.StandardSequences.MagneticField_cff") # for using TransientTrackBuilder
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff") # for using TransientTrackBuilder
#process.GlobalTag.globaltag = 'FT_53_V6_AN3::All' # for using TransientTrackBuilder
#process.GlobalTag.globaltag = 'START53_V19F::All'
#process.GlobalTag.globaltag = 'START53_V7C::All'
process.GlobalTag.globaltag = 'FT_R_53_V18::All' #Global tag for 2012B data
#process.GlobalTag.globaltag = 'FT53_V21A_AN6::All' #Global tag for 2012C data

process.load('Configuration/EventContent/EventContent_cff')
#
#  Load common sequences
#
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

####################################################################################
##################################good collisions############################################

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
######## adding tracks refitted with different mass
from RecoTracker.TrackProducer.TrackRefitters_cff import *
from TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi import *
#process.RungeKuttaTrackerPropagatorForMuons = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone( Mass = cms.double(0.10565837), ComponentName = cms.string('RungeKuttaTrackerPropagatorForMuons') )
#process.refittedGeneralTracksMuon = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForMuons" )
process.RungeKuttaTrackerPropagatorForPions = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone( Mass = cms.double(0.13957), ComponentName = cms.string('RungeKuttaTrackerPropagatorForPions') )
process.refittedGeneralTracksPion = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForPions" )
makeTrackCandidates( process,                                # patAODTrackCands
                     label = 'TrackCands',                   # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
                     tracks = cms.InputTag('generalTracks'), # input track collection
                     #tracks = cms.InputTag('refittedGeneralTracksMuon'), # input track collection               // AP changed from generalTracks
                     #tracks = cms.InputTag('refittedGeneralTracksPion'), # input track collection               // AP changed from generalTracks
                     #particleType = 'mu+',                   # particle type (for assigning a mass) # not working, everything is a pion
                     particleType = 'pi+',                   # particle type (for assigning a mass) # not working, everything is a pion
                     preselection = 'pt > 0.4',              # preselection cut on candidates. Only methods of 'reco::Candidate' are available
                     #selection = 'pt > 0.35',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     #selection = 'p > 0.5',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     selection = 'pt > 0.4 && p > 0.5',     # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
		     isolation = {},                         # Isolations to use ('source':deltaR; set to {} for None)
                     isoDeposits = [],
                     mcAs = None                           # Replicate MC match as the one used for Muons
             );                                    # you can specify more than one collection for this

l1cands = getattr(process, 'patTrackCands')
l1cands.addGenMatch = False

######## adding tracks refitted with Kaon mass
#process.RungeKuttaTrackerPropagator.Mass = cms.double(0.493677)
process.RungeKuttaTrackerPropagatorForKaons = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone(
        Mass = cms.double(0.493677), ComponentName = cms.string('RungeKuttaTrackerPropagatorForKaons') )
process.refittedGeneralTracksKaon = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForKaons" )
###################################################
makeTrackCandidates( process,                                        # patAODTrackCands
                     label = 'TrackKaonCands',                       # output collection will be 'allLayer0TrackCands', 'allLayer1TrackCands', 'selectedLayer1TrackCands'
                     #tracks = cms.InputTag('refittedGeneralTracksKaon'), # input track collection               // AP changed from generalTracks
                     tracks = cms.InputTag('generalTracks'), # input track collection               // AP changed from generalTracks
                     particleType = 'K+',                            # particle type (for assigning a mass)  // AP changed from pi to K # not working, everything is a pion
                     #particleType = 'pi+',                            # particle type (for assigning a mass)  // AP changed from pi to K # not working, everything is a pion
                     #particleType = 'mu+',                            # particle type (for assigning a mass)  // AP changed from pi to K # not working, everything is a pion
                     preselection = 'pt > 0.4',                      # preselection cut on candidates. Only methods of 'reco::Candidate' are available
                     #selection = 'pt > 0.35',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     #selection = 'p > 0.5',                 # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
                     selection = 'pt > 0.4 && p > 0.5',     # Selection on PAT Layer 1 objects ('selectedLayer1TrackCands')
		     isolation = {},                                 # Isolations to use ('source':deltaR; set to {} for None)
                     isoDeposits = [],
                     #mcAs = 'muon'                                   # Replicate MC match as the one used for Muons # AP "=None"  ??
                     mcAs = None                                    # Replicate MC match as the one used for Muons
             );                                                      # you can specify more than one collection for this

l1cands = getattr(process, 'patTrackKaonCands')
l1cands.addGenMatch = False

process.load("RecoTracker.DeDx.dedxHarmonic2_cfi")
process.dedxHarmonic2Kaon = RecoTracker.DeDx.dedxHarmonic2_cfi.dedxHarmonic2.clone (
        tracks = 'refittedGeneralTracksKaon',
        trajectoryTrackAssociation = 'refittedGeneralTracksKaon'
)
# dE/dx hits
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
#process.load("RecoTracker.TrackProducer.TrackRefitters_cff") #already imported above
#process.TrackRefitter.src = 'generalTracks'
#process.TrackRefitter.src = 'refittedGeneralTracksPion'

#process.dedxHitInfo = cms.EDProducer("HSCPDeDxInfoProducer",
#                                     #tracks = cms.InputTag("TrackRefitter"),
#                                     #trajectoryTrackAssociation = cms.InputTag("TrackRefitter"),
#                                     tracks = cms.InputTag("refittedGeneralTracksPion"),
#                                     trajectoryTrackAssociation = cms.InputTag("refittedGeneralTracksPion"),
#
#                                     UseStrip  = cms.bool(True),
#                                     UsePixel  = cms.bool(True),
#                                     MeVperADCStrip = cms.double(3.61e-06*265),
#                                     MeVperADCPixel = cms.double(3.61e-06),
#
#                                     UseCalibration = cms.bool(False),
#                                     calibrationPath = cms.string("/afs/cern.ch/user/q/querten/workspace/public/dEdx/CMSSW_5_2_4/src/dEdx/ppGridProject/Gains.root"),
#                                     ShapeTest = cms.bool(True),
#                             )
#
#process.dedxHitInfoKaon = cms.EDProducer("HSCPDeDxInfoProducer",
#                                         tracks = cms.InputTag("refittedGeneralTracksKaon"),
#                                         trajectoryTrackAssociation = cms.InputTag("refittedGeneralTracksKaon"),
#
#                                         UseStrip  = cms.bool(True),
#                                         UsePixel  = cms.bool(True),
#                                         MeVperADCStrip = cms.double(3.61e-06*265),
#                                         MeVperADCPixel = cms.double(3.61e-06),
#
#                                         UseCalibration = cms.bool(False),
#                                         calibrationPath = cms.string("/afs/cern.ch/user/q/querten/workspace/public/dEdx/CMSSW_5_2_4/src/dEdx/ppGridProject/Gains.root"),
#                                         ShapeTest = cms.bool(True),
#                                 )


#process.PATfilter = cms.EDFilter("X3872FilterPAT")
process.PATfilter = cms.EDFilter("X4140FilterPAT")

process.mkcands = cms.EDAnalyzer("mumukk",
                                 HLTriggerResults = cms.untracked.InputTag("TriggerResults","","HLT"),
                                 inputGEN  = cms.untracked.InputTag("genParticles"),
                                 VtxSample   = cms.untracked.string('offlinePrimaryVertices'),
                                 SameSign = cms.untracked.bool(False),
                                 DoMonteCarloTree = cms.untracked.bool( MC ),
                                 MonteCarloParticleId = cms.untracked.int32(20443),
                                 MonteCarloExclusiveDecay = cms.untracked.bool( MCExclusiveDecay ),
                                 MonteCarloMotherId = cms.untracked.int32( MCMotherId ),
                                 MonteCarloDaughtersN = cms.untracked.int32( 3 ), # 3 for exclusive B0->psi'KPi
                                 #
                                 DoMuMuMassConstraint = cms.untracked.bool(True),
                                 #SkipJPsi = cms.untracked.bool(True),
                                 SkipJPsi = cms.untracked.bool(False),
                                 SkipPsi2S = cms.untracked.bool(False), # SEMRA ask closing or not
                                 MinNumMuPixHits = cms.untracked.int32(1),
                                 MinNumMuSiHits = cms.untracked.int32(8),
                                 MaxMuNormChi2 = cms.untracked.double(7),
                                 MaxMuD0 = cms.untracked.double(10.0),
                                 sharedFraction = cms.untracked.double(0.5),

                                 MinJPsiMass = cms.untracked.double(2.8), # SEMRA changed
                                 MaxJPsiMass = cms.untracked.double(3.4), # SEMRA changed
                				 MinPhiMass = cms.untracked.double (0.97), # SEMRA added
                 				 MaxPhiMass = cms.untracked.double (1.07), # SEMRA added
                				 MaxJPsiPhiXMass = cms.untracked.double (5.1), # SEMRA added
                				 MinJPsiPhiB0Mass = cms.untracked.double (0.0), # SEMRA added
                				 MaxJPsiPhiB0Mass = cms.untracked.double (0.0), # SEMRA added

                                 MinNumTrSiHits = cms.untracked.int32(4),
                                 MinTrPt = cms.untracked.double(0.350),
                                 Chi2NDF_Track =  cms.untracked.double(7.0),
                				 # Delta R
                				 MaxMuMuTrackDR = cms.untracked.double(1.5),
                                 MaxXCandTrackDR = cms.untracked.double(1.5),
                                 UseXDr = cms.untracked.bool(True),

                                 resolvePileUpAmbiguity = cms.untracked.bool(False),
                                 addMuMulessPrimaryVertex = cms.untracked.bool(True),
                                 #addMuMulessPrimaryVertex = cms.untracked.bool(False),
                                 addXlessPrimaryVertex = cms.untracked.bool(True),
                                 Debug_Output = cms.untracked.bool(False), # true
                                 ##
                                 ##  use the correct trigger path
                                 ##
                                 TriggersForMatching = cms.untracked.vstring(
                                         #2012 displaced J/psi = Alessandra
                                         "HLT_DoubleMu4_Jpsi_Displaced_v9",
                    					 "HLT_DoubleMu4_Jpsi_Displaced_v10",
                    					 "HLT_DoubleMu4_Jpsi_Displaced_v11",
                                         "HLT_DoubleMu4_Jpsi_Displaced_v12",
                    					 "HLT_Dimuon8_Jpsi_v4",
                    					 "HLT_Dimuon8_Jpsi_v5",
                                         "HLT_Dimuon8_Jpsi_v6",
                                         "HLT_Dimuon8_Jpsi_v7",
                                 ),
				 FiltersForMatching = cms.untracked.vstring(
                                         "hltDisplacedmumuFilterDoubleMu4Jpsi",
                                         "hltDisplacedmumuFilterDoubleMu4Jpsi",
                                         "hltDisplacedmumuFilterDoubleMu4Jpsi",
                                         "hltDisplacedmumuFilterDoubleMu4Jpsi"
					                     "hltVertexmumuFilterDimuon8Jpsi",
                                         "hltVertexmumuFilterDimuon8Jpsi",
                                         "hltVertexmumuFilterDimuon8Jpsi",
                                         "hltVertexmumuFilterDimuon8Jpsi",
                                         "hltVertexmumuFilterDimuon8Jpsi",
                                )


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
