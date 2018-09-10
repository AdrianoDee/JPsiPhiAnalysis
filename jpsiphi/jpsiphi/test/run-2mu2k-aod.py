import FWCore.ParameterSet.Config as cms
process = cms.Process('PSIKK')

gen_file = "file:32B83273-030F-E811-9105-E0071B7AF7C0.root"
input_file = "file:006425F0-6DED-E711-850C-0025904C66E8.root"
mc_file = "file:py8_JPsiMM_EvtGen_13TeV_TuneCP5_cfi.root"
mc_file = "file:02CA3723-CEF3-E711-B1CC-4C79BA1810EF.root"
input_file = mc_file #gen_file

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
process.load('Configuration/EventContent/EventContent_cff')

process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
from PhysicsTools.PatAlgos.tools.trackTools import *
from RecoTracker.TrackProducer.TrackRefitters_cff import *
from TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi import *
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v1')
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v2') #F
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v10')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(input_file)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('rootuple-2017-dimuonditrak.root'),
)

kaonmass = 0.493677
pionmass = 0.13957061

pvSelection = cms.PSet(
        minNdof = cms.double( 4. )
        , maxZ    = cms.double( 24. )
        , maxRho  = cms.double( 2. )
)

process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
patMuons.embedTrack = cms.bool(True)
patMuons.embedPickyMuon = cms.bool(False)
patMuons.embedTpfmsMuon = cms.bool(False)

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import  addMCinfo, useExistingPATMuons, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution, addDiMuonTriggers
    # with some customization

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

###################################################################################################
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

#################################################################
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

####################################################################################################
##Tracks - adding tracks refitted with different mass

#process.RungeKuttaTrackerPropagatorForMuons = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone( Mass = cms.double(0.10565837), ComponentName = cms.string('RungeKuttaTrackerPropagatorForMuons') )
#process.refittedGeneralTracksMuon = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForMuons" )
process.RungeKuttaTrackerPropagatorForPions = TrackingTools.MaterialEffects.RungeKuttaTrackerPropagator_cfi.RungeKuttaTrackerPropagator.clone( Mass = cms.double(0.13957), ComponentName = cms.string('RungeKuttaTrackerPropagatorForPions') )
process.refittedGeneralTracksPion = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone( Propagator = "RungeKuttaTrackerPropagatorForPions" )

process.load("jpsiphi.jpsiphi.slimmedMuonsTriggerMatcher2017_cfi")
# process.load("jpsiphi.jpsiphi.slimmedTracksTriggerMatcher2017_cfi")

charmoniumHLT = [
#Phi
'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi',
#JPsi
'HLT_DoubleMu4_JpsiTrkTrk_Displaced',
'HLT_DoubleMu4_JpsiTrk_Displaced',
'HLT_DoubleMu4_Jpsi_Displaced',
'HLT_DoubleMu4_3_Jpsi_Displaced',
'HLT_Dimuon20_Jpsi_Barrel_Seagulls',
'HLT_Dimuon25_Jpsi',
]

hltList = charmoniumHLT #muoniaHLT

hltpaths = cms.vstring(hltList)

hltpathsV = cms.vstring([h + '_v*' for h in hltList])

filters = cms.vstring(
                                #HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi
                                #'hltDoubleMu2JpsiDoubleTrkL3Filtered',
                                'hltDoubleTrkmumuFilterDoubleMu2Jpsi',
                                'hltJpsiTkTkVertexFilterPhiDoubleTrk1v2',
                                #HLT_DoubleMu4_JpsiTrkTrk_Displaced_v4
                                'hltDoubleMu4JpsiDisplacedL3Filtered'
                                'hltDisplacedmumuFilterDoubleMu4Jpsi',
                                'hltJpsiTkTkVertexFilterPhiKstar',
                                #HLT_DoubleMu4_JpsiTrk_Displaced_v12
                                #'hltDoubleMu4JpsiDisplacedL3Filtered',
                                'hltDisplacedmumuFilterDoubleMu4Jpsi',
                                #'hltJpsiTkVertexProducer',
                                #'hltJpsiTkVertexFilter',
                                #HLT_DoubleMu4_Jpsi_Displaced
                                #'hltDoubleMu4JpsiDisplacedL3Filtered',
                                #'hltDisplacedmumuVtxProducerDoubleMu4Jpsi',
                                'hltDisplacedmumuFilterDoubleMu4Jpsi',
                                #HLT_DoubleMu4_3_Jpsi_Displaced
                                #'hltDoubleMu43JpsiDisplacedL3Filtered',
                                'hltDisplacedmumuFilterDoubleMu43Jpsi',
                                #HLT_Dimuon20_Jpsi_Barrel_Seagulls
                                #'hltDimuon20JpsiBarrelnoCowL3Filtered',
                                'hltDisplacedmumuFilterDimuon20JpsiBarrelnoCow',
                                #HLT_Dimuon25_Jpsi
                                'hltDisplacedmumuFilterDimuon25Jpsis'
                                )

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring(hltpathsV),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.unpackPatTriggers = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
  patTriggerObjectsStandAlone = cms.InputTag( 'slimmedstPatTrigger' ), #selectedPatTrigger for MC
  triggerResults              = cms.InputTag( 'TriggerResults::HLT' ),
  unpackFilterLabels          = cms.bool( True )
)

# process.softMuons = cms.EDFilter('PATMuonSelector',
#    src = cms.InputTag('patMuonsWithTrigger'),
#    cut = cms.string('muonID(\"TMOneStationTight\")'
#                     ' && abs(innerTrack.dxy) < 0.3'
#                     ' && abs(innerTrack.dz)  < 20.'
#                     ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
#                     ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
#                     ' && innerTrack.quality(\"highPurity\")'
#    ),
#    filter = cms.bool(True)
# )

process.JPsi2MuMuPAT = cms.EDProducer('DiMuonProducerPAT',
        muons                       = cms.InputTag('patMuonsWithTrigger'),
        primaryVertexTag            = cms.InputTag('goodOfflinePrimaryVertices'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("2.95 < mass && mass < 3.25 && charge==0"),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(True),
        resolvePileUpAmbiguity      = cms.bool(True),
        HLTFilters                  = filters
)

process.JPsi2MuMuFilter = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("JPsi2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.95 < mass && mass < 3.25 && userFloat('vProb') > 0.01 && pt > 2.0"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters
)

process.PsiPhiProducer = cms.EDProducer('DiMuonDiTrakProducerAOD',
    DiMuon = cms.InputTag('JPsi2MuMuPAT'),
    PFCandidates = cms.InputTag('cleanPatTrackKaonCands'),
    #TriggerInput            = cms.InputTag("unpackPatTriggers"),
    TriggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    DiMuonMassCuts = cms.vdouble(2.95,3.25),      # J/psi mass window 3.096916 +/- 0.150
    TrakTrakMassCuts = cms.vdouble(1.0,1.04),  # phi mass window 1.019461 +/- .015
    DiMuonDiTrakMassCuts = cms.vdouble(4.0,5.8),            # b-hadron mass window
    MassTraks = cms.vdouble(kaonmass,kaonmass),         # traks masses
    OnlyBest  = cms.bool(False),
    Product = cms.string("DiMuonDiTrakCandidates"),
    Filters = filters,
    IsMC = cms.bool(True),
)

process.PsiPhiFitter = cms.EDProducer('DiMuonDiTrakKinematicFit',
    DiMuonDiTrak              = cms.InputTag('PsiPhiProducer','DiMuonDiTrakCandidates'),
    DiMuonMass                = cms.double(3.096916),              # J/psi mass in GeV
    DiMuonTrakTrakMassCuts    = cms.vdouble(4.0,5.8),            # b-hadron mass window
    MassTraks                 = cms.vdouble(kaonmass,kaonmass),         # traks masses
    Product                   = cms.string('DiMuonDiTrakCandidatesRef')
)

process.rootuple = cms.EDAnalyzer('DiMuonDiTrakRootupler',
    dimuonditrk_cand = cms.InputTag('PsiPhiProducer','DiMuonDiTrakCandidates'),
    dimuonditrk_rf_cand = cms.InputTag("PsiPhiFitter","DiMuonDiTrakCandidatesRef"),
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    primaryVertices = cms.InputTag("goodOfflinePrimaryVertices"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    isMC = cms.bool(True),
    OnlyBest = cms.bool(False),
    OnlyGen = cms.bool(False),
    Mother_pdg = cms.uint32(531),
    JPsi_pdg = cms.uint32(443),
    Phi_pdg = cms.uint32(333),
    HLTs = hltpaths,
    Filters = filters,
    TreeName = cms.string('JPsi Phi Tree')
)

process.rootupleMuMu = cms.EDAnalyzer('DiMuonRootupler',
                          dimuons = cms.InputTag("JPsi2MuMuFilter"),
                          muons = cms.InputTag("replaceme"),
                          primaryVertices = cms.InputTag("goodOfflinePrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          dimuon_pdgid = cms.uint32(443),
                          dimuon_mass_cuts = cms.vdouble(2.5,3.5),
                          isMC = cms.bool(True),
                          OnlyBest = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          HLTs = hltpaths
                          )

process.p = cms.Path(process.triggerSelection *
                     #process.slimmedMuonsWithTriggerSequence *
                     #process.unpackPatTriggers *
                     #process.softMuons *
                     process.patDefaultSequence *
                     process.patMuonsWithTriggerSequence *
                     process.JPsi2MuMuPAT *
                     process.JPsi2MuMuFilter*
                     process.PsiPhiProducer *
                     process.PsiPhiFitter *
                     process.rootuple *
                     process.rootupleMuMu)# * process.Phi2KKPAT * process.patSelectedTracks *process.rootupleKK)
