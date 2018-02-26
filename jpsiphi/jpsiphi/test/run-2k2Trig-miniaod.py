import FWCore.ParameterSet.Config as cms
process = cms.Process('PSIKK')

input_file = "file:049F2D32-26F2-E711-A162-00266CFFC664.root"

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v1')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v2') #F

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

process.load("jpsiphi.jpsiphi.slimmedMuonsTriggerMatcher2017_cfi")
# process.load("jpsiphi.jpsiphi.slimmedTracksTriggerMatcher2017_cfi")

charmoniumHLT = [
#Phi
'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi',
#JPsi
# 'HLT_DoubleMu4_JpsiTrkTrk_Displaced',
# 'HLT_DoubleMu4_JpsiTrk_Displaced',
# 'HLT_DoubleMu4_Jpsi_Displaced',
# 'HLT_DoubleMu4_3_Jpsi_Displaced',
# 'HLT_Dimuon20_Jpsi_Barrel_Seagulls',
# 'HLT_Dimuon25_Jpsi',
]

hltList = charmoniumHLT #muoniaHLT

hltpaths = cms.vstring(hltList)

hltpathsV = cms.vstring([h + '_v*' for h in hltList])

filters = cms.vstring(
                                #HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi
                                'hltDoubleMu2JpsiDoubleTrkL3Filtered',
                                'hltDoubleTrkmumuFilterDoubleMu2Jpsi',
                                'hltJpsiTkTkVertexFilterPhiDoubleTrk1v2',
                                )

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring(hltpathsV),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )


process.softMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuonsWithTrigger'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
   ),
   filter = cms.bool(True)
)


process.JPsi2MuMuPAT = cms.EDProducer('DiMuonProducerPAT',
        muons                       = cms.InputTag('softMuons'),
        primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("2.95 < mass && mass < 3.25 && charge==0"),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(False),
        resolvePileUpAmbiguity      = cms.bool(True),
        HLTFilters                  = filters
)


process.JPsi2MuMuFilter = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("JPsi2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.95 < mass && mass < 3.25 && userFloat('vProb') > 0.01 && pt > 2.0"),
      do_trigger_match    = cms.bool(True),
      HLTFilters          = filters
)

process.DiMuonCounterJPsi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("JPsi2MuMuFilter"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)


process.DiTrakHLTProducer = cms.EDProducer('DiTrackHLTProducer',
    PFCandidates        = cms.InputTag("packedPFCandidates"),
    TriggerInput        = cms.InputTag("slimmedPatTrigger"),
    TrakTrakMassCuts    = cms.vdouble(0.6,1.3),
    MassTraks           = cms.vdouble(kaonmass,kaonmass),
    OnlyBest            = cms.bool(False),
    TTCandidate_name    = cms.string("DiTrakCandidate"),
    TTTrigger_name      = cms.string("DiTrigCandidate"),
    HLTFilters          = filters,
)

process.rootuple = cms.EDAnalyzer('DiTrakHLTRootupler',
    ditraks             = cms.InputTag('PsiPhiProducer','DiMuonDiTrakCandidates'),
    ditrigs             = cms.InputTag("PsiPhiFitter","DiMuonDiTrakCandidatesRef"),
    primaryVertices     = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults      = cms.InputTag("TriggerResults", "", "HLT"),
    TrakTrakMassCuts    = cms.vdouble(0.6,1.3),
    isMC                = cms.bool(False),
    OnlyBest            = cms.bool(False),
    HLTs                = hltpaths,
    Filters             = filters,
)
#
# process.rootupleMuMu = cms.EDAnalyzer('DiMuonRootupler',
#                           dimuons = cms.InputTag("JPsi2MuMuFilter"),
#                           muons = cms.InputTag("replaceme"),
#                           primaryVertices = cms.InputTag("offlinePrimaryVertices"),
#                           TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
#                           dimuon_pdgid = cms.uint32(443),
#                           dimuon_mass_cuts = cms.vdouble(2.5,3.5),
#                           isMC = cms.bool(False),
#                           OnlyBest = cms.bool(False),
#                           OnlyGen = cms.bool(False),
#                           HLTs = hltpaths
#                           )

process.p = cms.Path(process.triggerSelection *
                     process.slimmedMuonsWithTriggerSequence *
                     # process.slimmedPFCandsWithTriggerSequence *
                     process.softMuons *
                     process.JPsi2MuMuPAT *
                     process.JPsi2MuMuFilter*
                     process.DiMuonCounterJPsi*
                     # process.PsiPhiProducer *
                     # process.PsiPhiFitter *
                     process.DiTrakHLTProducer *
                     process.rootuple)
                     # process.rootupleMuMu)# * process.Phi2KKPAT * process.patSelectedTracks *process.rootupleKK)
