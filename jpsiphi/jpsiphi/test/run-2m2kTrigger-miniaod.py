import FWCore.ParameterSet.Config as cms
process = cms.Process('2mu2trkTrigger')

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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('rootuple-2017-ditraktrigger.root'),
)

kaonmass = 0.493677
pionmass = 0.13957061
muonmass = 0.1056583715;

process.load("jpsiphi.jpsiphi.slimmedMuonsTriggerMatcher2017_cfi")
# process.load("jpsiphi.jpsiphi.slimmedTracksTriggerMatcher2017_cfi")

charmoniumHLT = [
#Phi
'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi',
#JPsi
'HLT_DoubleMu4_JpsiTrkTrk_Displaced',
'HLT_DoubleMu4_JpsiTrk_Displaced'
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
                                'hltJpsiTkTkVertexFilterPhiDoubleTrk1v2',
                                'hltJpsiTkTkVertexFilterPhiKstar',
                                'hltJpsiTkVertexFilter'
                                )

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
    triggerConditions = cms.vstring(hltpathsV),
    hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
    l1tResults = cms.InputTag( "" ),
    throw = cms.bool(False)
)

process.unpackPatTriggers = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
      patTriggerObjectsStandAlone = cms.InputTag( 'slimmedPatTrigger' ),
      triggerResults              = cms.InputTag( 'TriggerResults::HLT' ),
      unpackFilterLabels          = cms.bool( True )
)

process.softMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuons'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
   ),
   filter = cms.bool(True)
)

process.DiMuonPAT = cms.EDProducer('DiMuonPAT',
    Muons                       = cms.InputTag('softMuons'),
    BeamSpot                     = cms.InputTag('offlineBeamSpot'),
    PrimaryVertex               = cms.InputTag('offlineSlimmedPrimaryVertices'),
    DiMuonCuts                  = cms.string("2.9 < mass && mass < 3.3 && charge==0"),
)


process.DiMuonCounterJPsi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("DiMuonPAT"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

process.DiTrakPAT = cms.EDProducer('DiTrakPAT',
    Traks                       = cms.InputTag('packedPFCandidates'),
    BeamSpot                    = cms.InputTag('offlineBeamSpot'),
    PrimaryVertex               = cms.InputTag('offlineSlimmedPrimaryVertices'),
    DiTrakCuts                  = cms.string("0.93 < mass && mass < 1.32 && charge==0"),
    TraksMasses                 = cms.vdouble(kaonmass,kaonmass)
)


process.DiTrakCounterPhi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("DiTrakPAT"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

process.DiMuonDiTrakPAT = cms.EDProducer('DiMuonDiTrakPAT',
    DiMuons             = cms.InputTag('DiMuonPAT'),
    DiTraks             = cms.InputTag('DiTrakPAT'),
    BeamSpot            = cms.InputTag('offlineBeamSpot'),
    PrimaryVertex       = cms.InputTag('offlineSlimmedPrimaryVertices'),
    DiMuonCuts          = cms.vdouble(2.9,3.3),
    DiTrakCuts          = cms.vdouble(0.9,1.3),
    DiMuonDiTrakCuts    = cms.vdouble(4.0,6.0),
    CandsMasses         = cms.vdouble(muonmass,muonmass,kaonmass,kaonmass)
)

process.Triggers  = cms.EDAnalyzer('TriggersByFilters',
    TriggerInput     = cms.InputTag("unpackPatTriggers"),
    TriggerResults   = cms.InputTag("TriggerResults", "", "HLT"),
    PrimaryVertex    = cms.InputTag("offlineSlimmedPrimaryVertices"),
    Filters          = filters,
    HLTs             = hltpaths
)

process.DiMuonRootuplerHLT = cms.EDAnalyzer('DiMuonRootupler',
    dimuons             = cms.InputTag("DiMuonPAT"),
    TriggerInput     = cms.InputTag("unpackPatTriggers"),
    primaryVertices     = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults      = cms.InputTag("TriggerResults", "", "HLT"),
    OnlyBest            = cms.bool(False),
    HLTs                = hltpaths,
    Filters          = filters
 )

process.DiTrakRootupler = cms.EDAnalyzer('DiTrakRootupler',
    ditraks             = cms.InputTag("DiTrakPAT"),
    TriggerInput     = cms.InputTag("unpackPatTriggers"),
    primaryVertices     = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults      = cms.InputTag("TriggerResults", "", "HLT"),
    OnlyBest            = cms.bool(False),
    HLTs                = hltpaths,
    Filters          = filters,
 )

process.DiMuonDiTrakRootupler = cms.EDAnalyzer('DiMuonDiTrakRootupler',
     dimuonditrks = cms.InputTag("DiMuonDiTrakPAT"),
     TriggerInput     = cms.InputTag("unpackPatTriggers"),
     primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
     TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
     OnlyBest = cms.bool(False),
     HLTs = hltpaths,
     Filters          = filters
  )

process.p = cms.Path(
                    process.triggerSelection *
                    process.unpackPatTriggers *
                    process.softMuons *
                    process.DiMuonPAT *
                    process.DiMuonCounterJPsi *
                    process.DiTrakPAT *
                    process.DiTrakCounterPhi *
                    process.DiMuonDiTrakPAT *
                    process.Triggers *
                    process.DiMuonRootuplerHLT *
                    process.DiTrakRootupler *
                    process.DiMuonDiTrakRootupler
                    )
