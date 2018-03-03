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

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('rootuple-2017-ditraktrigger.root'),
)

kaonmass = 0.493677
pionmass = 0.13957061

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


process.JPsi2MuMuPAT = cms.EDProducer('DiMuonPAT',
        muons                       = cms.InputTag('slimmedMuons'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        DiMuonCuts                  = cms.string("2.95 < mass && mass < 3.25 && charge==0"),
)


process.DiMuonCounterJPsi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("DiMuonPAT"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

process.Phi2KKPAT = cms.EDProducer('DiTrakPAT',
        muons                       = cms.InputTag('slimmedMuons'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        DiMuonCuts                  = cms.string("2.95 < mass && mass < 3.25 && charge==0"),
)


process.DiMuonCounterJPsi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("DiMuonPAT"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

process.DiTrakHLT  = cms.EDAnalyzer('DiTrakHLT',
    PFCandidates        = cms.InputTag("packedPFCandidates"),
    TriggerInput        = cms.InputTag("unpackPatTriggers"),
    primaryVertices     = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults      = cms.InputTag("TriggerResults", "", "HLT"),
    TrakTrakMassCuts    = cms.vdouble(0.5,1.5),
    MassTraks           = cms.vdouble(kaonmass,kaonmass),
    isMC                = cms.bool(False),
    OnlyBest            = cms.bool(False),
    HLTs                = hltpaths,
    Filters             = filters
)

process.DiMuonDiTrakProducerHLT = cms.EDProducer('DiMuonDiTrakProducerHLT',
    DiMuon = cms.InputTag('JPsi2MuMuPAT'),
    PFCandidates = cms.InputTag('packedPFCandidates'),
    TriggerInput        = cms.InputTag("unpackPatTriggers"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    DiMuonMassCuts = cms.vdouble(2.95,3.25),      # J/psi mass window 3.096916 +/- 0.150
    TrakTrakMassCuts = cms.vdouble(1.0,1.04),  # phi mass window 1.019461 +/- .015
    DiMuonDiTrakMassCuts = cms.vdouble(4.0,5.8),            # b-hadron mass window
    MassTraks = cms.vdouble(kaonmass,kaonmass),         # traks masses
    MaxDeltaRPt = cms.vdouble(0.01,2.0),
    OnlyBest  = cms.bool(False),
    Product = cms.string("DiMuonDiTrakCandidatesHLT"),
    HLTs                = hltpaths,
    Filters             = filters
)

process.DiMuonDiTrakKinematicFit = cms.EDProducer('DiMuonDiTrakKinematicFit',
    DiMuonDiTrak        = cms.InputTag('DiMuonDiTrakProducerHLT','DiMuonDiTrakCandidatesHLT'),
    DiMuonMass          = cms.double(3.096916),              # J/psi mass in GeV
    DiMuonTrakTrakMassCuts    = cms.vdouble(4.0,5.8),            # b-hadron mass window
    MassTraks           = cms.vdouble(kaonmass,kaonmass),         # traks masses
    Product             = cms.string('DiMuonDiTrakCandidatesRef')
)

process.DiMuonDiTrakRootuplerHLT = cms.EDAnalyzer('DiMuonDiTrakRootuplerHLT',
    dimuonditrk_cand = cms.InputTag('DiMuonDiTrakProducerHLT','DiMuonDiTrakCandidatesHLT'),
    dimuonditrk_rf_cand = cms.InputTag("DiMuonDiTrakKinematicFit","DiMuonDiTrakCandidatesRef"),
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    isMC = cms.bool(False),
    OnlyBest = cms.bool(False),
    HLTs = hltpaths,
    Filters = filters,
    TreeName = cms.string('JPsi Phi Tree')
)

process.DiMuonRootuplerHLT = cms.EDAnalyzer('DiMuonRootuplerHLT',
    dimuons = cms.InputTag("JPsi2MuMuFilter"),
    muons = cms.InputTag("slimmedMuonsWithTrigger"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    dimuon_pdgid = cms.uint32(443),
    dimuon_mass_cuts = cms.vdouble(2.8,3.3),
    isMC = cms.bool(False),
    OnlyBest = cms.bool(False),
    OnlyGen = cms.bool(False),
    HLTs = hltpaths
 )

process.p = cms.Path(process.triggerSelection *
                     process.slimmedMuonsWithTriggerSequence *
                     process.unpackPatTriggers *
                     # process.slimmedPFCandsWithTriggerSequence *
                     process.softMuons *
                     process.JPsi2MuMuPAT *
                     process.JPsi2MuMuFilter*
                     process.DiMuonCounterJPsi*
                     process.DiMuonDiTrakProducerHLT *
                     process.DiMuonDiTrakKinematicFit *
                     process.DiMuonDiTrakRootuplerHLT *
                     process.DiTrakHLT*
                     process.DiMuonRootuplerHLT)
                     # process.rootupleMuMu)# * process.Phi2KKPAT * process.patSelectedTracks *process.rootupleKK)
