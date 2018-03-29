import FWCore.ParameterSet.Config as cms
process = cms.Process('PSIKK')

input_file = "file:049F2D32-26F2-E711-A162-00266CFFC664.root"
input_file = "root://cms-xrd-global.cern.ch//store/data/Run2017F/Charmonium/RECO/PromptReco-v1/000/305/040/00000/0A96ADE7-E8B1-E711-8730-02163E01A45C.root"

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
        fileName = cms.string('rootuple-2017-ml.root'),
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

process.ditrakdimuon = cms.EDAnalyzer("DiMuonDiTrakMLAnalyzer",
                                        Muons           = cms.InputTag( "muons" ),
                                        Tracks          = cms.InputTag( "generalTracks" ),
                                        BeamSpot        = cms.InputTag("offlineBeamSpot"),
                                        PrimaryVertex   = cms.InputTag("offlinePrimaryVertices"),
                                        DiMuonCuts      = cms.vdouble(2.9,3.3),
                                        DiTrakCuts      = cms.vdouble(1.0,1.05),
                                        TriggerResults     = cms.InputTag("TriggerResults", "", "HLT"),
                                        DiMuonMass = cms.double(3.0916)
                                        )



process.p = cms.Path(
                     process.triggerSelection *
                     process.ditrakdimuon
                    )
