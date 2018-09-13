import FWCore.ParameterSet.Config as cms
process = cms.Process('PSIKK')

gen_file = "file:32B83273-030F-E811-9105-E0071B7AF7C0.root"
input_file = "file:006425F0-6DED-E711-850C-0025904C66E8.root"
mc_file = "file:py8_JPsiMM_EvtGen_13TeV_TuneCP5_cfi.root"

input_file = mc_file #gen_file

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
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
        fileName = cms.string('rootuple-2017-dimuonditrak_gen.root'),
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
                                        hltResults = cms.InputTag( "TriggerResults", "", "GEN" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )


process.rootupleMuMu = cms.EDAnalyzer('DiMuonRootuplerGEN',
                          dimuons = cms.InputTag("nothing"),
                          muons = cms.InputTag("replaceme"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          dimuon_pdgid = cms.uint32(443),
                          mother_pdgid = cms.uint32(531),
                          dimuon_mass_cuts = cms.vdouble(2.5,3.5),
                          isMC = cms.bool(True),
                          OnlyBest = cms.bool(True),
                          OnlyGen = cms.bool(True),
                          HLTs = hltpaths
                          )

process.p = cms.Path(process.rootupleMuMu)# * process.Phi2KKPAT * process.patSelectedTracks *process.rootupleKK)
