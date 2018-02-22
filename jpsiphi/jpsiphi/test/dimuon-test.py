#input_filename = '/store/data/Run2016B/MuOnia/MINIAOD/PromptReco-v1/000/297/723/00000/9040368C-DE5E-E711-ACFF-02163E0134FF.root'
ouput_filename = 'rootuple-2012-doubledimuon.root'
input_filename = 'file::004C60F9-188F-E211-B43F-001E6849D384.root'

import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration/EventContent/EventContent_cff')
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskAlgoTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
#process.GlobalTag.globaltag = 'FT_R_53_V18::All' #ABC
process.GlobalTag.globaltag = 'FT_R_53_V21::All'
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v16')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 20000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_filename))
process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

process.load("jpsiphi.jpsiphi.patMuonsTriggerMatcher2012_cfi")

process.isMC = cms.bool(False)

charmoniumHLT = [
#Chamonium
'HLT_DoubleMu4_Jpsi_Displaced',
'HLT_Dimuon8_Jpsi_v4'
]

hltList = charmoniumHLT

hltpaths = cms.vstring(hltList)

hltpathsV = cms.vstring([h + '_v*' for h in hltList])

filters = cms.vstring(
                                "hltDisplacedmumuFilterDoubleMu4Jpsi",
                                "hltVertexmumuFilterDimuon8Jpsi",
                                )


process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring(hltpathsV),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.JPsi2MuMuPAT = cms.EDProducer('DiMuonProducerPAT',
        muons                       = cms.InputTag('patMuonsWithTrigger'),
        primaryVertexTag            = cms.InputTag('offlinePrimaryVertices'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("2.9 < mass && mass < 3.3 && charge==0 "),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(False),
        resolvePileUpAmbiguity      = cms.bool(True),
        HLTFilters                  = filters
)

process.rootupleJPsi = cms.EDAnalyzer('DiMuonRootupler',
                          dimuons = cms.InputTag("JPsi2MuMuPAT"),
                          muons = cms.InputTag("patMuonsWithTrigger"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          onia_pdgid = cms.uint32(443),
                          onia_mass_cuts = cms.vdouble(2.8,3.3),
                          isMC = cms.bool(False),
                          OnlyBest = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          HLTs = hltpaths
                          )

process.sequence = cms.Sequence(
                    process.triggerSelection *
                    process.JPsi2MuMuPAT *
                    process.rootupleJPsi
)

process.p = cms.Path(process.sequence)
