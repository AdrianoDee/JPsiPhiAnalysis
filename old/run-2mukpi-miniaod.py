import FWCore.ParameterSet.Config as cms
process = cms.Process('PSIKK')

input_file = "file:006425F0-6DED-E711-850C-0025904C66E8.root"

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
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters
)

process.PsiPiKProducer = cms.EDProducer('DiMuonDiTrakProducer',
    DiMuon = cms.InputTag('JPsi2MuMuPAT'),
    PFCandidates = cms.InputTag('packedPFCandidates'),
    DiMuonMassCuts = cms.vdouble(2.95,3.25),      # J/psi mass window 3.096916 +/- 0.150
    TrakTrakMassCuts = cms.vdouble(1.0,1.04),  # phi mass window 1.019461 +/- .015
    DiMuonDiTrakMassCuts = cms.vdouble(4.0,5.8),            # b-hadron mass window
    MassTraks = cms.vdouble(pionmass,kaonmass),         # traks masses
    OnlyBest  = cms.bool(False),
    Product   = cms.string('DiMuonPiKCandidates')
)

process.PsiKPiProducer = cms.EDProducer('DiMuonDiTrakProducer',
    DiMuon = cms.InputTag('JPsi2MuMuPAT'),
    PFCandidates = cms.InputTag('packedPFCandidates'),
    DiMuonMassCuts = cms.vdouble(2.95,3.25),      # J/psi mass window 3.096916 +/- 0.150
    TrakTrakMassCuts = cms.vdouble(1.0,1.04),  # phi mass window 1.019461 +/- .015
    DiMuonDiTrakMassCuts = cms.vdouble(4.0,5.8),            # b-hadron mass window
    MassTraks = cms.vdouble(kaonmass,pionmass),         # traks masses
    OnlyBest  = cms.bool(False),
    Product   = cms.string('DiMuonKPiCandidates')
)

process.PsiPiKFitter = cms.EDProducer('DiMuonDiTrakKinematicFit',
    DiMuonDiTrak        = cms.InputTag('PsiPiKProducer','DiMuonPiKCandidates'),
    DiMuonMass          = cms.double(3.096916),              # J/psi mass in GeV
    DiMuonTrakTrakMassCuts    = cms.vdouble(4.0,5.8),            # b-hadron mass window
    MassTraks           = cms.vdouble(pionmass,kaonmass),         # traks masses
    Product             = cms.string('DiMuonPiKCandidatesRef')
)

process.PsiKPiFitter = cms.EDProducer('DiMuonDiTrakKinematicFit',
    DiMuonDiTrak        = cms.InputTag('PsiKPiProducer','DiMuonKPiCandidates'),
    DiMuonMass          = cms.double(3.096916),              # J/psi mass in GeV
    DiMuonTrakTrakMassCuts    = cms.vdouble(4.0,5.8),            # b-hadron mass window
    MassTraks           = cms.vdouble(kaonmass,pionmass),         # traks masses
    Product             = cms.string('DiMuonKPiCandidatesRef')
)


process.rootupleKPi = cms.EDAnalyzer('DiMuonDiTrakRootupler',
    dimuonditrk_cand = cms.InputTag('PsiKPiProducer','DiMuonKPiCandidates'),
    dimuonditrk_rf_cand = cms.InputTag("PsiKPiFitter","DiMuonKPiCandidatesRef"),
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    isMC = cms.bool(False),
    OnlyBest = cms.bool(False),
    HLTs = hltpaths,
    Filters = filters,
    TreeName = cms.string('DiMuon KPi Tree')
)

process.rootuplePiK = cms.EDAnalyzer('DiMuonDiTrakRootupler',
    dimuonditrk_cand = cms.InputTag('PsiPiKProducer','DiMuonPiKCandidates'),
    dimuonditrk_rf_cand = cms.InputTag("PsiPiKFitter","DiMuonPiKCandidatesRef"),
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    isMC = cms.bool(False),
    OnlyBest = cms.bool(False),
    HLTs = hltpaths,
    Filters = filters,
    TreeName = cms.string('DiMuon PiK Tree')
)

process.p = cms.Path(process.triggerSelection *
                     process.slimmedMuonsWithTriggerSequence *
                     # process.slimmedPFCandsWithTriggerSequence *
                     process.softMuons *
                     process.JPsi2MuMuPAT *
                     process.JPsi2MuMuFilter*
                     process.PsiKPiProducer *
                     process.PsiPiKProducer *
                     process.PsiKPiFitter *
                     process.PsiPiKFitter *
                     process.rootuplePiK *
                     process.rootupleKPi)# * process.Phi2KKPAT * process.patSelectedTracks *process.rootupleKK)
