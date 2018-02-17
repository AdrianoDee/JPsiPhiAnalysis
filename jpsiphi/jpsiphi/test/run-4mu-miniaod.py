#input_filename = '/store/data/Run2017B/MuOnia/MINIAOD/PromptReco-v1/000/297/723/00000/9040368C-DE5E-E711-ACFF-02163E0134FF.root'
ouput_filename = 'rootuple.root'
input_filename = 'file:FABC2662-9AC8-E711-BF94-02163E019BB9.root'

import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v11', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 20000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_filename))
process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

process.load("mmkk.mmkk.slimmedMuonsTriggerMatcher2017_cfi")
process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")

hltList = [
#Phi
'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi',
#'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi',
# 'HLT_Mu20_TkMu0_Phi',
#'HLT_Dimuon14_Phi_Barrel_Seagulls',
# 'HLT_Mu25_TkMu0_Phi',
# 'HLT_Dimuon24_Phi_noCorrL1',
#JPsi
'HLT_DoubleMu4_JpsiTrkTrk_Displaced',
'HLT_DoubleMu4_JpsiTrk_Displaced',
'HLT_DoubleMu4_Jpsi_Displaced',
'HLT_DoubleMu4_3_Jpsi_Displaced',
'HLT_Dimuon20_Jpsi_Barrel_Seagulls',
'HLT_Dimuon25_Jpsi',
'HLT_Dimuon0_Jpsi3p5_Muon2'
# 'HLT_Dimuon0_Jpsi'
]

hltpaths = cms.vstring(hltList)

hltpathsV = cms.vstring([h + '_v*' for h in hltList])

filters = cms.vstring(
                                #PHI TRIGGERS FILTER

                                #HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi
                                'hltDiMuonGlbOrTrkFiltered0v2', #Phi
                                ##'hltDiMuonGlbOrTrk0zFiltered0p2v2',
                                'hltDoubleMu2JpsiL3Filtered', ##JPsi
                                ##'hltMumuVtxProducerDoubleMu2Jpsi',
                                'hltMumuFilterDoubleMu2Jpsi',

                                #HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi
                                ##'hltDoubleMu2JpsiDoubleTrkL3Filtered',
                                ##'hltDoubleTrkmumuVtxProducerDoubleMu2Jpsi',
                                #'hltDoubleTrkmumuFilterDoubleMu2Jpsi',
                                ##'hltJpsiTkAllConeTracksIterDoubleTrk',
                                ##'hltJpsiTrkTrkVertexProducerPhiDoubleTrk1v2',
                                #'hltJpsiTkTkVertexFilterPhiDoubleTrk1v2',

                                #HLT_Mu20_TkMu0_Phi
                                ## 'hltL3fL1sMu16orMu18erorMu20L1f0L2f0L3Filtered20',
                                ## 'hltDiMuonGlbFiltered20TrkFiltered0',
                                ## 'hltDiMuonGlb20Trk0DzFiltered0p2',

                                #HLT_Dimuon14_Phi_Barrel_Seagulls
                                ##'hltDimuon14PhiBarrelnoCowL3Filtered',
                                ##'hltDisplacedmumuVtxProducerDimuon14PhiBarrelnoCow',
                                ##'hltDisplacedmumuFilterDimuon14PhiBarrelnoCow',

                                #HLT_Mu25_TkMu0_Phi
                                ##'hltL3fL1sMu16orMu18erorMu20L1f0L2f0L3Filtered20',
                                ##'hltDiMuonGlbFiltered25TrkFiltered0',
                                ##'hltDiMuonGlb25Trk0DzFiltered0p2',
                                #HLT_Dimuon24_Phi_noCorrL1
                                ##'hltDisplacedmumuFilterDimuon24PhiBarrelNoCorrL1',
                                ##'hltDisplacedmumuVtxProducerDimuon24PhiNoCorrL1',
                                ##'hltDimuon24PhiNoCorrL1L3fL3Filtered',

                                #JPSI Trigger Filters
                                #HLT_DoubleMu4_JpsiTrkTrk_Displaced_v4
                                ##'hltDoubleMu4JpsiDisplacedL3Filtered'
                                ##'hltDisplacedmumuVtxProducerDoubleMu4Jpsi',
                                'hltDisplacedmumuFilterDoubleMu4Jpsi',
                                ##'hltJpsiTkAllConeTracksIter',
                                ##'hltJpsiTrkTrkVertexProducerPhiKstar',
                                ##'hltJpsiTkTkVertexFilterPhiKstar',

                                #HLT_DoubleMu4_JpsiTrk_Displaced_v12
                                ##'hltDoubleMu4JpsiDisplacedL3Filtered',
                                'hltDisplacedmumuFilterDoubleMu4Jpsi',
                                ##'hltJpsiTkVertexProducer',
                                ##'hltJpsiTkVertexFilter',

                                #HLT_DoubleMu4_Jpsi_Displaced
                                ##'hltDoubleMu4JpsiDisplacedL3Filtered',
                                ##'hltDisplacedmumuVtxProducerDoubleMu4Jpsi',
                                'hltDisplacedmumuFilterDoubleMu4Jpsi',

                                #HLT_DoubleMu4_3_Jpsi_Displaced
                                ##'hltDoubleMu43JpsiDisplacedL3Filtered',
                                ##'hltDisplacedmumuVtxProducerDoubleMu43Jpsi',
                                'hltDisplacedmumuFilterDoubleMu43Jpsi',

                                #HLT_Dimuon20_Jpsi_Barrel_Seagulls
                                #'hltDimuon20JpsiBarrelnoCowL3Filtered',
                                #'hltDisplacedmumuVtxProducerDimuon20JpsiBarrelnoCow',
                                'hltDisplacedmumuFilterDimuon20JpsiBarrelnoCow',

                                #HLT_Dimuon25_Jpsi
                                'hltDisplacedmumuFilterDimuon25Jpsis',

                                #HLT_Dimuon0_Jpsi
                                #'hltDimuon0JpsiL3Filtered',
                                #'hltDisplacedmumuVtxProducerDimuon0Jpsi',
                                #'hltDisplacedmumuFilterDimuon0Jpsi'

                                #HLT_Dimuon0_Jpsi3p5_Muon2
                                'hltJpsiMuonL3Filtered3p5'
                                )


process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring(hltpathsV),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

# process.JPsi2MuMuPAT = cms.EDProducer('Onia2MuMuPAT',
#   muons = cms.InputTag("slimmedMuons"),
#   beamSpotTag = cms.InputTag("offlineBeamSpot"),
#   primaryVertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
#   higherPuritySelection = cms.string(""), ## At least one muon must pass this selection
#   lowerPuritySelection  = cms.string(""), ## BOTH muons must pass this selection
#   dimuonSelection  = cms.string("2.8 < mass < 3.3"), ## The dimuon must pass this selection before vertexing
#   addCommonVertex = cms.bool(True), ## Embed the full reco::Vertex out of the common vertex fit
#   addMuonlessPrimaryVertex = cms.bool(False), ## Embed the primary vertex re-made from all the tracks except the two muons
#   addMCTruth = cms.bool(False),      ## Add the common MC mother of the two muons, if any
#   resolvePileUpAmbiguity = cms.bool(True)   ## Order PVs by their vicinity to the J/psi vertex, not by sumPt
# )
#
# process.Phi2MuMuPAT = cms.EDProducer('Onia2MuMuPAT',
#   muons = cms.InputTag("slimmedMuons"),
#   beamSpotTag = cms.InputTag("offlineBeamSpot"),
#   primaryVertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
#   higherPuritySelection = cms.string(""), ## At least one muon must pass this selection
#   lowerPuritySelection  = cms.string(""), ## BOTH muons must pass this selection
#   dimuonSelection  = cms.string("0.55 < mass < 1.25"), ## The dimuon must pass this selection before vertexing
#   addCommonVertex = cms.bool(True), ## Embed the full reco::Vertex out of the common vertex fit
#   addMuonlessPrimaryVertex = cms.bool(False), ## Embed the primary vertex re-made from all the tracks except the two muons
#   addMCTruth = cms.bool(False),      ## Add the common MC mother of the two muons, if any
#   resolvePileUpAmbiguity = cms.bool(True)   ## Order PVs by their vicinity to the J/psi vertex, not by sumPt
# )


process.Phi2MuMuPAT = cms.EDProducer('FourOnia2MuMuPAT',
        muons                       = cms.InputTag('slimmedMuons'),
        primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        beamSpotTag                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("0.6 < mass && mass < 1.2 && charge==0 "),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(False),
        resolvePileUpAmbiguity      = cms.bool(True),
        HLTFilters                  = filters
)

process.JPsi2MuMuPAT = cms.EDProducer('FourOnia2MuMuPAT',
        muons                       = cms.InputTag('slimmedMuons'),
        primaryVertexTag            = cms.InputTag('offlineSlimmedPrimaryVertices'),
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

process.Onia2MuMuFilteredJpsi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("JPsi2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.9 < mass && mass < 3.3"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters
)

process.Onia2MuMuFilteredPhi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("Phi2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("0.9 < mass && mass < 1.2"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters

)

process.DiMuonCounterJPsi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("Onia2MuMuFilteredJpsi"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

process.DiMuonCounterPhi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("Onia2MuMuFilteredPhi"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

process.PsiPhiProducer = cms.EDProducer('PsiPhiFourMuonsProducer',
    JPsiCollection = cms.InputTag('JPsi2MuMuPAT'),
    PhiCollection = cms.InputTag('Phi2MuMuPAT'),
    JPsiMassCuts = cms.vdouble(2.9,3.2),      # J/psi mass window 3.096916 +/- 0.150
    PhiMassCuts = cms.vdouble(0.9,1.15),  # phi mass window 1.019461 +/- .015
    FourOniaMassCuts = cms.vdouble(4.0,6.0),            # b-hadron mass window
)

process.PsiPhiFitter = cms.EDProducer('PsiPhiFourMuKinematicFit',
    JPsiPhiCollection     = cms.InputTag('PsiPhiProducer','PsiPhiFourMuonsCandidates'),
    PhiConstraint = cms.double(1.019461),              # J/psi mass in GeV
    JPsiConstraint = cms.double(3.096916),
    FourOniaMassCuts = cms.vdouble(4.0,6.0),            # b-hadron mass window
    product_name    = cms.string('PsiPhiCandidatesRefit')
)

process.rootuple = cms.EDAnalyzer('PsiPhiFourMuonsRootupler',
    jpsiphi_cand = cms.InputTag('PsiPhiProducer','PsiPhiFourMuonsCandidates'),
    jpsiphi_rf_cand = cms.InputTag("PsiPhiFitter","PsiPhiCandidatesRefit"),
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    isMC = cms.bool(False),
    OnlyBest = cms.bool(False),
    HLTs = hltpaths,
    filters = filters
)

process.xCandSequence = cms.Sequence(
                   process.triggerSelection *
                   process.JPsi2MuMuPAT *
                   process.Phi2MuMuPAT *
                   process.PsiPhiProducer *
                   process.PsiPhiFitter *
                   process.rootuple
				   )

process.rootupleJPsi = cms.EDAnalyzer('Onia2MuMuRootupler',
                          dimuons = cms.InputTag("JPsi2MuMuPAT"),
                          muons = cms.InputTag("slimmedMuons"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          onia_pdgid = cms.uint32(443),
                          onia_mass_cuts = cms.vdouble(2.8,3.3),
                          isMC = cms.bool(False),
                          OnlyBest = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          HLTs = hltpaths
                          )

process.rootuplePhi = cms.EDAnalyzer('Onia2MuMuRootupler',
                          dimuons = cms.InputTag("Phi2MuMuPAT"),
                          muons = cms.InputTag("slimmedMuons"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          onia_pdgid = cms.uint32(331),
                          onia_mass_cuts = cms.vdouble(0.55,1.25),
                          isMC = cms.bool(False),
                          OnlyBest = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          HLTs = hltpaths
                          )

process.mumuSequence = cms.Sequence(
                    process.triggerSelection *
                    process.JPsi2MuMuPAT *
                    process.Phi2MuMuPAT *
                    process.rootupleJPsi *
                    process.rootuplePhi
)

process.p = cms.Path(process.xCandSequence * process.mumuSequence)
