import FWCore.ParameterSet.Config as cms
process = cms.Process('2mu2k')

from FWCore.ParameterSet.VarParsing import VarParsing

import sys
sys.path.append("mclists/")

par = VarParsing ('analysis')

par.register ('gtag',
                                  "101X_dataRun2_Prompt_v11",
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.string,
                                  "Global Tag")
par.register ('mc',
                                  "y4506_lhcb",
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.string,
                                  "MC Dataset")

par.register ('filein',
                                  "file:1401AF4A-447C-E811-8EEB-FA163E35DF95.root",
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.string,
                                  "Inputfile")

par.register ('dataset',
                                  "Av1",
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.string,
                                  "Dataset")

par.register ('ss',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Do Same Sign")

par.register ('isMC',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Is MC?")

par.register ('isLocal',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Is local?")

par.register ('isDebug',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Debug for you,sir?")

par.register ('kMass',
                                  0.493677,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "KMass")

par.register ('isGen',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Is gen counting only?")

par.register ('i',
                                  0,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.int,
                                  "i")

par.register ('n',
                                  50000,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.int,
                                  "n")

par.register ('isSix',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Also Five Tracks for you, sir?")

par.register ('isFive',
                                  False,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Also Six Tracks for you, sir?")


par.parseArguments()

i = par.i
IsMC = par.isMC
doss = par.ss

gen_file = "file:32B83273-030F-E811-9105-E0071B7AF7C0.root"
input_file = "file:006425F0-6DED-E711-850C-0025904C66E8.root"
mc_file = "file:py8_JPsiMM_EvtGen_13TeV_TuneCP5_cfi.root"
mc_file = "file:02CA3723-CEF3-E711-B1CC-4C79BA1810EF.root"
mc_file = "file:FCD01A2E-A6F5-E711-ACA1-003048F5ADF6.root"
runb2018 = "file:1401AF4A-447C-E811-8EEB-FA163E35DF95.root"
input_file = par.filein #runb2018 #gen_file

filename = "data" + par.dataset

if par.isLocal:

    from bujpsiphi_filelist import *
    from bsjpsiphi_filelist import *
    from bbbar_filelist import *
    from bbbar_soft_filelist import *
    from y4140_lhcb_filelist import *
    from y4140_zero_filelist import *
    from y4273_lhcb_filelist import *
    from y4273_zero_filelist import *
    from y4704_zero_filelist import *
    from y4704_lhcb_filelist import *
    from y4506_lhcb_filelist import *
    from y4506_zero_filelist import *
    from y4273_spin_filelist import *
    from y4506_spin_filelist import *
    from qcd_ml_filelist import *
    from bbbar_ml_filelist import *
    from bstojpsiphi_softqcd_filelist import *
    from bbhook_filelist import *
    from bbhook_filelist_v3 import *
    from BBbar_Hook_Samet import *
    from BBbar_Hook_v5 import *
    from bu_jpsiphi_k import *
    from bu_jpsiphi_k_2 import *
    from bu_jpsiphi_k_3 import *
    from bu_jpsiphi_k_4 import *
    from bs_psiphi import *
    from bs_psiphi_2 import *
    #from y4140_official import *
    y4140_official = 0
    filename = par.mc

    fileLists = {"bbhook_samet": BBbar_Hook_Samet , "qcd_ml" : qcd_ml_filelist,"bbbar_hard" : bbbar_file_list, "bbbar_soft" : bbbar_soft_filelist,
                 "y4273_zero" : y4273_zero_filelist, "y4273_lhcb" : y4273_lhcb_filelist , "BBbar_Hook_v5" :  BBbar_Hook_v5, "bs_psiphi_2" : bs_psiphi_2,
                 "y4140_lhcb" : y4140_lhcb_filelist, "y4140_zero" : y4140_zero_filelist, "bu_jpsiphi_k" : bu_jpsiphi_k,
                 "y4506_lhcb" : y4506_lhcb_filelist, "y4506_zero" : y4506_zero_filelist, "bu_jpsiphi_k_2" : bu_jpsiphi_k_2,
                 "y4704_lhcb" : y4704_lhcb_filelist, "y4704_zero" : y4704_zero_filelist, "bu_jpsiphi_k_3" : bu_jpsiphi_k_3,
                 "y4273_spin" : y4273_spin_filelist, "y4506_spin" : y4506_spin_filelist, "bu_jpsiphi_k_4" : bu_jpsiphi_k_4,
                 "bstojpsiphi_softqcd" : bstojpsiphi_softqcd_file_list, "bsjpsiphi" : bsjpsiphi_filelist, "bs_psiphi" : bs_psiphi,
                 "bujpsiphi" : bujpsiphi_filelist, "bbbar_hook": bbhook_filelist,"bbhook_filelist_v3" : bbhook_filelist_v3,
                 "y4140_official" : y4140_official}

    gtags = {"qcd_ml" : "100X_upgrade2018_realistic_v10", "bbhook_samet" : "100X_upgrade2018_realistic_v10", "bbbar_hard" : "100X_upgrade2018_realistic_v10",
                 "bbbar_soft" : "100X_upgrade2018_realistic_v10", "bbbar_hook" : "100X_upgrade2018_realistic_v10", "BBbar_Hook_v5" : "100X_upgrade2018_realistic_v10",
                 "y4273_zero" : "100X_upgrade2018_realistic_v10",  "y4273_lhcb" : "100X_upgrade2018_realistic_v10"  , "bu_jpsiphi_k" : "100X_upgrade2018_realistic_v10",
                 "y4140_lhcb" : "100X_upgrade2018_realistic_v10",  "y4140_zero" : "100X_upgrade2018_realistic_v10", "bu_jpsiphi_k_2" : "100X_upgrade2018_realistic_v10",
                 "y4506_lhcb" : "100X_upgrade2018_realistic_v10",  "y4506_zero" : "100X_upgrade2018_realistic_v10", "bu_jpsiphi_k_3" : "100X_upgrade2018_realistic_v10",
                 "y4704_lhcb" : "100X_upgrade2018_realistic_v10",  "y4704_zero" : "100X_upgrade2018_realistic_v10", "bu_jpsiphi_k_4" : "100X_upgrade2018_realistic_v10",
                 "y4273_spin" : "100X_upgrade2018_realistic_v10",  "y4506_spin" : "100X_upgrade2018_realistic_v10", "bs_psiphi" : "100X_upgrade2018_realistic_v10",
                 "bujpsiphi" : "100X_upgrade2018_realistic_v10", "bsjpsiphi" : "100X_upgrade2018_realistic_v10", "bs_psiphi_2" : "100X_upgrade2018_realistic_v10",
                 "bstojpsiphi_softqcd" : "94X_mc2017_realistic_v10", "bbhook_filelist_v3" : "100X_upgrade2018_realistic_v10",
                 "y4140_official" : "102X_upgrade2018_realistic_v15"}

    par.gtag = gtags[filename]
    n= par.n

    filelist = fileLists[filename] #bbbar_soft_list#bbbar_file_list
    size = (len(filelist) + n) / n
    input_file = filelist[min(i*size,len(filelist)):min((i+1)*size,len(filelist))]
    print min((i+1)*size,len(filelist))

if par.isGen:
    filename = filename + "_genOnly_"

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v1')
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v2') #F
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v11')
process.GlobalTag = GlobalTag(process.GlobalTag, par.gtag)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(input_file),
    #eventsToProcess = cms.untracked.VEventRange(eventsToProcess),
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('rootuple-2018-dimuondiTrack_five'+ filename  + '_' + str(i) +'.root'),
)

kaonmass = 0.493677
pionmass = 0.13957061

process.load("jpsiphi.jpsiphi.slimmedMuonsTriggerMatcher2017_cfi")
# process.load("jpsiphi.jpsiphi.slimmedTracksTriggerMatcher2017_cfi")

charmoniumHLT = [
#Phi
#'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi', #2017
'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05', #2018
#'HLT_Dimuon25_Jpsi'
#JPsi
'HLT_DoubleMu4_JpsiTrkTrk_Displaced',
'HLT_DoubleMu4_JpsiTrk_Displaced',
#'HLT_DoubleMu4_Jpsi_Displaced',
#'HLT_DoubleMu4_3_Jpsi_Displaced',
#'HLT_Dimuon20_Jpsi_Barrel_Seagulls',

]

year = "2018"

if "2017" in par.dataset:
    year = "2017"
if "2016" in par.dataset:
    year = "2016"

hlts = {}

hlts["2018"] = charmoniumHLT

hlts["2017"] =["HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi",
    "HLT_Dimuon25_Jpsi",
    "HLT_Dimuon20_Jpsi_Barrel_Seagulls",
    "HLT_DoubleMu4_JpsiTrk_Displaced",
    "HLT_DoubleMu4_JpsiTrkTrk_Displaced"]

hlts["2016"]=["HLT_Dimuon20_Jpsi",
    "HLT_Dimuon16_Jpsi",
    "HLT_Dimuon10_Jpsi_Barrel",
    "HLT_DoubleMu4_JpsiTrk_Displaced"]

filts = {}

filts["2018"] = ['hltDoubleMu2JpsiDoubleTrkL3Filtered',
'hltDoubleTrkmumuFilterDoubleMu2Jpsi',
'hltJpsiTkTkVertexFilterPhiDoubleTrk1v1',
'hltJpsiTkTkVertexFilterPhiDoubleTrk1v2',
'hltJpsiTkTkVertexFilterPhiDoubleTrk1v3',
'hltJpsiTkTkVertexFilterPhiDoubleTrk1v4',
'hltJpsiTkTkVertexFilterPhiDoubleTrk1v5',
'hltJpsiTkTkVertexFilterPhiDoubleTrk1v6',

'hltJpsiTrkTrkVertexProducerPhiDoubleTrk1v1',
'hltJpsiTrkTrkVertexProducerPhiDoubleTrk1v2',
'hltJpsiTrkTrkVertexProducerPhiDoubleTrk1v3',
'hltJpsiTrkTrkVertexProducerPhiDoubleTrk1v4',
'hltJpsiTrkTrkVertexProducerPhiDoubleTrk1v5',
'hltJpsiTrkTrkVertexProducerPhiDoubleTrk1v6',

'hltDimuon25JpsiL3fL3Filtered'
]

filts["2017"] =["hltDimuon25JpsiL3fL3Filtered",
"hltDimuon20JpsiBarrelnoCowL3Filtered",
"hltDisplacedmumuFilterDoubleMu4Jpsi",
"hltDoubleMu4JpsiDisplacedL3Filtered",
"hltJpsiTkTkVertexFilterPhiKstar"]

filts["2016"] =["hltDimuon20JpsiL3Filtered",
"hltDimuon16JpsiL3Filtered",
"hltDisplacedmumuFilterDimuon10JpsiBarrel",
"hltDisplacedmumuFilterDoubleMu4Jpsi",
"hltDoubleMu4JpsiDisplacedL3Filtered",
"hltJpsiTkVertexFilter"]


hltpaths = cms.vstring(hlts[year] )
hltpathsV = cms.vstring([h + '_v*' for h in hlts[year] ])

filters = cms.vstring(filts[year])

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring(hltpathsV),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.unpackPatTriggers = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
  patTriggerObjectsStandAlone = cms.InputTag( 'slimmedPatTrigger' ), #selectedPatTrigger for MC
  triggerResults              = cms.InputTag( 'TriggerResults::HLT' ),
  unpackFilterLabels          = cms.bool( True )
)
# process.softMuons = cms.EDFilter('PATMuonSelector',
#    src = cms.InputTag('slimmedMuonsWithTrigger'),
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
        muons                       = cms.InputTag('slimmedMuonsWithTrigger'),
        MuonPtCut                   = cms.double(0.9),
        PrimaryVertex            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        BeamSpot                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("2.95 < mass && mass < 3.25 && charge==0"),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(IsMC),
        resolvePileUpAmbiguity      = cms.bool(True),
        HLTFilters                  = filters,
        TriggerResults              = cms.InputTag("TriggerResults", "", "HLT"),
        TriggerInput                = cms.InputTag("unpackPatTriggers")
)

process.JPsi2MuMuFilter = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("JPsi2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.93 < mass && mass < 3.28 && userFloat('vProb') > 0.001 && pt > 2.0"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters
)

process.muonMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src     = cms.InputTag("slimmedMuonsWithTrigger"), # RECO objects to match
    matched = cms.InputTag("prunedGenParticles"),   # mc-truth particle collection
    mcPdgId     = cms.vint32(13), # one or more PDG ID (13 = muon); absolute values (see below)
    checkCharge = cms.bool(True), # True = require RECO and MC objects to have the same charge
    mcStatus = cms.vint32(1,3,91),     # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR = cms.double(0.5),  # Minimum deltaR for the match
    maxDPtRel = cms.double(0.5),  # Minimum deltaPt/Pt for the match
    resolveAmbiguities = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True), # False = just match input in order; True = pick lowest deltaR pair first
)

process.trackMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src     = cms.InputTag("packedPFCandidates"), # RECO objects to match
    matched = cms.InputTag("prunedGenParticles"),   # mc-truth particle collection
    mcPdgId     = cms.vint32(321,211,13,2212), # one or more PDG ID (13 = muon); absolute values (see below)
    checkCharge = cms.bool(True), # True = require RECO and MC objects to have the same charge
    mcStatus = cms.vint32(1,3,91),     # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR = cms.double(0.5),  # Minimum deltaR for the match
    maxDPtRel = cms.double(0.75),  # Minimum deltaPt/Pt for the match
    resolveAmbiguities = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True), # False = just match input in order; True = pick lowest deltaR pair first
)


process.PsiPhiProducer = cms.EDProducer('DiMuonDiTrackProducer',
    DiMuon              = cms.InputTag('JPsi2MuMuPAT'),
    TrackMatcher        = cms.InputTag("trackMatch"),
    PFCandidates        = cms.InputTag('packedPFCandidates'),
    TrackPtCut            = cms.double(0.7),
    BeamSpot             = cms.InputTag('offlineBeamSpot'),
    PrimaryVertex        = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerInput         = cms.InputTag("unpackPatTriggers"),
    TriggerResults       = cms.InputTag("TriggerResults", "", "HLT"),
    DiMuonMassCuts       = cms.vdouble(2.95,3.25),      # J/psi mass window 3.096916 +/- 0.150
    TrackTrackMassCuts     = cms.vdouble(0.95,1.05),  # phi mass window 1.019461 +/- .015
    DiMuonDiTrackMassCuts = cms.vdouble(4.0,6.0),            # b-hadron mass window
    MassTracks = cms.vdouble(kaonmass,kaonmass),         # Tracks masses
    JPsiMass = cms.double(3.096916),
    PhiMass  = cms.double(1.019461),
    OnlyBest  = cms.bool(False),
    Product = cms.string("DiMuonDiTrackCandidates"),
    Filters = filters,
    IsMC = cms.bool(IsMC),
    AddMCTruth = cms.bool(IsMC),
    AddSS    = cms.bool(doss),
    PionRefit = cms.bool(True)
)

process.FiveTracksProducer  = cms.EDProducer('FiveTracksProducer',

    DiMuoDiTrack             = cms.InputTag('PsiPhiProducer','DiMuonDiTrackCandidates'),
    PFCandidates            = cms.InputTag('packedPFCandidates'),
    TriggerInput            = cms.InputTag("unpackPatTriggers"),
    TrackPtCut               = cms.double(0.7),

    TrackMatcher            = cms.InputTag("trackMatch"),

    BeamSpot                = cms.InputTag('offlineBeamSpot'),
    PrimaryVertex           = cms.InputTag("offlineSlimmedPrimaryVertices"),

    TriggerResults          = cms.InputTag("TriggerResults", "", "HLT"),      # b-hadron mass window
    FiveTrackCuts            = cms.vdouble(4.5,6.5),

    Filters                 = filters,

    IsMC                    = cms.bool(IsMC),
)

process.SixTracksProducer  = cms.EDProducer('SixTracksProducer',

    FiveCollection          = cms.InputTag('FiveTracksProducer','FiveTracks'),
    PFCandidates            = cms.InputTag('packedPFCandidates'),
    TriggerInput            = cms.InputTag("unpackPatTriggers"),
    TrackPtCut               = cms.double(0.7),

    TrackMatcher            = cms.InputTag("trackMatch"),

    BeamSpot                = cms.InputTag('offlineBeamSpot'),
    PrimaryVertex           = cms.InputTag("offlineSlimmedPrimaryVertices"),

    TriggerResults          = cms.InputTag("TriggerResults", "", "HLT"),      # b-hadron mass window
    SixTrackCuts            = cms.vdouble(4.5,6.5),

    AddSS 		    = cms.bool(doss),
    Filters                 = filters,

    IsMC                    = cms.bool(IsMC),
)


# process.PsiPhiFitter = cms.EDProducer('DiMuonDiTrackKinematicFit',
#     DiMuonDiTrack              = cms.InputTag('PsiPhiProducer','DiMuonDiTrackCandidates'),
#     DiMuonMass                = cms.double(3.096916),              # J/psi mass in GeV
#     DiMuonTrackTrackMassCuts    = YMassCut,            # b-hadron mass window
#     MassTracks                 = cms.vdouble(kaonmass,kaonmass),         # Tracks masses
#     Product                   = cms.string('DiMuonDiTrackCandidatesRef')
# )

# process.PsiPhiFitter = cms.EDProducer('DiMuonDiTrackFits',
#     DiMuonDiTrack              = cms.InputTag('PsiPhiProducer','DiMuonDiTrackCandidates'),
#     JPsiMass                  = cms.double(3.096916),
#     PhiMass                   = cms.double(1.019461),              # J/psi mass in GeV
#     DiMuonTrackTrackMassCuts    = YMassCut,            # b-hadron mass window
#     MassTracks                 = cms.vdouble(kaonmass,kaonmass),         # Tracks masses
#     Product                   = cms.string('DiMuonDiTrackCandidatesRef')
# )


# process.rootuple = cms.EDAnalyzer('DiMuonDiTrackRootuplerFit',
#     dimuonditrk_cand = cms.InputTag('PsiPhiProducer','DiMuonDiTrackCandidates'),
#     BeamSpot = cms.InputTag("offlineBeamSpot"),
#     PrimaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
#     TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
#     IsMC = cms.bool(False),
#     OnlyBest = cms.bool(False),
#     OnlyGen = cms.bool(False),
#     Mother_pdg = cms.uint32(20443), #20443 #10441
#     JPsi_pdg = cms.uint32(443),
#     Phi_pdg = cms.uint32(333),
#     HLTs = hltpaths,
#     Filters = filters,
#     TreeName = cms.string('JPsiPhiTree')
# )

process.rootuple = cms.EDAnalyzer('DiMuonDiTrackRootupler',
    DiMuoDiTrack = cms.InputTag('PsiPhiProducer','DiMuonDiTrackCandidates'),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    PrimaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    IsMC = cms.bool(IsMC),
    OnlyGen = cms.bool(False),
    HLTs = hltpaths,
    Filters = filters,
    TreeName = cms.string('JPsiPhiTree')
)

process.rootupleFive = cms.EDAnalyzer('FiveTracksRootupler',
    FiveTracksCand = cms.InputTag('FiveTracksProducer','FiveTracks'),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    PrimaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    IsMC = cms.bool(IsMC),
    OnlyGen = cms.bool(False),
    HLTs = hltpaths,
    Filters = filters,
    TreeName = cms.string('FiveTracksTree')
)

process.rootupleSix = cms.EDAnalyzer('SixTracksRootupler',
    SixTracksCand = cms.InputTag('SixTracksProducer','SixTracks'),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    PrimaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    IsMC = cms.bool(IsMC),
    OnlyGen = cms.bool(False),
    HLTs = hltpaths,
    Filters = filters,
    TreeName = cms.string('SixTracksTree')
)


process.rootupleMuMu = cms.EDAnalyzer('DiMuonRootupler',
      dimuons = cms.InputTag("JPsi2MuMuFilter"),
      muons = cms.InputTag("replaceme"),
      PrimaryVertex = cms.InputTag("offlinePrimaryVertices"),
      TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
      dimuon_pdgid = cms.uint32(443),
      dimuon_mass_cuts = cms.vdouble(2.5,3.5),
      IsMC = cms.bool(IsMC),
      OnlyBest = cms.bool(False),
      OnlyGen = cms.bool(False),
      HLTs = hltpaths
      )

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.genstep = cms.EDAnalyzer('GenMCRootupler',
                      PdgIds          = cms.vint32(20443,10441,531),
                      GoodDaughters   = cms.vint32(443,333),
                      GoodGDaughters  = cms.vint32(13,-13,321,-321),
                      MaxDaughters    = cms.uint32(4),
                      TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
                      PrimaryVertex = cms.InputTag("offlinePrimaryVertices"),
                      HLTs = hltpaths
                       )

genparting = process.genstep
triggering = process.triggerSelection * process.slimmedMuonsWithTriggerSequence * process.unpackPatTriggers
mcmatching = process.trackMatch * process.muonMatch
jpsiing    = process.JPsi2MuMuPAT * process.JPsi2MuMuFilter
tracking   = process.PsiPhiProducer
rootupling = process.rootuple * process.rootupleMuMu

if par.isFive:
    tracking   = process.PsiPhiProducer * process.FiveTracksProducer
    rootupling = process.rootupleFive * process.rootuple * process.rootupleMuMu
if par.isSix:
    tracking   = process.PsiPhiProducer * process.FiveTracksProducer * process.SixTracksProducer
    rootupling = process.rootupleFive * process.rootuple * process.rootupleMuMu * process.rootupleSix

if IsMC:
    allsteps = genparting * triggering * mcmatching * jpsiing * tracking * rootupling
else:
    allsteps = triggering * jpsiing * tracking * rootupling

if par.isGen:
    allsteps = genparting

if par.isDebug:
    allsteps = allsteps * process.dump

process.p = cms.Path(allsteps)
