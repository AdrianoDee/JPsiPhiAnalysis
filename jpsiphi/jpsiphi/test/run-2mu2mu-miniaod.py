import FWCore.ParameterSet.Config as cms
process = cms.Process('2mu2mu')

from FWCore.ParameterSet.VarParsing import VarParsing
from Y_MC_Files import *

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
                                  True,
                                  VarParsing.multiplicity.singleton,
                                  VarParsing.varType.bool,
                                  "Do Same Sign")

par.register ('IsMC',
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


par.parseArguments()

i = par.i
IsMC = par.IsMC
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

    from bbbar_filelist import *
    from bbbar_soft_list import *
    from y4140_lhcb_filelist import *
    from y4140_zero_filelist import *
    from y4273_lhcb_filelist import *
    from y4273_zero_filelist import *
    from y4704_zero_filelist import *
    from y4704_lhcb_filelist import *
    from y4506_lhcb_filelist import *
    from y4506_zero_filelist import *
    from qcd_ml_filelist import *
    from bbbar_ml_filelist import *
    from y4140_official import *
    filename = par.mc

    fileLists = {"qcd_ml" : qcd_ml_filelist,"bbbar_hard" : bbbar_file_list,
                 "y4273_zero" : y4273_zero_filelist, "y4273_lhcb" : y4273_lhcb_filelist ,
                 "y4140_lhcb" : y4140_lhcb_filelist, "y4140_zero" : y4140_zero_filelist,
                 "y4506_lhcb" : y4506_lhcb_filelist, "y4506_zero" : y4506_zero_filelist,
                 "y4704_lhcb" : y4704_lhcb_filelist, "y4704_zero" : y4704_zero_filelist,
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

    n= par.n
    par.gtag = gtags[filename]

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
    fileNames = cms.untracked.vstring(input_file)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('rootuple-2018-fourmuons_'+ filename  + '_' + str(i) +'.root'),
)

kaonmass = 0.493677
pionmass = 0.13957061

process.load("jpsiphi.jpsiphi.slimmedMuonsTriggerMatcher2017_cfi")
# process.load("jpsiphi.jpsiphi.slimmedTracksTriggerMatcher2017_cfi")

year = "2018"

if "2017" in par.dataset:
    year = "2017"
if "2016" in par.dataset:
    year = "2016"

charmoniumHLT = ['HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi',
#'HLT_DoubleMu4_Jpsi_Displaced',
#'HLT_DoubleMu4_3_Jpsi_Displaced',
#'HLT_Dimuon20_Jpsi_Barrel_Seagulls',
#'HLT_Dimuon25_Jpsi',
#'HLT_Dimuon0_Jpsi3p5_Muon2'
]

hlts = {}

hlts["2018"] = charmoniumHLT

hlts["2017"] =["HLT_Dimuon0_Jpsi3p5_Muon2"]

hlts["2016"]=["HLT_Dimuon0_Jpsi_Muon"]

filts = {}

filts["2018"] = ['hltDoubleMu2JpsiDoubleTrkL3Filtered',
'hltDiMuonGlbOrTrk0zFiltered0p2v1',
'hltDiMuonGlbOrTrk0zFiltered0p2v2',
'hltDiMuonGlbOrTrk0zFiltered0p2v3',
'hltDiMuonGlbOrTrk0zFiltered0p2v4',
'hltDiMuonGlbOrTrk0zFiltered0p2v5',
'hltDiMuonGlbOrTrk0zFiltered0p2v6',

'hltDiMuonGlbOrTrkFiltered0v1',
'hltDiMuonGlbOrTrkFiltered0v2',
'hltDiMuonGlbOrTrkFiltered0v3',
'hltDiMuonGlbOrTrkFiltered0v4',
'hltDiMuonGlbOrTrkFiltered0v5',
'hltDiMuonGlbOrTrkFiltered0v6'

'hltJpsiMuonL3Filtered3p5',
'hltVertexmumuFilterJpsiMuon3p5'
]

filts["2016"] = ["hltJpsiMuonL3Filtered","hltVertexmumuFilterJpsiMuon"]

filts["2017"] = ["hltJpsiMuonL3Filtered3p5","hltVertexmumuFilterJpsiMuon3p5"]

hltpaths = cms.vstring(hlts[year])
hltpathsV = cms.vstring([h + '_v*' for h in hlts[year]])

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
process.Phi2MuMuPAT = cms.EDProducer('DiMuonProducerPAT',
        muons                       = cms.InputTag('slimmedMuonsWithTrigger'),
		MuonPtCut                   = cms.double(0.0),
		TriggerInput                = cms.InputTag("unpackPatTriggers"),
		TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
        PrimaryVertex            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        BeamSpot                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("0.6 < mass && mass < 1.2 && charge==0 "),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(IsMC),
        resolvePileUpAmbiguity      = cms.bool(True),
        HLTFilters                  = filters
)

process.JPsi2MuMuPAT = cms.EDProducer('DiMuonProducerPAT',
        muons                       = cms.InputTag('slimmedMuonsWithTrigger'),
		MuonPtCut                   = cms.double(0.5),
		TriggerInput                = cms.InputTag("unpackPatTriggers"),
		TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
        PrimaryVertex            = cms.InputTag('offlineSlimmedPrimaryVertices'),
        BeamSpot                 = cms.InputTag('offlineBeamSpot'),
        higherPuritySelection       = cms.string(""),
        lowerPuritySelection        = cms.string(""),
        dimuonSelection             = cms.string("2.9 < mass && mass < 3.3 && charge==0 "),
        addCommonVertex             = cms.bool(True),
        addMuonlessPrimaryVertex    = cms.bool(False),
        addMCTruth                  = cms.bool(IsMC),
        resolvePileUpAmbiguity      = cms.bool(True),
        HLTFilters                  = filters
)

process.DiMuonFilteredJpsi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("JPsi2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.9 < mass && mass < 3.3 && userFloat('vProb') > 0.005"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters
)

process.DiMuonFilteredPhi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("Phi2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("0.6 < mass && mass < 1.2 && userFloat('vProb') > 0.005 "),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters

)

process.DiMuonCounterJPsi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("DiMuonFilteredJpsi"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

process.DiMuonCounterPhi = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("DiMuonFilteredPhi"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

# process.PsiPhiProducer = cms.EDProducer('DiMuonDiTrakProducer',
#     DiMuon = cms.InputTag('JPsi2MuMuPAT'),
#     PFCandidates = cms.InputTag('packedPFCandidates'),
#     TriggerInput            = cms.InputTag("unpackPatTriggers"),
#     TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
#     DiMuonMassCuts = cms.vdouble(2.95,3.25),      # J/psi mass window 3.096916 +/- 0.150
#     TrakTrakMassCuts = cms.vdouble(1.0,1.04),  # phi mass window 1.019461 +/- .015
#     DiMuonDiTrakMassCuts = cms.vdouble(4.0,5.8),            # b-hadron mass window
#     MassTraks = cms.vdouble(kaonmass,kaonmass),         # traks masses
#     OnlyBest  = cms.bool(False),
#     Product = cms.string("DiMuonDiTrakCandidates"),
#     Filters = filters,
#     IsMC = cms.bool(True),
# )

# process.PsiPhiFitter = cms.EDProducer('DiMuonDiTrakKinematicFit',
#     DiMuonDiTrak              = cms.InputTag('PsiPhiProducer','DiMuonDiTrakCandidates'),
#     DiMuonMass                = cms.double(3.096916),              # J/psi mass in GeV
#     DiMuonTrakTrakMassCuts    = cms.vdouble(4.1,5.5),            # b-hadron mass window
#     MassTraks                 = cms.vdouble(kaonmass,kaonmass),         # traks masses
#     Product                   = cms.string('DiMuonDiTrakCandidatesRef')
# )

# process.rootuple = cms.EDAnalyzer('DiMuonDiTrakRootupler',
#     dimuonditrk_cand = cms.InputTag('PsiPhiProducer','DiMuonDiTrakCandidates'),
#     dimuonditrk_rf_cand = cms.InputTag("PsiPhiFitter","DiMuonDiTrakCandidatesRef"),
#     BeamSpot = cms.InputTag("offlineBeamSpot"),
#     PrimaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
#     TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
#     IsMC = cms.bool(True),
#     OnlyBest = cms.bool(False),
#     OnlyGen = cms.bool(False),
#     Mother_pdg = cms.uint32(10441), #20443 #10441
#     JPsi_pdg = cms.uint32(443),
#     Phi_pdg = cms.uint32(333),
#     HLTs = hltpaths,
#     Filters = filters,
#     TreeName = cms.string('JPsi Phi Tree')
# )


process.PsiPhiProducer = cms.EDProducer('DoubleDiMuonProducer',
    HighDiMuonCollection  = cms.InputTag('JPsi2MuMuPAT'),
    LowDiMuonCollection  = cms.InputTag('Phi2MuMuPAT'),
	BeamSpot                 = cms.InputTag('offlineBeamSpot'),
	PrimaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerInput            = cms.InputTag("unpackPatTriggers"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    HighDiMuonMassCuts	= cms.vdouble(2.95,3.25),      # J/psi mass window 3.096916 +/- 0.150
    LowDiMuonMassCuts = cms.vdouble(0.6,1.2),  # phi mass window 1.019461 +/- .015
    DoubleDiMuonMassCuts = cms.vdouble(4.0,6.0),            # b-hadron mass window
	JPsiMass = cms.double(3.096916),
	PhiMass  = cms.double(1.019461),
    OnlyBest  = cms.bool(False),
    Product = cms.string("FourMuonCandidates"),
    Filters = filters,
    IsMC = cms.bool(IsMC),
    AddMCTruth = cms.bool(IsMC),
	DoDouble = cms.bool(False),
	AddSS    = cms.bool(doss),
)

process.rootuple = cms.EDAnalyzer('DoubleDiMuonRootupler',
    FourMuons = cms.InputTag('PsiPhiProducer','FourMuonCandidates'),
	FiveTrakPos = cms.InputTag('FiveTracksProducer','FiveTracksPos'),
	FiveTrakNeg = cms.InputTag('FiveTracksProducer','FiveTracksNeg'),
	FiveTrakNeu = cms.InputTag('FiveTracksProducer','FiveTracksNeu'),
    PrimaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
	AddMC = cms.bool(IsMC),
    IsMC = cms.bool(IsMC),
    OnlyBest = cms.bool(False),
    OnlyGen = cms.bool(False),
    Mother_pdg = cms.uint32(20443), #20443 #20443
    JPsi_pdg = cms.uint32(443),
    Phi_pdg = cms.uint32(333),
    HLTs = hltpaths,
    Filters = filters,
    TreeName = cms.string('JPsiPhiTree')
)

process.genmc = cms.EDAnalyzer('GenMCRootupler',
    PdgIds = cms.vint32(20443,20443,531),
)

process.rootupleJPsi = cms.EDAnalyzer('DiMuonRootupler',
                          dimuons = cms.InputTag("JPsi2MuMuPAT"),
                          muons = cms.InputTag("replaceme"),
                          PrimaryVertex = cms.InputTag("offlinePrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          dimuon_pdgid = cms.uint32(443),
                          dimuon_mass_cuts = cms.vdouble(2.5,3.5),
                          IsMC = cms.bool(False),
                          OnlyBest = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          HLTs = hltpaths
                          )

process.rootuplePhi = cms.EDAnalyzer('DiMuonRootupler',
                          dimuons = cms.InputTag("Phi2MuMuPAT"),
                          muons = cms.InputTag("replaceme"),
                          PrimaryVertex = cms.InputTag("offlinePrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          dimuon_pdgid = cms.uint32(333),
                          dimuon_mass_cuts = cms.vdouble(0.6,1.2),
                          IsMC = cms.bool(False),
                          OnlyBest = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          HLTs = hltpaths
                          )
process.genstep = cms.EDAnalyzer('GenMCRootupler',
                      PdgIds          = cms.vint32(20443,10441,531),
                      GoodDaughters   = cms.vint32(443,333),
                      GoodGDaughters  = cms.vint32(13,-13,13,-13),
                      MaxDaughters    = cms.uint32(4),
                      TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
                      PrimaryVertex = cms.InputTag("offlinePrimaryVertices"),
                      HLTs = hltpaths
                       )

process.dump=cms.EDAnalyzer('EventContentAnalyzer')


genparting = process.genstep
triggering = process.triggerSelection * process.slimmedMuonsWithTriggerSequence * process.unpackPatTriggers
##mcmatching = process.trackMatch * process.muonMatch
dimunoing  = process.JPsi2MuMuPAT * process.Phi2MuMuPAT
tracking   = process.PsiPhiProducer
rootupling = process.rootuple * process.rootuplePhi * process.rootupleJPsi
debugging  = process.dump

if IsMC:
    allsteps = genparting * triggering * dimunoing * tracking * rootupling
else:
    allsteps = triggering * dimunoing * tracking * rootupling

if par.isGen:
    allsteps = genparting

if par.isDebug:
    allsteps = allsteps * process.dump

process.p = cms.Path(allsteps)
