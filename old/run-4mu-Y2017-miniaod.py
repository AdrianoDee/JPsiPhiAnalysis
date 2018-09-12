import FWCore.ParameterSet.Config as cms
process = cms.Process("4MU")
from FWCore.ParameterSet.VarParsing import VarParsing
from Y_MC_Files import *

options = VarParsing ('analysis')

options.register ('yMass',
				  4700,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "MC sample Y mass ")

options.register ('trigger',
				  True,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.bool,
				  "Adding triggers")

options.parseArguments()

#input_filename = '/store/data/Run2016B/MuOnia/MINIAOD/PromptReco-v1/000/297/723/00000/9040368C-DE5E-E711-ACFF-02163E0134FF.root'

input_file = 'file:00000113-58FF-E711-AB19-002590E7E02E.root'#Muonia
input_file = "file:006425F0-6DED-E711-850C-0025904C66E8.root" #Charmonium

Y = str(options.yMass)

input_files = {"4700" : input_files_4700,"4500" : input_files_4500,"4300" : input_files_4300,"4100" : input_files_4100}
YMassCuts   = {"4700" : cms.vdouble(4.0,5.4),"4500" : cms.vdouble(3.9,5.1),"4300" : cms.vdouble(3.7,4.9),"4100" : cms.vdouble(3.5,4.7)}
motherPdgs  = {"4700" : 10441, "B" : 532 , "4100" : 20443, "4300": 20443, "4500": 10441}

motherPdg = motherPdgs[Y]
input_file = input_files[Y]
YMassCut = YMassCuts[Y]

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v1')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v11')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 20000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_file))

filename = "rootuple-2017_MC_Y_" + str(Y) + "_doubledimuon.root"

process.TFileService = cms.Service("TFileService",fileName = cms.string(filename))
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False))

process.load("jpsiphi.jpsiphi.slimmedMuonsTriggerMatcher2017_cfi")


charmoniumHLT = [
'HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi',
#'HLT_DoubleMu4_Jpsi_Displaced',
#'HLT_DoubleMu4_3_Jpsi_Displaced',
#'HLT_Dimuon20_Jpsi_Barrel_Seagulls',
#'HLT_Dimuon25_Jpsi',
#'HLT_Dimuon0_Jpsi3p5_Muon2'
]

muoniaHLT = [
#'HLT_Dimuon14_Phi_Barrel_Seagulls',
]

hltList = charmoniumHLT

hltpaths = cms.vstring(hltList)

hltpathsV = cms.vstring([h + '_v*' for h in hltList])

charmoniumFilters = cms.vstring(
                                "hltDoubleMu2JpsiDoubleTrkL3Filtered",
                                "hltDoubleTrkmumuFilterDoubleMu2Jpsi",
                                "hltJpsiTkTkVertexFilterPhiDoubleTrk1v2",

                                'hltDisplacedmumuFilterDoubleMu4Jpsi',

                                'hltDisplacedmumuFilterDoubleMu43Jpsi',

                                'hltDisplacedmumuFilterDimuon20JpsiBarrelnoCow',

                                'hltDisplacedmumuFilterDimuon25Jpsis',

                                'hltJpsiMuonL3Filtered3p5',
                                )

muoniaFilters = cms.vstring("hltDisplacedmumuFilterDimuon14PhiBarrelnoCow")

filters = charmoniumFilters

process.yfilter = cms.EDFilter("GenFilter",
								PdgId = cms.uint32(motherPdg)
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

process.Phi2MuMuPAT = cms.EDProducer('DiMuonProducerPAT',
        muons                       = cms.InputTag('slimmedMuonsWithTrigger'),
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

process.JPsi2MuMuPAT = cms.EDProducer('DiMuonProducerPAT',
        muons                       = cms.InputTag('slimmedMuonsWithTrigger'),
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

process.DiMuonFilteredJpsi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("JPsi2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.95 < mass && mass < 3.25 && userFloat('vProb') > 0.01"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = filters
)

process.DiMuonFilteredPhi = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("Phi2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("0.6 < mass && mass < 1.12 && userFloat('vProb') > 0.01 "),
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

process.DoubleDiMuonProducer = cms.EDProducer('DoubleDiMuonProducer',
    HighDiMuonCollection    = cms.InputTag('JPsi2MuMuPAT'),
    LowDiMuonCollection     = cms.InputTag('Phi2MuMuPAT'),
    HighDiMuonMassCuts      = cms.vdouble(2.9,3.2),      # J/psi mass window 3.096916 +/- 0.150
    LowDiMuonMassCuts       = cms.vdouble(0.9,1.15),  # phi mass window 1.019461 +/- .015
    DoubleDiMuonMassCuts    = cms.vdouble(4.0,6.0),            # b-hadron mass window
	AddMCTruth				= cms.bool(True)
)

process.DoubleDiMuonFitter = cms.EDProducer('DoubleDiMuonKinematicFit',
    DoubleDiMuonCollection  = cms.InputTag('DoubleDiMuonProducer','DoubleDiMuonCandidates'),
    LowMassConstraint       = cms.double(1.019461),              # J/psi mass in GeV
    HighMassConstraint      = cms.double(3.096916),
    DoubleDiMuonMassCuts    = YMassCut,            # b-hadron mass window
    product_name            = cms.string('DoubleDiMuonCandidatesRefit')
)

process.rootuplefourmu = cms.EDAnalyzer('DoubleDiMuonRootupler',
    doubledimuon_cand       = cms.InputTag('DoubleDiMuonProducer','DoubleDiMuonCandidates'),
    doubledimuon_rf_cand    = cms.InputTag("DoubleDiMuonFitter","DoubleDiMuonCandidatesRefit"),
    beamSpotTag             = cms.InputTag("offlineBeamSpot"),
    primaryVertices         = cms.InputTag("offlineSlimmedPrimaryVertices"),
    TriggerResults          = cms.InputTag("TriggerResults", "", "HLT"),
    Mother_pdg              = cms.uint32(motherPdg),
    isMC                    = cms.bool(True),
    OnlyBest                = cms.bool(False),
    OnlyGen                 = cms.bool(False),
    HLTs                    = hltpaths,
    filters                 = filters
)


process.rootupleJPsi = cms.EDAnalyzer('DiMuonRootupler',
                          dimuons = cms.InputTag("JPsi2MuMuPAT"),
                          muons = cms.InputTag("softMuons"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          dimuon_pdgid = cms.uint32(443),
                          dimuon_mass_cuts = cms.vdouble(2.8,3.3),
                          isMC = cms.bool(False),
                          OnlyBest = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          HLTs = hltpaths
                          )

process.rootuplePhi = cms.EDAnalyzer('DiMuonRootupler',
                          dimuons = cms.InputTag("Phi2MuMuPAT"),
                          muons = cms.InputTag("softMuons"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          TriggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                          dimuon_pdgid = cms.uint32(331),
                          dimuon_mass_cuts = cms.vdouble(0.55,1.25),
                          isMC = cms.bool(False),
                          OnlyBest = cms.bool(False),
                          OnlyGen = cms.bool(False),
                          HLTs = hltpaths
                          )

if options.trigger:
	process.sequence = cms.Sequence(
						process.yfilter *
	                    process.triggerSelection *
	                    process.slimmedMuonsWithTriggerSequence *
	                    #process.softMuons *
	                    process.JPsi2MuMuPAT *
	                    process.Phi2MuMuPAT *
	                    process.DiMuonFilteredJpsi *
	                    process.DiMuonFilteredPhi *
	                    process.DoubleDiMuonProducer *
	                    process.DoubleDiMuonFitter *
	                    process.rootuplefourmu
	                    #process.rootupleJPsi *
	                    #process.rootuplePhi
	)
else:
		process.sequence = cms.Sequence(
							process.yfilter *
		                    process.slimmedMuonsWithTriggerSequence *
		                    #process.softMuons *
		                    process.JPsi2MuMuPAT *
		                    process.Phi2MuMuPAT *
		                    process.DiMuonFilteredJpsi *
		                    process.DiMuonFilteredPhi *
		                    process.DoubleDiMuonProducer *
		                    process.DoubleDiMuonFitter *
		                    process.rootuplefourmu
		                    #process.rootupleJPsi *
		                    #process.rootuplePhi
		)

process.p = cms.Path(process.sequence)
