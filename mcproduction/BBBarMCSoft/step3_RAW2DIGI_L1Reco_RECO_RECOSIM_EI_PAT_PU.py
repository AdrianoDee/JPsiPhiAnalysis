# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step3 --conditions auto:phase1_2018_realistic --pileup_input das:/RelValMinBias_13/CMSSW_10_2_1-102X_upgrade2018_realistic_v9_gcc7-v1/GEN-SIM-DIGI-RAW -n -1 --era Run2_2018 --pileup AVE_45_BX_25ns --geometry DB:Extended --datatier GEN-SIM-DIGI-RAW --filein file:step2_DIGI_L1_DIGI2RAW_HLT_RAW2DIGI_L1Reco_PU.root --datatier GEN-SIM-RECO,MINIAODSIM --eventcontent RECOSIM,MINIAODSIM -s RAW2DIGI,L1Reco,RECO,RECOSIM,EI,PAT --runUnscheduled
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:BBbar_JpsiFilter_SoftQCD_DIGI_HLT_RAW_L1_PU40.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('BBbar_JpsiFilter_SoftQCD_RECOSIM_PU40.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('BBbar_JpsiFilter_SoftQCD_MINIAODSIM_PU40.root'),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ),
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.input.nbPileupEvents.averageNumber = cms.double(45.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = cms.untracked.vstring(['/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/02B5104E-5294-E811-AAB5-0CC47A4D7658.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/1ADB395A-5294-E811-BCEF-0CC47A7C34E6.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/E2F2E8AD-5194-E811-BA97-0CC47A78A418.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/461E8C94-5194-E811-B16B-0CC47A4C8F10.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/DEECC36D-5294-E811-9D8F-0CC47A7C34B0.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/3ED1A853-5294-E811-9FAA-0CC47A4D7638.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/1C26F353-5294-E811-BE64-0CC47A7C3404.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/60B0CF4E-5294-E811-A652-0CC47A78A3EE.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/88CE8462-5294-E811-B103-0CC47A78A456.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/9C595FAD-5194-E811-B442-0CC47A78A3E8.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/642C6312-5394-E811-8D59-0CC47A4D767E.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/5AEB7564-5294-E811-9583-0025905B856C.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/A2FF2C32-5394-E811-8633-0025905A60B8.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/BA69222E-5394-E811-8EBA-0025905A607A.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/7C81DC1C-5394-E811-BD15-0025905A60B4.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/8E3C7633-5394-E811-904C-0025905B85BA.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/24D90787-5494-E811-8598-0025905B860E.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/EAA22B88-5994-E811-B480-0025905A60B4.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/CE9255CD-5294-E811-99F8-0CC47A4D761A.root'])
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,
                                process.eventinterpretaion_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,
                                process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,
                                process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,
                                process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,
                                process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,
                                process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,
                                process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,
                                process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,
                                process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,
                                process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,
                                process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,
                                process.Flag_METFilters,process.endjob_step,process.RECOSIMoutput_step,process.MINIAODSIMoutput_step)
process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
