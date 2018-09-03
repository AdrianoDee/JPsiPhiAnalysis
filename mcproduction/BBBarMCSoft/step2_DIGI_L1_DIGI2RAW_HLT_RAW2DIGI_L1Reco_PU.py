# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step2 --pileup_input das:/RelValMinBias_13/CMSSW_10_2_1-102X_upgrade2018_realistic_v9_gcc7-v1/GEN-SIM-DIGI-RAW -n -1 --era Run2_2018 --pileup AVE_45_BX_25ns --geometry DB:Extended --datatier GEN-SIM-DIGI-RAW --filein file:BBbar_Jpsi_Phi_Filter_HardQCD_10_cfg_py_GEN_SIM.root --eventcontent FEVTDEBUGHLT -s DIGI:pdigi_valid,L1,DIGI2RAW,HLT:@relval2018,RAW2DIGI,L1Reco --conditions auto:phase1_2018_realistic
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('HLT',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_2018v32_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('file:BBbar_JpsiFilter_SoftQCD_GEN_SIM.root'),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_genParticles_*_*',
        'drop *_genParticlesForJets_*_*',
        'drop *_kt4GenJets_*_*',
        'drop *_kt6GenJets_*_*',
        'drop *_iterativeCone5GenJets_*_*',
        'drop *_ak4GenJets_*_*',
        'drop *_ak7GenJets_*_*',
        'drop *_ak8GenJets_*_*',
        'drop *_ak4GenJetsNoNu_*_*',
        'drop *_ak8GenJetsNoNu_*_*',
        'drop *_genCandidatesForMET_*_*',
        'drop *_genParticlesForMETAllVisible_*_*',
        'drop *_genMetCalo_*_*',
        'drop *_genMetCaloAndNonPrompt_*_*',
        'drop *_genMetTrue_*_*',
        'drop *_genMetIC5GenJs_*_*'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:BBbar_JpsiFilter_SoftQCD_DIGI_HLT_RAW_L1_PU40.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.input.nbPileupEvents.averageNumber = cms.double(45.000000)
process.mix.bunchspace = cms.int32(25)
process.mix.minBunch = cms.int32(-12)
process.mix.maxBunch = cms.int32(3)
process.mix.input.fileNames = cms.untracked.vstring(['/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/02B5104E-5294-E811-AAB5-0CC47A4D7658.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/1ADB395A-5294-E811-BCEF-0CC47A7C34E6.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/E2F2E8AD-5194-E811-BA97-0CC47A78A418.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/461E8C94-5194-E811-B16B-0CC47A4C8F10.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/DEECC36D-5294-E811-9D8F-0CC47A7C34B0.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/3ED1A853-5294-E811-9FAA-0CC47A4D7638.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/1C26F353-5294-E811-BE64-0CC47A7C3404.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/60B0CF4E-5294-E811-A652-0CC47A78A3EE.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/88CE8462-5294-E811-B103-0CC47A78A456.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/9C595FAD-5194-E811-B442-0CC47A78A3E8.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/642C6312-5394-E811-8D59-0CC47A4D767E.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/5AEB7564-5294-E811-9583-0025905B856C.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/A2FF2C32-5394-E811-8633-0025905A60B8.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/BA69222E-5394-E811-8EBA-0025905A607A.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/7C81DC1C-5394-E811-BD15-0025905A60B4.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/8E3C7633-5394-E811-904C-0025905B85BA.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/24D90787-5494-E811-8598-0025905B860E.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/EAA22B88-5994-E811-B480-0025905A60B4.root', '/store/relval/CMSSW_10_2_1/RelValMinBias_13/GEN-SIM-DIGI-RAW/102X_upgrade2018_realistic_v9_gcc7-v1/10000/CE9255CD-5294-E811-99F8-0CC47A4D761A.root'])
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.endjob_step,process.FEVTDEBUGHLToutput_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
