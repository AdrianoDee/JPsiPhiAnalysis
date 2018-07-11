import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from GeneratorInterface.EvtGenInterface.EvtGenSetting_cff import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.0),
    maxEventsToPrint = cms.untracked.int32(0),
    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            operates_on_particles = cms.vint32(),
            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_2014.pdl'),
            user_decay_file = cms.vstring('GeneratorInterface/EvtGenInterface/data/Onia_mumu.dec'),
            list_forced_decays = cms.vstring('MyPhi'),
            convertPythiaCodes = cms.untracked.bool(False)
        ),
        parameterSets = cms.vstring('EvtGen130')
    ),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        processParameters = cms.vstring(
        'HardQCD:all = on',
            'PhaseSpace:pTHatMin = 8.',

       ),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            'pythia8CUEP8M1Settings',
            'processParameters',
        )
    )
)

generator.PythiaParameters.processParameters.extend(EvtGenExtraParticles)

bfilter = cms.EDFilter("PythiaFilter",
                       ParticleID = cms.untracked.int32(5)
                       )

kkgenfilter = cms.EDFilter(
    "PythiaDauVFilter",
    verbose         = cms.untracked.int32(0),
    NumberDaughters = cms.untracked.int32(2),
    ParticleID      = cms.untracked.int32(333),
    DaughterIDs     = cms.untracked.vint32(321,-321),
    MinPt           = cms.untracked.vdouble(1.0, 1.0),
    MinEta          = cms.untracked.vdouble(-2.5, -2.5),
    MaxEta          = cms.untracked.vdouble( 2.5,  2.5)
    )

ProductionFilterSequence = cms.Sequence(generator*bfilter*kkgenfilter)
