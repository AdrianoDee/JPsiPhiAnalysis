import FWCore.ParameterSet.Config as cms

### ==== This is our version of the patPFCandsWithTrigger using MINIAOD

### unpack them

unpackedPatTriggerPF = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
  patTriggerObjectsStandAlone = cms.InputTag( 'slimmedPatTrigger' ),
  triggerResults              = cms.InputTag( 'TriggerResults::HLT' ),
  unpackFilterLabels          = cms.bool( True )
)

### then perform a match for all HLT triggers of interest
PATPFCandsTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "packedPFCandidates" ),
    matched = cms.InputTag( "unpackedPatTriggerPF" ),
    matchedCuts = cms.string(""),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)

slimmedPFCandsTriggerMatchers = cms.Sequence(
      PATPFCandsTriggerMatchHLT
)

slimmedPFCandsTriggerMatchersTags = [
    cms.InputTag('PATPFCandsTriggerMatchHLT'),
]

### Embed them
slimmedPFCandsWithTrigger = cms.EDProducer( "PATTriggerMatchPFEmbedder",
    src     = cms.InputTag(  "packedPFCandidates" ),
    matches = cms.VInputTag()
)
slimmedPFCandsWithTrigger.matches += slimmedPFCandsTriggerMatchersTags

### Run the whole trigger Sequence
slimmedPFCandsWithTriggerSequence = cms.Sequence(
    unpackedPatTriggerPF *
    slimmedPFCandsTriggerMatchers *
    slimmedPFCandsWithTrigger
)
