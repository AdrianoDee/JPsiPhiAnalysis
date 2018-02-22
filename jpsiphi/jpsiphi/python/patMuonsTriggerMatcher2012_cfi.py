import FWCore.ParameterSet.Config as cms

### ==== This is our version of the patMuonsWithTrigger using MINIAOD

### unpack them
unpackedPatTrigger = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
  patTriggerObjectsStandAlone = cms.InputTag( 'selectedPatTrigger' ),
  triggerResults              = cms.InputTag( 'TriggerResults::HLT' ),
  unpackFilterLabels          = cms.bool( True )
)

### then perform a match for all HLT triggers of interest
PATmuonTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRDPtLessByR",
    src     = cms.InputTag( "patMuons" ),
    matched = cms.InputTag( "unpackedPatTrigger" ),
    matchedCuts = cms.string(""),
    maxDPtRel = cms.double( 0.5 ),
    maxDeltaR = cms.double( 0.5 ),
    resolveAmbiguities    = cms.bool( True ),
    resolveByMatchQuality = cms.bool( True )
)

PATmuonMatchHLTL2   = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltL2MuonCandidates")'), maxDeltaR = 0.3, maxDPtRel = 10.0)           #maxDeltaR Changed accordingly to Zoltan tuning. It was: 1.2
PATmuonMatchHLTL3   = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltIterL3FromL2MuonCandidates")'), maxDeltaR = 0.1, maxDPtRel = 10.0) #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5
PATmuonMatchHLTL3v2 = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltIterL3MuonCandidates")'), maxDeltaR = 0.1, maxDPtRel = 10.0)       #needed because name changed
PATmuonMatchHLTL3T  = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltGlbTrkMuonCands")'), maxDeltaR = 0.1, maxDPtRel = 10.0)            #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5
PATmuonMatchHLTTkMu = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('coll("hltHighPtTkMuonCands")'), maxDeltaR = 0.1, maxDPtRel = 10.0)          #maxDeltaR Changed accordingly to Zoltan tuning. It was: 0.5
PATmuonMatchHLTMuon = PATmuonTriggerMatchHLT.clone(matchedCuts = cms.string('type( "TriggerL1Mu" ) || type( "TriggerMuon" )'), maxDeltaR = 0.1, maxDPtRel = 10.0)           #maxDeltaR Changed accordingly to Zoltan tuning. It was: 1.2

patMuonsTriggerMatchers1Mu = cms.Sequence(
      PATmuonMatchHLTL2 +
      PATmuonMatchHLTL3 +
      PATmuonMatchHLTL3v2 +
      PATmuonMatchHLTL3T +
      PATmuonMatchHLTTkMu +
      PATmuonMatchHLTMuon
)

patMuonsTriggerMatchers1MuInputTags = [
    cms.InputTag('PATmuonMatchHLTL2'),
    cms.InputTag('PATmuonMatchHLTL3'),
    cms.InputTag('PATmuonMatchHLTL3v2'),
    cms.InputTag('PATmuonMatchHLTL3T'),
    cms.InputTag('PATmuonMatchHLTTkMu'),
    cms.InputTag('PATmuonMatchHLTMuon'),
]

### Embed them
patMuonsWithTrigger = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
    src     = cms.InputTag(  "patMuons" ),
    matches = cms.VInputTag()
)
patMuonsWithTrigger.matches += patMuonsTriggerMatchers1MuInputTags

### Run the whole trigger Sequence
patMuonsWithTriggerSequence = cms.Sequence(
    unpackedPatTrigger *
    patMuonsTriggerMatchers1Mu *
    patMuonsWithTrigger
)



from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import  addMCinfo, useExistingPATMuons, useL1MatchingWindowForSinglets, changeTriggerProcessName, switchOffAmbiguityResolution, addDiMuonTriggers


# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")

process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)
patMuons.embedTrack = cms.bool(True)
patMuons.embedPickyMuon = cms.bool(False)
patMuons.embedTpfmsMuon = cms.bool(False)


if process.isMC:
        addMCinfo(process)
        # since we match inner tracks, keep the matching tight and make it one-to-one
        process.muonMatch.maxDeltaR = 0.05
        process.muonMatch.resolveByMatchQuality = True

addDiMuonTriggers(process)
useExistingPATMuons(process,'cleanPatMuons',addL1Info=False)
changeTriggerProcessName(process, 'HLT')
switchOffAmbiguityResolution(process) # Switch off ambiguity resolution: allow multiple reco muons to match to the same trigger muon
useL1MatchingWindowForSinglets(process)

process.muonL1Info.maxDeltaR     = 0.3
process.muonL1Info.fallbackToME1 = True
process.muonMatchHLTL1.maxDeltaR     = 0.3
process.muonMatchHLTL1.fallbackToME1 = True
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL2.maxDPtRel = 10.0
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0
process.muonMatchHLTCtfTrack.maxDeltaR = 0.1
process.muonMatchHLTCtfTrack.maxDPtRel = 10.0
process.muonMatchHLTTrackMu.maxDeltaR = 0.1
process.muonMatchHLTTrackMu.maxDPtRel = 10.0
