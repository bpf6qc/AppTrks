import FWCore.ParameterSet.Config as cms
import copy

genMinimal = cms.PSet(
    name = cms.string("genMinimal"),
    triggers = cms.vstring(),
    cuts = cms.VPSet(
        cms.PSet(
            inputCollection = cms.vstring("mcparticles"),
            cutString = cms.string("abs ( pdgId ) == 1000024"),
            numberRequired = cms.string(">= 1")
        ),
        cms.PSet(
            inputCollection = cms.vstring("mcparticles"),
            cutString = cms.string("pt > 10"),
            numberRequired = cms.string(">= 1")
        ),
    )
)

muMinimal = cms.PSet(
    name = cms.string("muMinimal"),
    triggers = cms.vstring(),
    cuts = cms.VPSet(
        # good pv
        cms.PSet(
            inputCollection = cms.vstring("primaryvertexs"),
            cutString = cms.string("isValid > 0 && ndof >= 4"),
            numberRequired = cms.string(">= 1")
        ),
        # muon pt
        cms.PSet(
            inputCollection = cms.vstring("muons"),
            cutString = cms.string("pt > 10"),
            numberRequired = cms.string(">= 1")
        ),
        # muon eta
        cms.PSet(
            inputCollection = cms.vstring("muons"),
            cutString = cms.string("abs(eta) < 2.4"),
            numberRequired = cms.string(">= 1")
        ),
        # ask for innerTrack to exist
        cms.PSet(
            inputCollection = cms.vstring("muons"),
            cutString = cms.string("innerTrack.hitPattern_.trackerLayersWithMeasurement >= 0"),
            numberRequired = cms.string(">= 1")
        ),
        
    )
)
