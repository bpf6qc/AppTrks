import sys
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process ("Demo")

process.load ("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet ( input = cms.untracked.int32 (-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source ("PoolSource",
    fileNames = cms.untracked.vstring (
        'file:../../DisplacedMuonStudies/test/step3_250gev_300cm.root',
        )
)

process.source.duplicateCheckMode = cms.untracked.string ('noDuplicateCheck')

process.demo = cms.EDAnalyzer ('AnalyzeDecays',
    genParticleTag = cms.InputTag ("genParticlePlusGeant", ""),
    recoParticleTag = cms.InputTag ("standAloneMuons", "", "RECO"),
    isParticleGun = cms.untracked.bool (False),
    MaxEta = cms.untracked.double (1.0e12),
    quiet = cms.untracked.bool (True),
)

process.TFileService = cms.Service ("TFileService",
    fileName = cms.string ('file:output.root')
)

process.p = cms.Path (process.allMuonsGenParticlesMatch*process.demo)

