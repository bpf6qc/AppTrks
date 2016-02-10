import sys
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process ("Demo")

process.load ("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet ( input = cms.untracked.int32 (-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source ("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring (
        'file:../../StandardAnalysis/test/toymodel_250gev_300cm.root',
        )
)

process.source.duplicateCheckMode = cms.untracked.string ('noDuplicateCheck')

process.demo = cms.EDAnalyzer ('AnalyzeDecays',
    genParticleTag = cms.InputTag ("genParticlePlusGeant", ""),
    isParticleGun = cms.untracked.bool (False),
    MaxEta = cms.untracked.double (1.0e12),
    quiet = cms.untracked.bool (True),
)


process.TFileService = cms.Service ("TFileService",
    fileName = cms.string ('file:output.root')
)

process.p = cms.Path (process.demo)

