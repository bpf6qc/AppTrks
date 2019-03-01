import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import math
import os
import sys

process = cms.Process('ZPeakAnalyzer')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'file:step3_250gev_300cm.root',
                                )
    )

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('output_zpeak.root')
    )

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.maxEvents = cms.untracked.PSet (
    input = cms.untracked.int32 (-1)
    )

process.ZPeakAnalyzer = cms.EDAnalyzer('ZPeakAnalyzer',
                                        triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
                                        trackCollectionLabel = cms.InputTag("displacedGlobalMuons", "", "RECO"),
                                        genCollectionLabel = cms.InputTag("genParticlePlusGeant", "")
    )
process.myPath = cms.Path (process.ZPeakAnalyzer)
