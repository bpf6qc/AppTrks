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

                                )
    )

process.TFileService = cms.Service('TFileService',
                                   fileName = cms.string('ZPeak.root')
    )

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.maxEvents = cms.untracked.PSet (
    input = cms.untracked.int32 (-1)
    )

process.ZPeakPlotter = cms.EDAnalyzer('ZPeakPlotter',
                                      triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),
                                      recoTrack = cms.InputTag('displacedGlobalMuons', '', 'RECO'),
                                      recoTrack2 = cms.InputTag('displacedGlobalMuons', '', 'RECO'),
                                      #globalMuons, displacedStandAloneMuons
    )
process.myPath = cms.Path (process.ZPeakPlotter)
