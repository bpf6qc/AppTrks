import FWCore.ParameterSet.Config as cms

###########################################################
##### Set up process #####
###########################################################

process = cms.Process ('OSUAnalysis')
process.load ('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.source = cms.Source ('PoolSource',
    fileNames = cms.untracked.vstring (
#        'file:charginoPartGun_GEN_SIM.root', 
	'file:step1/step1_250gev_100cm.root',
    )
)

process.maxEvents = cms.untracked.PSet (
    input = cms.untracked.int32 (1)
)

###########################################################
##### Set up the analyzer #####
###########################################################

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.particleListDrawer = cms.EDAnalyzer ('ParticleListDrawer',
    maxEventsToPrint = cms.untracked.int32(-1),
    printVertex = cms.untracked.bool(True),
#    src = cms.InputTag("genParticles")
    src = cms.InputTag("genParticlePlusGeant")
)

process.myPath = cms.Path (process.particleListDrawer)
