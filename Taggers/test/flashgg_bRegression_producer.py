import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggMicroAOD")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v13')
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag = GlobalTag(process.GlobalTag,'80X_mcRun2_asymptotic_2016_miniAODv2')
else:
raise Exception,"The default setup for microAODstd.py does not support releases other than 76X and 80X"


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/group/phys_higgs/cmshgg/micheli/flashgg/RunIISummer16-2_4_2-25ns_Moriond17/2_4_2/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16-2_4_2-25ns_Moriond17-2_4_2-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180118_055550/0000/myMicroAODOutputFile_1.root"))

process.flashggbRegression = cms.EDProducer('FlashggbRegressionProducer',
                                                JetTag="flashggJets"
                                                )

p = cms.Path(process.flashggbRegression)
