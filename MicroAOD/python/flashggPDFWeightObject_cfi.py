import FWCore.ParameterSet.Config as cms

flashggPDFWeightObject = cms.EDProducer('FlashggPDFWeightProducer',
		LHEEventTag = cms.InputTag('externalLHEProducer'),
                GenTag      = cms.InputTag('generator'),
		tag = cms.untracked.string("initrwgt"),
                isStandardSample = cms.untracked.bool(False), # set to True for centrally produced Hgg samples, set to False for HH2B2G samples
                isThqSample = cms.untracked.bool(False), #set to False for centrally produced Hgg samples, set to True for THq samples
                pdfset = cms.untracked.string("NNPDF31_nnlo_as_0118_nf_4.LHgrid"), # for HH2B2G samples set here the pdf set name
		doAlphasWeights = cms.untracked.bool(True), # set to True for centrally produced Hgg samples, set to False for HH2B2G samples 
                doScaleWeights  = cms.untracked.bool(True),
		nPdfEigWeights = cms.uint32(60),
		mc2hessianCSV = cms.FileInPath('PhysicsTools/HepMCCandAlgos/data/NNPDF30_nlo_as_0118_nf_4_hessian_60.csv'),
                LHERunLabel = cms.string("externalLHEProducer"),
                Debug = cms.bool(True) 
	)
