import FWCore.ParameterSet.Config as cms

flashggTagTester = cms.EDAnalyzer('FlashggTagTestAnalyzer',
                                  TagSorter = cms.untracked.InputTag('flashggTagSorter'),
                                  )
