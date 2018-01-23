import FWCore.ParameterSet.Config as cms

from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
from flashgg.Taggers.flashggTags_cff import flashggUnpackedJets

#recoJetCollections = UnpackedJetCollectionVInputTag

#print recoJetCollections

#for icoll,coll in enumerate(recoJetCollections):
flashggbRegressionProducer= cms.EDProducer('flashggbRegressionProducer',
#                                               JetTag=coll,
                                               JetTag=cms.InputTag("flashggUnpackedJets","0"),
                                               )

