import FWCore.ParameterSet.Config as cms

def variablesToDump(customize):
    var_workspace = [
             "Mjj := dijet().M()"
    ]
    variables = [ "leadingJet_bDis := leadJet().bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",#FIXME make the btag type configurable?
             "subleadingJet_bDis := subleadJet().bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",
             "absCosThetaStar_CS := abs(getCosThetaStar_CS(6500))",#FIXME get energy from somewhere?
             "absCosTheta_bb := abs(CosThetaAngles()[1])",
             "absCosTheta_gg := abs(CosThetaAngles()[0])",
             "diphotonCandidatePtOverdiHiggsM := diphotonPtOverM()",
             "dijetCandidatePtOverdiHiggsM := dijetPtOverM()",
             "customLeadingPhotonIDMVA := diPhoton.leadingView.phoIdMvaWrtChosenVtx",
             "customSubLeadingPhotonIDMVA := diPhoton.subLeadingView.phoIdMvaWrtChosenVtx",
             "leadingPhotonSigOverE := diPhoton.leadingPhoton.sigEOverE",
             "subleadingPhotonSigOverE := diPhoton.subLeadingPhoton.sigEOverE",
             "sigmaMOverM := sqrt(0.5*(diPhoton.leadingPhoton.sigEOverE*diPhoton.leadingPhoton.sigEOverE + diPhoton.subLeadingPhoton.sigEOverE*diPhoton.subLeadingPhoton.sigEOverE))",
             "sigmaMOverMDecorr := getSigmaMDecorr()",
             "PhoJetMinDr := getPhoJetMinDr()",#up to here input variables to MVA
             "HHbbggMVA := MVA()",
             "MX := MX()",
             "Mjj := dijet().M()",
             "dijet_pt := dijet().pt",
             "dijet_eta := dijet().eta",
             "dijet_phi := dijet().phi",
             "diphoton_pt := diPhoton.pt",
             "diphoton_eta := diPhoton.eta",
             "diphoton_phi := diPhoton.phi",
 
             "diHiggs_pt := getdiHiggsP4().pt()",
             "diHiggs_mass := getdiHiggsP4().M()",
             "diHiggs_eta :=  getdiHiggsP4().eta()",
             "diHiggs_phi := getdiHiggsP4().phi()",
             "category := categoryNumber()",

             "leadingPhoton_pt := diPhoton.leadingPhoton.pt",
             "leadingPhoton_eta := diPhoton.leadingPhoton.eta",
             "leadingPhoton_phi := diPhoton.leadingPhoton.phi",
             "subleadingPhoton_pt := diPhoton.subLeadingPhoton.pt",
             "subleadingPhoton_eta := diPhoton.subLeadingPhoton.eta",
             "subleadingPhoton_phi := diPhoton.subLeadingPhoton.phi",

             "leadingJet_pt := leadJet().pt",
             "leadingJet_eta := leadJet().eta",
             "leadingJet_phi := leadJet().phi",
             "leadingJet_mass := leadJet().p4().M()",

             "subleadingJet_pt := subleadJet().pt",
             "subleadingJet_eta := subleadJet().eta",
             "subleadingJet_phi := subleadJet().phi",
             "subleadingJet_mass := subleadJet().p4().M()",

                  
             "HHttHMVA_sumEt := getttHMVAVariableValues()[0]",#FIXME add proper names of variables
             "HHttHMVA_MET := getttHMVAVariableValues()[1]",
             "HHttHMVA_dPhi1 := getttHMVAVariableValues()[2]",
             "HHttHMVA_dPhi2 := getttHMVAVariableValues()[3]",
             "HHttHMVA_PhoJetMinDr := getttHMVAVariableValues()[4]",
             "HHttHMVA_njets8 := getttHMVAVariableValues()[5]",
             "HHttHMVA_Xtt0 := getttHMVAVariableValues()[6]",
             "HHttHMVA_Xtt1 := getttHMVAVariableValues()[7]",
             "HHttHMVA_pte1 := getttHMVAVariableValues()[8]",
             "HHttHMVA_pte2 := getttHMVAVariableValues()[9]",
             "HHttHMVA_ptmu1 := getttHMVAVariableValues()[10]",
             "HHttHMVA_ptmu2 := getttHMVAVariableValues()[11]",
             "HHttHMVA_fabs_CosThetaStar_CS := getttHMVAVariableValues()[12]",
             "HHttHMVA_fabs_CosTheta_bb := getttHMVAVariableValues()[13]",

             "HHttHMVA_mva := getttHMVA()"  


    ]
    if customize.doBJetRegression : variables +=[
             "leadingJet_bRegNNCorr := leadJet().userFloat('bRegNNCorr')",
             "leadingJet_bRegNNResolution := leadJet().userFloat('bRegNNResolution')",
             "subleadingJet_bRegNNCorr := subleadJet().userFloat('bRegNNCorr')",
             "subleadingJet_bRegNNResolution := subleadJet().userFloat('bRegNNResolution')",
             "sigmaMJets := getSigmaMOverMJets()"
    ]
    if customize.dumpWorkspace == False : return variables
    else : return var_workspace

def tagList(customize,process):
    return [ ["DoubleHTag",12] ]#12 is the number of categories?


def customizeTagSequence(customize,process):
    process.load("flashgg.Taggers.flashggDoubleHTag_cff")
    from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
    if customize.doBJetRegression : process.flashggDoubleHTag.JetTags = cms.VInputTag( ["bRegProducer%d" % icoll for icoll,coll in enumerate(UnpackedJetCollectionVInputTag) ] )
    ## customize here (regression, kin-fit, MVA...)
    ## process.flashggTagSequence += process.flashggDoubleHTagSequence
#    pos = process.flashggTagSequence.index(process.flashggTagSorter) - 1 
#    pos = process.flashggTagSequence.index(process.flashggTagSorter) -1  
#    process.flashggTagSequence.insert( pos, process.flashggDoubleHTagSequence )
    ## for step in process.flashggDoubleHTagSequence:
    ##     process.flashggTagSequence.append(step)

    if customize.doubleHTagsOnly:
        process.flashggTagSequence.remove(process.flashggVBFTag)
        process.flashggTagSequence.remove(process.flashggTTHLeptonicTag)
        process.flashggTagSequence.remove(process.flashggTTHHadronicTag)
        process.flashggTagSequence.remove(process.flashggVHEtTag)
        process.flashggTagSequence.remove(process.flashggVHLooseTag)
        process.flashggTagSequence.remove(process.flashggVHTightTag)
        process.flashggTagSequence.remove(process.flashggVHMetTag)
        process.flashggTagSequence.remove(process.flashggWHLeptonicTag)
        process.flashggTagSequence.remove(process.flashggZHLeptonicTag)
        process.flashggTagSequence.remove(process.flashggVHLeptonicLooseTag)
        process.flashggTagSequence.remove(process.flashggVHHadronicTag)
        process.flashggTagSequence.remove(process.flashggVBFMVA)
        process.flashggTagSequence.remove(process.flashggVBFDiPhoDiJetMVA)


        process.flashggTagSequence.replace(process.flashggUntagged, process.flashggDoubleHTagSequence)   



