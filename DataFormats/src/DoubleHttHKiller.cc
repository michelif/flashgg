#include "flashgg/DataFormats/interface/DoubleHttHKiller.h"
#include <cstdlib>

using namespace flashgg;
using namespace std;
using namespace edm;

double mW = 80.3;
double mT = 173.5;

DoubleHttHKiller::DoubleHttHKiller()
{  
  TMVAReady=0;
  //  mva_=2.;
}

DoubleHttHKiller::~DoubleHttHKiller() {}

void DoubleHttHKiller::initializeSelectionThresholds(){//FIXME do this configurable


    looseLeptonPtThreshold              =       10.;
    muEtaThreshold               =       2.4;
    muPFIsoSumRelThreshold         =       0.25;
    
    dRPhoLeptonThreshold      =   0.3;//FIXME we could do 0.4 for both of them
    dRJetLeptonThreshold      =   0.5;

    useElecMVARecipe            =   false;
    useElecLooseId              =   true;

    elecEtaThresholds.push_back(1.4442);elecEtaThresholds.push_back(1.566),elecEtaThresholds.push_back(2.5);
}

void DoubleHttHKiller::dummyPrint(){
  std::cout<<"mva is:::::"<<mva_<<std::endl;
}


void DoubleHttHKiller::initializeMVAVariables(std::string weightFile, std::vector<string> varvec){ 
  std::vector<string> defaultVariables = {"sumEt", "MET", "dPhi1", "dPhi2", "PhoJetMinDr", "njets>8", "Xtt0", "Xtt1", "pte1", "pte2", "ptmu1", "ptmu2", "fabs_CosThetaStar_CS", "fabs_CosTheta_bb"};
  if (varvec.size()<1){//FIXME check if this works
    varvec=defaultVariables;
  }
  
  for(unsigned int i=0;i<varvec.size();++i){
    ttHMVAVars_.push_back(varvec[i]);
  }
  
  if (weightFile.compare("") != 0){
    ttHMVAWeights_.assign(weightFile);
  }else{
    //ttHMVAWeights_.assign(std::string(getenv("$CMSSW_BASE"))+std::string("src/flashgg/MetaData/data/HHTagger/ttH_BDT.weights.xml"));
    ttHMVAWeights_ = std::string(getenv("CMSSW_BASE"))+std::string("/src/flashgg/MetaData/data/HHTagger/ttH_BDT.weights.xml");
  }

}

void DoubleHttHKiller::setupMVA(std::string File, std::vector<std::string> inVars){
  TMVAReady = 1;
  Reader_ = new TMVA::Reader();
  for (unsigned int iv = 0; iv < inVars.size(); iv++)
    mvaVars[inVars[iv]] = -10;

  orderedVars = inVars;

  std::cout << "[MVA::SetupMVA] MVA set with the following variables (attention, the order is important): " << std::endl;
  //    for ( std::map<std::string,float>::iterator it = mvaVars.begin(); it != mvaVars.end(); it++ ) {
  for ( unsigned int vv = 0; vv < orderedVars.size(); vv++) {
    std::map<std::string,float>::iterator it = mvaVars.find(orderedVars[vv]);
    Reader_->AddVariable(it->first, &it->second);
    std::cout << "\t" << it->first << std::endl;
  }
  Reader_->BookMVA("BDT", File); 

}

float DoubleHttHKiller::mvaDiscriminants(std::map<std::string,float> params)
{

  float mvaDis = -10.;

  if(TMVAReady == 0){
    std::cout << "[bbggNonResMVA::mvaDiscriminants] You haven't yet setup the NonResMVA reader!" << std::endl;
    return mvaDis;
  }

  if( params.size() != mvaVars.size()) {
    std::cout << "[bbggNonResMVA::mvaDiscriminants] Your input list does not match the variables list in the MVA method!" << std::endl;
    return mvaDis;
  }

  //  for ( std::map<std::string,float>::iterator it = mvaVars.begin(); it != mvaVars.end(); it++ ) {
  for ( unsigned int vv = 0; vv < orderedVars.size(); vv++) {
    std::map<std::string,float>::iterator it = params.find(orderedVars[vv]);
    if ( it == mvaVars.end() ) {
      std::cout << "[bbggNonResMVA::mvaDiscriminants] Variable that is in the training ## " << it->first << " not found in the parameters map provided!!" <<std::endl;
      return mvaDis;
    } else {
      mvaVars[it->first] = params[it->first];
      std::cout<<"before eval------"<<it->first<<":"<< mvaVars[it->first]<<" "<<params[it->first]<<std::endl;//FIXME remove
      varValues.push_back(mvaVars[it->first]);
    }
  }

  float dis = Reader_->EvaluateMVA("BDT");
   
  return dis;
}



std::vector<float> DoubleHttHKiller::XttCalculation(std::vector<flashgg::Jet> jetsCol, std::vector<flashgg::Jet> DiJet)
{

  std::vector<float> Xtt;
  Xtt.push_back(1000);  Xtt.push_back(0);   Xtt.push_back(0);  Xtt.push_back(1000); Xtt.push_back(0);  Xtt.push_back(0); 
  unsigned int WJetIndex00 = 1000,  WJetIndex01 = 1000, bJetIndex0 = 1000;

  
  for( unsigned int jetIndex0 = 0; jetIndex0 < jetsCol.size(); jetIndex0++ ) {
    //flashgg::Jet jet = jetsCol[jetIndex];
    const flashgg::Jet *jet0 = &(jetsCol[jetIndex0]);
    for( unsigned int jetIndex1 = jetIndex0+1; jetIndex1 < jetsCol.size(); jetIndex1++ ) {
      const flashgg::Jet *jet1 = &(jetsCol[jetIndex1]);
         
      LorentzVector Wcand = jet0->p4() + jet1->p4();
      LorentzVector tcand0 = Wcand + DiJet[0].p4();
      LorentzVector tcand1 = Wcand + DiJet[1].p4();
       
      float XW = TMath::Power((Wcand.M() - mW)/(0.1*mW),2);
      float XT0 = XW + TMath::Power((tcand0.M() - mT)/(0.1*mT),2);
      float XT1 = XW + TMath::Power((tcand1.M() - mT)/(0.1*mT),2);

      if (Xtt[0] > XT0) {
	Xtt[0] = XT0; Xtt[1] = Wcand.M(); Xtt[2] = tcand0.M(); 
	WJetIndex00 = jetIndex0, WJetIndex01 = jetIndex1, bJetIndex0 = 0;
      }
      if (Xtt[0] > XT1) {
	Xtt[0] = XT1; Xtt[1] = Wcand.M(); Xtt[2] = tcand1.M(); 
	WJetIndex00 = jetIndex0, WJetIndex01 = jetIndex1, bJetIndex0 = 1;
      }
    }
  }

  if (jetsCol.size() < 4) return Xtt;
  
  int  bJetIndex = 0;
  if ( bJetIndex0 == 0) bJetIndex = 1;
   
  for( unsigned int jetIndex0 = 0; jetIndex0 < jetsCol.size(); jetIndex0++ ) {
    //flashgg::Jet jet = jetsCol[jetIndex];
    const flashgg::Jet *jet0 = &(jetsCol[jetIndex0]);
    
    for( unsigned int jetIndex1 = jetIndex0+1; jetIndex1 < jetsCol.size(); jetIndex1++ ) {
        
      if( ( WJetIndex00 == jetIndex0 && WJetIndex01 == jetIndex1 ) || 
	  ( WJetIndex01 == jetIndex0 && WJetIndex00 == jetIndex1 )) continue;

      const flashgg::Jet *jet1 = &(jetsCol[jetIndex1]);

      LorentzVector Wcand = jet0->p4() + jet1->p4();
      LorentzVector tcand = Wcand + DiJet[bJetIndex].p4();

      float XW = TMath::Power((Wcand.M() - mW)/(0.1*mW),2);
      float XT = XW + TMath::Power((tcand.M() - mT)/(0.1*mT),2);

      if (Xtt[3] > XT) {
	Xtt[3] = XT;  Xtt[4] = Wcand.M(); Xtt[5] = tcand.M(); 
      }

    }
  }
  
  
  return Xtt;
}

std::vector<edm::Ptr<flashgg::Electron> > DoubleHttHKiller::filterElectrons( const std::vector<edm::Ptr<flashgg::Electron> > &electronPointers, 
								      flashgg::DiPhotonCandidate dipho, 
								      LorentzVector leadingJet, LorentzVector subleadingJet,
								      double dRPhoLeptonThreshold, double dRJetLeptonThreshold){
  std::vector<edm::Ptr<flashgg::Electron> > goodElectrons;

  double Phi_Pho1 = dipho.leadingPhoton()->superCluster()->phi();
  double Phi_Pho2 = dipho.subLeadingPhoton()->superCluster()->phi();
  double Eta_Pho1 = dipho.leadingPhoton()->superCluster()->eta();
  double Eta_Pho2 = dipho.subLeadingPhoton()->superCluster()->eta();

  double Phi_Jet1 =  leadingJet.Phi();
  double Phi_Jet2 =  subleadingJet.Phi();
  double Eta_Jet1 =  leadingJet.Eta();
  double Eta_Jet2 =  subleadingJet.Eta();

  for( unsigned int electronIndex = 0; electronIndex < electronPointers.size(); electronIndex++ ) {
    edm::Ptr<flashgg::Electron> electron = electronPointers[electronIndex];
    double Eta_Lepton = electron->superCluster()->eta(), Phi_Lepton = electron->superCluster()->phi();

    float dRLeadPhoLepton      = deltaR( Eta_Lepton, Phi_Lepton, Eta_Pho1, Phi_Pho1);
    float dRSubLeadPhoLepton   = deltaR( Eta_Lepton, Phi_Lepton, Eta_Pho2, Phi_Pho2);
    
    float dRLeadJetLepton    = deltaR( Eta_Lepton, Phi_Lepton, Eta_Jet1, Phi_Jet1);
    float dRSubLeadJetLepton = deltaR( Eta_Lepton, Phi_Lepton, Eta_Jet2, Phi_Jet2);

    
    if( dRLeadPhoLepton < dRPhoLeptonThreshold || dRSubLeadPhoLepton < dRPhoLeptonThreshold || 
	dRLeadJetLepton < dRJetLeptonThreshold || dRSubLeadJetLepton < dRJetLeptonThreshold) continue; 
    
    goodElectrons.push_back( electron );
  }
  return goodElectrons;
  
}

std::vector<edm::Ptr<flashgg::Muon> > DoubleHttHKiller::filterMuons( const std::vector<edm::Ptr<flashgg::Muon> > &muonPointers, 
							      flashgg::DiPhotonCandidate dipho, 
							      LorentzVector leadingJet, LorentzVector subleadingJet,
							      double dRPhoLeptonThreshold, double dRJetLeptonThreshold){

  std::vector<edm::Ptr<flashgg::Muon> > goodMuons;

  double Phi_Pho1 = dipho.leadingPhoton()->superCluster()->phi();
  double Phi_Pho2 = dipho.subLeadingPhoton()->superCluster()->phi();
  double Eta_Pho1 = dipho.leadingPhoton()->superCluster()->eta();
  double Eta_Pho2 = dipho.subLeadingPhoton()->superCluster()->eta();

  double Phi_Jet1 =  leadingJet.Phi();
  double Phi_Jet2 =  subleadingJet.Phi();
  double Eta_Jet1 =  leadingJet.Eta();
  double Eta_Jet2 =  subleadingJet.Eta();

  for( unsigned int muonIndex = 0; muonIndex < muonPointers.size(); muonIndex++ ) {
    edm::Ptr<flashgg::Muon> muon = muonPointers[muonIndex];
    double Eta_Lepton = muon->eta(), Phi_Lepton = muon->phi();

    float dRLeadPhoLepton      = deltaR( Eta_Lepton, Phi_Lepton, Eta_Pho1, Phi_Pho1);
    float dRSubLeadPhoLepton   = deltaR( Eta_Lepton, Phi_Lepton, Eta_Pho2, Phi_Pho2);
    
    float dRLeadJetLepton    = deltaR( Eta_Lepton, Phi_Lepton, Eta_Jet1, Phi_Jet1);
    float dRSubLeadJetLepton = deltaR( Eta_Lepton, Phi_Lepton, Eta_Jet2, Phi_Jet2);

    
    if( dRLeadPhoLepton < dRPhoLeptonThreshold || dRSubLeadPhoLepton < dRPhoLeptonThreshold || 
	dRLeadJetLepton < dRJetLeptonThreshold || dRSubLeadJetLepton < dRJetLeptonThreshold) continue; 
    
    goodMuons.push_back( muon );
  }
  return goodMuons;
  
}


std::vector<float> DoubleHttHKiller::getVariableValues() const{
  return varValues;
}
