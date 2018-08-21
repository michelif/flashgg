#include "flashgg/DataFormats/interface/DoubleHttHKiller.h"

using namespace flashgg;
using namespace std;
using namespace edm;


DoubleHttHKiller::DoubleHttHKiller()
{  
  TMVAReady=0;
  //  mva_=2.;
}

DoubleHttHKiller::~DoubleHttHKiller() {}

void DoubleHttHKiller::initializeSelectionThresholds(){


    muPtThreshold              =       20;
    muEtaThreshold               =       2.4;
    muPFIsoSumRelThreshold         =       0.25;

    dRPhoLeptonThreshold      =   0.3;
    dRJetLeptonThreshold      =   0.5;

    elecPtThreshold             =       20;
    useElecMVARecipe            =   false;
    useElecLooseId              =   true;

    elecEtaThresholds.push_back(1.4442);elecEtaThresholds.push_back(1.566),elecEtaThresholds.push_back(2.5);//FIXME this is not needed    
    // Xtt0                MW0                 Mt0                Xtt1                 MW1               Mt1
    Xtt.push_back(1000); Xtt.push_back(0); Xtt.push_back(0); Xtt.push_back(1000); Xtt.push_back(0); Xtt.push_back(0); 
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
    ttHMVAWeights_.assign("MetaData/data/HHTagger/ttH_BDT.weights.xml");
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

  for ( std::map<std::string,float>::iterator it = mvaVars.begin(); it != mvaVars.end(); it++ ) {
    std::map<std::string,float>::iterator toFind = params.find(it->first);
    if ( toFind == mvaVars.end() ) {
      std::cout << "[bbggNonResMVA::mvaDiscriminants] Variable that is in the training ## " << it->first << " not found in the parameters map provided!!" <<std::endl;
      return mvaDis;
    } else {
      mvaVars[it->first] = params[it->first];
    }
  }

  float dis = Reader_->EvaluateMVA("BDT");
   
  return dis;
}
