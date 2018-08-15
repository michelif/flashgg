#include "flashgg/DataFormats/interface/DoubleHttHKiller.h"

using namespace flashgg;

DoubleHttHKiller::DoubleHttHKiller()
{  
  mva_=2.;
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
