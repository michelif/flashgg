#ifndef flashgg_DoubleHttHKiller
#define flashgg_DoubleHttHKiller

#include "TLorentzVector.h"

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"

#include "flashgg/Taggers/interface/FunctionHelpers.h"
#include "TMVA/Reader.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"

namespace flashgg {

    class DoubleHttHKiller: public DiPhotonTagBase, public reco::LeafCandidate
    {
    public:
        DoubleHttHKiller();
        ~DoubleHttHKiller();

        //variables
        double muEtaThreshold, muPFIsoSumRelThreshold;// deltaRMuonPhoThreshold;
        double dRPhoLeptonThreshold, dRJetLeptonThreshold;

        double looseLeptonPtThreshold;
        bool useElecMVARecipe, useElecLooseId;
        std::vector<double> elecEtaThresholds;

        int njets, nelecs, nmus,  nelecs_loose, nmus_loose;
        vector<float> Xtt;
        float Xtt0, Xtt1, MjjW0, MjjW1, Mjjbt0, Mjjbt1;


        TMVA::Reader *Reader_;
        bool TMVAReady;
        //MVA variables
        std::map<std::string, float> mvaVars;
        std::vector<std::string> orderedVars;


        void initializeSelectionThresholds();
        void dummyPrint();
        void initializeMVAVariables(std::string weightFile="",std::vector<string> varvec={});  

        std::string ttHMVAWeights_;
        std::vector<std::string> ttHMVAVars_;
        void setupMVA(std::string File, std::vector<std::string> inVars);
        float mvaDiscriminants(std::map<std::string,float>) ;
        std::vector<float> XttCalculation(std::vector<flashgg::Jet> allJetsCol, std::vector<flashgg::Jet> DiJet);
        std::vector<edm::Ptr<flashgg::Electron> > filterElectrons( const std::vector<edm::Ptr<flashgg::Electron> > &electronPointers, 
                                                                                     flashgg::DiPhotonCandidate dipho, 
                                                                                     LorentzVector leadingJet, LorentzVector subleadingJet,
                                                                                     double dRPhoLeptonThreshold, double dRJetLeptonThreshold);
        std::vector<edm::Ptr<flashgg::Muon> > filterMuons( const std::vector<edm::Ptr<flashgg::Muon> > &muonPointers, 
                                                                      flashgg::DiPhotonCandidate dipho, 
                                                                      LorentzVector leadingJet, LorentzVector subleadingJet,
                                                                      double dRPhoLeptonThreshold, double dRJetLeptonThreshold);
        std::vector<float> varValues;
        std::vector<float> getVariableValues() const;

    private:
        double mva_;
    };
}


#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

