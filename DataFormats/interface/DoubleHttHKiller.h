#ifndef flashgg_DoubleHttHKiller
#define flashgg_DoubleHttHKiller

#include "TLorentzVector.h"

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "flashgg/Taggers/interface/FunctionHelpers.h"
#include "flashgg/Taggers/interface/genericMVA.h"

namespace flashgg {

    class DoubleHttHKiller: public DiPhotonTagBase, public reco::LeafCandidate
    {
    public:
        DoubleHttHKiller();
        ~DoubleHttHKiller();

        /// ttH tokens move it to producer
        //        edm::EDGetTokenT<double> rhoToken_;
        //        edm::EDGetTokenT<edm::View<flashgg::Electron> > electronToken_;
        //        edm::EDGetTokenT<edm::View<flashgg::Muon> > muonToken_;
        //variables
        double muPtThreshold, muEtaThreshold, muPFIsoSumRelThreshold;// deltaRMuonPhoThreshold;
        double dRPhoLeptonThreshold, dRJetLeptonThreshold;

        double elecPtThreshold;
        bool useElecMVARecipe, useElecLooseId;
        std::vector<double> elecEtaThresholds;

        int njets, nelecs, nmus,  nelecs_loose, nmus_loose;
        vector<float> Xtt;
        float Xtt0, Xtt1, MjjW0, MjjW1, Mjjbt0, Mjjbt1;


        genericMVA ttHMVA_;

        void initializeSelectionThresholds();
        void dummyPrint();

    private:
        double mva_;
        //        ttHMVAWeights_;
        //        ttHMVAVars_;
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

