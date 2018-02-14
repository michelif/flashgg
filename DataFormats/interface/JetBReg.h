#ifndef FLASHgg_JetBReg_h
#define FLASHgg_JetBReg_h

#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"
//#include "flashgg/DataFormats/interface/WeightedObject.h"

namespace flashgg {
    
    class JetBReg : public flashgg::Jet
    {
        
    public:
        JetBReg();
        JetBReg( const flashgg::Jet & );
        ~JetBReg();
        void setPuJetId( const edm::Ptr<reco::Vertex> vtx, const PileupJetIdentifier & );
        bool hasPuJetId( const edm::Ptr<reco::Vertex> vtx ) const;
        bool passesPuJetId( const edm::Ptr<reco::Vertex> vtx, PileupJetIdentifier::Id level = PileupJetIdentifier::kLoose ) const;
        void setSimpleRMS( float theRMS )  { simpleRMS_ = theRMS; }
        void setSimpleMVA( float theMVA )  { simpleMVA_ = theMVA; }
        float rms() const { return simpleRMS_; }
        float rms( const edm::Ptr<reco::Vertex> vtx ) const;
        float betaStar( const edm::Ptr<reco::Vertex> vtx ) const;
        bool passesPuJetId( const edm::Ptr<DiPhotonCandidate> dipho, PileupJetIdentifier::Id level = PileupJetIdentifier::kLoose )const;
        float rms( const edm::Ptr<DiPhotonCandidate> dipho ) const;
        float betaStar( const edm::Ptr<DiPhotonCandidate> dipho ) const;
        float puJetIdMVA() const { return simpleMVA_; }
        JetBReg *clone() const { return ( new JetBReg( *this ) ); }
        
        void  setQGL(const float qglikelihood=-99) {qglikelihood_ = qglikelihood;}
        float QGL () const {return qglikelihood_;}
        
        bool passesJetID( JetIDLevel level = Loose ) const; 

        const bool hasGenMatch() const { return (genJet() != 0); }
        
        const std::vector<float>& chEnergies() const { return chEnergies_; }
        const std::vector<float>& emEnergies() const { return emEnergies_; } 
        const std::vector<float>& neEnergies() const { return neEnergies_; } 
        const std::vector<float>& muEnergies() const { return muEnergies_; }

        void setChEnergies(std::vector<float> val) { chEnergies_ = val; }
        void setEmEnergies(std::vector<float> val) { emEnergies_ = val; } 
        void setNeEnergies(std::vector<float> val) { neEnergies_ = val; } 
        void setMuEnergies(std::vector<float> val) { muEnergies_ = val; }
        void setbRegMVA(float val) {bRegMVA_=val; }

        const float getBRegMVA() const { return bRegMVA_; }

    private:
        std::map<edm::Ptr<reco::Vertex>, MinimalPileupJetIdentifier> puJetId_;
        float qglikelihood_;
        float simpleRMS_; // simpler storage for PFCHS where this is not vertex-dependent
        float simpleMVA_;
        std::vector<float> chEnergies_, emEnergies_, neEnergies_, muEnergies_;
        float bRegMVA_;
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

