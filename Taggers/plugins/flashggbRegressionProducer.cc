#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/JetBReg.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/Math/interface/deltaR.h"


using namespace std;
using namespace edm;

namespace flashgg {

    class bRegressionProducer : public EDProducer
    {

    public:
        bRegressionProducer( const ParameterSet & );
        ~bRegressionProducer(){};
    private:
        void produce( Event &, const EventSetup & ) override;
        //        std::vector<edm::InputTag> inputTagJets_;
        edm::InputTag inputTagJets_;
        EDGetTokenT<View<flashgg::Jet> > jetToken_;
        edm::EDGetTokenT<double> rhoToken_;        

        unique_ptr<TMVA::Reader>bRegressionReader_;
        FileInPath bRegressionWeightfile_;

        //add vector of mva for eache jet

        //mva variables
        float Jet_pt ;
        float Jet_eta ;
        float rho ;
        float Jet_mt ;
        float Jet_leadTrackPt ;
        float Jet_leptonPtRel ;
        float Jet_leptonDeltaR ;
        float Jet_neHEF ;
        float Jet_neEmEF ;
        float Jet_vtxPt ;
        float Jet_vtxMass ;
        float Jet_vtx3dL ;
        float Jet_vtxNtrk ;
        float Jet_vtx3deL ;
        float Jet_energyRing_dR0_em_Jet_e ;
        float Jet_energyRing_dR1_em_Jet_e ;
        float Jet_energyRing_dR2_em_Jet_e ;
        float Jet_energyRing_dR3_em_Jet_e ;
        float Jet_energyRing_dR4_em_Jet_e ;
        float Jet_energyRing_dR0_neut_Jet_e ;
        float Jet_energyRing_dR1_neut_Jet_e ;
        float Jet_energyRing_dR2_neut_Jet_e ;
        float Jet_energyRing_dR3_neut_Jet_e ;
        float Jet_energyRing_dR4_neut_Jet_e ;
        float Jet_energyRing_dR0_ch_Jet_e ;
        float Jet_energyRing_dR1_ch_Jet_e ;
        float Jet_energyRing_dR2_ch_Jet_e ;
        float Jet_energyRing_dR3_ch_Jet_e ;
        float Jet_energyRing_dR4_ch_Jet_e ;
        float Jet_energyRing_dR0_mu_Jet_e ;
        float Jet_energyRing_dR1_mu_Jet_e ;
        float Jet_energyRing_dR2_mu_Jet_e ;
        float Jet_energyRing_dR3_mu_Jet_e ;
        float Jet_energyRing_dR4_mu_Jet_e ;
        float Jet_numDaughters_pt03 ;

        
    };


    bRegressionProducer::bRegressionProducer( const ParameterSet &iConfig ) :
        //     inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "JetTag" )) {
        inputTagJets_( iConfig.getParameter<edm::InputTag>( "JetTag" )) ,
        rhoToken_( consumes<double>(iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),
        bRegressionWeightfile_(iConfig.getParameter<edm::FileInPath>("bRegressionWeightfile"))
    {
        jetToken_= consumes<View<flashgg::Jet> >(inputTagJets_);
        bRegressionReader_.reset( new TMVA::Reader( "!Color:Silent" ) );

        //for variables for breg check this PR https://github.com/cms-analysis/flashgg/pull/968
        Jet_pt = 0.;
        Jet_eta = 0.;
        rho = 0.;
        Jet_mt = 0.;
        Jet_leadTrackPt = 0.;
        Jet_leptonPtRel = 0.;
        Jet_leptonDeltaR = 0.;
        Jet_neHEF = 0.;
        Jet_neEmEF = 0.;
        Jet_vtxPt = 0.;
        Jet_vtxMass = 0.;
        Jet_vtx3dL = 0.;
        Jet_vtxNtrk = 0.;
        Jet_vtx3deL = 0.;
        Jet_energyRing_dR0_em_Jet_e = 0.;
        Jet_energyRing_dR1_em_Jet_e = 0.;
        Jet_energyRing_dR2_em_Jet_e = 0.;
        Jet_energyRing_dR3_em_Jet_e = 0.;
        Jet_energyRing_dR4_em_Jet_e = 0.;
        Jet_energyRing_dR0_neut_Jet_e = 0.;
        Jet_energyRing_dR1_neut_Jet_e = 0.;
        Jet_energyRing_dR2_neut_Jet_e = 0.;
        Jet_energyRing_dR3_neut_Jet_e = 0.;
        Jet_energyRing_dR4_neut_Jet_e = 0.;
        Jet_energyRing_dR0_ch_Jet_e = 0.;
        Jet_energyRing_dR1_ch_Jet_e = 0.;
        Jet_energyRing_dR2_ch_Jet_e = 0.;
        Jet_energyRing_dR3_ch_Jet_e = 0.;
        Jet_energyRing_dR4_ch_Jet_e = 0.;
        Jet_energyRing_dR0_mu_Jet_e = 0.;
        Jet_energyRing_dR1_mu_Jet_e = 0.;
        Jet_energyRing_dR2_mu_Jet_e = 0.;
        Jet_energyRing_dR3_mu_Jet_e = 0.;
        Jet_energyRing_dR4_mu_Jet_e = 0.;
        Jet_numDaughters_pt03 = 0;



        bRegressionReader_->AddVariable( "Jet_pt", &Jet_pt );
        bRegressionReader_->AddVariable( "Jet_eta", &Jet_eta );
        bRegressionReader_->AddVariable( "rho", &rho );
        bRegressionReader_->AddVariable( "Jet_mt", &Jet_mt );
        bRegressionReader_->AddVariable( "Jet_leadTrackPt", &Jet_leadTrackPt );
        bRegressionReader_->AddVariable( "Jet_leptonPtRel", &Jet_leptonPtRel );
        bRegressionReader_->AddVariable( "Jet_leptonDeltaR", &Jet_leptonDeltaR );
        bRegressionReader_->AddVariable( "Jet_neHEF", &Jet_neHEF );
        bRegressionReader_->AddVariable( "Jet_neEmEF", &Jet_neEmEF );
        bRegressionReader_->AddVariable( "Jet_vtxPt", &Jet_vtxPt );
        bRegressionReader_->AddVariable( "Jet_vtxMass", &Jet_vtxMass );
        bRegressionReader_->AddVariable( "Jet_vtx3dL", &Jet_vtx3dL );
        bRegressionReader_->AddVariable( "Jet_vtxNtrk", &Jet_vtxNtrk );
        bRegressionReader_->AddVariable( "Jet_vtx3deL", &Jet_vtx3deL );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR0_em_Jet_e", &Jet_energyRing_dR0_em_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR1_em_Jet_e", &Jet_energyRing_dR1_em_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR2_em_Jet_e", &Jet_energyRing_dR2_em_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR3_em_Jet_e", &Jet_energyRing_dR3_em_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR4_em_Jet_e", &Jet_energyRing_dR4_em_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR0_neut_Jet_e", &Jet_energyRing_dR0_neut_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR1_neut_Jet_e", &Jet_energyRing_dR1_neut_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR2_neut_Jet_e", &Jet_energyRing_dR2_neut_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR3_neut_Jet_e", &Jet_energyRing_dR3_neut_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR4_neut_Jet_e", &Jet_energyRing_dR4_neut_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR0_ch_Jet_e", &Jet_energyRing_dR0_ch_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR1_ch_Jet_e", &Jet_energyRing_dR1_ch_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR2_ch_Jet_e", &Jet_energyRing_dR2_ch_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR3_ch_Jet_e", &Jet_energyRing_dR3_ch_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR4_ch_Jet_e", &Jet_energyRing_dR4_ch_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR0_mu_Jet_e", &Jet_energyRing_dR0_mu_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR1_mu_Jet_e", &Jet_energyRing_dR1_mu_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR2_mu_Jet_e", &Jet_energyRing_dR2_mu_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR3_mu_Jet_e", &Jet_energyRing_dR3_mu_Jet_e );
        bRegressionReader_->AddVariable( "Jet_energyRing_dR4_mu_Jet_e", &Jet_energyRing_dR4_mu_Jet_e );
        bRegressionReader_->AddVariable( "Jet_numDaughters_pt03", &Jet_numDaughters_pt03 );

        bRegressionReader_->BookMVA( "BDT" , bRegressionWeightfile_.fullPath() );


        produces<vector<flashgg::JetBReg> > ();
        
    }



    void bRegressionProducer::produce( Event &evt, const EventSetup & )
    {
        
        // input jets
        Handle<View<flashgg::Jet> > jets;
        evt.getByToken( jetToken_, jets );//just to try get the first one
        std::cout<<"jet size: "<<jets->size()<<std::endl;
        
        unique_ptr<vector<flashgg::JetBReg> > jetColl( new vector<flashgg::JetBReg> );
        for( unsigned int i = 0 ; i < jets->size() ; i++ ) {

            
            Ptr<flashgg::Jet> pjet = jets->ptrAt( i );
            //            flashgg::Jet fjet = flashgg::Jet( *pjet );
            flashgg::JetBReg fjet = flashgg::JetBReg( *pjet );

            //variables needed for regression
            Jet_pt = fjet.pt();
            Jet_eta = fjet.eta() ;
            Jet_leadTrackPt = fjet.userFloat("leadTrackPt");
            edm::Handle<double> rhoHandle;
            evt.getByToken( rhoToken_, rhoHandle );
            const double rhoFixedGrd = *( rhoHandle.product() );
            rho = rhoFixedGrd;
            Jet_mt = sqrt(fjet.energy()*fjet.energy()-fjet.pz()*fjet.pz());
            Jet_leptonPtRel = fjet.userFloat("softLepPtRel");
            Jet_leptonDeltaR = fjet.userFloat("softLepDr");
            Jet_neHEF = fjet.neutralHadronEnergyFraction();
            Jet_neEmEF = fjet.chargedEmEnergyFraction();
            if(fjet.userFloat("nSecVertices")>0){
                float vertexX=fjet.userFloat("vtxPosX")-fjet.userFloat("vtxPx");//check if it's correct
                float vertexY=fjet.userFloat("vtxPosY")-fjet.userFloat("vtxPy");                
                Jet_vtxPt = sqrt(vertexX*vertexX+vertexY*vertexY);
                Jet_vtxMass = fjet.userFloat("vtxMass");
                Jet_vtx3dL = fjet.userFloat("vtx3DVal");
                Jet_vtxNtrk = fjet.userFloat("vtxNTracks");
                Jet_vtx3deL = fjet.userFloat("vtx3DSig");
            }
            if(fjet.emEnergies().size()>0){
                Jet_energyRing_dR0_em_Jet_e = fjet.emEnergies()[0]/fjet.energy();//remember to divide by jet energy
                Jet_energyRing_dR1_em_Jet_e = fjet.emEnergies()[1]/fjet.energy();
                Jet_energyRing_dR2_em_Jet_e = fjet.emEnergies()[2]/fjet.energy();
                Jet_energyRing_dR3_em_Jet_e = fjet.emEnergies()[3]/fjet.energy();
                Jet_energyRing_dR4_em_Jet_e = fjet.emEnergies()[4]/fjet.energy();
            }
            if(fjet.neEnergies().size()>0){
                Jet_energyRing_dR0_neut_Jet_e = fjet.neEnergies()[0]/fjet.energy();
                Jet_energyRing_dR1_neut_Jet_e = fjet.neEnergies()[1]/fjet.energy();
                Jet_energyRing_dR2_neut_Jet_e = fjet.neEnergies()[2]/fjet.energy();
                Jet_energyRing_dR3_neut_Jet_e = fjet.neEnergies()[3]/fjet.energy();
                Jet_energyRing_dR4_neut_Jet_e = fjet.neEnergies()[4]/fjet.energy();
            }
            if(fjet.chEnergies().size()>0){
                Jet_energyRing_dR0_ch_Jet_e = fjet.chEnergies()[0]/fjet.energy();
                Jet_energyRing_dR1_ch_Jet_e = fjet.chEnergies()[1]/fjet.energy();
                Jet_energyRing_dR2_ch_Jet_e = fjet.chEnergies()[2]/fjet.energy();
                Jet_energyRing_dR3_ch_Jet_e = fjet.chEnergies()[3]/fjet.energy();
                Jet_energyRing_dR4_ch_Jet_e = fjet.chEnergies()[4]/fjet.energy();
            }
            if(fjet.muEnergies().size()>0){
                Jet_energyRing_dR0_mu_Jet_e = fjet.muEnergies()[0]/fjet.energy();
                Jet_energyRing_dR1_mu_Jet_e = fjet.muEnergies()[1]/fjet.energy();
                Jet_energyRing_dR2_mu_Jet_e = fjet.muEnergies()[2]/fjet.energy();
                Jet_energyRing_dR3_mu_Jet_e = fjet.muEnergies()[3]/fjet.energy();
                Jet_energyRing_dR4_mu_Jet_e = fjet.muEnergies()[4]/fjet.energy();
            }
            Jet_numDaughters_pt03 = fjet.userInt("numDaug03");
            
            float bRegMVA=-999;
            

            int debug=0;
            if(debug){
                cout<<"Jet_pt :"<<Jet_pt <<endl;
                cout<<"Jet_eta :"<<Jet_eta <<endl;
                cout<<"rho :"<<rho <<endl;
                cout<<"Jet_mt :"<<Jet_mt <<endl;
                cout<<"Jet_leadTrackPt :"<<Jet_leadTrackPt <<endl;
                cout<<"Jet_leptonPtRel :"<<Jet_leptonPtRel <<endl;
                cout<<"Jet_leptonDeltaR :"<<Jet_leptonDeltaR <<endl;
                cout<<"Jet_neHEF :"<<Jet_neHEF <<endl;
                cout<<"Jet_neEmEF :"<<Jet_neEmEF <<endl;
                cout<<"Jet_vtxPt :"<<Jet_vtxPt <<endl;
                cout<<"Jet_vtxMass :"<<Jet_vtxMass <<endl;
                cout<<"Jet_vtx3dL :"<<Jet_vtx3dL <<endl;
                cout<<"Jet_vtxNtrk :"<<Jet_vtxNtrk <<endl;
                cout<<"Jet_vtx3deL :"<<Jet_vtx3deL <<endl;
                cout<<"Jet_energyRing_dR0_em_Jet_e :"<<Jet_energyRing_dR0_em_Jet_e <<endl;
                cout<<"Jet_energyRing_dR1_em_Jet_e :"<<Jet_energyRing_dR1_em_Jet_e <<endl;
                cout<<"Jet_energyRing_dR2_em_Jet_e :"<<Jet_energyRing_dR2_em_Jet_e <<endl;
                cout<<"Jet_energyRing_dR3_em_Jet_e :"<<Jet_energyRing_dR3_em_Jet_e <<endl;
                cout<<"Jet_energyRing_dR4_em_Jet_e :"<<Jet_energyRing_dR4_em_Jet_e <<endl;
                cout<<"Jet_energyRing_dR0_neut_Jet_e :"<<Jet_energyRing_dR0_neut_Jet_e <<endl;
                cout<<"Jet_energyRing_dR1_neut_Jet_e :"<<Jet_energyRing_dR1_neut_Jet_e <<endl;
                cout<<"Jet_energyRing_dR2_neut_Jet_e :"<<Jet_energyRing_dR2_neut_Jet_e <<endl;
                cout<<"Jet_energyRing_dR3_neut_Jet_e :"<<Jet_energyRing_dR3_neut_Jet_e <<endl;
                cout<<"Jet_energyRing_dR4_neut_Jet_e :"<<Jet_energyRing_dR4_neut_Jet_e <<endl;
                cout<<"Jet_energyRing_dR0_ch_Jet_e :"<<Jet_energyRing_dR0_ch_Jet_e <<endl;
                cout<<"Jet_energyRing_dR1_ch_Jet_e :"<<Jet_energyRing_dR1_ch_Jet_e <<endl;
                cout<<"Jet_energyRing_dR2_ch_Jet_e :"<<Jet_energyRing_dR2_ch_Jet_e <<endl;
                cout<<"Jet_energyRing_dR3_ch_Jet_e :"<<Jet_energyRing_dR3_ch_Jet_e <<endl;
                cout<<"Jet_energyRing_dR4_ch_Jet_e :"<<Jet_energyRing_dR4_ch_Jet_e <<endl;
                cout<<"Jet_energyRing_dR0_mu_Jet_e :"<<Jet_energyRing_dR0_mu_Jet_e <<endl;
                cout<<"Jet_energyRing_dR1_mu_Jet_e :"<<Jet_energyRing_dR1_mu_Jet_e <<endl;
                cout<<"Jet_energyRing_dR2_mu_Jet_e :"<<Jet_energyRing_dR2_mu_Jet_e <<endl;
                cout<<"Jet_energyRing_dR3_mu_Jet_e :"<<Jet_energyRing_dR3_mu_Jet_e <<endl;
                cout<<"Jet_energyRing_dR4_mu_Jet_e :"<<Jet_energyRing_dR4_mu_Jet_e <<endl;
                cout<<"Jet_numDaughters_pt03 :"<<Jet_numDaughters_pt03 <<endl;



            }

            bRegMVA = bRegressionReader_->EvaluateMVA("BDT");
            std::cout<<bRegMVA<<std::endl;

            std::cout<<"Jet index: "<<i<<" Pt:"<<fjet.pt()<<std::endl;

            //..... gen jets info                                                                                   
            int cflav = 0; //~correct flavour definition
            if ( !evt.isRealData() ) {
                int hflav = fjet.hadronFlavour();//4 if c, 5 if b, 0 if light jets
                int pflav = fjet.partonFlavour();

                if( hflav != 0 ) {
                    cflav = hflav;
                } else { //not a heavy jet                                              
                    cflav = std::abs(pflav) == 4 || std::abs(pflav) == 5 ? 0 : pflav;
                }
                std::cout<<cflav<<std::endl;
                if (cflav != 5) continue;//i want only bjets
            }
            if (fjet.pt()<15. || fabs(fjet.eta())>2.5) continue;
            std::cout<<"found a b-jet of pt"<<fjet.pt()<<" eta:"<<fjet.eta()<<std::endl;

            //            fjet.addUserFloat("bRegMVA", bRegMVA);
            fjet.setbRegMVA(bRegMVA);

            std::cout<<fjet.getBRegMVA()<<std::endl;
            jetColl->push_back( fjet );

        }
        evt.put( std::move( jetColl ) );
    }
    


}

typedef flashgg::bRegressionProducer flashggbRegressionProducer;
DEFINE_FWK_MODULE( flashggbRegressionProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4


