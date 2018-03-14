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
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/Math/interface/deltaR.h"


#include <iostream>
#include <string>
#include <vector>

#include "DNN/Tensorflow/interface/Graph.h"
#include "DNN/Tensorflow/interface/Tensor.h"

#define doBDTAnalysis 0 
#define debug 0

using namespace std;
using namespace edm;

namespace flashgg {

    class bRegressionProducer : public EDProducer
    {

    public:
        bRegressionProducer( const ParameterSet & );
        ~bRegressionProducer(){};
        void InitJet();
        void SetNNVectorVar();
        std::vector<float> EvaluateNN();
    private:
        void produce( Event &, const EventSetup & ) override;
        //        std::vector<edm::InputTag> inputTagJets_;
        edm::InputTag inputTagJets_;
        EDGetTokenT<View<flashgg::Jet> > jetToken_;
        edm::EDGetTokenT<double> rhoToken_;        

        unique_ptr<TMVA::Reader>bRegressionReader_;
        FileInPath bRegressionWeightfile_;
        dnn::tf::Graph NNgraph_;
        std::vector<float> NNvectorVar_; 
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


        NNgraph_ = *(new dnn::tf::Graph("/afs/cern.ch/work/m/micheli/CMSSW_8_0_28/src/flashgg/MetaData/data/DNN_models/model-09"); //FIXME make this configurable
        //        NNgraph_ = *(new dnn::tf::Graph("/afs/cern.ch/work/m/micheli/CMSSW_8_0_28/src/flashgg/MetaData/data/DNN_models/model-09",dnn::LogLevel::ALL)); //FIXME make this configurable

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

        if (doBDTAnalysis){ 
        bRegressionReader_.reset( new TMVA::Reader( "!Color:Silent" ) );
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

        }
        //FIXME        produces<vector<flashgg::JetBReg> > ();
        produces<vector<flashgg::Jet> > ();
    }



    void bRegressionProducer::produce( Event &evt, const EventSetup & )
    {
        InitJet();
        // input jets
        Handle<View<flashgg::Jet> > jets;
        evt.getByToken( jetToken_, jets );//just to try get the first one
       unique_ptr<vector<flashgg::Jet> > jetColl( new vector<flashgg::Jet> );
        for( unsigned int i = 0 ; i < jets->size() ; i++ ) {

            
            Ptr<flashgg::Jet> pjet = jets->ptrAt( i );
            flashgg::Jet fjet = flashgg::Jet( *pjet );

            if (fjet.pt()<15. || fabs(fjet.eta())>2.5) continue;


            //variables needed for regression
            Jet_pt = fjet.pt();
            Jet_eta = fjet.eta() ;
            Jet_leadTrackPt = fjet.userFloat("leadTrackPt");
            edm::Handle<double> rhoHandle;
            evt.getByToken( rhoToken_, rhoHandle );
            const double rhoFixedGrd = *( rhoHandle.product() );
            rho = rhoFixedGrd;
            Jet_mt = sqrt(fjet.energy()*fjet.energy()-fjet.pz()*fjet.pz());//seems correct but check again

            //this max probably not needed, it's just heppy
            Jet_leptonPtRel = std::max(float(0.),fjet.userFloat("softLepPtRel"));
            Jet_leptonDeltaR = std::max(float(0.),fjet.userFloat("softLepDr"));
            Jet_neHEF = fjet.neutralHadronEnergyFraction();
            Jet_neEmEF = fjet.neutralEmEnergyFraction();
            if(fjet.userFloat("nSecVertices")>0){
//                float vertexX=fjet.userFloat("vtxPosX")-fjet.userFloat("vtxPx");//check if it's correct
//                float vertexY=fjet.userFloat("vtxPosY")-fjet.userFloat("vtxPy");                
//                Jet_vtxPt = sqrt(vertexX*vertexX+vertexY*vertexY);
                Jet_vtxPt=sqrt(fjet.userFloat("vtxPx")*fjet.userFloat("vtxPx")+fjet.userFloat("vtxPy")*fjet.userFloat("vtxPy"));
                Jet_vtxMass = std::max(float(0.),fjet.userFloat("vtxMass"));
                Jet_vtx3dL = std::max(float(0.),fjet.userFloat("vtx3DVal"));
                Jet_vtxNtrk = std::max(float(0.),fjet.userFloat("vtxNTracks"));
                Jet_vtx3deL = std::max(float(0.),fjet.userFloat("vtx3DSig"));
            }
            Jet_energyRing_dR0_em_Jet_e = fjet.emEnergies()[0]/fjet.energy();//remember to divide by jet energy
            Jet_energyRing_dR1_em_Jet_e = fjet.emEnergies()[1]/fjet.energy();
            Jet_energyRing_dR2_em_Jet_e = fjet.emEnergies()[2]/fjet.energy();
            Jet_energyRing_dR3_em_Jet_e = fjet.emEnergies()[3]/fjet.energy();
            Jet_energyRing_dR4_em_Jet_e = fjet.emEnergies()[4]/fjet.energy();
            Jet_energyRing_dR0_neut_Jet_e = fjet.neEnergies()[0]/fjet.energy();
            Jet_energyRing_dR1_neut_Jet_e = fjet.neEnergies()[1]/fjet.energy();
            Jet_energyRing_dR2_neut_Jet_e = fjet.neEnergies()[2]/fjet.energy();
            Jet_energyRing_dR3_neut_Jet_e = fjet.neEnergies()[3]/fjet.energy();
            Jet_energyRing_dR4_neut_Jet_e = fjet.neEnergies()[4]/fjet.energy();
            Jet_energyRing_dR0_ch_Jet_e = fjet.chEnergies()[0]/fjet.energy();
            Jet_energyRing_dR1_ch_Jet_e = fjet.chEnergies()[1]/fjet.energy();
            Jet_energyRing_dR2_ch_Jet_e = fjet.chEnergies()[2]/fjet.energy();
            Jet_energyRing_dR3_ch_Jet_e = fjet.chEnergies()[3]/fjet.energy();
            Jet_energyRing_dR4_ch_Jet_e = fjet.chEnergies()[4]/fjet.energy();
            Jet_energyRing_dR0_mu_Jet_e = fjet.muEnergies()[0]/fjet.energy();
            Jet_energyRing_dR1_mu_Jet_e = fjet.muEnergies()[1]/fjet.energy();
            Jet_energyRing_dR2_mu_Jet_e = fjet.muEnergies()[2]/fjet.energy();
            Jet_energyRing_dR3_mu_Jet_e = fjet.muEnergies()[3]/fjet.energy();
            Jet_energyRing_dR4_mu_Jet_e = fjet.muEnergies()[4]/fjet.energy();

            Jet_numDaughters_pt03 = fjet.userInt("numDaug03");
            
            float bRegMVA=-999;
            std::vector<float> bRegNN(3,-999);
            

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
                //                std::cout<<cflav<<std::endl;
                if (cflav != 5) continue;//i want only bjets
            }
            //            std::cout<<"found a b-jet of pt"<<fjet.pt()<<" eta:"<<fjet.eta()<<bRegMVA<<std::endl;


            if(doBDTAnalysis) bRegMVA = bRegressionReader_->EvaluateMVA("BDT");
            SetNNVectorVar();
            bRegNN = EvaluateNN();
            NNvectorVar_.clear();

            //FIXME read through file, config is here /afs/cern.ch/user/n/nchernya/public/100M_2018-03-01_job23_rawJetsJECtarget/config.json
            float y_mean= 1.0454729795455933;
            float y_std = 0.3162831664085388;

            fjet.addUserFloat("bRegMVA", bRegMVA);
            fjet.addUserFloat("bRegNNCorr", bRegNN[0]*y_std+y_mean);
            fjet.addUserFloat("bRegNNResolution",0.5*(bRegNN[2]-bRegNN[1])*y_std);
            fjet.addUserFloat("energyRing_dR0_em_Jet_e", Jet_energyRing_dR0_em_Jet_e) ;//remember to divide by jet energy
            fjet.addUserFloat("energyRing_dR1_em_Jet_e", Jet_energyRing_dR1_em_Jet_e) ;
            fjet.addUserFloat("energyRing_dR2_em_Jet_e", Jet_energyRing_dR2_em_Jet_e) ;
            fjet.addUserFloat("energyRing_dR3_em_Jet_e", Jet_energyRing_dR3_em_Jet_e) ;
            fjet.addUserFloat("energyRing_dR4_em_Jet_e", Jet_energyRing_dR4_em_Jet_e) ;
            fjet.addUserFloat("energyRing_dR0_neut_Jet_e", Jet_energyRing_dR0_neut_Jet_e) ;
            fjet.addUserFloat("energyRing_dR1_neut_Jet_e", Jet_energyRing_dR1_neut_Jet_e) ;
            fjet.addUserFloat("energyRing_dR2_neut_Jet_e", Jet_energyRing_dR2_neut_Jet_e) ;
            fjet.addUserFloat("energyRing_dR3_neut_Jet_e", Jet_energyRing_dR3_neut_Jet_e) ;
            fjet.addUserFloat("energyRing_dR4_neut_Jet_e", Jet_energyRing_dR4_neut_Jet_e) ;
            fjet.addUserFloat("energyRing_dR0_ch_Jet_e", Jet_energyRing_dR0_ch_Jet_e) ;
            fjet.addUserFloat("energyRing_dR1_ch_Jet_e", Jet_energyRing_dR1_ch_Jet_e) ;
            fjet.addUserFloat("energyRing_dR2_ch_Jet_e", Jet_energyRing_dR2_ch_Jet_e) ;
            fjet.addUserFloat("energyRing_dR3_ch_Jet_e", Jet_energyRing_dR3_ch_Jet_e) ;
            fjet.addUserFloat("energyRing_dR4_ch_Jet_e", Jet_energyRing_dR4_ch_Jet_e) ;
            fjet.addUserFloat("energyRing_dR0_mu_Jet_e", Jet_energyRing_dR0_mu_Jet_e) ;
            fjet.addUserFloat("energyRing_dR1_mu_Jet_e", Jet_energyRing_dR1_mu_Jet_e) ;
            fjet.addUserFloat("energyRing_dR2_mu_Jet_e", Jet_energyRing_dR2_mu_Jet_e) ;
            fjet.addUserFloat("energyRing_dR3_mu_Jet_e", Jet_energyRing_dR3_mu_Jet_e) ;
            fjet.addUserFloat("energyRing_dR4_mu_Jet_e", Jet_energyRing_dR4_mu_Jet_e ) ;
            fjet.addUserFloat("numDaughters_pt03", fjet.userInt("numDaug03"));


            jetColl->push_back( fjet );

            

        }
        evt.put( std::move( jetColl ) );
    }
    
    void bRegressionProducer::InitJet(){
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

    }//end InitJet

    void bRegressionProducer::SetNNVectorVar(){

        NNvectorVar_.push_back(Jet_pt) ;//0
        NNvectorVar_.push_back(Jet_eta) ;
        NNvectorVar_.push_back(rho) ;
        NNvectorVar_.push_back(Jet_mt) ;
        NNvectorVar_.push_back(Jet_leadTrackPt) ;
        NNvectorVar_.push_back(Jet_leptonPtRel) ;//5
        NNvectorVar_.push_back(Jet_leptonDeltaR) ;
        NNvectorVar_.push_back(Jet_neHEF) ;
        NNvectorVar_.push_back(Jet_neEmEF) ;
        NNvectorVar_.push_back(Jet_vtxPt) ;
        NNvectorVar_.push_back(Jet_vtxMass) ;//10
        NNvectorVar_.push_back(Jet_vtx3dL) ;
        NNvectorVar_.push_back(Jet_vtxNtrk) ;
        NNvectorVar_.push_back(Jet_vtx3deL) ;
        NNvectorVar_.push_back(Jet_numDaughters_pt03) ;//this variable has changed order, in bdt it was last, check why
        NNvectorVar_.push_back(Jet_energyRing_dR0_em_Jet_e) ;//15
        NNvectorVar_.push_back(Jet_energyRing_dR1_em_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR2_em_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR3_em_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR4_em_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR0_neut_Jet_e) ;//20
        NNvectorVar_.push_back(Jet_energyRing_dR1_neut_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR2_neut_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR3_neut_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR4_neut_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR0_ch_Jet_e) ;//25
        NNvectorVar_.push_back(Jet_energyRing_dR1_ch_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR2_ch_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR3_ch_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR4_ch_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR0_mu_Jet_e) ;//30
        NNvectorVar_.push_back(Jet_energyRing_dR1_mu_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR2_mu_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR3_mu_Jet_e) ;
        NNvectorVar_.push_back(Jet_energyRing_dR4_mu_Jet_e) ;

    }
    
    std::vector<float> bRegressionProducer::EvaluateNN(){
        dnn::tf::Shape xShape[] = { 1, 35 };

        dnn::tf::Tensor* x = NNgraph_.defineInput(new dnn::tf::Tensor("ffwd_inp:0", 2, xShape));
        dnn::tf::Tensor* y = NNgraph_.defineOutput(new dnn::tf::Tensor("ffwd_out/BiasAdd:0"));
        for (int i = 0; i < x->getShape(1); i++){
            //            std::cout<<"i:"<<i<<" x:"<<NNvectorVar_[i]<<std::endl;
            x->setValue<float>(0, i, NNvectorVar_[i]);
        }
        NNgraph_.eval();
        std::vector<float> correction(3);//3 outputs, first value is mean and then other 2 quantiles
        correction[0] = y->getValue<float>(0, 0);
        correction[1] = y->getValue<float>(0, 1);            
        correction[2] = y->getValue<float>(0, 2);            
        return correction;
    }//end EvaluateNN
    
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


