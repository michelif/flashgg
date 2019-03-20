#ifndef flashgg_MVAComputer_h
#define flashgg_MVAComputer_h

#include <tuple>
#include <vector>
#include <string>

#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "flashgg/Taggers/interface/StringHelpers.h"

#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "TMVA/Reader.h"
#include "flashgg/Taggers/interface/GlobalVariablesDumper.h"

#include "XGBoostCMSSW/XGBoostInterface/interface/XGBComputer.h"

namespace flashgg {

    template<class ObjectT, class FunctorT = StringObjectFunction<ObjectT, true> >
    class MVAComputer
    {
    public:
        typedef ObjectT object_type;
        typedef FunctorT functor_type;
        
        MVAComputer( const edm::ParameterSet &cfg, GlobalVariablesComputer *global = 0 );
        ~MVAComputer();

        float operator()( const object_type &obj ) const;
        std::vector<float> predict_prob( const object_type &obj ) const;
        void getMVAVar();
        XGBComputer xgbComputer_;
        XGBComputer::mva_variables  xgbVars_;
        void setupXGB();
        std::vector<float> predict_probXGB();
        
    private:
        void bookMVA() const;
        
        mutable TMVA::Reader *reader_;
        GlobalVariablesComputer *global_;

        bool regression_;
        bool multiclass_;
        std::string classifier_, weights_;
        std::vector<std::tuple<std::string, int> > variables_;
        mutable std::vector<float> values_;
        std::vector<functor_type> functors_;
    };

    template<class F, class O>
    MVAComputer<F, O>::MVAComputer( const edm::ParameterSet &cfg, GlobalVariablesComputer *global ) :
        reader_( 0 ),
        global_( global ),
        regression_( cfg.exists("regression") ? cfg.getParameter<bool>("regression") : false ),
        multiclass_( cfg.exists("multiclass") ? cfg.getParameter<bool>("multiclass") : false ),
        classifier_( cfg.getParameter<std::string>( "classifier" ) )
    {
        using namespace std;

        weights_    = cfg.getParameter<edm::FileInPath>( "weights" ).fullPath();

        auto variables = cfg.getParameter<vector<edm::ParameterSet> >( "variables" );
        for( auto &var : variables ) {
            auto expr = var.getParameter<string>( "expr" );
            auto name = var.getUntrackedParameter<string>( "name", expr );
            auto pos = expr.find( "global." );
            if( pos == 0 ) {
                assert( global != 0 );
                variables_.push_back( std::make_tuple( name + "::" + expr.substr( 7 ), -1 ) );
            } else {
                functors_.push_back( functor_type( expr ) );
                variables_.push_back( std::make_tuple( name, functors_.size() - 1 ) );
            }
        }
    }

    template<class F, class O>
    MVAComputer<F, O>::~MVAComputer()
    {
        if( reader_ ) {
            delete reader_;
        }
    }

    template<class F, class O>
    void MVAComputer<F, O>::bookMVA() const
    {
        if( reader_ ) { return; }
        reader_ = new TMVA::Reader( "!Color" );
        values_.resize( functors_.size(), 0. );
        for( auto &var : variables_ ) {
            auto &name = std::get<0>( var );
            auto ivar = std::get<1>( var );
            if( ivar >= 0 ) {
                reader_->AddVariable( name, &values_[ivar] );
            } else {
                assert( global_ != 0 );
                auto pos = name.find( "::" );
                auto vname = name.substr( 0, pos );
                auto gname = name.substr( pos + 2 );
                reader_->AddVariable( vname, global_->addressOf( gname ) );
            }
        }
        reader_->BookMVA( classifier_, weights_ );

    }

    template<class F, class O>
    float MVAComputer<F, O>::operator()( const object_type &obj ) const
    {
        assert (multiclass_==false);
        if( ! reader_ ) { bookMVA(); }
        for( size_t ivar = 0; ivar < functors_.size(); ++ivar ) {
            values_[ivar] = functors_[ivar]( obj );
        }
        return ( regression_ ? reader_->EvaluateRegression(0, classifier_.c_str() ) : reader_->EvaluateMVA( classifier_.c_str() ) );
    }
    template<class F, class O>
    std::vector<float> MVAComputer<F, O>::predict_prob( const object_type &obj ) const
    {
        assert (multiclass_==true);
        if( ! reader_ ) { bookMVA(); }
        for( size_t ivar = 0; ivar < functors_.size(); ++ivar ) {
            values_[ivar] = functors_[ivar]( obj );
        }
        return reader_->EvaluateMulticlass(classifier_.c_str());
    }

    template<class F, class O>
    void MVAComputer<F, O>::getMVAVar() 
    {
        for( size_t ivar = 0; ivar < values_.size(); ++ivar ) {
            xgbVars_.push_back(std::make_tuple(get<0>(variables_[ivar]),values_[ivar]));
            //            std::cout<<vars[0]<<std::endl;
        }
    }

    template<class F, class O>
    void MVAComputer<F, O>::setupXGB()
    {
        std::string modelfile="Taggers/data/HHTagger/training_with_20190201_test_2.pkl";
        getMVAVar();
        xgbComputer_=XGBComputer(&xgbVars_,modelfile);
    }

    template<class F, class O>
    std::vector<float> MVAComputer<F, O>::predict_probXGB() 
    {
        return xgbComputer_();
    }

}


#endif  // flashgg_MVAComputer_h
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

