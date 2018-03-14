#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/UtilAlgos/interface/EDAnalyzerWrapper.h"

#include "flashgg/Taggers/interface/bRegressionDumper.h"

typedef edm::AnalyzerWrapper<flashgg::bRegressionDumper> bRegressionDumper;
typedef edm::AnalyzerWrapper<flashgg::CutBasedbRegressionDumper> CutBasedbRegressionDumper;

DEFINE_FWK_MODULE( bRegressionDumper );
DEFINE_FWK_MODULE( CutBasedbRegressionDumper );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
