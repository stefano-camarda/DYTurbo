#ifndef DYTurbo_unittest_CXX
#define DYTurbo_unittest_CXX
/**
 * @file DYTurbo_unittest.cxx
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-16
 */


#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "config.h"
#include "src/dyturbo.h"
#include "src/settings.h"
#include "src/cubacall.h"
#include "src/clock_real.h"


#include "src/resintegr.h"
#include "src/ctintegr.h"
#include "src/finintegr.h"
#include "src/bornintegr.h"
#include "mcfm/mcfm_interface.h"

#include "histo/HistoHandler.h"

#include <iostream>

::testing::AssertionResult ApproxDoubles(double exp, double equalto) {
    double significance = exp - equalto;
    if (significance==0) return ::testing::AssertionSuccess();
    significance/=exp;
    significance=fabs(significance);
    char exp_hex[40];
    char eqt_hex[40];
    sprintf(exp_hex, "%a", exp);
    sprintf(eqt_hex, "%a", equalto);
    if (fabs(significance) > 5e-12 ) {
        return ::testing::AssertionFailure() << endl <<
            "Expected       : " << exp     << "    " << exp_hex << std::endl <<
            "To be equal to : " << equalto << "    " << eqt_hex << std::endl <<
            "Significance   : " << significance << std::endl;
    }
    return ::testing::AssertionSuccess();
}

TEST(DYTurbo,Initialization){
    int argc = 4;
    char *argv[] =  {
        ((char *) "DYTurbo_unittest"),
        ((char *) "../input/test.in"),
        ((char *) "--qtbins"),
        ((char *) "10,0,100"),
    };
    DYTurbo::Init(argc,argv);
    //test histogram binning inherits from CLI
    VecDbl mbins;
    bins.GetBins("qt", mbins);
    ASSERT_EQ(11, mbins.size());
}

TEST(DYTurbo,TermIteration){
    DYTurbo::WarmUp();
    std::vector <void (*)(std::vector<double>&, double&)> term_funs;
    term_funs.push_back(bornintegrMC6d);
    if (opts.makecuts) term_funs.push_back(vjlointegr7d);
    else term_funs.push_back(vjintegr3d);
    ASSERT_EQ(term_funs.size(), DYTurbo::ActiveTerms.size() )
        << "Incorrect number of terms. Check your input file." ;
    for ( DYTurbo::TermIterator it_term; !it_term.IsEnd(); ++it_term) {
        ASSERT_EQ( term_funs[it_term.icurrent], (*it_term).integrate )
             << "Wrong term in active list index: " << it_term.icurrent << endl
             << "name: " << (*it_term).name << endl
             << "decr: " << (*it_term).description << endl
             << "Check your inputfile." << endl;
    }
}

TEST(DYTurbo,BoundIteration){
    size_t expected;
    size_t counter;
    // simple binner
    bins.qtbins = {0.,30.};
    bins.ybins  = {-3.,3.};
    bins.mbins  = {50.,500. };
    expected = (bins.qtbins.size() - 1 ) * (bins.ybins.size() - 1 ) * (bins.mbins.size() - 1 ) ;
    counter = 0 ;
    DYTurbo::WarmUp();
    for ( DYTurbo::BoundIterator it_bound; !it_bound.IsEnd(); ++it_bound) {
        counter ++;
    }
    ASSERT_EQ ( expected , counter  ) << " Incorrect number of bins from boundaries.";
    opts.force_binsampling=true;
    //if (opts.makecuts) return; // in case of leton cut it automatically changes to Vegas
    // every step
    bins.qtbins = {0., 10., 20., 30.};
    bins.ybins  = {-3., 0. ,3.};
    expected = (bins.qtbins.size() - 1 ) * (bins.ybins.size() - 1 ) * (bins.mbins.size() - 1 ) ;
    counter = 0 ;
    DYTurbo::WarmUp();
    for ( DYTurbo::BoundIterator it_bound; !it_bound.IsEnd(); ++it_bound) {
        counter ++;
    }
    ASSERT_EQ ( expected , counter  ) << " Incorrect number of bins from boundaries.";
    opts.force_binsampling=false;
}


/**
 * @brief Save results to file.
 *
 * We have to save results as binary values. For this you need first
 * calibration run to create files. This happens when you set
 * `doSaveResults=true`
 *
 * For actual testing you have to keep it `false`. Then values will be compared
 * to calibration run.
 */
bool doSaveResults=false;
FILE *F = 0;
void writefloat(double v) {
  fwrite((void*)(&v), sizeof(v), 1, F);
}
double readfloat() {
  double v=666.;
  if (!doSaveResults) {
      int i = fread((void*)(&v), sizeof(v), 1, F);
  }
  return v;
}


TEST(DYTurbo,MainLoop){
    DYTurbo::isDryRun = false;
    DYTurbo::WarmUp();
    DYTurbo::PrintTable::Header();
    F= fopen( opts.makecuts ? "integration_cuts.res" : "integration.res", (!DYTurbo::isDryRun && doSaveResults) ? "w" : "r");
    for ( DYTurbo::BoundIterator bounds; !bounds.IsEnd(); ++bounds) {
        DYTurbo::SetBounds(bounds);
        DYTurbo::PrintTable::Bounds();
        for (DYTurbo::TermIterator term; !term.IsEnd(); ++term){
            (*term).RunIntegration();
            DYTurbo::PrintTable::Result((*term));
            if (!DYTurbo::isDryRun){
                double val = (*term).last_int[0];
                double exp = readfloat();
                if (doSaveResults) writefloat(val);
                //else ASSERT_TRUE(ApproxDoubles(exp,val));
                else ASSERT_DOUBLE_EQ(exp,val);
            }
        }
        DYTurbo::PrintTable::ResultSubTotal();
    }
    DYTurbo::PrintTable::ResultGrandTotal();
    fclose(F);
}


// Helper interface to DYTurbo for testing all terms.
namespace DYTurbo{ 
    extern bool TestAllTerms; 
    void init_params();
}



//! RunIntegrand for cubature like integrands.
void RunIntegrand( int (* (*fun)(const int&, const double*, const int&, double*))(const int*, const double*, const int*, double*, void*),
        int & dim, double *point,double * result
        ){
    int ncomp = 1;
    fun(dim,point,ncomp,result);
}
//! RunIntegrand for Vegas like integrands.
void RunIntegrand( int (* (*fun)(const int&, const double*, const int&, double*, void*, const int&, const int&, double&, const int&))(const int*, const double*, const int*, double*, void*),
        int & dim, double *point,double * result
        ){
    int ncomp = 1;
    void* userdata=0;
    const int nvec=0;
    const int core=0;
    double weight=1.0;
    const int iter=1;
    fun(dim,point,ncomp,result,userdata,nvec,core,weight,iter);
}
//! Run and check integrand output.
template<typename IntFun>
::testing::AssertionResult CheckIntegrand(int &ord, const char *name, IntFun fun, int dim){
    /// @todo For finite born level set random number corresponding to pt to 0
    // point: Fake random point for integrands, set to lower bound of integration.
    VecDbl point  (dim, 0.8);
    VecDbl result (opts.totpdf,0.);
    double time_ellapsed = clock_real();
    RunIntegrand(fun, dim,&point[0],&result[0]);
    time_ellapsed = clock_real()-time_ellapsed;
    if (doSaveResults){
        printf(" saving to file: %d %s : %.10e ( %fs) \n", dim, name, result[0], time_ellapsed);
        writefloat(result[0]);
        return ::testing::AssertionSuccess();
    } else {
        double expected = readfloat();
        return ApproxDoubles(result[0],expected) << "Wrong value of integrand " <<  name << " order="<<ord <<endl;
    }
}

/**
 * @brief Testing the output of every integrand function separatelly
 */
TEST(DYTurbo,CheckIntegrandFunctions){
    // turn off dryrun
    DYTurbo::isDryRun = true;
    // turn on all terms 
    DYTurbo::TestAllTerms = true;
    // set boundaries
    bins.qtbins = {0., 10. , 80., 100.  };
    bins.ybins  = {0., 1. , 3.   };
    bins.mbins  = {80. , 100. };
    opts.costhmin=-1;
    opts.costhmax=+1;
    // warm up
    DYTurbo::WarmUp();
    // for all orders and for all terms
    F = fopen( opts.makecuts ? "terms_cut.res" : "terms.res", doSaveResults ? "w" : "r");
    for (int ord = 0; ord < 3; ++ord) {
        // NOTE: If you are adding new term, dont forget to run first with
        // `doSaveResults=true;` to update values in file.
        opts.order=ord;
        // Resummed part
        opts.fixedorder=false;
        DYTurbo::init_params();
        for (DYTurbo::BoundIterator bound; !bound.IsEnd(); ++bound){
            DYTurbo::SetBounds(bound);
            ASSERT_TRUE(CheckIntegrand ( ord , "Resintegr 2D"      , resintegrand2d , 2 ) ) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "Resintegr 3D"      , resintegrand3d , 3 ) ) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "Resintegr MC"      , resintegrandMC , 6 ) ) ;

            ASSERT_TRUE(CheckIntegrand ( ord , "Counter Resum 2D"  , ctintegrand2d  , 2 ) ) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "Counter Resum 3D"  , ctintegrand3d  , 3 ) ) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "Counter Resum MC6" , ctintegrandMC  , 6 ) ) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "Counter Resum MC8" , ctintegrand    , 8 ) ) ;
        }
        // Fixed part
        opts.fixedorder=true;
        DYTurbo::init_params();
        for ( DYTurbo::BoundIterator bound ; !bound.IsEnd(); ++bound){
            DYTurbo::SetBounds(bound);
            ASSERT_TRUE(CheckIntegrand ( ord , "Born Fixed 2D"     , lointegrand2d       , 2  ))  ;
            ASSERT_TRUE(CheckIntegrand ( ord , "Born Fixed MC4"    , lointegrandMC       , 4  )) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "Born DYNNLO MC6"   , doublevirtintegrand , 6  )) ;

            ASSERT_TRUE(CheckIntegrand ( ord , "Counter Fixed 2D"  , ctintegrand2d       , 2  )) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "Counter Fixed 3D"  , ctintegrand3d       , 3  )) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "Counter Fixed MC6" , ctintegrandMC       , 6  )) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "Counter Fixed MC8" , ctintegrand         , 8  )) ;

            ASSERT_TRUE(CheckIntegrand ( ord , "VJ LO 3D"          , vjintegrand         , 3  )) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "VJ LO 5D"          , vjlointegrand       , 5  )) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "VJ LO MC7"         , vjlointegrandMC     , 7  )) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "VJ LO MCFM7"       , lowintegrand        , 7  )) ;

            ASSERT_TRUE(CheckIntegrand ( ord , "VJ REAL MCFM"      , realintegrand       , 10 )) ;
            ASSERT_TRUE(CheckIntegrand ( ord , "VJ VIRT MCFM"      , virtintegrand       , 8  )) ;

            ASSERT_TRUE(CheckIntegrand ( ord , "VJJ MC10"          , v2jintegrand        , 10 )) ;
        }

    }
    fclose(F);
    /**
     *  @note To see values saved in file run
     *
     *         hexdump -v -e '"%010_ad :" 7/8 " %e " "\n"' test/terms.res
     *
     *    - one line is seven numbers, there are 6 kinematic points
     *    - first 6 lines is for resummed part
     *    - next 12 lines is for fixed order part
     *    - This repeated 3 times (per each order)
     */
}

#ifdef USEROOT
#include "TFile.h"
#include "TH1.h"
#endif

void CheckResultFile(string fname = ""){
    if (fname.size()==0){
        fname =  HistoHandler::result_filename;
        fname += HistoHandler::file_suffix;
    }
    VecStr names = { "s_qt" };
    F= fopen( opts.makecuts ? "filecheck_cuts.res" : "filecheck.res", (doSaveResults) ? "w" : "r");
#ifdef USEROOT
    TFile *f = TFile::Open(fname.c_str());
    ASSERT_NE(f,NULL) << "Wrong file name " << fname;
    for (String hname: names){
        TH1* hist = (TH1*) f->Get(hname.c_str());
        ASSERT_NE(hist,NULL) << "Wrong hist name " << hname << " in filename " << fname;
        double val = hist->Integral();
        double exp = readfloat();
        if (doSaveResults) writefloat(val);
        //else ASSERT_TRUE(ApproxDoubles(val,exp))
        else ASSERT_DOUBLE_EQ(val,exp)
            << "Incorrect number of entries in saved histogram " << hname
            << " inside file" << fname;
    }
#else // STL
#endif
    fclose(F);
}


TEST(DYTurbo,Termination){
    DYTurbo::Terminate();
    CheckResultFile("../src/results.root");
    CheckResultFile();
}

#endif /* ifndef DYTurbo_unittest_CXX */
