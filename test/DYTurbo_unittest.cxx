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

#include "src/dyturbo.h"
#include "src/settings.h"
#include "src/cubacall.h"

#include "src/resintegr.h"
#include "src/ctintegr.h"
#include "src/finintegr.h"
#include "mcfm/mcfm_interface.h"

typedef std::vector<string> VecStr;
typedef std::vector<double> VecDbl;
typedef std::stringstream SStream;



TEST(DYTurbo,Initialization){
    int argc = 2;
    char *argv[] =  {
        ((char *) "DYTurbo_unittest"),
        ((char *) "../input/test.in"),
    };
    DYTurbo::Init(argc,argv);
}

TEST(DYTurbo,TermIteration){
    DYTurbo::WarmUp();
    std::vector <void (*)(std::vector<double>&, double&)> term_funs;
    term_funs.push_back(resintegr2d);
    term_funs.push_back(ctintegr2d);
    term_funs.push_back(vjrealintegr);
    term_funs.push_back(vjvirtintegr);
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
}

TEST(DYTurbo,DryLooping){
    DYTurbo::isDryRun = true;
    DYTurbo::WarmUp();
    DYTurbo::PrintTable::Header();
    for ( DYTurbo::BoundIterator bounds; !bounds.IsEnd(); ++bounds) {
        DYTurbo::SetBounds(bounds);
        DYTurbo::PrintTable::Bounds();
        for (DYTurbo::TermIterator term; !term.IsEnd(); ++term){
            (*term).RunIntegration();
            DYTurbo::PrintTable::Result((*term));
        }
        DYTurbo::PrintTable::ResultSubTotal();
    }
    DYTurbo::PrintTable::ResultGrandTotal();
}

// Testing term
namespace DYTurbo{ 
    extern bool TestAllTerms; 
    void init_params();
}

bool onlyPrintResults=false;
FILE *F = 0;
void writefloat(double v) {
  fwrite((void*)(&v), sizeof(v), 1, F);
}
double readfloat() {
  double v=666.;
  if (!onlyPrintResults) {
      int i = fread((void*)(&v), sizeof(v), 1, F);
  }
  return v;
}



void RunIntegrand( int (* (*fun)(const int&, const double*, const int&, double*))(const int*, const double*, const int*, double*, void*),
        int & dim, double *point,double * result
        ){
    int ncomp = 1;
    fun(dim,point,ncomp,result);
}
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
template<typename IntFun>
void CheckIntegrand(int &ord, const char *name, IntFun fun, int dim){
    VecDbl point  (dim, 0.0); // fake random point for integrands, set to lower bound of integration
    VecDbl result (opts.totpdf,0.);
    RunIntegrand(fun, dim,&point[0],&result[0]);
    if (onlyPrintResults){
        printf(" saving to file: %d %s : %.10e \n", dim, name, result[0]);
        writefloat(result[0]);
    } else {
        double expected = readfloat();
        ASSERT_EQ(result[0],expected) << "Wrong value of integrand " <<  name << " order="<<ord;
    }
}

TEST(DYTurbo,CheckIntegrandFunctions){
    // turn off dryrun
    DYTurbo::isDryRun = true;
    // turn on all terms 
    DYTurbo::TestAllTerms = true;
    // set boundaries
    opts.nproc=3;
    nproc_.nproc_ = opts.nproc;
    bins.qtbins = {0., 10. , 50., 100.  };
    bins.ybins  = {0., 1. , 3.   };
    bins.mbins  = {80. , 100. };
    opts.costhmin=-1;
    opts.costhmax=+1;
    // warm up
    DYTurbo::WarmUp();
    // for all orders and for all terms
    if (onlyPrintResults){
        F = fopen("terms.res","w");
    } else {
        F = fopen("terms.res","r");
    }
        for (int ord = 0; ord < 3; ++ord) {
            // NOTE: If you are adding new term, dont forget to run first with
            // `onlyPrintResults=true;` to update values in file.
            opts.order=ord;
            // Resummed part
            opts.fixedorder=false;
            DYTurbo::init_params();
            for (DYTurbo::BoundIterator bound; !bound.IsEnd(); ++bound){
                DYTurbo::SetBounds(bound);
                CheckIntegrand ( ord , "Resintegr 2D"      , resintegrand2d , 2 ) ;
                CheckIntegrand ( ord , "Resintegr 3D"      , resintegrand3d , 3 ) ;
                CheckIntegrand ( ord , "Resintegr MC"      , resintegrandMC , 6 ) ;

                CheckIntegrand ( ord , "Counter Resum 2D"  , ctintegrand2d  , 2 ) ;
                CheckIntegrand ( ord , "Counter Resum 3D"  , ctintegrand3d  , 3 ) ;
                CheckIntegrand ( ord , "Counter Resum MC6" , ctintegrandMC  , 6 ) ;
                CheckIntegrand ( ord , "Counter Resum MC8" , ctintegrand    , 8 ) ;
            }
            // Fixed part
            opts.fixedorder=true;
            DYTurbo::init_params();
            for ( DYTurbo::BoundIterator bound ; !bound.IsEnd(); ++bound){
                DYTurbo::SetBounds(bound);
                CheckIntegrand ( ord , "Born Fixed 2D"     , lointegrand2d       , 2   ) ;
                CheckIntegrand ( ord , "Born Fixed MC4"    , lointegrandMC       , 4   ) ;

                CheckIntegrand ( ord , "Counter Fixed 2D"  , ctintegrand2d       , 2   ) ;
                CheckIntegrand ( ord , "Counter Fixed 3D"  , ctintegrand3d       , 3   ) ;
                CheckIntegrand ( ord , "Counter Fixed MC6" , ctintegrandMC       , 6   ) ;
                CheckIntegrand ( ord , "Counter Fixed MC8" , ctintegrand         , 8   ) ;

                CheckIntegrand ( ord , "VJ LO 3D"          , vjintegrand         , 3   ) ;
                CheckIntegrand ( ord , "VJ LO 5D"          , vjlointegrand       , 5   ) ;
                CheckIntegrand ( ord , "VJ LO MC7"         , vjlointegrandMC     , 7   ) ;

                CheckIntegrand ( ord , "VJ LO MC7"         , lowintegrand        , 7   ) ;
                CheckIntegrand ( ord , "VJ REAL MC"        , realintegrand       , 10  ) ;
                CheckIntegrand ( ord , "VJ VIRT MC"        , virtintegrand       , 8   ) ;
                CheckIntegrand ( ord , "VJ VV MC6"         , doublevirtintegrand , 6   ) ;
                CheckIntegrand ( ord , "VJ VV MC10"        , v2jintegrand        , 10  ) ;
            }

        }
    fclose(F);
    // To see values saved in file run
    //  `hexdump -v -e '"%010_ad :" 7/8 " %e " "\n"' test/terms.res`
    //   - one line is seven numbers, there are 6 kinematic points
    //   - first 6 lines is for resummed part 
    //   - next 12 lines is for fixed order part
    //   - This repeated 3 times (per each order)
}

TEST(DYTurbo,Termination){
    DYTurbo::PrintTable::Footer();
    DYTurbo::Terminate();
}

#endif /* ifndef DYTurbo_unittest_CXX */
