#ifndef TurboHist_unittest_CXX
#define TurboHist_unittest_CXX
/**
 * @file TurboHist_unittest.cxx
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-05
 */

#include "config.h"
#include "TurboHist_File.h"
#include "TurboHist_H1.h"
#include "TurboHist_H2.h"

#include <algorithm>
using std::lower_bound;
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <map>
using std::map;

typedef vector<string> VecStr;
typedef vector<double> VecDbl;
typedef vector<vector<double>> VecVecDbl;

#include "gtest/gtest.h"
#include "gmock/gmock.h"



TEST(TurboHist, Counter) {
    TurboHist::Counter c1;
    ASSERT_DOUBLE_EQ ( 0 , c1.sum_w  ) << "Construct ";
    ASSERT_DOUBLE_EQ ( 0 , c1.sum_w2 ) << "Construct ";
    // add double
    c1+=1;
    ASSERT_DOUBLE_EQ ( 1 , c1.sum_w  ) << "Add double";
    ASSERT_DOUBLE_EQ ( 1 , c1.sum_w2 ) << "Add double";
    // mult double
    c1*=2;
    ASSERT_DOUBLE_EQ ( 2 , c1.sum_w  ) << "Mult double";
    ASSERT_DOUBLE_EQ ( 4 , c1.sum_w2 ) << "Mult double";
    // 
    TurboHist::Counter c2;
    c2+=1;
    ASSERT_DOUBLE_EQ ( 1 , c2.sum_w  ) << "Contruct2";
    ASSERT_DOUBLE_EQ ( 1 , c2.sum_w2 ) << "Contruct2";
    // mult
    ASSERT_DOUBLE_EQ ( (c2*2).sum_w  , c1.sum_w  ) << "num*Counter";
    ASSERT_DOUBLE_EQ ( (c2*2).sum_w2 , c1.sum_w2 ) << "num*Counter";
    // add counter
    c2+=c1;
    ASSERT_DOUBLE_EQ ( 3 , c2.sum_w  ) << "Add Counter";
    ASSERT_DOUBLE_EQ ( 5 , c2.sum_w2 ) << "Add Counter";
}

template <class HType>
void CheckSpecifiers(HType &histo, const int &dim, const char &type, const size_t &maxbins){
    ASSERT_EQ ( dim     , histo.GetDim    ( ) ) << "Wrong dimmension.";
    ASSERT_EQ ( type    , histo.GetType   ( ) ) << "Wrong type.";
    ASSERT_EQ ( maxbins , histo.GetMaxBin ( ) ) << "Wrong number of bins.";
}

template <class HType>
void CheckBinContent(HType &histo, VecDbl binval){
    for (size_t ibin = 0; ibin < binval.size(); ++ibin) {
        ASSERT_DOUBLE_EQ(binval[ibin], histo.GetBinContent(ibin)) << "Wrong bin CONTENT in bin " << ibin;
    }
}

template <class HType>
void CheckBinError(HType &histo, VecDbl binval){
    for (size_t ibin = 0; ibin < binval.size(); ++ibin) {
        ASSERT_DOUBLE_EQ(binval[ibin], histo.GetBinError(ibin)) << "Wrong bin ERROR in bin " << ibin;
    }
}

//===============================
// Init, binning and filling
//===============================
TEST(TurboHist, InitH1) {
    TurboHist::H1 h;
    CheckSpecifiers(h, 1, 'h', 3);
    CheckBinContent(h, {0.,0.,0.});
    CheckBinError(h, {0.,0.,0.});
}

TEST(TurboHist, InitH2) {
    TurboHist::H2 h;
    CheckSpecifiers(h, 2, 'h', 9);
    CheckBinContent(h, {0.,0.,0., 0.,0.,0., 0.,0.,0. });
    CheckBinError  (h, {0.,0.,0., 0.,0.,0., 0.,0.,0. });
}

// Helper functions
// Array comparison : http://stackoverflow.com/a/10062016
// You don't need to add a dependency on googlemock if you don't want, you
// could write your own simple function that returns a
// testing::AssertionResult, e.g.
template<typename T> ::testing::AssertionResult ArraysMatch(const T &expected, const T &actual){
    if ( expected.size() != actual.size() )
        return ::testing::AssertionFailure() << "incorrect length " <<  actual.size()
           << "instead of "  << expected.size();
    for (size_t i(0); i < expected.size(); ++i){
        if (expected[i] != actual[i]){
            return ::testing::AssertionFailure() << "array[" << i
                << "] (" << actual[i] << ") != expected[" << i
                << "] (" << expected[i] << ")";
        }
    }
    return ::testing::AssertionSuccess();
}
// Then in your test, call:
//EXPECT_TRUE(ArraysMatch(two_sorted, two));

TEST(TurboHist, CreateEquidistantVector) {
    VecDbl bins = {0.,2.,4.,6.};
    VecDbl equid = TurboHist::Binning::GetEquidistantVector(3,0.,6.);
    ASSERT_TRUE(ArraysMatch(bins,equid)) << "Wrong implementation of equidistant binning";
}



TEST(TurboHist, Binning) {
    VecDbl bins = {0.,2.,4.,6.};
    //
    TurboHist::H1 turbo;
    // equidistant
    turbo.SetBins(3,0.,6.);
    ASSERT_TRUE(turbo.IsEquidistant()) << "From range: Wrong test of equidistancy." << turbo.Print();
    // Find
    ASSERT_EQ( 0, turbo.FindBin(-100)) << "From min: Incorrect bin found.";
    ASSERT_EQ( 4, turbo.FindBin(100)) << "From max: Incorrect bin found.";
    ASSERT_EQ( 2, turbo.FindBin(3.33)) << "From range: Incorrect bin found.";
    //ASSERT_EQ( 2, turbo.GetBin(2)) << "Incorrect bin calculated.";
    // equidistant
    turbo.SetBins(bins);
    ASSERT_TRUE(turbo.IsEquidistant()) << "From array: Wrong test of equidistancy.";
    // Find
    ASSERT_EQ( 2, turbo.FindBin(3.33)) << "From equid array: Incorrect bin found.";
    //ASSERT_EQ( 2, turbo.GetBin(2)) << "Incorrect bin calculated.";
    // non-equidistant
    bins = {0.,3.,4.,6.};
    turbo.SetBins(bins);
    ASSERT_FALSE(turbo.IsEquidistant()) << "Wrong test of equidistancy.";
    // Find
    ASSERT_EQ( 2, turbo.FindBin(3.33)) << "From array: Incorrect bin found.";
    //ASSERT_EQ( 2, turbo.GetBin(2)) << "Incorrect bin calculated.";
    // uniq
    bins = {0., 2.,2.,4.,6.};
    turbo.SetBins(bins);
    EXPECT_EQ( 3, turbo.FindBin(4.33)) << "Repeating bin edges: Incorrect bin found.";
    ASSERT_TRUE(turbo.IsEquidistant()) << "Repeating bin edges: Wrong test of equidistancy.";
}


TEST(TurboHist, Filling) {
    TurboHist::H1 turbo;
    //
    VecDbl bins = {0.,2.,4.,6.};
    turbo.SetBins(bins);
    //
    turbo.Fill(1.5);
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinContent(0)) << turbo.Print();
    ASSERT_DOUBLE_EQ(1.,turbo.GetBinContent(1)) << turbo.Print();
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinContent(2)) << turbo.Print();
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinContent(3)) << turbo.Print();
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinContent(4)) << turbo.Print();
    //
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinError (0)) << turbo.Print();
    ASSERT_DOUBLE_EQ(1.,turbo.GetBinError (1)) << turbo.Print();
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinError (2)) << turbo.Print();
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinError (3)) << turbo.Print();
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinError (4)) << turbo.Print();
    //
    turbo.Fill(3.3,3.3);
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinContent(0))  << turbo.Print();
    ASSERT_DOUBLE_EQ(1.,turbo.GetBinContent(1))  << turbo.Print();
    ASSERT_DOUBLE_EQ(3.3,turbo.GetBinContent(2)) << turbo.Print();
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinContent(3))  << turbo.Print();
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinContent(4))  << turbo.Print();
    //
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinError  (0)) << turbo.Print();
    ASSERT_DOUBLE_EQ(1.,turbo.GetBinError  (1)) << turbo.Print();
    ASSERT_DOUBLE_EQ(3.3,turbo.GetBinError (2)) << turbo.Print();
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinError  (3)) << turbo.Print();
    ASSERT_DOUBLE_EQ(0.,turbo.GetBinError  (4)) << turbo.Print();
}

TEST(TurboHist, AddHist) {
    VecDbl bins = {0.,2.,4.,6.};
    //
    TurboHist::H1 turbo1;
    turbo1.SetBins(bins);
    turbo1.Fill(1.5);
    //
    ASSERT_DOUBLE_EQ(0.,turbo1.GetBinContent(0)) << turbo1.Print();
    ASSERT_DOUBLE_EQ(1.,turbo1.GetBinContent(1)) << turbo1.Print();
    ASSERT_DOUBLE_EQ(0.,turbo1.GetBinContent(2)) << turbo1.Print();
    ASSERT_DOUBLE_EQ(0.,turbo1.GetBinContent(3)) << turbo1.Print();
    ASSERT_DOUBLE_EQ(0.,turbo1.GetBinContent(4)) << turbo1.Print();
    //
    TurboHist::H1 turbo2;
    turbo2.SetBins(bins);
    //
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(0)) << turbo2.Print();
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(1)) << turbo2.Print();
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(2)) << turbo2.Print();
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(3)) << turbo2.Print();
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(4)) << turbo2.Print();
    //
    turbo2.Fill(2.5);
    //
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(0)) << turbo2.Print();
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(1)) << turbo2.Print();
    ASSERT_DOUBLE_EQ(1.,turbo2.GetBinContent(2)) << turbo2.Print();
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(3)) << turbo2.Print();
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(4)) << turbo2.Print();
    //
    turbo2.Add(turbo1);
    //
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(0)) << "After add :" << turbo2.Print();
    ASSERT_DOUBLE_EQ(1.,turbo2.GetBinContent(1)) << "After add :" << turbo2.Print();
    ASSERT_DOUBLE_EQ(1.,turbo2.GetBinContent(2)) << "After add :" << turbo2.Print();
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(3)) << "After add :" << turbo2.Print();
    ASSERT_DOUBLE_EQ(0.,turbo2.GetBinContent(4)) << "After add :" << turbo2.Print();
}

TEST(TurboHist, IOoperations) {
    // init and fill
    TurboHist::H1 turbo;
    VecDbl bins = {0.,2.,4.,6.};
    turbo.SetBins(bins);
    turbo.Fill(1.5);
    turbo.Fill(3.3,3.3);
    // open file
    TurboHist::File f;
    f.Open("turbo.dat","RECREATE");
    f << turbo;
    f.Close();
    //
    f.Open("turbo.dat"); // READ
    TurboHist::H1 turbo2;
    ASSERT_TRUE(turbo2.LoadFromFile(f,true)) << "Test the hist is written correctly";
    f >> turbo2;
    f.Close();
    //
    ASSERT_EQ(turbo2.FindBin(1.1), turbo.FindBin(1.1));
    ASSERT_DOUBLE_EQ(turbo2.GetBinContent(2), turbo.GetBinContent(2));
    ASSERT_DOUBLE_EQ(turbo2.GetBinError(2), turbo.GetBinError(2));
    //
    f.Open("turbo2.dat","RECREATE");
    f << turbo2;
    f.Close();
    turbo.Add(turbo2);
    //
    TurboHist::MergeFiles( "turbo3.dat", {"turbo.dat", "turbo2.dat"});
    //
    //f.Open("turbo3.data");
    //TurboHist::H1 turbo3;
    //f >> turbo3;
    //ASSERT_EQ(turbo3.GetBin(1.1), turbo.GetBin(1.1));
    //ASSERT_DOUBLE_EQ(turbo3.GetBinContent(2), turbo.GetBinContent(2));
    //ASSERT_DOUBLE_EQ(turbo3.GetBinError(2), turbo.GetBinError(2));
}

#ifdef USEROOT
#include "TH1D.h"
#include "TRandom3.h"

//===============================
// Versus ROOT
//===============================

// Time measurement
template <class T>
double tp_find(T* h, VecVecDbl & datastream ){
    clock_t startTime = clock();
    for (vector<double> &a : datastream) h->FindBin(a[0]);
    return double( clock() - startTime ) / (double)CLOCKS_PER_SEC;
}
template <class T>
double tp_find2(T* h, VecVecDbl & datastream ){
    clock_t startTime = clock();
    for (vector<double> &a : datastream) h->FindBin(a[0],a[1]);
    return double( clock() - startTime ) / (double)CLOCKS_PER_SEC;
}
template <class T>
double tp_fill(T* h, VecVecDbl & datastream ){
    clock_t startTime = clock();
    for (vector<double> &a : datastream) h->Fill(a[0]);
    return double( clock() - startTime ) / (double)CLOCKS_PER_SEC;
}
template <class T>
double tp_fill2(T* h, VecVecDbl & datastream ){
    clock_t startTime = clock();
    for (vector<double> &a : datastream) h->Fill(a[0],a[1]);
    return double( clock() - startTime ) / (double)CLOCKS_PER_SEC;
}

void generateData(VecVecDbl & vec, int N = 1e6){
    TRandom * rnd = new TRandom3();
    for (int i = 0; i < N; ++i) {
        VecDbl in;
        in.push_back(rnd->Uniform(0.,100.));
        in.push_back(rnd->Uniform(0.,100.));
        vec.push_back(in);
    }
}

struct VectoBin{
    VecDbl data;
    VecDbl::const_iterator beg;

    VectoBin(VecDbl& newbins){
        for (size_t ilo = 0; ilo < newbins.size(); ++ilo) {
            double edge = newbins[ilo];
            data.push_back(edge); // check for duplication
        }
        beg = data.cbegin();
    }

    inline const size_t FindBin(const double & val) const {
        return lower_bound(data.begin(), data.end(), val)-beg;
    }
};

struct ArrBin{
    VecDbl data;
    VecDbl::const_iterator beg;
    size_t N;
    const double* BEG;
    const double* END;

    ArrBin(VecDbl& newbins){
        for (size_t ilo = 0; ilo < newbins.size(); ++ilo) {
            double edge = newbins[ilo];
            data.push_back(edge); // check for duplication
        }
        N = data.size();
        BEG = &data[0];
        END = &data[N];
    }

    inline const size_t FindBin(const double & val) const {
        return std::lower_bound(BEG, BEG+N, val)-BEG;
    }
};

struct MapoBin{
    map<double,size_t> data;
    MapoBin(VecDbl& newbins){
        size_t ibin = 1;
        for (size_t ilo = 0; ilo < newbins.size(); ++ilo) {
            double edge = newbins[ilo];
            if(data.count(edge)==0) data[edge]=ibin++; // check for duplication
        }
    }

    inline const size_t FindBin(const double & val) const {
        return data.lower_bound(val)->second-1;
    }
};

// We have to be faster then ROOT
TEST(TurboHistVsROOT, FourTimesFasterThanROOT_FindBin) {
    // prepare data
    VecDbl bins ={ 0., 1., 2., 3., 4., 5., 6., 7., 8. , 9., 10., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54., 55., 56., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69., 70., 71., 72., 73., 74., 75., 76., 77., 78., 79., 80., 81., 82., 83., 84., 85., 86., 87., 88., 89., 90., 91., 92., 93., 94., 95., 96., 97., 98., 99., 100. };
    VecVecDbl testdata;
    generateData(testdata);
    //
    //VectoBin vb (bins);
    //MapoBin mb (bins);
    //ArrBin ab (bins);
    //double timeVec = tp_find(&vb,testdata);
    //timeVec = tp_find(&vb,testdata);
    //double timeMap =  tp_find(&mb,testdata);
    //double timeArr =  tp_find(&ab,testdata);
    //EXPECT_GT(timeMap,timeVec) << "time is in seconds -- BinarySearch Vec vs Map";
    //EXPECT_GT(timeArr,timeVec) << "time is in seconds -- BinarySearch Vec vs Arr";
    //
    // init
    TurboHist::H1 *turbo = new TurboHist::H1 ("name","title", 100, 0, 100);
    TH1D          *root  = new TH1D          ("name5","title", 100, 0, 100);
    // run time
    double timeTurbo = tp_find(turbo,testdata);
    double timeROOT =  tp_find(root,testdata);
    EXPECT_GT(timeROOT,timeTurbo) << "time is in seconds -- Equidistant";
    //EXPECT_GT(timeROOT,4*timeTurbo) << "time is in seconds -- Equidistant";
    // non-equidistant
    turbo->SetBins(bins);
    delete root;
    root  = new TH1D          ("name0","title", bins.size()-1, &bins[0]);
    // run time
    timeROOT =  tp_find(root,testdata);
    timeTurbo = tp_find(turbo,testdata);
    timeROOT =  tp_find(root,testdata);
    EXPECT_GT(timeROOT,timeTurbo) << "time is in seconds -- BinarySearch";
    EXPECT_GT(timeROOT,4*timeTurbo) << "time is in seconds -- BinarySearch";
}

// We have to be faster then ROOT
TEST(TurboHistVsROOT, FourTimesFasterThanROOT_Fill) {
    // prepare data
    VecDbl bins ={ 0., 1., 2., 3., 4., 5., 6., 7., 8. , 9., 10., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51., 52., 53., 54., 55., 56., 57., 58., 59., 60., 61., 62., 63., 64., 65., 66., 67., 68., 69., 70., 71., 72., 73., 74., 75., 76., 77., 78., 79., 80., 81., 82., 83., 84., 85., 86., 87., 88., 89., 90., 91., 92., 93., 94., 95., 96., 97., 98., 99., 100. };
    VecVecDbl testdata;
    generateData(testdata);
    //
    // init
    TurboHist::H1 *turbo = new TurboHist::H1 ("name","title", 100, 0, 100);
    TH1D          *root  = new TH1D          ("name","title", 100, 0, 100);
    // run time
    double timeTurbo = tp_fill(turbo,testdata);
    double timeROOT =  tp_fill(root,testdata);
    EXPECT_GT(timeROOT,timeTurbo) << "time is in seconds -- Equidistant";
    //EXPECT_GT(timeROOT,4*timeTurbo) << "time is in seconds -- Equidistant";
    // non-equidistant
    turbo->SetBins(bins);
    delete root;
    root  = new TH1D          ("name1","title", bins.size()-1, &bins[0]);
    // run time
    timeROOT =  tp_fill(root,testdata);
    timeTurbo = tp_fill(turbo,testdata);
    timeROOT =  tp_fill(root,testdata);
    EXPECT_GT(timeROOT,timeTurbo) << "time is in seconds -- BinarySearch";
    EXPECT_GT(timeROOT,4*timeTurbo) << "time is in seconds -- BinarySearch";
}


// Make sure that filling result is same
TEST(TurboHistVsROOT, SameNumbersAfterFill) {
    // init
    size_t N = 100;
    TurboHist::H1 *turbo = new TurboHist::H1 ("name","title", N, 0, 100);
    TH1D          *root  = new TH1D          ("name2","title", N, 0, 100);
    root->Sumw2();
    //
    for (size_t i = 0; i < N+2; ++i){
        ASSERT_DOUBLE_EQ(turbo->GetBinContent(i),root->GetBinContent(i)) << "Before fill, unequal bin content for i=" << i;
        ASSERT_DOUBLE_EQ(turbo->GetBinError  (i),root->GetBinError  (i)) << "Before fill, unequal bin error for i=" << i;
    }
    VecVecDbl testdata;
    generateData(testdata,400);
    for (auto a : testdata) {
        turbo->Fill(a[0],a[1]);
        root->Fill(a[0],a[1]);
    }
    //
    for (size_t i = 0; i < N+2; ++i){
        ASSERT_DOUBLE_EQ(turbo->GetBinContent(i),root->GetBinContent(i)) << "After fill, unequal bin content for i=" << i;
        ASSERT_DOUBLE_EQ(turbo->GetBinError  (i),root->GetBinError  (i)) << "After fill, unequal bin error for i=" << i;
    }
}
#endif /* ROOT */



#endif /* ifndef TurboHist_unittest_CXX */
