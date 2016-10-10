/**
 * @file HistoSvc.cxx
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-26
 */

#include <sys/wait.h>
#include <unistd.h>

#include "histo/HistoHandler.h"
#include "histo/HistoObjects.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <random>

#include "src/handy_typdefs.h"
using HistoHandler::histos;


class HistoServiceTest : public testing::Test{
    protected :

        virtual void SetUp(){
            gen.seed(12345);
            gausT.param( std::normal_distribution<>::param_type(0.,10.)  );
            gausZ.param( std::normal_distribution<>::param_type(0.,30.)  );
            gausW.param( std::normal_distribution<>::param_type(1.,0.1) );

            Nterm=2;
            Nevents=100;
            Nnodes=4;
            isParent=true;
            parent_pid=getpid();
        }

        virtual void TearDown(){
        }


        void random_kin(){

            l1[0]= gausT(gen);
            l1[1]= gausT(gen);
            l1[2]= gausZ(gen);
            l1[3]=sqrt( l1[0]*l1[0] + l1[1]*l1[1] + l1[2]*l1[2]);

            l2[0]= gausT(gen);
            l2[1]= gausT(gen);
            l2[2]= gausZ(gen);
            l2[3]=sqrt( l2[0]*l2[0] + l2[1]*l2[1] + l2[2]*l2[2]);

            wgt = gausW(gen);
        }

        void CheckEntries(double entries=0., double integr=0.){
            for (size_t ihist = 0; ihist < histos.size(); ++ihist) {
                double hentries = histos[ihist]->GetEntries();
                double exp = entries + (histos[ihist]->IsIntegrationSafe() ? integr : 0) ;
                ASSERT_DOUBLE_EQ(exp, histos[ihist]->GetEntries())
                    << "Histogram " << histos[ihist]->GetName()
                    << " has incorect number of entries: " << hentries 
                    << " instead of " << exp ;
            }
        }

        void Loop(){
            for (int iTerm = 0; iTerm < Nterm; ++iTerm) {
                for (int iEvent = 0; iEvent < Nevents; ++iEvent) {
                    random_kin();
                    HistoHandler::FillEvent(l1,l2,wgt);
                }
            }
        }
        void LoopIntegrator(){
            for (int iTerm = 0; iTerm < Nterm; ++iTerm) {
                // Set Middle point
                // Add to bin
                HistoHandler::FillResult(1.0,0.1);
            }
        }

        bool IsParent(){
            return parent_pid == getpid();
            // parent skip and fork again ( double condition just to be sure)
            //if ((child_pid!=0) || (this_pid==parent_pid)) continue;
        }

        void Fork(){
            if (IsParent()){
                child_pid = fork();
            }
        }

        void Wait(){
            // from: http://stackoverflow.com/a/279761
            //printf("parent is waiting...\n"); fflush(stdout);
            while (true) {
                int status;
                pid_t done = wait(&status);
                if (done == -1) {
                    if (errno == ECHILD) break; // no more child processes
                } else {
                    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
                        printf("pid %d failed!!!\n",done); fflush(stdout);
                        exit(1);
                    }
                }
            }
        }

        void InitWorker(){
            HistoHandler::Reset();
        }

        void ExitWorker(int inode){
            HistoHandler::Save(inode);
            exit(0);
        }

        double l1[4];
        double l2[4];
        double wgt;

        int Nterm;
        int Nevents;
        int Nnodes;
        bool isParent;
        static double total_entries;
        static double total_entries_integ;
        static VecStr hist_names;


        int child_pid;
        int parent_pid;

        std::mt19937 gen;
        std::normal_distribution<> gausT; // (0.,10.);
        std::normal_distribution<> gausZ; // (0.,30.);
        std::normal_distribution<> gausW; // (1.,0.1);
};

double HistoServiceTest::total_entries = 0.;
double HistoServiceTest::total_entries_integ = 0.;
VecStr HistoServiceTest::hist_names;

TEST_F(HistoServiceTest, HistCreateAndFill){
    HistoHandler::Histo1D<BosPT> htmp("qt");
    random_kin();
    Kinematics::SetKinematics(l1,l2,wgt);
    htmp.FillEvent();
    ASSERT_DOUBLE_EQ(1., htmp.GetEntries()) << "FillEvent failed.";
    //
    random_kin();
    Kinematics::SetKinematics(l1,l2,wgt);
    htmp.FillDipole();
    htmp.FillDipole();
    htmp.FillRealEvent();
    ASSERT_DOUBLE_EQ(2., htmp.GetEntries()) << "FillRealEvent failed.";
};

TEST_F(HistoServiceTest, Init){
    HistoHandler::Init();
    total_entries=0.;
    CheckEntries(total_entries);
    // store all names of all histograms
    for (auto h_it = HistoHandler::histos.begin();h_it!=HistoHandler::histos.end();h_it++){
        hist_names.push_back((*h_it)->GetName());
    }
};

TEST_F(HistoServiceTest, LoopNoThread){
    Loop();
    total_entries+=Nterm*Nevents;
    CheckEntries(total_entries);
};

TEST_F(HistoServiceTest, LoopNoThreadIntegrator){
    LoopIntegrator();
    total_entries_integ+=Nterm;
    CheckEntries(total_entries,total_entries_integ);
}

TEST_F(HistoServiceTest, LoopInsideThreads){
    for (int inodes = 0;  inodes < Nnodes; inodes++) {
        Fork();
        if(IsParent()) continue;
        // only child goes here
        InitWorker();
        Loop();
        CheckEntries(Nterm*Nevents);
        cout << "exit worker" << inodes << endl;
        ExitWorker(inodes); // turn off worker
    }
    cout << "Waiting" << endl;
    Wait(); // wait for all workers to finish
    cout << "all finished" << endl;
    Loop();
    total_entries+=Nterm*Nevents;
    CheckEntries(total_entries,total_entries_integ);
};

TEST_F(HistoServiceTest, Terminate){
    HistoHandler::Save();
    HistoHandler::Terminate();
};

bool IsIntegName(string hname){
    if ( hname == "s_qt"      ) return true;
    if ( hname == "s_qt_vs_y" ) return true;
    if ( hname == "s_qt_vs_y_vs_m" ) return true;
    if ( hname == "user_qt"   ) return true;
    return false;
}

#ifdef USEROOT
#include "TH1.h"
#include "TFile.h"


TEST_F(HistoServiceTest, CheckSavedFiles){
    // open root file or stl file
    String fname = HistoHandler::result_filename;
    fname+=HistoHandler::file_suffix;
    TFile *f = TFile::Open(fname.c_str());
    total_entries+=Nnodes*Nterm*Nevents;
    // for all histograms
    for (auto hname : hist_names){
        // check number of entries
        TH1 * hist = ((TH1*) f->Get(hname.c_str()));
        double hentries = hist->GetEntries();
        double exp = total_entries + (IsIntegName(hname)  ? total_entries_integ : 0) ;
        ASSERT_DOUBLE_EQ(exp, hentries)
            << "Saved histogram " << hname
            << " has incorect number of entries: " << hentries 
            << " instead of " << exp ;
    }
    f->Close();
};

#else  // STL
TEST_F(HistoServiceTest, CheckSavedFiles){
    total_entries+=Nnodes*Nterm*Nevents;
    ASSERT_TRUE(false) << "Not implemented for STL";
};
#endif

