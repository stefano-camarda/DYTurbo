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

#include "HistoHandler.h"
#include "HistoObjects.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <random>


class HistoServiceTest : public testing::Test{
    protected :

        virtual void SetUp(){
            gen.seed(12345);
            gausT.param( std::normal_distribution<>::param_type(0.,10.)  );
            gausZ.param( std::normal_distribution<>::param_type(0.,30.)  );
            gausW.param( std::normal_distribution<>::param_type(1.,0.1) );

            Nterm=1;
            Nevents=1;
            Nnodes=1;
            isParent=true;
            inode=0;
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

        void CheckEntries(double entries=0.){
            for (auto h_it = HistoHandler::hists.begin();h_it!=HistoHandler::hists.end();h_it++){
                double hentries = (*h_it)->GetEntries();
                ASSERT_DOUBLE_EQ(entries, (*h_it)->GetEntries())
                    << "Histogram " << (*h_it)->GetName()
                    << " has incorect number of entries: " << hentries 
                    << " instead of " << entries ;
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
            HistoHandler::Clear();
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
        int inode;
        static double total_entries;


        int child_pid;
        int parent_pid;

        std::mt19937 gen;
        std::normal_distribution<> gausT; // (0.,10.);
        std::normal_distribution<> gausZ; // (0.,30.);
        std::normal_distribution<> gausW; // (1.,0.1);
};

double HistoServiceTest::total_entries = 0.;

TEST_F(HistoServiceTest, HistCreateAndFill){
    HistoHandler::Histo1D<BosPT> htmp("dummy");
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
};

TEST_F(HistoServiceTest, LoopNoThread){
    Loop();
    total_entries+=Nterm*Nevents;
    CheckEntries(total_entries);
};

TEST_F(HistoServiceTest, LoopNoThreadIntegrator){
    LoopIntegrator();
    total_entries+=1;
    CheckEntries(total_entries);
}

TEST_F(HistoServiceTest, LoopInsideThreads){
    for (int inodes = 0;  inodes < Nnodes; inodes++) {
        Fork();
        if(IsParent()) continue;
        // only child goes here
        InitWorker();
        Loop();
        CheckEntries(Nterm*Nevents);
        ExitWorker(inode); // turn off worker
    }
    Wait(); // wait for all workers to finish
    Loop();
    total_entries+=Nterm*Nevents;
    CheckEntries(total_entries);
    HistoHandler::Save(); // parent always saves its instance and check for temporary files
};

TEST_F(HistoServiceTest, Terminate){
    HistoHandler::Save();
    HistoHandler::Terminate();
};


// TODO: fork code and merge back
// TODO: why 2 histograms

// testing HistoSvc functionality
/*
void Test_HistoSvc(){

    const int parent_pid=getpid();

    bool usePrallel=true;
    int this_pid=parent_pid;
    int child_pid=0;

    double l1[4] = {0.3,0.3,0.3,0.33};
    double l2[4] = {0.3,0.3,0.3,0.33};
    double wgt=1;
    double* pwgt=&wgt;
    printf ("\n\nPID: %d Testing Kinematics\n",this_pid); fflush(stdout);
    BosPX var_px1;
    BosPX var_px2;
    BosPT var_qt1;
    BosPT var_qt2;
    SetKinematics(l1,l2,wgt);
    printf( "PID: %d px1 %f is same as  px2 %f; pt1 %f is same as pt2 %f\n",this_pid,var_px1(),var_px2(),var_qt1(),var_qt2()); fflush(stdout);

    // init histograms
    Init();
    int ipdf = 2;

    // simulate different term calculation
    for (int iterm = 0; iterm < 2; ++iterm) {
        int Nnodes = 3;
        // from http://stackoverflow.com/a/16890759
        for (int inode = 0; inode < Nnodes; ++inode) {
            if (usePrallel){
                if (this_pid==parent_pid){
                    printf("forking %d \n",this_pid); fflush(stdout);
                    child_pid = fork();
                    this_pid = getpid();
                    printf("I am %d \n",this_pid); fflush(stdout);
                }
                // parent skip and fork again ( double condition just to be sure)
                if ((child_pid!=0) || (this_pid==parent_pid)) continue;
            }

            printf ("\n\nTesting Looping pid %d inode %d \n",this_pid, inode); fflush(stdout);
            for (size_t ievent=0; ievent < 3; ievent++){
                printf ("PID %d Event %zu\n",this_pid,ievent); fflush(stdout);
                random_kin(l1,l2,wgt);
                //
                printf ("PID %d SetPDF0 and fill %zu\n",this_pid,ievent); fflush(stdout);
                SetVariation(0);
                histo_fill(l1,l2,pwgt);
                //
                printf ("PID %d SetPDF2 and filldipole %zu\n",this_pid,ievent); fflush(stdout);
                histo_setpdf(&ipdf);
                histo_filldipole(l1,l2,pwgt);
                wgt+=0.3;
                printf ("PID %d filldipole2 %zu\n",this_pid,ievent); fflush(stdout);
                histo_filldipole(l1,l2,pwgt);
                histo_fillreal();
            }

            if (usePrallel){
                // save worker
                printf ("\n\nPID %d Testing Save worker\n",this_pid); fflush(stdout);
                Save(inode);
                printf("child %d is dying... bye bye\n",this_pid); fflush(stdout);
                return;  // only child should go here and then die
            }
        } // end of parallel
        if (usePrallel){
            // from: http://stackoverflow.com/a/279761
            printf("parent is waiting...\n"); fflush(stdout);
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
        // save tmp and parent thread
        Save();
    }
    // save main
    printf ("\n\nTesting Save\n"); fflush(stdout);
    SetVariation(0);
    histo_fill(l1,l2,pwgt);
    Save();
    Terminate();
    return;
}
*/

