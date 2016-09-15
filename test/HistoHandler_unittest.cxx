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

using namespace HistoHandler;
using namespace Kinematics;

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
        }

        void Loop(){
        }

        void Fork(){
        }

        void Wait(){
        }

        void InitWorker(){
        }

        void ExitWorker(){
        }

        double l1[4];
        double l2[4];
        double wgt;

        int Nterm;
        int Nevents;
        int Nnodes;

        std::mt19937 gen;
        std::normal_distribution<> gausT; // (0.,10.);
        std::normal_distribution<> gausZ; // (0.,30.);
        std::normal_distribution<> gausW; // (1.,0.1);

};

TEST_F(HistoServiceTest, Init){
    HistoHandler::Init();
};

TEST_F(HistoServiceTest, LoopNoThread){
    CheckEntries(0.);
    Loop();
    CheckEntries(Nterm*Nevents);
    HistoHandler::Save();
    CheckEntries(0.);
};

TEST_F(HistoServiceTest, LoopInsideThreads){
    CheckEntries(0.);
    for (int inodes = 0;  inodes < Nnodes; inodes++) {
        Fork();
        if(isParent) continue;
        // only child goes here
        InitWorker();
        Loop();
        CheckEntries(Nterm*Nevents);
        HistoHandler::Save(inode);
        CheckEntries(0.);
        ExitWorker();
    }
    Wait();
    HistoHandler::Save();
    Loop();
    HistoHandler::Save();
    CheckEntries(0.);
};

TEST_F(HistoServiceTest, Terminate){
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

