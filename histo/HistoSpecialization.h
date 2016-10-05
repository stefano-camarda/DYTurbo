#ifndef HistoSpecialization_H
#define HistoSpecialization_H
/**
 * @file HistoSpecialization.h
 * @brief A brief description
 *
 * Detailed description of file
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-10-05
 */

#include "HistoBase.h"

namespace HistoHandler {

    // specialization
    template<class TX> class Histo1D : public HistoImp<H1> {
        private:
            TX varX;

        public:
            Histo1D(const String &bin_name_X){
                // general init
                Init(bin_name_X);
                name = "s_"+name;
                title += ";#sigma[fb]";
                // create new
                current = New<H1>(this);
                isIntegrationSafe = varX.IsIntegrableObservable();
            }

            virtual void FillEvent(){
                current->Fill(varX(),Kinematics::event_weight);
            };

            virtual void FillDipole(){
                current_point.valX = varX();
                AddPoint();
            };

            virtual int CurrentBin(){
                return current->FindBin(varX());
            };

            void FillPoint(DipPt point){
                current->Fill(point.valX,point.weight);
            }
    };

    template<class TX, class TY> class Histo2D : public HistoImp<H2> {
        public :
            Histo2D(const String &bin_name_X, const String &bin_name_Y){
                // general init
                Init(bin_name_X,bin_name_Y);
                name = "s_"+name;
                title += ";#sigma[fb]";
                // create new
                current = New<H2>(this);
                isIntegrationSafe = varX.IsIntegrableObservable() && varY.IsIntegrableObservable();
            }

            virtual void FillEvent(){
                current->Fill(varX(),varY(),Kinematics::event_weight);
            };

            virtual void FillDipole(){
                current_point.valX = varX();
                current_point.valY = varY();
                AddPoint();
            };

            virtual int CurrentBin(){
                return current->FindBin(varX(),varY());
            };

            void FillPoint(DipPt point){
                current->Fill(point.valX, point.valY ,point.weight);
            }


        private:
            TX varX;
            TY varY;
    };

    template<class TX, class TY, class TZ> class Histo3D : public HistoImp<H3> {
        public :
            Histo3D(const String &bin_name_X, const String &bin_name_Y, const String &bin_name_Z){
                // general init
                Init(bin_name_X,bin_name_Y,bin_name_Z);
                name = "s_"+name;
                // create new
                current = New<H3>(this);
                isIntegrationSafe = varX.IsIntegrableObservable() && varY.IsIntegrableObservable() && varZ.IsIntegrableObservable();
            }

            virtual void FillEvent(){
                current->Fill(varX(),varY(), varZ(), Kinematics::event_weight);
            };

            virtual void FillDipole(){
                current_point.valX = varX();
                current_point.valY = varY();
                current_point.valZ = varZ();
                AddPoint();
            };

            virtual int CurrentBin(){
                return current->FindBin(varX(),varY(),varX());
            };

            void FillPoint(DipPt point){
                current->Fill(point.valX, point.valY, point.valZ, point.weight);
            }


        private:
            TX varX;
            TY varY;
            TZ varZ;
    };


    template<class TX, class TY> class HistoProfile : public HistoImp<P1> {
        public :
            HistoProfile(const String &bin_name_X,const String &bin_name_Y){
                // general init
                Init(bin_name_X);
                // change name, add axis title
                title+=";";
                title+=bin_name_Y;
                name=bin_name_Y+"_"+name;
                // create new
                current = New<P1>(this);
                isIntegrationSafe = varX.IsIntegrableObservable() && varY.IsIntegrableObservable();
            }

            virtual void FillEvent(){
                current->Fill(varX(),varY(),Kinematics::event_weight);
            };

            virtual void FillDipole(){
                // Since the Ai moments are actualy weighted mean we need to do
                // weighted mean per each dipole point. We started by storing
                // the profiled value times weight.
                current_point.valX = varX();
                current_point.valY = varY()*Kinematics::event_weight;
                AddPoint();
            };

            virtual int CurrentBin(){
                return current->FindBin(varX());
            };

            void FillPoint(DipPt point){
                if (point.weight!=0) current->Fill(point.valX, point.valY/point.weight ,point.weight);
            }


        private:
            TX varX;
            TY varY;
    };

    template<class TX, class TY, class TZ> class HistoProfile2D : public HistoImp<P2> {
        public :
            HistoProfile2D(const String &bin_name_X,const String &bin_name_Y, const String &bin_name_Z){
                // general init
                Init(bin_name_X, bin_name_Y);
                // change name, add axis title
                title+=";";
                title+=bin_name_Z;
                name=bin_name_Z+"_"+name;
                // create new
                current = New<P2>(this);
                isIntegrationSafe = varX.IsIntegrableObservable() && varY.IsIntegrableObservable() && varZ.IsIntegrableObservable();
            }

            virtual void FillEvent(){
                current->Fill(varX(),varY(),varZ(),Kinematics::event_weight);
            };

            virtual void FillDipole(){
                current_point.valX = varX();
                current_point.valY = varY();
                current_point.valZ = varZ()*Kinematics::event_weight;
                AddPoint();
            };

            virtual int CurrentBin(){
                return current->FindBin(varX(),varY());
            };

            void FillPoint(DipPt point){
                if (point.weight!=0) current->Fill(point.valX, point.valY, point.valZ/point.weight ,point.weight);
            }


        private:
            TX varX;
            TY varY;
            TZ varZ;
    };
}

#endif /* ifndef HistoSpecialization_H */
