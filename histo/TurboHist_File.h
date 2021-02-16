#ifndef TurboHist_File_H
#define TurboHist_File_H
/**
 * @file TurboHist_File.h
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */

#include "TurboHist.h"

#include <fstream>
using std::fstream;

namespace TurboHist{
    void MergeFiles(String outname,VecStr inNames);

    class File{
        public :
            void Open(String name, String method="READ");
            inline void Close(){fst.close();};
            String fname;

            File& operator << (const int           & rhs) { return SaveBase    (rhs); };
            File& operator << (const size_t        & rhs) { return SaveBase    (rhs); };
            File& operator << (const char          & rhs) { return SaveBase    (rhs); };
            File& operator << (const double        & rhs) { return SaveBase    (rhs); };
            File& operator << (const String        & rhs) { return SaveVector  (rhs); };
            File& operator << (const VecObs        & rhs) { return SaveVector  (rhs); };
            File& operator << (const VecStr        & rhs) { return SaveVector  (rhs); };
            File& operator << (const Vec<Counter>  & rhs) { return SaveVector  (rhs); };
            File& operator << (const Vec<Averager> & rhs) { return SaveVector  (rhs); };
            File& operator << (const Counter       & rhs) { return SaveReverse (rhs); };
            File& operator << (const Averager      & rhs) { return SaveReverse (rhs); };
            File& operator << (const Binning       & rhs) { return SaveReverse (rhs); };
            template<class H,class C> File& operator << (const HBase<H,C> & rhs) { return SaveReverse(rhs); };

            File& operator >> (int           & rhs) { return LoadBase    (rhs); };
            File& operator >> (size_t        & rhs) { return LoadBase    (rhs); };
            File& operator >> (char          & rhs) { return LoadBase    (rhs); };
            File& operator >> (double        & rhs) { return LoadBase    (rhs); };
            File& operator >> (String        & rhs) { return LoadVector  (rhs); };
            File& operator >> (VecObs        & rhs) { return LoadVector  (rhs); };
            File& operator >> (Vec<Counter>  & rhs) { return LoadVector  (rhs); };
            File& operator >> (Vec<Averager> & rhs) { return LoadVector  (rhs); };
            File& operator >> (VecStr        & rhs) { return LoadVector  (rhs); };
            File& operator >> (Counter       & rhs) { return LoadReverse (rhs); };
            File& operator >> (Averager      & rhs) { return LoadReverse (rhs); };
            //File& operator >> (Binning & rhs) { return ReadReverse (rhs); }; // use histogram to read binning
            template<class H, class C> File& operator >> (HBase<H,C> & rhs) { return LoadReverse(rhs); };

            inline size_t tellg () {return fst.tellg();}
            inline void   seekg (size_t pos) {fst.seekg(pos);}

            fstream fst;

        private:
            template<class OrdType> File & SaveBase(const OrdType & rhs ){
                fst << rhs;
                fst << ' ';
                return (*this);
            }

            template<class VecType> File & SaveVector(const VecType & rhs ){
                (*this) << rhs.size();
                for (const auto item : rhs ) (*this) << item;
                return (*this);
            }

            template<class MyClass> File & SaveReverse(const MyClass & rhs ){
                rhs.SaveToFile(*this);
                return (*this);
            }

            template<class OrdType> File & LoadBase(OrdType & rhs ){
                fst >> rhs;
                return (*this);
            }

            template<class ItemType> File & LoadVector(Vec<ItemType> & rhs ){
                size_t N;
                rhs.clear();
                (*this) >> N;
                for (size_t i = 0; i < N; ++i) {
                    ItemType item;
                    (*this) >> item;
                    rhs.push_back(item);
                }
                return (*this);
            }

            File & LoadVector(String & rhs ){
                size_t N;
                rhs.clear();
                (*this) >> N;
                for (size_t i = 0; i < N; ++i) {
                    char item;
                    (*this) >> item;
                    rhs.push_back(item);
                }
                return (*this);
            }

            template<class MyClass> File & LoadReverse(MyClass & rhs ){
                rhs.LoadFromFile(*this);
                return (*this);
            }
    };
}

#endif /* ifndef TurboHist_File_H */
