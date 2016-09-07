#ifndef TurboHist_H
#define TurboHist_H
/**
 * @file TurboHist.h
 * Simple histograming with STL
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-01
 */

#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>

using std::fstream;
using std::string;
using std::map;
using std::sort;
using std::vector;
using std::lower_bound;

#include <iostream>
using std::cout;
using std::endl;

namespace TurboHist {
    // OBS: observable ususaly double but I think float is fine
    typedef double OBS;
    // PRE: precission of data -- here we need doubles
    typedef double PRE;

    typedef vector<string> VecStr;
    typedef vector<OBS> VecObs;
    typedef VecObs::iterator VecObsItr;
    typedef VecObs::const_iterator VecObsCItr;
    typedef vector<PRE> VecPre;

    enum AxisName {
        X=0, Y=1, Z=2
    };

    template<class HTYPE> class HBase;
    class Counter;
    class Binning;
    typedef vector<Counter> Data;


    void MergeFiles(string outname,VecStr inNames);

    class File{
        public :
            void Open(string name, string method="READ");
            inline void Close(){fst.close();};
            string fname;

            File& operator << (const int     & rhs) { return SaveBase    (rhs); };
            File& operator << (const size_t  & rhs) { return SaveBase    (rhs); };
            File& operator << (const char    & rhs) { return SaveBase    (rhs); };
            File& operator << (const double  & rhs) { return SaveBase    (rhs); };
            File& operator << (const string  & rhs) { return SaveVector  (rhs); };
            File& operator << (const VecObs  & rhs) { return SaveVector  (rhs); };
            File& operator << (const Data    & rhs) { return SaveVector  (rhs); };
            File& operator << (const VecStr  & rhs) { return SaveVector  (rhs); };
            File& operator << (const Counter & rhs) { return SaveReverse (rhs); };
            File& operator << (const Binning & rhs) { return SaveReverse (rhs); };
            template<class HTYPE> File& operator << (const HBase<HTYPE> & rhs) { return SaveReverse(rhs); };

            File& operator >> (int     & rhs) { return LoadBase    (rhs); };
            File& operator >> (size_t  & rhs) { return LoadBase    (rhs); };
            File& operator >> (char    & rhs) { return LoadBase    (rhs); };
            File& operator >> (double  & rhs) { return LoadBase    (rhs); };
            File& operator >> (string  & rhs) { return LoadVector  (rhs); };
            File& operator >> (VecObs  & rhs) { return LoadVector  (rhs); };
            File& operator >> (Data    & rhs) { return LoadVector  (rhs); };
            File& operator >> (VecStr  & rhs) { return LoadVector  (rhs); };
            File& operator >> (Counter & rhs) { return LoadReverse (rhs); };
            //File& operator >> (Binning & rhs) { return ReadReverse (rhs); }; // use histogram to read binning
            template<class HTYPE> File& operator >> (HBase<HTYPE> & rhs) { return LoadReverse(rhs); };

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

            template<class ItemType> File & LoadVector(vector<ItemType> & rhs ){
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

            File & LoadVector(string & rhs ){
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



    struct Binning {
        VecObs data; ///< Mapping value to index.
        size_t N; ///< Number of bins, without underflow and overflow.
        OBS min; ///< Low edge of first bin
        OBS max; ///< High edge of last bin
        bool isEquidistant=false; ///<Is equidistant flag.
        VecObsCItr BEG;
        VecObsCItr END;
        const double *aBEG; ///< trick std::lower_bound to think its array not vector

        inline size_t     size()   const{ return data .size()   ;}
        inline VecObsItr  begin()  {      return data .begin()  ;}
        inline VecObsItr  end()    {      return data .end()    ;}
        inline VecObsCItr cbegin() const{ return data .cbegin() ;}
        inline VecObsCItr cend()   const{ return data .cend()   ;}

        bool operator!=(const Binning &rhs) const {
            // mismatch returns first non-equal pair if all items match it
            // returns pair ot iterators pointing to end() of containers
            return END != mismatch(BEG,END,rhs.BEG).first;
        }

        void Print() const {
            printf(" nbins: %zu min: %f max: %f isEquidistant: %s \n",N,min,max, isEquidistant ? "True" : "False");
            for (VecObsCItr itr=BEG;itr!=END;++itr) printf("%f : %zu\n", (*itr), itr-BEG);
            printf("\n");
        }

        void SetBins(VecObs newbins){
            data.clear();
            sort(newbins.begin(),newbins.end());
            // TODO: remove Repeat
            min = newbins.front();
            max = newbins.back();
            //
            for (size_t ilo = 0; ilo < newbins.size(); ++ilo) {
                OBS edge = newbins[ilo];
                if (count(begin(),end(),edge)==0) data.push_back(edge); // check for duplication
            }
            N = data.size()-1; // high edge
            BEG = cbegin();
            END = cend();
            aBEG = &data[0];
            // test equidistancy
            VecObs equid = GetEquidistantVector(N,min,max);
            isEquidistant = (END == mismatch(BEG,END,equid.begin()) .first);
        }

        inline size_t FindBin(const OBS &val) const {
            if (val < min) return 0;
            if (val >= max) return N+1;
            if (isEquidistant) return 1 +  size_t( N* (val-min)/(max-min) );
            else return lower_bound(BEG,END,val)-BEG;
        }


        static VecObs GetEquidistantVector(const size_t &N, const OBS &min, const OBS &max){
            VecObs newbins;
            OBS step = (max-min)/ OBS(N);
            for (OBS bine = min; bine <= max; bine+=step) newbins.push_back(bine);
            return newbins;
        }

        void SaveToFile(File &file) const {
            file << data;
        }

    };
    typedef vector<Binning> VecBins;

    struct Counter{
        PRE sum_w;
        PRE sum_w2;
        Counter() : sum_w(0.), sum_w2(0.) { };

        Counter(const Counter& rhs) : sum_w(rhs.sum_w), sum_w2(rhs.sum_w2) { };

        inline Counter operator = (const Counter &rhs){
            sum_w=rhs.sum_w;
            sum_w2=rhs.sum_w2;
            return *this;
        };

        inline Counter operator += (const PRE &weight){
            sum_w+=weight;
            sum_w2+=weight*weight;
            return *this;
        };

        inline Counter & operator += (const Counter &rhs){
            sum_w+=rhs.sum_w;
            sum_w2+=rhs.sum_w2;
            return *this;
        };

        inline Counter & operator *= (const double &weight){
            sum_w*=weight;
            sum_w2*=weight*weight;
            return *this;
        };

        void SaveToFile(File &file) const {
            file << sum_w;
            file << sum_w2;
        }

        void LoadFromFile(File &file){
            file >> sum_w;
            file >> sum_w2;
        }
    };
    Counter operator* ( const Counter &a , const double &b );
    Counter operator* ( const double &b, const Counter &a  );


    template<class HTYPE>
    class HBase {
        public :

            inline int GetDim() const { return dim;};
            inline char GetType() const { return type;};

            inline const bool IsEquidistant(size_t iaxis=X) const {
                return GetBinning(iaxis).isEquidistant;
            }

            inline size_t GetBin(size_t ibinX, size_t ibinY=0, size_t ibinZ=0) const {
                // column major order: hard-code per each case should speed up
                size_t ibin = 0;
                if (dim==1) {
                    return ibinX;
                }
                if (dim==2) {
                    ibin = ibinY;
                    return ibinX + (binsX.N+2)*ibin;
                }
                if (dim==3) {
                    ibin = ibinZ;
                    ibin = ibinY + (binsY.N+2)*ibin;
                    return ibinX + (binsX.N+2)*ibin;
                }
                return ibin;
            };

            inline PRE GetBinContent (size_t ibin) const { return d[ibin].sum_w; };
            inline PRE GetBinError   (size_t ibin) const { return sqrt(d[ibin].sum_w2); };


            void SaveToFile(File &file) const {
                file << type;
                file << dim;
                file << binsX;
                file << binsY;
                file << binsZ;
                file << data;
                file << name;
                file << title;
                file << '\n';
            }

            bool LoadFromFile(File &file, bool doOnlyCheck=false){
                size_t ptr = file.tellg();
                bool isGoodSave = true;
                //read and check type and dim
                char _type;
                int  _dim;
                VecObs _ax;
                VecObs _ay;
                VecObs _az;
                Data _data;
                string _name;
                VecStr _titles;
                // Check this should be same order as in operator write
                file >> _type;
                file >> _dim;
                file >> _ax;
                file >> _ay;
                file >> _az;
                file >> _data;
                file >> _name;
                file >> _titles;
                //
                // Tests
                try {
                    isGoodSave&= _dim==dim;
                    isGoodSave&= _type==type;
                    if (!isGoodSave) throw false;
                    if (dim>0) isGoodSave&= _ax.size()>0;
                    if (dim>1) isGoodSave&= _ay.size()>0;
                    if (dim>2) isGoodSave&= _az.size()>0;
                    if (!isGoodSave) throw false;
                    isGoodSave&= int(_titles.size())==dim;
                    if (!isGoodSave) throw false;
                    if (!doOnlyCheck){
                        //set bin points
                        if (dim >0) SetBinsAxis(X,_ax,_titles[X]);
                        if (dim >1) SetBinsAxis(Y,_ay,_titles[Y]);
                        if (dim >2) SetBinsAxis(Z,_az,_titles[Z]);
                        //set data
                        for (size_t i=0; i<GetMaxBin(); i++){
                            d[i]=_data[i];
                        }
                    }
                } catch (bool) { };
                // it was only check, put read pointer back where we started
                if (doOnlyCheck) file.seekg(ptr);
                return isGoodSave;
            }

            const char * Print(string opt="") const {
                printf("name '%s' tiltles : '", name.c_str());
                for (auto tit : title ) printf(";%s",tit.c_str());
                printf("'\n");
                if (dim>0) binsX.Print();
                if (dim>1) binsY.Print();
                if (dim>2) binsZ.Print();
                size_t ibin=0;
                for (auto value : data) printf ("%zu: %f +- %f \n",ibin++, value.sum_w , sqrt(value.sum_w2));
                return "\n";
            }

            void Add(const HBase<HTYPE> &rhs , double c=1. ){
                // check dim
                if (dim!=rhs.dim) return;
                if (type!=rhs.type) return;
                // check binning
                if(binsX!=rhs.binsX) return;
                if(binsY!=rhs.binsY) return;
                if(binsZ!=rhs.binsZ) return;
                // loop over all data
                for (size_t ibin = 0; ibin < data.size(); ++ibin) {
                    d[ibin] += rhs.d[ibin]*c;
                }
            }

            inline size_t GetMaxBin() const {
                size_t    ibin =  binsX.N+2;
                if(dim>1) ibin *= binsY.N+2;
                if(dim>2) ibin *= binsZ.N+2;
                return ibin;
            }


        protected :

            inline Binning &GetBinning(size_t iaxis){
                switch (iaxis) {
                    case X:
                        return binsX;
                    case Y:
                        return binsY;
                    case Z:
                        return binsZ;
                    default:
                        return binsX;
                }
            }

            inline const Binning &GetBinning(size_t iaxis) const{
                switch (iaxis) {
                    case X:
                        return binsX;
                    case Y:
                        return binsY;
                    case Z:
                        return binsZ;
                    default:
                        return binsX;
                }
            }

            void SetBinsAxis(AxisName iaxis, const VecObs &newbins, const string &tit){
                title .resize(dim);
                title[iaxis]=tit;
                //
                GetBinning(iaxis).SetBins(newbins);
                data.assign( GetMaxBin() , Counter());
                //data.assign( GetMaxBin() , 0.);
                d=&data[0];
            };

            static int dim;
            static char type; // new to write file
            string name;
            VecStr title; // size dim
            Binning binsX; // bin edges size: (nbinsX+2)
            Binning binsY; // bin edges size: (nbinsY+2)
            Binning binsZ; // bin edges size: (nbinsZ+2)
            Data data; // bin content (nbinX+2)*(nbinX+Y)*(nbinZ+2)
            Counter* d; // bin content -- faster access
    };


    // Specializations
    struct H1 : public  HBase<H1> {
        H1() {SetBins(1,0.,1.);};

        H1(const string &_name, const string &_title,const size_t &N, const OBS &min, const OBS &max){
            name = _name;
            SetBins(N,min,max,_title);
        };

        inline size_t FindBin(const OBS &val){
            return binsX.FindBin(val);
        };

        void SetBins(size_t N, OBS min, OBS max, const string &tit="X"){
            VecObs newbins=Binning::GetEquidistantVector(N,min,max);
            SetBinsAxis(X,newbins,tit);
        };

        void SetBins(const VecObs &newbins, const string &tit="X"){
            SetBinsAxis(X,newbins,tit);
        };

        void Fill (const OBS &val,const PRE &weight=1.0){
            d[FindBin(val)]+=weight;
        };
    };



    //struct H2;
    //struct H3;
    //struct Prof1;
    //struct Prof2;
};

#endif /* ifndef TurboHist_H */
