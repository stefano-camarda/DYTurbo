#ifndef TurboHist_Binning_H
#define TurboHist_Binning_H
/**
 * @file TurboHist_Binning.h
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */

#include "TurboHist.h"

#include <algorithm>
using std::lower_bound;
using std::sort;

namespace TurboHist {
    struct Binning {
        VecObs data; ///< Mapping value to index.
        size_t N= 0.; ///< Number of bins, without underflow and overflow.
        OBS min = 0.; ///< Low edge of first bin
        OBS max = 0.; ///< High edge of last bin
        bool isEquidistant=false; ///<Is equidistant flag.
        VecObsCItr BEG;
        VecObsCItr END;
        const double *aBEG = NULL; ///< trick std::lower_bound to think its array not vector

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
            /// @todo remove repeated bin edges
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
}

#endif /* ifndef TurboHist_Binning_H */
