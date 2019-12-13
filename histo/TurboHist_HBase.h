#ifndef TurboHist_HBase_H
#define TurboHist_HBase_H
/**
 * @file TurboHist_HBase.h
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */

#include "TurboHist.h"
#include "TurboHist_File.h"
#include "TurboHist_Binning.h"


#include <cmath>

namespace TurboHist {

    template<class HistoType, class CountType>
        class HBase : public OBase {
            public :
                typedef Vec<CountType> Data;

                inline int GetDim() const { return dim;};
                inline char GetType() const { return type;};

                inline const bool IsEquidistant(const size_t iaxis=X) const {
                    return GetBinning(iaxis).isEquidistant;
                }

                OBase * Clone(String name= "") const {
                    HBase *nh = new HBase(*this);
                    if (name.size()!=0) nh->SetName(name);
                    return nh;
                };

                void Reset(){
                    for (size_t i=0; i<GetMaxBin(); i++) data[i].reset(); 
                    entries=0;
                }

                inline size_t GetBin(size_t ibinX, size_t ibinY) const {
                    // column major order: hard-code per each case should speed up
                    size_t ibin = ibinY;
                    return ibinX + (binsX.N+2)*ibin;
                };
                inline size_t GetBin(size_t ibinX, size_t ibinY, size_t ibinZ) const {
                    // column major order: hard-code per each case should speed up
                    size_t ibin = ibinZ;
                    ibin = ibinY + (binsY.N+2)*ibin;
                    return ibinX + (binsX.N+2)*ibin;
                };

                inline PRE GetBinContent (size_t ibin) const { return data[ibin].value(); };
                inline PRE GetBinError   (size_t ibin) const { return data[ibin].error(); };
                inline size_t GetEntries () const { return entries; };
                inline const char * GetName () const { return name.c_str(); };

                inline void SetBinContent (size_t ibin, PRE val) { data[ibin].value(val); entries++; };
                inline void SetBinError   (size_t ibin, PRE err) { data[ibin].error(err); };


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
                    String _name;
                    String _title;
                    // Check this should be same order as in operator write
                    file >> _type;
                    file >> _dim;
                    file >> _ax;
                    file >> _ay;
                    file >> _az;
                    file >> _data;
                    file >> _name;
                    file >> _title;
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
                        if (!isGoodSave) throw false;
                        if (!doOnlyCheck){
                            //set bin points
                            if (dim >0) SetBinsAxis(X,_ax,_title);
                            if (dim >1) SetBinsAxis(Y,_ay,_title);
                            if (dim >2) SetBinsAxis(Z,_az,_title);
                            //set data
                            for (size_t i=0; i<GetMaxBin(); i++){
                                data[i]=_data[i];
                            }
                        }
                    } catch (bool) { };
                    // it was only check, put read pointer back where we started
                    if (doOnlyCheck) file.seekg(ptr);
                    return isGoodSave;
                }

                const char * Print(String opt="") const {
                    printf("name '%s' titles '%s' \n", name.c_str(), title.c_str());
                    if (dim>0) binsX.Print();
                    if (dim>1) binsY.Print();
                    if (dim>2) binsZ.Print();
                    size_t ibin=0;
                    for (auto value : data) printf ("%zu: %f +- %f \n",ibin++, value.sum_w , sqrt(value.sum_w2));
                    return "\n";
                }

                void Add(const HBase<HistoType,CountType> &rhs , double c=1. ){
                    // check dim
                    if (dim!=rhs.dim) return;
                    if (type!=rhs.type) return;
                    // check binning
                    if(binsX!=rhs.binsX) return;
                    if(binsY!=rhs.binsY) return;
                    if(binsZ!=rhs.binsZ) return;
                    // loop over all data
                    for (size_t ibin = 0; ibin < data.size(); ++ibin) {
                        data[ibin] += rhs.data[ibin]*c;
                    }
                }

                inline size_t GetMaxBin() const {
                    size_t    ibin =  binsX.N+2;
                    if(dim>1) ibin *= binsY.N+2;
                    if(dim>2) ibin *= binsZ.N+2;
                    return ibin;
                }

                inline void SetName(String _name){ name=_name; }

                void SetBinsAxis(AxisName iaxis, const VecObs &newbins, const String &tit){
                    title=tit;
                    //
                    GetBinning(iaxis).SetBins(newbins);
                    data.assign( GetMaxBin() , CountType());
                };


            protected :
                // Static variables will be shared per each template specification
                static int dim; // dimension (number of axis)
                static char type; // type of counter
                String name;
                String title; // title of axes divided by semicolon
                Binning binsX; // bin edges size: (nbinsX+2)
                Binning binsY; // bin edges size: (nbinsY+2)
                Binning binsZ; // bin edges size: (nbinsZ+2)
                Data data; // bin content (nbinX+2)*(nbinX+Y)*(nbinZ+2)
                size_t entries=0; // number off called fill or set bin content

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
        };
}

#endif /* ifndef TurboHist_HBase_H */
