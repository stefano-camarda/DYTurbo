#ifndef HistoHandler_H
#define HistoHandler_H

/**
 * @file HistoHandler.h
 * Description of this macro
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-27
 */

#include <vector>
#include <map>
#include <string>
#include <sstream>

// handy typedefs
typedef std::vector<double> VecDbl;
typedef std::string String;
typedef std::ostringstream SStream;

#include <sys/wait.h>
#include <unistd.h>



#define PARENT_PROC 32768

extern "C" {
    //! Interface to HistoHander: SetPDF
    void histo_setpdf(int *imember);
    //! Interface to HistoHander: FillEvent
    void histo_fill(double *l1,double *l2,double *wgt);
    //! Interface to HistoHander: FillDipole
    void histo_filldipole(double *l1,double *l2, double *wgt);
    //! Interface to HistoHander: FillRealEvent
    void histo_fillreal();
}

/** 
 * @brief Histogramming services
 */
namespace HistoHandler {
    // Typedef and forward declarations (implemented in HistoObjects.h)
    class HistoBase;
    //! Container of all histograms
    typedef std::vector<HistoBase*> HistoContainer;
    //! Container of all histograms
    extern HistoContainer histos;

    /**
     * @brief Switch between filling and integration mode.
     *
     *  * `true`  : _Filling mode_ is for Vegas integration.
     *  * `false` : _Integration mode_ is for quadrature integration.
     */
    extern bool isFillMode;

    /**
     * @brief Output file without suffix.
     *
     *  Suffix will be added according to used implementation (ROOT or STL).
     *
     *  @todo set from inputfile
     */
    extern String result_filename;

     /**
      * @brief Pid for main thread.
      *
      *  Saving the process ID for multi-thread calculation.
      */
    extern int parent_pid;

    /**
     * @defgroup HistoFun Main histogramming functions
     * @{
     */
    //! Initialization of histograming service.
    void Init();
    //! Remove data from histograms. Does not delete objects, just reset them.
    void Reset();
    //! Book general and user histograms histograms.
    void Book();
    //! Remove all histograms from container. It deletes all histograms.
    void DeleteHists();

    /**
     * @brief Change histogram to another variation.
     *
     *  Variations (like PDF) are treated by \ref HistoBase. Current
     *  implementation is matching index of PDF member to index of variation.
     *  In future there is possible to extend functionality to scale variations
     *  by using negative values or values larger than number of PDF members.
     *
     *  Suffix of histogram is stored in \ref variation_suffixes. This is
     *  container which maps integer value of variation (it can be negative) to
     *  suffix and non-negative index (used as vector index inside HistoBase).
     *
     * @param imember Index of variation.
     */
    void SetVariation(int imember);

    /**
     * @brief Standard filling of histograms with current kinematics.
     *
     * This is used by all Vegas calculations, except V+J NLO REAL term.
     *
     * @param p3 Four-momentum (x,y,z,E) of first lepton.
     * @param p4 Four-momentum (x,y,z,E) of second lepton.
     * @param wgt Event weight.
     */
    void FillEvent(double p3[4],double p4[4], double wgt);

    /**
     * @brief Fill result of integration.
     *
     * This is used in integration mode to store result of integration inside
     * appropriate bin. Bin is found according to phasespace boundaries.
     *
     * @param int_val The calculated value of integration.
     * @param int_err Uncertainty of integral.
     */
    void FillResult(double int_val,double int_err);

    /**
     * @brief Save dipole kinematics.
     *
     * This is used by Vegas integration of V+J NLO REAL term. The integration
     * produces per one random phasespace point several kinematic results with
     * different weights. To get proper uncertainty we need to treat weight
     * properly i.e.:
     *
     *  > The weights corresponding to same bin are summed up and filled as one
     *  > event. Contribution to different bins are treated as statistically
     *  > uncorrelated (summed in squares)
     *
     * @attention This function does not fill histogram, just store temporary
     * values. To really fill histograms it is necessary to call \ref
     * FillRealEvent after adding all dipoles.
     *
     * @param p3 Four-momentum (x,y,z,E) of first lepton.
     * @param p4 Four-momentum (x,y,z,E) of second lepton.
     * @param wgt Event weight.
     */
    void FillDipole(double p3[4],double p4[4], double wgt);

    /**
     * @brief Fill dipoles properly to histogram.
     *
     * This is called after all dipoles for current event are filled. This
     * functions decides which dipoles are filled as same event and which are
     * filled as independent event. For more info \ref FillDipole.
     *
     * After this all dipoles are removed and histograms are filled.
     */
    void FillRealEvent();


    /**
     * @brief File extension (suffix) for output file.
     *
     * Depends on implementation
     *  - _ROOT_ : `root`
     *  - _STL : `dat`
     */
    extern String file_suffix;

    /** 
     * @brief Save histograms to file.
     *
     *  Always saves as temporary file. To create final call Merge (via Terminate).
     */
    void Save(int iworker=PARENT_PROC);

    /** 
     * @brief Close and clean histogramming services.
     */
    void Terminate(int iworker=PARENT_PROC);


    //! Add histogram to the list.
    template<class T>
    void Add(T* newhist){
        histos.push_back(newhist);
    }
    //! @}

    //! Variation suffix struct.
    struct KeySuffix{
        String suffix; //!< Suffix added to histogram name
        size_t index; //!< index for HistoBase
    };
    //! Container matching integer index (can be negative) to variation suffix.
    typedef std::map<int,const KeySuffix> MapSuffix;

    //! Index of last added variation
    extern size_t last_index;
    //! All variation suffixes (so far only implemented for PDF)
    extern MapSuffix variation_suffixes;

}



#endif /* ifndef HistoHandler_H */
