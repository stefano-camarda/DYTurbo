#ifndef dyturbo_H
#define dyturbo_H
/**
 * @file dyturbo.h
 * Steering class for Drell-Yan calculation.
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-17
 */

#include <iomanip>
using std::setw;
using std::setprecision;

#include "handy_typdefs.h"

/**
 * @brief Interface to properly control calculation.
 *
 */
namespace DYTurbo {

    /**
     * @brief Flag for if we need Integration mode.
     *
     *  If we have only Vegas integration we don't need to run per each bin and
     *  boundaries can be simplified. Only one bin from lowest to highest
     *  value.
     */
    extern bool HasOnlyVegas;

    //! Debug flag for testing code.
    extern bool isDryRun;


    //forward
    struct Term;
    struct Boundaries;
    struct BoundIterator;

    /** @defgroup MainFunctions Main DYTurbo functions.
     * @{
     */

    /**
     * @brief Initialization of program.
     *
     * It parse command line and input file and set interfaced values. It's
     * also calling initialization of submodules.
     */
    void Init(int argc, char * argv[]);
    /**
     * @brief Initialization of submodules.
     */
    void init_params();

    /**
     * @brief Warming up integration.
     *
     * Some of terms needed to be pre-run before actual integration. This done
     * inside this function. It also prints out current integration settings.
     */
    void WarmUp();

    /**
     * @brief Set integration bounds.
     *
     * The values from Boundaries are read and provided to \ref phasespace.
     *
     * @attention This needs to be called before calling the \ref RunIntegration.
     */
    void SetBounds(BoundIterator);

    /**
     * @brief Clean up and save results.
     */
    void Terminate();

    namespace PrintTable {
        void Settings();
        void Header() ;
        void Bounds(bool use_full_bound=false) ;
        void Result(const Term &term,bool printGrandTotal=false) ;
        void ResultSubTotal(bool is_grandtotal=false) ;
        void ResultGrandTotal() ;
        void Footer() ;
    }
     /** @} */

    /**
     * @brief Class for keeping information about active term.
     *
     * This class contains name, description and function to be called for
     * selected term. Terms are then stored in TermList and looped in main
     * program.
     */
    struct Term {
        String name = "dummy"; //!< Name of term will be used in header of table.
        String description = ""; //!< Description is used while calling \ref Print function.
        bool isVegas = true; //!< If term is not Vegas fill result.

        /**
         * @brief Pointer to integration function. This par
         *
         * This will be called within \ref RunIntegration function.
         */
        void (*integrate)(VecDbl &val,double &err) = NULL;

        double total_int = 0; //!< Total cross-section summed from all bins.
        double total_err2 = 0; //!< Squared uncertainties summed from all bins.
        double total_time = 0; //!< Total time spend for all bins.
        VecDbl last_int = {}; //!< Here is last cross-section stored for all variations (PDF).
        double last_err2 =0; //!< Uncertainty of last integration.
        double last_time =0; //!< Time spend in last iteration.

        /**
         * @brief Main functionality of this class.
         *
         * This will call integration, measure time for processing, store and
         * cumulate the result.
         */
        void RunIntegration();

        /**
         * @brief Print settings of this term to screen.
         *
         * This is used while printing \ref IntegrationSettings and format
         * follow 4 column mode.
         */
        void Print();

        //! Overloading stream operator to set description.
        template<class Streamable> Term & operator<<(const Streamable &data);

        //! Reset to zero values of last_int, last_err2 and last_time
        void last_reset();
    };

    //! Helper instance for summing results from all terms.
    extern Term subtotal;

    //! Container class for Terms
    typedef std::vector<Term> TermList;

    /**
     * @brief Helper instance storing currenlty active terms.
     *
     * This is where DYTurbo stores active terms. They are looped by \ref
     * TermIterator.
     */
    extern TermList ActiveTerms;

    /**
     * @brief Helper struct for looping active terms.
     *
     * Term Iterator always loops instance \ref ActiveTerms from beginning.
     */
    struct TermIterator {
        TermIterator(); //!< By construction always starts at begining of ActiveTerms.
        size_t icurrent; //!< Storing position in vector.
        bool IsEnd(); //!< Test iterator is in the end.
        TermIterator& operator++(); //! Increase operator simple increases \ref icurrent variable.
        Term & operator*(); //!< Dereferencing operator uses `vector::operator[]`
    };



    //! Container of boundary values for one variable.
    typedef VecDbl BoundaryValues;
    //! Iterator of boundary values for one variable.
    typedef BoundaryValues::iterator BoundValueItr;

    /**
     * @brief Holder of boundaries for one variable
     */
    struct Boundaries{
        std::string name =""; //!< Name of variable for printing purpose.
        BoundaryValues data = {}; //!< Actual values of boundaries for this variable.

        inline BoundValueItr begin() { return data.begin();}; //!< Iterator interface to data.
        inline BoundValueItr end()   { return data.end(); }; //!< Iterator interface to data.
        inline double front() { return data.front();}; //!< Iterator interface to data.
        inline double back()  { return data.back(); }; //!< Iterator interface to data.
        inline size_t size()  { return data.size(); }; //!< Iterator interface to data.
    };

    /** 
     * @brief Define position of variable inside \ref ActiveBoundaries
     *
     * When adding new variable please also add filling in \ref AddBoundaries
     * and interface to phasespace boundaries in \ref SetBounds.
     *
     * @note Trick is to define one extra enum item in the end for counting number of
     * defined items.
     *
     * @attention The `N_boundaries` must be always the last item.
     *
     * @todo b_Phi, b_LepPt
     */
    enum BoundaryIndex {
        b_M=0, //!< Vector boson invariant mass
        b_Y, //!< Vector boson rapidity
        b_QT, //!< Vector boson transverse momentum
        b_CsTh, //!< Cosine of longitudinal angle in Collins-Soper frame
        N_boundaries //! This here to count number of boundaries.
    };

    /**
     * @brief Container of all variable boundaries.
     *
     * Position of variable is set by BoundaryIndex. This is filled inside function \ref AddBoundaries.
     */
    typedef std::vector<Boundaries> BoundariesList;

    /**
     * @brief Helper instance containing active boundaries.
     *
     * This is looped by BoundIterator.
     */
    extern BoundariesList ActiveBoundaries;

    /**
     * @brief Container of boundary value iterators 
     */
    typedef std::vector<BoundValueItr> BoundValIterList;


    /**
     * @brief Iterating over all variable boundaries (bins) for integration.
     */
    struct BoundIterator {
        BoundIterator(); //!< After construction always point to first boundary.
        bool IsEnd(); //!< Check we already looped through all boundaries.

        /**
         * @brief Go to next kinematic boundary.
         *
         * This function assures that all combinations of boundaries are taken.
         * For more detail see implementation.
         */
        BoundIterator& operator++();

        //! Print to screen. Used by `PrintTable::Boundary`.
        void Print();

        /**
         * @brief Current lower boundary.
         * @param ib Index of variable (\ref BoundaryIndex can be used as input)
         * @return Lower boundary of selected variable.
         */
        inline double loBound(size_t ib){return *current[ib]     ;}

        /**
         * @brief Current upper boundary.
         * @param ib Index of variable (\ref BoundaryIndex can be used as input)
         * @return Upper boundary of selected variable.
         */
        inline double hiBound(size_t ib){return *(current[ib]+1) ;}

        /**
         * @brief List of pointers to all variable boundary values.
         */
        BoundValIterList current;

        /**
         * @brief Control flag.
         *
         * Set to true in beginning and after first increment is changed to
         * false.
         */
        bool isFirst;
    };
 
    //! Helper bounds: remember last bounds for printing.
    extern BoundIterator last_bounds;


    namespace PrintTable{
        //! Helper structs for printing: Four column
        const int colw1 = 25;
        const int colw2 = 20;
        const int colw3 = 15;
        const int colw4 = 12;

        struct Col4 {
            String data;
            template<class S1, class S2, class S3, class S4 > Col4( S1 col1, S2 col2, S3 col3, S4 col4, int prec=6){
                SStream tmp;
                tmp << setprecision(prec) << setw(colw1) << col1;
                tmp << setprecision(prec) << setw(colw2) << col2;
                tmp << setprecision(prec) << setw(colw3) << col3;
                tmp << setprecision(prec) << setw(colw4) << col4;
                tmp << '\n';
                data = tmp.str();
            }
        };
        inline std::ostream & operator<< (std::ostream & strm, const Col4 &col){ strm << col.data; return strm; }

        //! Helper structs for printing: Four column, skip first
        struct Col3{
            String data;
            template<class S2, class S3, class S4 > Col3(S2 col2, S3 col3, S4 col4){
                SStream tmp;
                tmp << setw(colw2) << col2;
                tmp << setw(colw3) << col3;
                tmp << setw(colw4) << col4;
                tmp << '\n';
                data = tmp.str();
            }
        };
        inline std::ostream & operator<< (std::ostream & strm, const Col3 &col){ strm << col.data; return strm; }
    }

}

#endif /* ifndef dyturbo_H */
