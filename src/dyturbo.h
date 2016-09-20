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

#include <vector>
#include <string>
#include <sstream>

namespace DYTurbo {
    extern bool HasOnlyVegas;

    struct Term;
    struct TermIterator;
    typedef std::vector<Term> TermList;
    typedef std::vector<double> VecDbl;
    typedef std::string String;
    typedef std::stringstream SStream;
    struct Term {
        String name = "dummy";
        String description = "";
        void (*integrate)(VecDbl &val,double &err) = NULL;
        double total_int = 0;
        double total_err2 = 0;
        double total_time = 0;
        VecDbl last_int = {};
        double last_err2 =0;
        double last_time =0;
        void RunIntegration();
        void Print();
        template<class Streamable> Term & operator<<(const Streamable &data);
        void last_reset();
    };
    extern Term subtotal;
    struct TermIterator {
        size_t icurrent;
        TermIterator();
        bool IsEnd();
        TermIterator& operator++();
        Term & operator*();
    };
    extern TermList ActiveTerms;

    struct Boundaries;
    struct BoundIterator;
    typedef std::vector<Boundaries> BoundariesList;
    typedef VecDbl::iterator BoundariesListItr;
    typedef std::vector<BoundariesListItr> BoundIterList;

    struct BoundIterator {
        BoundIterList current;
        bool isFirst;
        BoundIterator();
        bool IsEnd();
        BoundIterator& operator++();
        void Print();
        inline double loBound(size_t ib){return *current[ib]     ;}
        inline double hiBound(size_t ib){return *(current[ib]+1) ;}
        //inline BoundIterator& operator=(const BoundIterator rhs){
            //this.isFirst = rhs.isFirst;
            //this.current = rhs.current;
            //return (*this);
        //}
    };

    void Init(int argc, char * argv[]);
    void SetBounds(BoundIterator);
    void WarmUp();
    void Terminate();

    namespace PrintTable {
        void IntegrationSettings();
        void Header() ;
        void Bounds(bool use_full_bound=false) ;
        void Result(const Term &term,bool printGrandTotal=false) ;
        void ResultSubTotal(bool is_grandtotal=false) ;
        void ResultGrandTotal() ;
        void Footer() ;
    }
}

#endif /* ifndef dyturbo_H */
