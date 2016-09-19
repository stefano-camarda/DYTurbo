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
        String name;
        String description;
        void (*integrate)(VecDbl &val,double &err);
        double total_int;
        double total_err2;
        double total_time;
        void RunIntegration();
        void Print();
        template<class Streamable> Term & operator<<(Streamable data);
    };
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
    };

    void Init(int argc, char * argv[]);
    void SetBounds(BoundIterator);
    void WarmUp();
    void Terminate();

    namespace PrintTable {
        void IntegrationSettings();
        void Header() ;
        void Bounds() ;
        void Result(const Term &term,bool printGrandTotal=false) ;
        void ResultSubTotal() ;
        void ResultGrandTotal() ;
        void Footer() ;
    }
}

#endif /* ifndef dyturbo_H */
