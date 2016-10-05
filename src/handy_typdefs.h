#ifndef handy_typdefs_H
#define handy_typdefs_H
/**
 * @file handy_typdefs.h
 * @brief A brief description
 *
 * Detailed description of file
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-10-05
 */


#include <vector>
#include <string>
#include <sstream>
//! Handy typedef: string.
typedef std::string String;
//! Handy typedef: stringstream.
typedef std::ostringstream SStream;
//! Handy typedef: vector of doubles.
typedef std::vector<double> VecDbl;
//! Handy typedef: vector of vectors of doubles.
typedef vector<vector<double>> VecVecDbl;
//! Handy typedef: string vector.
typedef std::vector<String> VecStr;

#endif /* ifndef handy_typdefs_H */
