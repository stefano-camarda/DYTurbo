#ifndef HistoBook_H
#define HistoBook_H
/**
 * @file HistoBook.h
 * Description of cpp file
 *
 * @brief A brief description
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-08-29
 */

#include "KinematicDefinitions.h"

#include "user_book.h"

using namespace Kinematics;
namespace HistoHandler{
    void Book() {
        DeleteHists();
        Add( new Histo1D   <BosPT         > ("qt"          ) );
        Add( new Histo2D   <BosPT,BosY    > ("qt","y"      ) );
        //Add( new HistoProfile   <BosPT,A1      > ("qt","a1"     ) );
        //Add( new HistoProfile2D <BosPT,BosY,A1 > ("qt","y","a1" ) );

        UserBook();
    }
}

#endif /* ifndef HistoBook_H */
