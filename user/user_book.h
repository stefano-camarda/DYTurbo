#ifndef user_book_H
#define user_book_H
/**
 * @file user_book.h
 * User definitions of histogram booking.
 *
 * @brief User definitions will be included in `histo/HistoBook.h`.
 * For example of definition check `histo/HistoBook.h` file.
 *
 * @author Jakub Cuth <Jakub.Cuth@cern.ch>
 * @date 2016-09-07
 */

using namespace Kinematics;

namespace HistoHandler {
    void UserBook(){
        // dont forget to define binning
        Add( new Histo1D <BigAnswer>("biganswer") );
        // For example you can change binning and set different name
        Add( new Histo1D <BosPT>("qt") )->SetName("user_qt");

        //  histograms filled by specific weight: fill(variable1, event_weight*variable2)
        //  hack create nominator of ai profiles :
        Add( new HistoWeighted <BosPT,A1>("qt","a1") )->SetName("wgt_a1_qt");
    }
}

#endif /* ifndef user_book_H */
