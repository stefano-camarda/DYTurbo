#ifndef TLVUtils_h
#define TLVUtils_h

#include <TLorentzVector.h>

using namespace std;

namespace Utils {

  // Functions to compute deltaRs

  // min phi  distance between a tlv and a collection of tlvs
  // returns the min distance
  // index will point to the position of the tlv that minimizes the phi distance
  double deltaPhimin(const TLorentzVector& tlv, const vector<TLorentzVector>& vtlv, int& index);

  // min distance between a tlv and a collection of tlvs
  // returns the min distance
  // index will point to the position of the tlv that minimizes the distance
  double deltaRmin(const TLorentzVector& tlv, const vector<TLorentzVector>& vtlv, int& index, const float min=-1);

  // min distance among a collection of tlvs
  // returns the min distance
  // index1 and index2 will point to the positions of the tlvs that minimizes the distance
  // NB: index2>index1
  double deltaRmin(const vector<TLorentzVector>& vtlv, int& index1, int& index2);

  // Return decay angles in the Collins-Soper Frame
  // Require tlv of both decay-product particles
  // Return cos(theta_{CS}) and phi_{CS} (passed by reference in argument)
  // Return Ai moment computation is pointer to vector given (passed by reference by pointer)
  void getCSFAngles(const TLorentzVector & lep1, const int &charge1, const TLorentzVector & lep2, double ebeam, double &costh, double &phi, std::vector<double> *aimom=0);

  // Return decay angles in the Vectorian-boson Helicity Frame 
  // Require tlv of both decay-product particles
  // Return cos(theta_{VH}) and phi_{VH} (passed by reference in argument)
  // Return Ai moment computation is pointer to vector given (passed by reference by pointer)
  void getVHFAngles(const TLorentzVector & lep1, const int &charge1, const TLorentzVector & lep2, double &costh, double &phi, std::vector<double> *aimom=0);

  // Return decay angles in the Beam Direction Frame
  // Require tlv of both decay-product particles
  // Return cos(theta_{BD}) and phi_{BD} (passed by reference in argument)
  void getBDFAngles(const TLorentzVector& tlv1, const int &charge1, const TLorentzVector& tlv2, double &costh, double &phi);

  // remove from values all elements that are not present in indexes.
  void remove_if_not_in(vector<int>& values, const vector<int>& indexes);

  // Generic functions to sort vector<TLorentzVector>
  vector<int> rank     ( const vector<float> & vecfloat, bool decreasing=true );
  vector<int> rankPt   ( const vector<TLorentzVector> & vectlv, bool decreasing=true );
  vector<int> rankEta  ( const vector<TLorentzVector> & vectlv, bool decreasing=true );
  vector<int> rankMass ( const vector<TLorentzVector> & vectlv, bool decreasing=true );

  struct larger {
    bool operator()(const pair<int,float>& p1, const pair<int,float>& p2) {
      return p1.second > p2.second;
    }
  };

  struct smaller {
    bool operator()(const pair<int,float>& p1, const pair<int,float>& p2) {
      return p2.second > p1.second;
    }
  };

}

#endif
