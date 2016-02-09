#ifndef TLVUtils_cxx
#define TLVUtils_cxx

#include "TLVUtils.h"
#include <iostream>

namespace Utils 
{

  double deltaPhimin(const TLorentzVector& tlv, const vector<TLorentzVector>& vtlv, int& index) 
  {
    double dphimin=666;
    for(unsigned int i=0; i<vtlv.size(); i++) {
      double dphi = tlv.DeltaPhi(vtlv[i]);
      if(dphi<dphimin) {
        dphimin=dphi;
        index=i;
      }
    }
    return dphimin;
  }


  double deltaRmin(const TLorentzVector& tlv, const vector<TLorentzVector>& vtlv, int& index, const float min) 
  {
    double drmin=666;
    for(unsigned int i=0; i<vtlv.size(); i++) {
      double dr = tlv.DeltaR(vtlv[i]);
      if(dr<drmin && dr>min) {
        drmin=dr;
        index=i;
      }
    }
    return drmin;
  }

  double deltaRmin(const vector<TLorentzVector>& vtlv, int& index1, int& index2) 
  {
    double drmin=666;
    for(unsigned int i=0; i<vtlv.size(); i++) {
      for(unsigned int j=i+1; j<vtlv.size(); j++) {
        double dr = vtlv[i].DeltaR(vtlv[j]);
        if(dr<drmin) {
          drmin=dr;
          index1=i;
          index2=j;
        }
      }
    }
    return drmin;
  }

  void remove_if_not_in(vector<int>& values, const vector<int>& indexes) 
  {
    vector<int>::iterator idx = values.begin();
    while( idx!=values.end() ) {
      vector<int>::const_iterator pos = find(indexes.begin(), indexes.end(), *idx);
      if(pos==indexes.end()) // index not in indexes => has to be removed
        idx = values.erase(idx); // then idx points to the element after the erased one
      else
        ++idx;
    }
  }


  // Generic functions to sort vector<TLorentzVector>

  vector<int> rank     ( const vector<float> & vecfloat, bool decreasing ) 
  {

    vector<pair<int,float> > vecpair;
    for(unsigned int i=0; i<vecfloat.size(); i++)
      vecpair.push_back( pair<int,float>( i, vecfloat[i]) );

    if ( decreasing )
      sort( vecpair.begin(), vecpair.end(), larger() );
    else
      sort( vecpair.begin(), vecpair.end(), smaller() );

    vector<int> vecint;
    for(unsigned int i=0; i<vecfloat.size(); i++)
      vecint.push_back( vecpair[i].first );

    return vecint;
  }

  vector<int> rankPt   ( const vector<TLorentzVector> & vectlv, bool decreasing ) 
  {

    vector<float> vecfloat;
    for(unsigned int i=0; i<vectlv.size(); i++) {
      vecfloat.push_back( vectlv[i].Pt() );
    }

    return rank( vecfloat, decreasing );
  }

  vector<int> rankEta  ( const vector<TLorentzVector> & vectlv, bool decreasing ) 
  {

    vector<float> vecfloat;
    for(unsigned int i=0; i<vectlv.size(); i++)
      vecfloat.push_back( vectlv[i].Eta() );

    return rank( vecfloat, decreasing );
  }

  vector<int> rankMass ( const vector<TLorentzVector> & vectlv, bool decreasing ) 
  {

    vector<float> vecfloat;
    for(unsigned int i=0; i<vectlv.size(); i++)
      vecfloat.push_back( vectlv[i].M() );

    return rank( vecfloat, decreasing );
  }

  void getBDFAngles(const TLorentzVector& lep1, const int &charge1, const TLorentzVector& lep2, double &costh, double &phi)
  {
    const TLorentzVector dilep = lep1 + lep2;
    
    double cosThetap = (charge1 > 0) ? TMath::TanH( (lep2.Eta() - lep1.Eta())/2.) : TMath::TanH( (lep1.Eta() - lep2.Eta())/2.);
    costh = ( dilep.Rapidity() < 0 ) ? -cosThetap : cosThetap;
    
    double sintheta = TMath::Sqrt(1. - costh*costh);
    double phiacop = TMath::Pi() - TVector2::Phi_mpi_pi(lep1.Phi() - lep2.Phi());
    phi = TMath::Tan(phiacop/2)*sintheta;
  }
  
  void getVHFAngles(const TLorentzVector & lep1, const int &charge1, const TLorentzVector & lep2, double &costh, double &phi, std::vector<double>  *aimom)
  {
    TLorentzVector boson = lep1+lep2;
    TVector3 boostV = boson.BoostVector();
    TLorentzVector lep1_boosted = (charge1 > 0) ? lep2 : lep1;
    lep1_boosted.Boost(-boostV);
    phi = TVector2::Phi_mpi_pi(lep1_boosted.Phi()-boson.Phi());
    double theta = lep1_boosted.Angle(boson.Vect());

    if( charge1 == 0 )  std::cout<<"[93mCharge1 est nulle "<<charge1<<" AH AH AH AH[0m"<<std::endl;

    costh = TMath::Cos(theta);

    if( aimom != 0)
      {
	aimom->resize(8);
	double sintheta = TMath::Sin(theta);
	double sin2theta = TMath::Sin(2.0*theta);

	double cosphi = TMath::Cos(phi);
	double cos2phi = TMath::Cos(2.0*phi);
	double sinphi = TMath::Sin(phi);
	double sin2phi = TMath::Sin(2.0*phi);

	aimom->at(0) = 0.5 - 1.5*costh*costh;
	aimom->at(1) = sin2theta*cosphi;
	aimom->at(2) = sintheta*sintheta*cos2phi;
	aimom->at(3) = sintheta*cosphi;
	aimom->at(4) = costh;
	aimom->at(5) = sintheta*sintheta*sin2phi;
	aimom->at(6) = sin2theta*sinphi;
	aimom->at(7) = sintheta*sinphi;
      }

  }

  void getCSFAngles(const TLorentzVector & lep1, const int &charge1, const TLorentzVector & lep2, double ebeam, double &costh, double &phi, std::vector<double> *aimom)
  {
    TLorentzVector boson = lep1+lep2;
    double Lplus  = (charge1 < 0) ? lep1.E()+lep1.Pz() : lep2.E()+lep2.Pz();
    double Lminus = (charge1 < 0) ? lep1.E()-lep1.Pz() : lep2.E()-lep2.Pz();
    double Pplus  = (charge1 < 0) ? lep2.E()+lep2.Pz() : lep1.E()+lep1.Pz();
    double Pminus = (charge1 < 0) ? lep2.E()-lep2.Pz() : lep1.E()-lep1.Pz();
    
    costh  = (Lplus*Pminus - Lminus*Pplus);
    costh *= TMath::Abs(boson.Pz());
    costh /= (boson.Mag()*boson.Pz());
    costh /= TMath::Sqrt(boson.Mag2() + boson.Pt()*boson.Pt());

    TVector3 boostV = -boson.BoostVector();
    TLorentzVector lep1_boosted = (charge1 < 0 ) ? lep1 : lep2;
    lep1_boosted.Boost(boostV);

    TVector3 CSAxis, xAxis, yAxis;
    TLorentzVector p1, p2;
    double sign = +1.;
    if (boson.Z() < 0)
      sign = -1.;
    p1.SetXYZM(0., 0., sign*ebeam, 0.938); // quark (?)
    p2.SetXYZM(0., 0., -sign*ebeam, 0.938); // antiquark (?)

    p1.Boost(boostV);
    p2.Boost(boostV);
    CSAxis = (p1.Vect().Unit()-p2.Vect().Unit()).Unit();
    yAxis = (p1.Vect().Unit()).Cross(p2.Vect().Unit());
    yAxis = yAxis.Unit();
    xAxis = yAxis.Cross(CSAxis);
    xAxis = xAxis.Unit();

    phi = TMath::ATan2(lep1_boosted.Vect()*yAxis, lep1_boosted.Vect()*xAxis);
    if( aimom != 0)
      {
	aimom->resize(8);
	double theta = TMath::ACos(costh);

	double sintheta = TMath::Sin(theta);
	double sin2theta = TMath::Sin(2.0*theta);

	double cosphi = TMath::Cos(phi);
	double cos2phi = TMath::Cos(2.0*phi);
	double sinphi = TMath::Sin(phi);
	double sin2phi = TMath::Sin(2.0*phi);

	aimom->at(0) = 0.5 - 1.5*costh*costh;
	aimom->at(1) = sin2theta*cosphi;
	aimom->at(2) = sintheta*sintheta*cos2phi;
	aimom->at(3) = sintheta*cosphi;
	aimom->at(4) = costh;
	aimom->at(5) = sintheta*sintheta*sin2phi;
	aimom->at(6) = sin2theta*sinphi;
	aimom->at(7) = sintheta*sinphi;
      }
  }

}

#endif
