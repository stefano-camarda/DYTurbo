
#ifndef AIMOMENTS_CXX
#define AIMOMENTS_CXX

#include "AiMoments.h"

AiMoments::AiMoments( string AlgoName ) {

  m_algoname = AlgoName;
  pt_low = -1.;
  pt_hgh = 1e9;

}

AiMoments::~AiMoments() {}

// all kinematics Before FSR
// pt1, pt2, ebeam in GeV
// type1, type2 = PDG type of lepton

bool AiMoments::Execute(double pt1, double eta1, double phi1, double m1, int type1, 
			double pt2, double eta2, double phi2, double m2, int type2, 
			double ebeam, double weight) {

  TLorentzVector lep1, lep2, boson;
  lep1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
  lep2.SetPtEtaPhiM(pt2, eta2, phi2, m2);
  boson = lep1 + lep2;
  //printf(" ai maarten lep1 "); lep1.Print();
  //printf(" ai maarten lep2 "); lep2.Print();
  //printf(" ai maarten VB   "); boson.Print();

  //if (boson.Pt() <  pt_low) return false;
  //if (boson.Pt() >= pt_hgh) return false;

  //Collins Soper frame (CS)

  double costhcs, phics; vector<double> aimomcs;  
  Utils::getCSFAngles(lep1, ((type1>0) ? -1 : 1), lep2, ebeam, costhcs, phics, &aimomcs);
  //printf(" ai maarten costh %f phi %f \n", costhcs, phics);

  //Z Direction (ZD)

  double costhzd, phizd;; vector<double> aimomzd;
  Utils::getVHFAngles(lep1, ((type1>0) ? -1 : 1), lep2, costhzd, phizd, &aimomzd);  

  for(int ii=0; ii<8; ii++) {

    A0shift = (ii==0) ? 2./3. : 0;

    aimomTruthCS_vs_pt.at(ii)->Fill(boson.Pt(), m_Corr.at(ii)*aimomcs[ii]+A0shift, weight);
    aimomTruthZD_vs_pt.at(ii)->Fill(boson.Pt(), m_Corr.at(ii)*aimomzd[ii]+A0shift, weight);

    aimomTruthCS_vs_y.at(ii)->Fill(fabs(boson.Rapidity()), m_Corr.at(ii)*aimomcs[ii]+A0shift, weight);
    aimomTruthZD_vs_y.at(ii)->Fill(fabs(boson.Rapidity()), m_Corr.at(ii)*aimomzd[ii]+A0shift, weight);
    
    aimomTruthCS_vs_pty.at(ii)->Fill(fabs(boson.Rapidity()), boson.Pt(), m_Corr.at(ii)*aimomcs[ii]+A0shift, weight);
    aimomTruthZD_vs_pty.at(ii)->Fill(fabs(boson.Rapidity()), boson.Pt(), m_Corr.at(ii)*aimomzd[ii]+A0shift, weight);
    
  }
  
  m_pt->Fill(boson.Pt(), weight);
  m_yvspt->Fill(boson.Pt(), fabs(boson.Rapidity()), weight);
  m_yvsptS->Fill(boson.Pt(), fabs(boson.Rapidity()), weight);
  m_yvsptvsm->Fill(boson.M(), boson.Pt(), fabs(boson.Rapidity()), weight);

  if( type1+type2==0 && pt1>20. && pt2>20. && fabs(eta1)<2.4 && fabs(eta2)<2.4 && boson.M()>66. && boson.M()<116. ) {

    if(fabs(boson.Rapidity())<2.4)
      m_ptfid_inc->Fill(boson.Pt(), weight);

    if(fabs(boson.Rapidity())<1.0)
      m_ptfid_yb1->Fill(boson.Pt(), weight);
    else if(fabs(boson.Rapidity())<2.0)
      m_ptfid_yb2->Fill(boson.Pt(), weight);
    else if(fabs(boson.Rapidity())<2.4)
      m_ptfid_yb3->Fill(boson.Pt(), weight);

  }

  m_cosHvsy->Fill(fabs(boson.Rapidity()), costhzd, weight);
  m_cosHvspt->Fill(boson.Pt(), costhzd, weight);
  m_cosCSvsy->Fill(fabs(boson.Rapidity()), costhcs, weight);
  m_cosCSvspt->Fill(boson.Pt(), costhcs, weight);
  m_cosCSvsyvspt->Fill(boson.Pt(), fabs(boson.Rapidity()), costhcs, weight);

  m_phiHvsy->Fill(fabs(boson.Rapidity()), phizd, weight);
  m_phiHvspt->Fill(boson.Pt(), phizd, weight);
  m_phiCSvsy->Fill(fabs(boson.Rapidity()), phics, weight);
  m_phiCSvspt->Fill(boson.Pt(), phics, weight);
  m_phiCSvsyvspt->Fill(boson.Pt(), fabs(boson.Rapidity()), phics, weight);

  m_phiHvscosH->Fill(costhzd, phizd, weight);
  m_phiCSvscosCS->Fill(costhcs, phics, weight);
  
  return true;
}


bool AiMoments::Initialize() {

  fOut = new TFile((m_algoname+".root").c_str(),"recreate");

  // Ai moment coefficients

  float corr[] = { 20./3., 5., 10., 4., 4., 5., 5., 4. };
  for(unsigned int ii=0; ii<8; ii++) 
    m_Corr.push_back(corr[ii]);

  // Initializing truth plots
  stringstream sstr, tit; 

  for(unsigned int ii=0; ii<8; ii++) {

    sstr.str(""); tit.str(""); 
    sstr << "a" << ii << "TruthCS_vs_pt"; 
    tit << "; p_{T}^{W,Z} [GeV]; A_{" << ii << "}";
    aimomTruthCS_vs_pt.push_back( new TProfile( sstr.str().c_str(), tit.str().c_str(), 50, 0., 100.) ); 
    aimomTruthCS_vs_pt.back()->Sumw2();

    sstr.str(""); tit.str(""); 
    sstr << "a" << ii << "TruthCS_vs_y"; 
    tit << "; y_{W,Z}; A_{" << ii << "}" ;
    aimomTruthCS_vs_y.push_back( new TProfile( sstr.str().c_str(), tit.str().c_str(), 50, 0., 5.) ); 
    aimomTruthCS_vs_y.back()->Sumw2();

    sstr.str(""); tit.str(""); 
    sstr << "a" << ii << "TruthCS_vs_pty"; 
    tit << "; y_{W,Z}; p_{T}^{W,Z} [GeV]; A_{" << ii << "}" ;
    aimomTruthCS_vs_pty.push_back( new TProfile2D( sstr.str().c_str(), tit.str().c_str(), 50, 0., 5., 50, 0., 100.) ); 
    aimomTruthCS_vs_pty.back()->Sumw2();
    aimomTruthCS_vs_pty.back()->SetOption("colZ");

    sstr.str(""); tit.str(""); 
    sstr << "a" << ii << "TruthZD_vs_pt"; 
    tit << "; p_{T}^{W,Z} [GeV]; A_{" << ii << "}" ;
    aimomTruthZD_vs_pt.push_back( new TProfile( sstr.str().c_str(), tit.str().c_str(), 50, 0., 100.) ); 
    aimomTruthZD_vs_pt.back()->Sumw2();

    sstr.str(""); tit.str(""); 
    sstr << "a" << ii << "TruthZD_vs_y"; 
    tit << "; y_{W,Z}; A_{" << ii << "}" ;
    aimomTruthZD_vs_y.push_back( new TProfile( sstr.str().c_str(), tit.str().c_str(), 50, 0., 5.) ); 
    aimomTruthZD_vs_y.back()->Sumw2();
    
    sstr.str(""); tit.str(""); 
    sstr << "a" << ii << "TruthZD_vs_pty"; 
    tit << "; y_{W,Z}; p_{T}^{W,Z} [GeV]; A_{" << ii << "}" ;
    aimomTruthZD_vs_pty.push_back( new TProfile2D( sstr.str().c_str(), tit.str().c_str(), 50, 0., 5., 50, 0., 100.) ); 
    aimomTruthZD_vs_pty.back()->Sumw2();
    aimomTruthZD_vs_pty.back()->SetOption("colZ");
    
  }

  double Mass[9] = {0.,40.,60.,75.,85.,100.,120.,160.,1000.};

  double Pt[92] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0,
		   16.5,17.0,17.5,18.0,18.5,19.0,19.5,20.0,20.5,21.0,21.5,22.0,22.5,23.0,23.5,24.0,24.5,25.0,25.5,26.0,26.5,27.0,27.5,28.0,28.5,29.0,29.5,30.0,31.,32.,
		   33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,47.,49.,51.,53.,55.,58.,61.,64.,67.,70.,75.,80.,85.,90.,95.,100.};  

  double Yv[20] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.7,3.0,3.3,3.6,4.0,4.5,5.0};

  double PtFid[27] = {0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 22., 26., 30., 34., 38., 42., 46., 50., 54., 60., 70., 80., 100., 150., 200., 300., 800.};  

  m_pt           = new TH1D("pt",      "pt; p_{T}^{W,Z} [GeV]",                              200, 0, 100);                m_pt->Sumw2();
  m_yvspt        = new TH2D("yvspt",   "yvspt; p_{T}^{W,Z} [GeV]; y_{W,Z}",                  200, 0, 100, 25, 0, 5);      m_yvspt->SetOption("colZ"); m_yvspt->Sumw2();
  m_yvsptS       = new TH2D("yvsptS",  "yvsptS; p_{T}^{W,Z} [GeV]; y_{W,Z}",                 400, 0, 100, 50, 0, 5);      m_yvsptS->SetOption("colZ"); m_yvsptS->Sumw2();
  m_yvsptvsm     = new TH3D("yvsptvsm","yvsptvsm; m_{W,Z}[GeV]; p_{T}^{W,Z} [GeV]; y_{W,Z}",   8, Mass, 91, Pt, 19, Yv);  m_yvsptvsm->SetOption("colZ"); m_yvsptvsm->Sumw2();

  m_ptfid_inc    = new TH1D("ptfid_inc", "ptfid_inc; p_{T}^{Z} [GeV]", 26, PtFid);                                        m_ptfid_inc->Sumw2();
  m_ptfid_yb1    = new TH1D("ptfid_yb1", "ptfid_yb1; p_{T}^{Z} [GeV]", 26, PtFid);                                        m_ptfid_yb1->Sumw2();
  m_ptfid_yb2    = new TH1D("ptfid_yb2", "ptfid_yb2; p_{T}^{Z} [GeV]", 26, PtFid);                                        m_ptfid_yb2->Sumw2();
  m_ptfid_yb3    = new TH1D("ptfid_yb3", "ptfid_yb3; p_{T}^{Z} [GeV]", 26, PtFid);                                        m_ptfid_yb3->Sumw2();

  m_cosHvsy      = new TH2D("cosHvsy",      "cosHvsy",      25, 0,   5, 50, -1, 1);                                       m_cosHvsy->SetOption("colZ"); m_cosHvsy->Sumw2();
  m_cosHvspt     = new TH2D("cosHvspt",     "cosHvspt",     20, 0, 100, 50, -1, 1);                                       m_cosHvspt->SetOption("colZ"); m_cosHvspt->Sumw2();
  m_cosCSvsy     = new TH2D("cosCSvsy",     "cosCSvsy",     25, 0,   5, 50, -1, 1);                                       m_cosCSvsy->SetOption("colZ"); m_cosCSvsy->Sumw2();
  m_cosCSvspt    = new TH2D("cosCSvspt",    "cosCSvspt",    20, 0, 100, 50, -1, 1);                                       m_cosCSvspt->SetOption("colZ"); m_cosCSvspt->Sumw2();
  m_cosCSvsyvspt = new TH3D("cosCSvsyvspt", "cosCSvsyvspt", 10, 0, 100, 18, 0, 4.5, 50, -1, 1);                           m_cosCSvsyvspt->SetOption("colZ"); m_cosCSvsyvspt->Sumw2();

  m_phiHvsy      = new TH2D("phiHvsy",      "phiHvsy",      25, 0,   5, 64, -3.2, 3.2);                                   m_phiHvsy->SetOption("colZ"); m_phiHvsy->Sumw2();
  m_phiHvspt     = new TH2D("phiHvspt",     "phiHvspt",     20, 0, 100, 64, -3.2, 3.2);                                   m_phiHvspt->SetOption("colZ"); m_phiHvspt->Sumw2();
  m_phiCSvsy     = new TH2D("phiCSvsy",     "phiCSvsy",     25, 0,   5, 64, -3.2, 3.2);                                   m_phiCSvsy->SetOption("colZ"); m_phiCSvsy->Sumw2();
  m_phiCSvspt    = new TH2D("phiCSvspt",    "phiCSvspt",    20, 0, 100, 64, -3.2, 3.2);                                   m_phiCSvspt->SetOption("colZ"); m_phiCSvspt->Sumw2();
  m_phiCSvsyvspt = new TH3D("phiCSvsyvspt", "phiCSvsyvspt", 10, 0, 100, 18, 0, 4.5, 64, -3.2, 3.2);                       m_phiCSvsyvspt->SetOption("colZ"); m_phiCSvsyvspt->Sumw2();

  m_phiHvscosH   = new TH2D("phiHvscosH", "phiHvscosH",     50, -1, 1, 64, -3.2, 3.2);                                    m_phiHvscosH->SetOption("colZ"); m_phiHvscosH->Sumw2();
  m_phiCSvscosCS = new TH2D("phiCSvscosCS", "phiCSvscosCS", 50, -1, 1, 64, -3.2, 3.2);                                    m_phiCSvscosCS->SetOption("colZ"); m_phiCSvscosCS->Sumw2();

  return true;

}

bool AiMoments::Finalize()
{

  fOut->cd();

  m_pt->Write();
  m_yvspt->Write();
  m_yvsptS->Write();
  m_yvsptvsm->Write();

  m_ptfid_inc->Write();
  m_ptfid_yb1->Write();
  m_ptfid_yb2->Write();
  m_ptfid_yb3->Write();

  m_cosHvsy->Write();
  m_cosHvspt->Write();
  m_cosCSvsy->Write();
  m_cosCSvspt->Write();
  m_cosCSvsyvspt->Write();

  m_phiHvsy->Write();
  m_phiHvspt->Write();
  m_phiCSvsy->Write();
  m_phiCSvspt->Write();
  m_phiCSvsyvspt->Write();

  m_phiHvscosH->Write();
  m_phiCSvscosCS->Write();

  for(int ii=0; ii<8; ii++) {

    aimomTruthCS_vs_pt.at(ii)->Write(); 
    aimomTruthCS_vs_y.at(ii)->Write(); 
    aimomTruthCS_vs_pty.at(ii)->Write(); 

    aimomTruthZD_vs_pt.at(ii)->Write(); 
    aimomTruthZD_vs_y.at(ii)->Write();
    aimomTruthZD_vs_pty.at(ii)->Write();
    
  }

  fOut->Close();

  return true;

}



#endif
