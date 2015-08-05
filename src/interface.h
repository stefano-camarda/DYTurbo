#ifndef interface_h
#define interface_h



extern "C" {
  double resumm_(double &costh, double &mm, double &qtt, double &yy, int& mode);
  void setup_();
  void dyinit_();
  void pdfini_();
  bool cuts_(double p[4][12], int &njet);

  void breitw_(double& x1, double& mminsq, double& mmaxsq, double& rmass, double& rwidth, double& msq, double& wt);
  void boost_(double& mass, double p1[],double p_in[], double p_out[]);

  
  void rapintegrals_(double &ymin,double &ymax, double& mass, int& nocuts);

  extern struct {
      double md_;
      double mu_;
      double ms_;
      double mc_;
      double mb_;
      double mt_;
      double mel_;
      double mmu_;
      double mtau_;
      double hmass_;
      double hwidth_;
      double wmass_;
      double wwidth_;
      double zmass_;
      double zwidth_;
      double twidth_;
      double mtausq_;
      double mcsq_;
      double mbsq_;
  } masses_;

  extern struct {
    int approxpdf_;
  } opts_;

  extern struct {
    double amz_;
  } couple_;

  extern struct {
    int nlooprun_;
  } nlooprun_;

  extern struct {
    int lhapdfs_;
  } lhapdfs_;


  // input file related variables
  extern struct {
      double sroot_;
  } energy_;

  extern struct {
      int ih1_;
      int ih2_;
  } density_;

  extern struct {
      int nproc_;
  }  nproc_;

  extern struct {
      double scale_;
      double musq_;
  } scale_;

  extern struct {
      double facscale_;
  } facscale_;

  extern struct {
      double a_param_;
      double b0p_;
  } a_param_;

  extern struct {
      double g_param_;
  } g_param_;

  extern struct {
      int order_;
  } nnlo_;

  extern struct {
      char part_[4];
  } part_;

  extern struct {
      int zerowidth_;
  } zerowidth_;

  extern struct {
      double Mwmin_;
      double Mwmax_;
  } mwminmax_;

  extern struct {
      int itmx1_;
      int ncall1_;
      int itmx2_;
      int ncall2_;
  } iterat_;

  extern struct {
      int rseed_;
  } rseed_;

  extern struct {
      int iset_;
  } pdfiset_;

  extern struct {
      double nset_;
      char prefix_[50];
  } prefix_;

  extern struct {
      char PDFname_[30];
  } lhapdf_char_;

  extern struct {
      int PDFmember_;
  } lhapdf_int_;

  extern struct {
      char runstring_[30];
  } runstring_;

  extern struct {
      int pr_;
  } pr_;

  extern struct {
    int rtsmin_;
  } rtsmin_;

  extern struct {
    int xqtcut_;
    } qtcut_;


    //  double resumm_(double &costh, double &mm, double &qtt, double &yy);
  double realvirt2_(double r[22], double &wgt);
  double countint_(double r[22], double &wgt);
  double realint_(double r[22], double &wgt);
  double virtint_(double r[22], double &wgt);
  
  double countterm_(double &costh, double &mm, double &qtt, double &yy, double &alfa, double &beta, double &cthmom0, double &cthmom1, double &cthmom2);
}

#endif
