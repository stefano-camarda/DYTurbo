#ifndef interface_h
#define interface_h

// constants.f
#define NF 5
#define FN -5
#define NLOOP 2
#define mxpart 12

extern "C" {

    // rewrited functions
    double resumm_(double &costh, double &mm, double &qtt, double &yy, int& mode);
    void setup_();
    void dyinit_();
    void pdfini_();
    void gaussinit_();
    double dyalphas_(double &q, double &amz, int &nloop);
    bool cuts_(double p[4][12], int &njet);

    void breitw_(double& x1, double& mminsq, double& mmaxsq, double& rmass, double& rwidth, double& msq, double& wt);
    void boost_(double& mass, double p1[],double p_in[], double p_out[]);

    void rapintegrals_(double &ymin,double &ymax, double& mass, int& nocuts);
    void cacheyrapint_(double &ymin,double &ymax);

    void ctqtint_(double &m, double &y, double &qtmin, double &qtmax);

    void initmoments_();
    // fortran common spaces

    // z coupling
    extern struct {
        double l_[NF];
        double r_[NF];
        double q1_;
        double l1_;
        double r1_;
        double q2_;
        double l2_;
        double r2_;
        double le_;
        double ln_;
        double re_;
        double rn_;
        double sin2w_;
    } zcouple_;

    // masses
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
    } dymasses_;

    // ewinput
    extern struct {
        double Gf_inp_;
        double aemmz_inp_;
        double xw_inp_;
        double wmass_inp_;
        double zmass_inp_;
    } ewinput_;

    extern struct {
        int ewscheme_;
    } ewscheme_;

    // CKM
    extern struct {
        double Vud_;
        double Vus_;
        double Vub_;
        double Vcd_;
        double Vcs_;
        double Vcb_;
    } cabib_;

    //QCD coupling
    extern struct {
      double gsq_;
      double as_;
      double ason2pi_;
      double ason4pi_;
    } qcdcouple_;

    // H+b mb msbar value
    extern struct {
        double mb_msbar_;
    } mb_msbar_;

    // dimensional regularization parameters
    extern struct {
        double epinv_;
    } epinv_;

    extern struct {
        double epinv2_;
    } epinv2_;

    // controls
    extern struct {
      int approxpdf_;
      int pdfintervals_;
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
    double rtsmin_;
  } rtsmin_;

  extern struct {
    double xqtcut_;
    } qtcut_;

  extern struct {
    int doFill_;
    } dofill_;


  double realvirt2_(double r[22], double &wgt);

  double lowint_(double r[22], double &wgt);
  double realint_(double r[22], double &wgt);
  double virtint_(double r[22], double &wgt);
  double countint_(double r[22], double &wgt);
  
  double countterm_(double &costh, double &mm, double &qtt, double &yy, int &mode);

  int binner_(double p3[4], double p4[4]);
  void hists_fill_(double p3[4], double p4[4], double weight);
}

#endif
