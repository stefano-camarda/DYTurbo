#include "pmom.h"
#include "anomalous.h"
#include "mesq.h"
#include "pegasus.h"
#include "settings.h"

#include <iostream>

complex <double> *pmom::gamma1qq;
complex <double> *pmom::gamma1qqb;
complex <double> *pmom::gamma1qqp;
complex <double> *pmom::gamma1qqbp;
complex <double> *pmom::gamma1qg;
complex <double> *pmom::gamma1gq;
complex <double> *pmom::gamma1gg;

complex <double> *pmom::gamma2qq;
complex <double> *pmom::gamma2qqb;
complex <double> *pmom::gamma2qqp;
complex <double> *pmom::gamma2qqbp;
complex <double> *pmom::gamma2qg;
complex <double> *pmom::gamma2gq;
complex <double> *pmom::gamma2gg;

complex <double> *pmom::gamma3qq;
complex <double> *pmom::gamma3qqb;
complex <double> *pmom::gamma3qqp;
complex <double> *pmom::gamma3qqbp;
complex <double> *pmom::gamma3qg;
complex <double> *pmom::gamma3gq;
complex <double> *pmom::gamma3gg;

complex <double> *pmom::gamma1qq_1;
complex <double> *pmom::gamma1qqb_1;
complex <double> *pmom::gamma1qqp_1;
complex <double> *pmom::gamma1qqbp_1;
complex <double> *pmom::gamma1qg_1;
complex <double> *pmom::gamma1gq_1;
complex <double> *pmom::gamma1gg_1;
complex <double> *pmom::gamma2qq_1;
complex <double> *pmom::gamma2qqb_1;
complex <double> *pmom::gamma2qqp_1;
complex <double> *pmom::gamma2qqbp_1;
complex <double> *pmom::gamma2qg_1;
complex <double> *pmom::gamma2gq_1;
complex <double> *pmom::gamma2gg_1;
complex <double> *pmom::gamma3qq_1;
complex <double> *pmom::gamma3qqb_1;
complex <double> *pmom::gamma3qqp_1;
complex <double> *pmom::gamma3qqbp_1;
complex <double> *pmom::gamma3qg_1;
complex <double> *pmom::gamma3gq_1;
complex <double> *pmom::gamma3gg_1;

complex <double> *pmom::gamma1qq_2;
complex <double> *pmom::gamma1qqb_2;
complex <double> *pmom::gamma1qqp_2;
complex <double> *pmom::gamma1qqbp_2;
complex <double> *pmom::gamma1qg_2;
complex <double> *pmom::gamma1gq_2;
complex <double> *pmom::gamma1gg_2;
complex <double> *pmom::gamma2qq_2;
complex <double> *pmom::gamma2qqb_2;
complex <double> *pmom::gamma2qqp_2;
complex <double> *pmom::gamma2qqbp_2;
complex <double> *pmom::gamma2qg_2;
complex <double> *pmom::gamma2gq_2;
complex <double> *pmom::gamma2gg_2;
complex <double> *pmom::gamma3qq_2;
complex <double> *pmom::gamma3qqb_2;
complex <double> *pmom::gamma3qqp_2;
complex <double> *pmom::gamma3qqbp_2;
complex <double> *pmom::gamma3qg_2;
complex <double> *pmom::gamma3gq_2;
complex <double> *pmom::gamma3gg_2;

void pmom::allocate()
{
  if (opts.mellin1d)
    {
      gamma1qq = new complex <double> [mellinint::mdim*2];
      gamma1qqb = new complex <double> [mellinint::mdim*2];
      gamma1qqp = new complex <double> [mellinint::mdim*2];
      gamma1qqbp = new complex <double> [mellinint::mdim*2];
      gamma1qg = new complex <double> [mellinint::mdim*2];
      gamma1gq = new complex <double> [mellinint::mdim*2];
      gamma1gg = new complex <double> [mellinint::mdim*2];

      gamma2qq = new complex <double> [mellinint::mdim*2];
      gamma2qqb = new complex <double> [mellinint::mdim*2];
      gamma2qqp = new complex <double> [mellinint::mdim*2];
      gamma2qqbp = new complex <double> [mellinint::mdim*2];
      gamma2qg = new complex <double> [mellinint::mdim*2];
      gamma2gq = new complex <double> [mellinint::mdim*2];
      gamma2gg = new complex <double> [mellinint::mdim*2];

      gamma3qq = new complex <double> [mellinint::mdim*2];
      gamma3qqb = new complex <double> [mellinint::mdim*2];
      gamma3qqp = new complex <double> [mellinint::mdim*2];
      gamma3qqbp = new complex <double> [mellinint::mdim*2];
      gamma3qg = new complex <double> [mellinint::mdim*2];
      gamma3gq = new complex <double> [mellinint::mdim*2];
      gamma3gg = new complex <double> [mellinint::mdim*2];
    }
  else
    {
      gamma1qq_1 = new complex <double> [mellinint::mdim*2];
      gamma1qqb_1 = new complex <double> [mellinint::mdim*2];
      gamma1qqp_1 = new complex <double> [mellinint::mdim*2];
      gamma1qqbp_1 = new complex <double> [mellinint::mdim*2];
      gamma1qg_1 = new complex <double> [mellinint::mdim*2];
      gamma1gq_1 = new complex <double> [mellinint::mdim*2];
      gamma1gg_1 = new complex <double> [mellinint::mdim*2];
      gamma2qq_1 = new complex <double> [mellinint::mdim*2];
      gamma2qqb_1 = new complex <double> [mellinint::mdim*2];
      gamma2qqp_1 = new complex <double> [mellinint::mdim*2];
      gamma2qqbp_1 = new complex <double> [mellinint::mdim*2];
      gamma2qg_1 = new complex <double> [mellinint::mdim*2];
      gamma2gq_1 = new complex <double> [mellinint::mdim*2];
      gamma2gg_1 = new complex <double> [mellinint::mdim*2];
      gamma3qq_1 = new complex <double> [mellinint::mdim*2];
      gamma3qqb_1 = new complex <double> [mellinint::mdim*2];
      gamma3qqp_1 = new complex <double> [mellinint::mdim*2];
      gamma3qqbp_1 = new complex <double> [mellinint::mdim*2];
      gamma3qg_1 = new complex <double> [mellinint::mdim*2];
      gamma3gq_1 = new complex <double> [mellinint::mdim*2];
      gamma3gg_1 = new complex <double> [mellinint::mdim*2];

      gamma1qq_2 = new complex <double> [mellinint::mdim*2];
      gamma1qqb_2 = new complex <double> [mellinint::mdim*2];
      gamma1qqp_2 = new complex <double> [mellinint::mdim*2];
      gamma1qqbp_2 = new complex <double> [mellinint::mdim*2];
      gamma1qg_2 = new complex <double> [mellinint::mdim*2];
      gamma1gq_2 = new complex <double> [mellinint::mdim*2];
      gamma1gg_2 = new complex <double> [mellinint::mdim*2];
      gamma2qq_2 = new complex <double> [mellinint::mdim*2];
      gamma2qqb_2 = new complex <double> [mellinint::mdim*2];
      gamma2qqp_2 = new complex <double> [mellinint::mdim*2];
      gamma2qqbp_2 = new complex <double> [mellinint::mdim*2];
      gamma2qg_2 = new complex <double> [mellinint::mdim*2];
      gamma2gq_2 = new complex <double> [mellinint::mdim*2];
      gamma2gg_2 = new complex <double> [mellinint::mdim*2];
      gamma3qq_2 = new complex <double> [mellinint::mdim*2];
      gamma3qqb_2 = new complex <double> [mellinint::mdim*2];
      gamma3qqp_2 = new complex <double> [mellinint::mdim*2];
      gamma3qqbp_2 = new complex <double> [mellinint::mdim*2];
      gamma3qg_2 = new complex <double> [mellinint::mdim*2];
      gamma3gq_2 = new complex <double> [mellinint::mdim*2];
      gamma3gg_2 = new complex <double> [mellinint::mdim*2];
    }
}

void pmom::calc()
{
  //Evaluate Mellin moments of the QCD splitting functions at LO, NLO, and NNLO
  complex <double> psg_qq,psg_qg,psg_gq,psg_gg,pns_p,pns_m,pns_v;
  complex <double> pqqS_p_pqqbS,pqqS_m_pqqbS,pqqS,pqqbS,pqqV,pqqbV;
  double norm;
  double nf = 5.;
  int nfl = 5-3;
  int sign = mesq::positive;
  if (opts.mellin1d)
    for (int i = 0; i < mellinint::mdim; i++)
      {
	//Flavour decompostion from https://arxiv.org/pdf/hep-ph/0408244.pdf Eqs. (2.14), (2.15)	

	//gamma1 (LO)
	psg_qq = cx(psg0_.p0sg_[0][0][nfl][i]);
	psg_qg = cx(psg0_.p0sg_[1][0][nfl][i]);
	psg_gq = cx(psg0_.p0sg_[0][1][nfl][i]);
	psg_gg = cx(psg0_.p0sg_[1][1][nfl][i]);
	pns_p = cx(pns0_.p0ns_[nfl][i]);
	pns_m = cx(pns0_.p0ns_[nfl][i]);
	pns_v = cx(pns0_.p0ns_[nfl][i]);

	norm = 4.;
	pqqS_p_pqqbS = (psg_qq-pns_p)/nf/norm;
	pqqS_m_pqqbS = (pns_v-pns_m)/nf/norm;
	pqqS         = (pqqS_p_pqqbS+pqqS_m_pqqbS)/2.;
	pqqbS        = (pqqS_p_pqqbS-pqqS_m_pqqbS)/2.;
	pqqV         = (pns_p+pns_m)/2./norm;
	pqqbV        = (pns_p-pns_m)/2./norm;

	gamma1qq[index(i,sign)]   = pqqV +pqqS;
	gamma1qqb[index(i,sign)]  = 0.;//pqqbV+pqqbS;
	gamma1qqp[index(i,sign)]  = 0.;//pqqS;
	gamma1qqbp[index(i,sign)] = 0.;//pqqbS;
	gamma1qg[index(i,sign)]   = psg_qg/2./nf/norm;
	gamma1gq[index(i,sign)]   = psg_gq/norm;
	gamma1gg[index(i,sign)]   = psg_gg/norm;

	//gamma2 (NLO)
	psg_qq = cx(psg1_.p1sg_[0][0][nfl][i]);
	psg_qg = cx(psg1_.p1sg_[1][0][nfl][i]);
	psg_gq = cx(psg1_.p1sg_[0][1][nfl][i]);
	psg_gg = cx(psg1_.p1sg_[1][1][nfl][i]);
	pns_p = cx(pns1_.p1ns_[0][nfl][i]);
	pns_m = cx(pns1_.p1ns_[1][nfl][i]);
	pns_v = cx(pns1_.p1ns_[2][nfl][i]);
	
	norm = pow(4.,2);
	pqqS_p_pqqbS = (psg_qq-pns_p)/nf/norm;
	pqqS_m_pqqbS = (pns_v-pns_m)/nf/norm;
	pqqS         = (pqqS_p_pqqbS+pqqS_m_pqqbS)/2.;
	pqqbS        = (pqqS_p_pqqbS-pqqS_m_pqqbS)/2.;
	pqqV         = (pns_p+pns_m)/2./norm;
	pqqbV        = (pns_p-pns_m)/2./norm;

	gamma2qq[index(i,sign)]   = pqqV +pqqS;
	gamma2qqb[index(i,sign)]  = pqqbV+pqqbS;
	gamma2qqp[index(i,sign)]  = pqqS;
	gamma2qqbp[index(i,sign)] = pqqbS;
	gamma2qg[index(i,sign)]   = psg_qg/2./nf/norm;
	gamma2gq[index(i,sign)]   = psg_gq/norm;
	gamma2gg[index(i,sign)]   = psg_gg/norm;
	  
	//gamma3 (NNLO)
	psg_qq = cx(psg2_.p2sg_[0][0][nfl][i]);
	psg_qg = cx(psg2_.p2sg_[1][0][nfl][i]);
	psg_gq = cx(psg2_.p2sg_[0][1][nfl][i]);
	psg_gg = cx(psg2_.p2sg_[1][1][nfl][i]);
	pns_p = cx(pns2_.p2ns_[0][nfl][i]);
	pns_m = cx(pns2_.p2ns_[1][nfl][i]);
	pns_v = cx(pns2_.p2ns_[2][nfl][i]);
	
	norm = pow(4.,3);
	pqqS_p_pqqbS = (psg_qq-pns_p)/nf/norm;
	pqqS_m_pqqbS = (pns_v-pns_m)/nf/norm;
	pqqS         = (pqqS_p_pqqbS+pqqS_m_pqqbS)/2.;
	pqqbS        = (pqqS_p_pqqbS-pqqS_m_pqqbS)/2.;
	pqqV         = (pns_p+pns_m)/2./norm;
	pqqbV        = (pns_p-pns_m)/2./norm;

	gamma3qq[index(i,sign)]   = pqqV +pqqS;
	gamma3qqb[index(i,sign)]  = pqqbV+pqqbS;
	gamma3qqp[index(i,sign)]  = pqqS;
	gamma3qqbp[index(i,sign)] = pqqbS;
	gamma3qg[index(i,sign)]   = psg_qg/2./nf/norm;
	gamma3gq[index(i,sign)]   = psg_gq/norm;
	gamma3gg[index(i,sign)]   = psg_gg/norm;
      }
  else
    for (int i = 0; i < 2*mellinint::mdim; i++)
      {
	//Flavour decompostion from https://arxiv.org/pdf/hep-ph/0408244.pdf Eqs. (2.14), (2.15)	

	//gamma1 (LO)
	psg_qq = cx(psg0_.p0sg_[0][0][nfl][i]);
	psg_qg = cx(psg0_.p0sg_[1][0][nfl][i]);
	psg_gq = cx(psg0_.p0sg_[0][1][nfl][i]);
	psg_gg = cx(psg0_.p0sg_[1][1][nfl][i]);
	pns_p = cx(pns0_.p0ns_[nfl][i]);
	pns_m = cx(pns0_.p0ns_[nfl][i]);
	pns_v = cx(pns0_.p0ns_[nfl][i]);

	norm = 4.;
	pqqS_p_pqqbS = (psg_qq-pns_p)/nf/norm;
	pqqS_m_pqqbS = (pns_v-pns_m)/nf/norm;
	pqqS         = (pqqS_p_pqqbS+pqqS_m_pqqbS)/2.;
	pqqbS        = (pqqS_p_pqqbS-pqqS_m_pqqbS)/2.;
	pqqV         = (pns_p+pns_m)/2./norm;
	pqqbV        = (pns_p-pns_m)/2./norm;

	if (i < mellinint::mdim)
	  {
	    gamma1qq_1[index(i,sign)]   = pqqV +pqqS;
	    gamma1qqb_1[index(i,sign)]  = 0.;//pqqbV+pqqbS;
	    gamma1qqp_1[index(i,sign)]  = 0.;//pqqS;
	    gamma1qqbp_1[index(i,sign)] = 0.;//pqqbS;
	    gamma1qg_1[index(i,sign)]   = psg_qg/2./nf/norm;
	    gamma1gq_1[index(i,sign)]   = psg_gq/norm;
	    gamma1gg_1[index(i,sign)]   = psg_gg/norm;
	    //cout << i << " Np1 " << mellinint::Np_1[i] << "  " << gamma1qg_1[index(i,sign)] << endl;
	  }
	else
	  {
	    gamma1qq_2[index(i-mellinint::mdim,sign)]   = pqqV +pqqS;
	    gamma1qqb_2[index(i-mellinint::mdim,sign)]  = 0.;//pqqbV+pqqbS;
	    gamma1qqp_2[index(i-mellinint::mdim,sign)]  = 0.;//pqqS;
	    gamma1qqbp_2[index(i-mellinint::mdim,sign)] = 0.;//pqqbS;
	    gamma1qg_2[index(i-mellinint::mdim,sign)]   = psg_qg/2./nf/norm;
	    gamma1gq_2[index(i-mellinint::mdim,sign)]   = psg_gq/norm;
	    gamma1gg_2[index(i-mellinint::mdim,sign)]   = psg_gg/norm;
	    //cout << i-mellinint::mdim << " Np2  " << mellinint::Np_2[i-mellinint::mdim] << "  " << gamma1qg_2[index(i-mellinint::mdim,sign)] << endl;
	  }

	//gamma2 (NLO)
	psg_qq = cx(psg1_.p1sg_[0][0][nfl][i]);
	psg_qg = cx(psg1_.p1sg_[1][0][nfl][i]);
	psg_gq = cx(psg1_.p1sg_[0][1][nfl][i]);
	psg_gg = cx(psg1_.p1sg_[1][1][nfl][i]);
	pns_p = cx(pns1_.p1ns_[0][nfl][i]);
	pns_m = cx(pns1_.p1ns_[1][nfl][i]);
	pns_v = cx(pns1_.p1ns_[2][nfl][i]);

	norm = pow(4.,2);
	pqqS_p_pqqbS = (psg_qq-pns_p)/nf/norm;
	pqqS_m_pqqbS = (pns_v-pns_m)/nf/norm;
	pqqS         = (pqqS_p_pqqbS+pqqS_m_pqqbS)/2.;
	pqqbS        = (pqqS_p_pqqbS-pqqS_m_pqqbS)/2.;
	pqqV         = (pns_p+pns_m)/2./norm;
	pqqbV        = (pns_p-pns_m)/2./norm;

	if (i < mellinint::mdim)
	  {
	    gamma2qq_1[index(i,sign)]   = pqqV +pqqS;
	    gamma2qqb_1[index(i,sign)]  = pqqbV+pqqbS;
	    gamma2qqp_1[index(i,sign)]  = pqqS;
	    gamma2qqbp_1[index(i,sign)] = pqqbS;
	    gamma2qg_1[index(i,sign)]   = psg_qg/2./nf/norm;
	    gamma2gq_1[index(i,sign)]   = psg_gq/norm;
	    gamma2gg_1[index(i,sign)]   = psg_gg/norm;
	  }
	else
	  {
	    gamma2qq_2[index(i-mellinint::mdim,sign)]   = pqqV +pqqS;
	    gamma2qqb_2[index(i-mellinint::mdim,sign)]  = pqqbV+pqqbS;
	    gamma2qqp_2[index(i-mellinint::mdim,sign)]  = pqqS;
	    gamma2qqbp_2[index(i-mellinint::mdim,sign)] = pqqbS;
	    gamma2qg_2[index(i-mellinint::mdim,sign)]   = psg_qg/2./nf/norm;
	    gamma2gq_2[index(i-mellinint::mdim,sign)]   = psg_gq/norm;
	    gamma2gg_2[index(i-mellinint::mdim,sign)]   = psg_gg/norm;
	  }
	
	//gamma3 (NNLO)
	psg_qq = cx(psg2_.p2sg_[0][0][nfl][i]);
	psg_qg = cx(psg2_.p2sg_[1][0][nfl][i]);
	psg_gq = cx(psg2_.p2sg_[0][1][nfl][i]);
	psg_gg = cx(psg2_.p2sg_[1][1][nfl][i]);
	pns_p = cx(pns2_.p2ns_[0][nfl][i]);
	pns_m = cx(pns2_.p2ns_[1][nfl][i]);
	pns_v = cx(pns2_.p2ns_[2][nfl][i]);
	
	norm = pow(4.,3);
	pqqS_p_pqqbS = (psg_qq-pns_p)/nf/norm;
	pqqS_m_pqqbS = (pns_v-pns_m)/nf/norm;
	pqqS         = (pqqS_p_pqqbS+pqqS_m_pqqbS)/2.;
	pqqbS        = (pqqS_p_pqqbS-pqqS_m_pqqbS)/2.;
	pqqV         = (pns_p+pns_m)/2./norm;
	pqqbV        = (pns_p-pns_m)/2./norm;

	if (i < mellinint::mdim)
	  {
	    gamma3qq_1[index(i,sign)]   = pqqV +pqqS;
	    gamma3qqb_1[index(i,sign)]  = pqqbV+pqqbS;
	    gamma3qqp_1[index(i,sign)]  = pqqS;
	    gamma3qqbp_1[index(i,sign)] = pqqbS;
	    gamma3qg_1[index(i,sign)]   = psg_qg/2./nf/norm;
	    gamma3gq_1[index(i,sign)]   = psg_gq/norm;
	    gamma3gg_1[index(i,sign)]   = psg_gg/norm;
	  }
	else
	  {
	    gamma3qq_2[index(i-mellinint::mdim,sign)]   = pqqV +pqqS;
	    gamma3qqb_2[index(i-mellinint::mdim,sign)]  = pqqbV+pqqbS;
	    gamma3qqp_2[index(i-mellinint::mdim,sign)]  = pqqS;
	    gamma3qqbp_2[index(i-mellinint::mdim,sign)] = pqqbS;
	    gamma3qg_2[index(i-mellinint::mdim,sign)]   = psg_qg/2./nf/norm;
	    gamma3gq_2[index(i-mellinint::mdim,sign)]   = psg_gq/norm;
	    gamma3gg_2[index(i-mellinint::mdim,sign)]   = psg_gg/norm;
	  }
      }


  //Compute negative branch
  if (opts.mellin1d)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	int idxp = index(m,mesq::positive);
	int idxm = index(m,mesq::negative);

	gamma1qq[idxm]   = conj(gamma1qq[idxp]);
	gamma1qqb[idxm]  = conj(gamma1qqb[idxp]);
	gamma1qqp[idxm]  = conj(gamma1qqp[idxp]);
	gamma1qqbp[idxm] = conj(gamma1qqbp[idxp]);
	gamma1qg[idxm]   = conj(gamma1qg[idxp]);
	gamma1gq[idxm]   = conj(gamma1gq[idxp]);
	gamma1gg[idxm]   = conj(gamma1gg[idxp]);

	gamma2qq[idxm]   = conj(gamma2qq[idxp]);
	gamma2qqb[idxm]  = conj(gamma2qqb[idxp]);
	gamma2qqp[idxm]  = conj(gamma2qqp[idxp]);
	gamma2qqbp[idxm] = conj(gamma2qqbp[idxp]);
	gamma2qg[idxm]   = conj(gamma2qg[idxp]);
	gamma2gq[idxm]   = conj(gamma2gq[idxp]);
	gamma2gg[idxm]   = conj(gamma2gg[idxp]);

	gamma3qq[idxm]   = conj(gamma3qq[idxp]);
	gamma3qqb[idxm]  = conj(gamma3qqb[idxp]);
	gamma3qqp[idxm]  = conj(gamma3qqp[idxp]);
	gamma3qqbp[idxm] = conj(gamma3qqbp[idxp]);
	gamma3qg[idxm]   = conj(gamma3qg[idxp]);
	gamma3gq[idxm]   = conj(gamma3gq[idxp]);
	gamma3gg[idxm]   = conj(gamma3gg[idxp]);
      }
  else
    {
    for (int m = 0; m < mellinint::mdim; m++)
      {
	int idxp = index(m,mesq::positive);
	int idxm = index(m,mesq::negative);

	gamma1qq_1[idxm]   = conj(gamma1qq_1[idxp]);
	gamma1qqb_1[idxm]  = conj(gamma1qqb_1[idxp]);
	gamma1qqp_1[idxm]  = conj(gamma1qqp_1[idxp]);
	gamma1qqbp_1[idxm] = conj(gamma1qqbp_1[idxp]);
	gamma1qg_1[idxm]   = conj(gamma1qg_1[idxp]);
	gamma1gq_1[idxm]   = conj(gamma1gq_1[idxp]);
	gamma1gg_1[idxm]   = conj(gamma1gg_1[idxp]);
	gamma2qq_1[idxm]   = conj(gamma2qq_1[idxp]);
	gamma2qqb_1[idxm]  = conj(gamma2qqb_1[idxp]);
	gamma2qqp_1[idxm]  = conj(gamma2qqp_1[idxp]);
	gamma2qqbp_1[idxm] = conj(gamma2qqbp_1[idxp]);
	gamma2qg_1[idxm]   = conj(gamma2qg_1[idxp]);
	gamma2gq_1[idxm]   = conj(gamma2gq_1[idxp]);
	gamma2gg_1[idxm]   = conj(gamma2gg_1[idxp]);
	gamma3qq_1[idxm]   = conj(gamma3qq_1[idxp]);
	gamma3qqb_1[idxm]  = conj(gamma3qqb_1[idxp]);
	gamma3qqp_1[idxm]  = conj(gamma3qqp_1[idxp]);
	gamma3qqbp_1[idxm] = conj(gamma3qqbp_1[idxp]);
	gamma3qg_1[idxm]   = conj(gamma3qg_1[idxp]);
	gamma3gq_1[idxm]   = conj(gamma3gq_1[idxp]);
	gamma3gg_1[idxm]   = conj(gamma3gg_1[idxp]);

	gamma1qq_2[idxm]   = conj(gamma1qq_2[idxp]);
	gamma1qqb_2[idxm]  = conj(gamma1qqb_2[idxp]);
	gamma1qqp_2[idxm]  = conj(gamma1qqp_2[idxp]);
	gamma1qqbp_2[idxm] = conj(gamma1qqbp_2[idxp]);
	gamma1qg_2[idxm]   = conj(gamma1qg_2[idxp]);
	gamma1gq_2[idxm]   = conj(gamma1gq_2[idxp]);
	gamma1gg_2[idxm]   = conj(gamma1gg_2[idxp]);
	gamma2qq_2[idxm]   = conj(gamma2qq_2[idxp]);
	gamma2qqb_2[idxm]  = conj(gamma2qqb_2[idxp]);
	gamma2qqp_2[idxm]  = conj(gamma2qqp_2[idxp]);
	gamma2qqbp_2[idxm] = conj(gamma2qqbp_2[idxp]);
	gamma2qg_2[idxm]   = conj(gamma2qg_2[idxp]);
	gamma2gq_2[idxm]   = conj(gamma2gq_2[idxp]);
	gamma2gg_2[idxm]   = conj(gamma2gg_2[idxp]);
	gamma3qq_2[idxm]   = conj(gamma3qq_2[idxp]);
	gamma3qqb_2[idxm]  = conj(gamma3qqb_2[idxp]);
	gamma3qqp_2[idxm]  = conj(gamma3qqp_2[idxp]);
	gamma3qqbp_2[idxm] = conj(gamma3qqbp_2[idxp]);
	gamma3qg_2[idxm]   = conj(gamma3qg_2[idxp]);
	gamma3gq_2[idxm]   = conj(gamma3gq_2[idxp]);
	gamma3gg_2[idxm]   = conj(gamma3gg_2[idxp]);
      }
    }
  
//  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
//    for (int i = 0; i < mellinint::mdim; i++)
//      {
//	cout << endl;
//	cout << i << endl;
//	cout << "gamma 1 " << endl;
//	cout << "qq        " << anomalous::gamma1qq[anomalous::index(i,sign)]/2. << endl;
//	cout << "qg        " << anomalous::gamma1qg[anomalous::index(i,sign)]/2. << endl;
//	cout << "gq        " << anomalous::gamma1gq[anomalous::index(i,sign)]/2. << endl;
//	cout << "gg        " << anomalous::gamma1gg[anomalous::index(i,sign)]/2. << endl;
//	cout << endl;
//	cout << "qq        " << gamma1qq[index(i,sign)] << endl;
//	cout << "qqb       " << gamma1qqb[index(i,sign)] << endl;
//	cout << "qqp       " << gamma1qqp[index(i,sign)] << endl;
//	cout << "qqbp      " << gamma1qqbp[index(i,sign)] << endl;
//	cout << "qg        " << gamma1qg[index(i,sign)] << endl;
//	cout << "gq        " << gamma1gq[index(i,sign)] << endl;
//	cout << "gg        " << gamma1gg[index(i,sign)] << endl;
//	cout << endl;
//
//	cout << "gamma 2 " << endl;
//	cout << "qq        " << (anomalous::gamma2qqV[anomalous::index(i,sign)]+anomalous::gamma2qqS[anomalous::index(i,sign)])/4. << endl;
//	cout << "qqb       " << (anomalous::gamma2qqbV[anomalous::index(i,sign)]+anomalous::gamma2qqbS[anomalous::index(i,sign)])/4. << endl;
//	cout << "qqp       " << anomalous::gamma2qqbS[anomalous::index(i,sign)]/4. << endl;
//	cout << "qqbp      " << anomalous::gamma2qqbS[anomalous::index(i,sign)]/4. << endl;
//	cout << "qg        " << anomalous::gamma2qg[anomalous::index(i,sign)]/4. << endl;
//	cout << "gq        " << anomalous::gamma2gq[anomalous::index(i,sign)]/4. << endl;
//	cout << "gg        " << anomalous::gamma2gg[anomalous::index(i,sign)]/4. << endl;
//	cout << endl;
//	cout << "qq        " << gamma2qq[index(i,sign)] << endl;
//	cout << "qqb       " << gamma2qqb[index(i,sign)] << endl;
//	cout << "qqp       " << gamma2qqp[index(i,sign)] << endl;
//	cout << "qqbp      " << gamma2qqbp[index(i,sign)] << endl;
//	cout << "qg        " << gamma2qg[index(i,sign)] << endl;
//	cout << "gq        " << gamma2gq[index(i,sign)] << endl;
//	cout << "gg        " << gamma2gg[index(i,sign)] << endl;
//	cout << endl;
//
//	cout << "gamma 3 " << endl;
//	cout << "qq        " << gamma3qq[index(i,sign)] << endl;
//	cout << "qqb       " << gamma3qqb[index(i,sign)] << endl;
//	cout << "qqp       " << gamma3qqp[index(i,sign)] << endl;
//	cout << "qqbp      " << gamma3qqbp[index(i,sign)] << endl;
//	cout << "qg        " << gamma3qg[index(i,sign)] << endl;
//	cout << "gq        " << gamma3gq[index(i,sign)] << endl;
//	cout << "gg        " << gamma3gg[index(i,sign)] << endl;
//	cout << endl;
//      }
/*
  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
    for (int i = 0; i < mellinint::mdim; i++)
      {
	cout << endl;
	cout << i << endl;
	cout << "gamma 1 " << endl;
	cout << "qq        " << gamma1qq_1[index(i,sign)] << endl;
	cout << "qqb       " << gamma1qqb_1[index(i,sign)] << endl;
	cout << "qqp       " << gamma1qqp_1[index(i,sign)] << endl;
	cout << "qqbp      " << gamma1qqbp_1[index(i,sign)] << endl;
	cout << "qg        " << gamma1qg_1[index(i,sign)] << endl;
	cout << "gq        " << gamma1gq_1[index(i,sign)] << endl;
	cout << "gg        " << gamma1gg_1[index(i,sign)] << endl;
	cout << endl;

	cout << "gamma 2 " << endl;
	cout << "qq        " << gamma2qq_1[index(i,sign)] << endl;
	cout << "qqb       " << gamma2qqb_1[index(i,sign)] << endl;
	cout << "qqp       " << gamma2qqp_1[index(i,sign)] << endl;
	cout << "qqbp      " << gamma2qqbp_1[index(i,sign)] << endl;
	cout << "qg        " << gamma2qg_1[index(i,sign)] << endl;
	cout << "gq        " << gamma2gq_1[index(i,sign)] << endl;
	cout << "gg        " << gamma2gg_1[index(i,sign)] << endl;
	cout << endl;

	cout << "gamma 3 " << endl;
	cout << "qq        " << gamma3qq_1[index(i,sign)] << endl;
	cout << "qqb       " << gamma3qqb_1[index(i,sign)] << endl;
	cout << "qqp       " << gamma3qqp_1[index(i,sign)] << endl;
	cout << "qqbp      " << gamma3qqbp_1[index(i,sign)] << endl;
	cout << "qg        " << gamma3qg_1[index(i,sign)] << endl;
	cout << "gq        " << gamma3gq_1[index(i,sign)] << endl;
	cout << "gg        " << gamma3gg_1[index(i,sign)] << endl;
	cout << endl;
      }
*/
  
}

void pmom::init()
{
  if (opts.melup <= 1)
    {
      allocate();
      calc();
    }
}

void pmom::release()
{
  if (opts.melup <= 1)
    pmom::free();
}

void pmom::free()
{
  if (opts.mellin1d)
    {
      delete[] gamma1qq;
      delete[] gamma1qqb;
      delete[] gamma1qqp;
      delete[] gamma1qqbp;
      delete[] gamma1qg;
      delete[] gamma1gq;
      delete[] gamma1gg;

      delete[] gamma2qq;
      delete[] gamma2qqb;
      delete[] gamma2qqp;
      delete[] gamma2qqbp;
      delete[] gamma2qg;
      delete[] gamma2gq;
      delete[] gamma2gg;

      delete[] gamma3qq;
      delete[] gamma3qqb;
      delete[] gamma3qqp;
      delete[] gamma3qqbp;
      delete[] gamma3qg;
      delete[] gamma3gq;
      delete[] gamma3gg;
    }
  else
    {
      delete[] gamma1qq_1;
      delete[] gamma1qqb_1;
      delete[] gamma1qqp_1;
      delete[] gamma1qqbp_1;
      delete[] gamma1qg_1;
      delete[] gamma1gq_1;
      delete[] gamma1gg_1;
      delete[] gamma2qq_1;
      delete[] gamma2qqb_1;
      delete[] gamma2qqp_1;
      delete[] gamma2qqbp_1;
      delete[] gamma2qg_1;
      delete[] gamma2gq_1;
      delete[] gamma2gg_1;
      delete[] gamma3qq_1;
      delete[] gamma3qqb_1;
      delete[] gamma3qqp_1;
      delete[] gamma3qqbp_1;
      delete[] gamma3qg_1;
      delete[] gamma3gq_1;
      delete[] gamma3gg_1;

      delete[] gamma1qq_2;
      delete[] gamma1qqb_2;
      delete[] gamma1qqp_2;
      delete[] gamma1qqbp_2;
      delete[] gamma1qg_2;
      delete[] gamma1gq_2;
      delete[] gamma1gg_2;
      delete[] gamma2qq_2;
      delete[] gamma2qqb_2;
      delete[] gamma2qqp_2;
      delete[] gamma2qqbp_2;
      delete[] gamma2qg_2;
      delete[] gamma2gq_2;
      delete[] gamma2gg_2;
      delete[] gamma3qq_2;
      delete[] gamma3qqb_2;
      delete[] gamma3qqp_2;
      delete[] gamma3qqbp_2;
      delete[] gamma3qg_2;
      delete[] gamma3gq_2;
      delete[] gamma3gg_2;
    }
}
