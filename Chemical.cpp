#include "stdafx.h"
#include "Chemical.h"

void Chemical::operator=(const Chemical &other)
{
  m_mw = other.m_mw;
  m_K_ow = other.m_K_ow;
  m_pKa = other.m_pKa;
  m_frac_non_ion = other.m_frac_non_ion;
  m_frac_unbound = other.m_frac_unbound;
  m_acid_base = other.m_acid_base;
  m_r_s = other.m_r_s;
}

void Chemical::Init(double mw, double K_ow, double pKa, double frac_non_ion, double frac_unbound, char acid_base)
{
  m_mw = mw;   // 119.12; // Da, i.e. g/mol
  m_K_ow = K_ow;
  m_pKa = pKa;

  m_frac_non_ion = frac_non_ion;
  if ( frac_non_ion < 0 )
    calcIon();

  m_frac_unbound = frac_unbound;
  if ( frac_unbound < 0 )
    calcBinding();

  m_acid_base = acid_base;
  m_r_s = pow( 0.9087 * mw * 3/4/M_PI, 1.0/3 )*1e-10; // from A to meter
}

/* calculate the fraction of solute non-ionised at pH 7.4 (m_frac_non_ion) 
  Refs: Florence AT, Attwood D (2006). Physicochemical Principles of Pharmacy, Pharmaceutical Press, London, p. 77. */
void Chemical::calcIon()
{
  switch (m_acid_base) {
  case 'A' : // weak acid
    m_frac_non_ion = 1 / ( 1 + pow(10, 7.4-m_pKa) );
    break;
  case 'B' : // weak base
    m_frac_non_ion = 1 / ( 1 + pow(10, m_pKa-7.4) );
    break;
  default :
    SayBye ("Needs to provide whether it's acid or base");
  }
}

/* calculate the fraction of unbound in a 2.7% albumin solution at pH 7.4 (m_frac_unbound)
   Refs:  Yamazaki K, Kanaoka M (2004). Journal of Pharmaceutical Sciences, 93: 1480. */
void Chemical::calcBinding()
{
  switch (m_acid_base) {
    case 'A' : // weak acid
      m_frac_unbound = 1 - ( 0.7936 * exp(log10(m_K_ow)) + 0.2239 ) / ( 0.7936 * exp(log10(m_K_ow)) + 1.2239 );
      break;
    case 'B' : // weak base
      m_frac_unbound = 1 - ( 0.5578 * exp(log10(m_K_ow)) + 0.0188 ) / ( 0.5578 * exp(log10(m_K_ow)) + 1.0188 );
      break;
    default :
      SayBye ("Needs to provide whether it's acid or base");
  }
}
