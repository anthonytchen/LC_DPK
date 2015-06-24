#include "stdafx.h"
#include "Chemical.h"

void Chemical::Init(double mw, double K_ow, double pKa, double frac_non_ion, double frac_unbound, char acid_base)
{
  m_mw = mw;   // 119.12; // Da, i.e. g/mol
  m_K_ow = K_ow;
  m_pKa = pKa;
  m_frac_non_ion = frac_non_ion;
  m_frac_unbound = frac_unbound;

  m_acid_base = acid_base;
}
