#include "stdafx.h"
#include "Chemical.h"

void Chemical::Init(double mw, double K_ow, double pKa)
{
  m_mw = mw;   // 119.12; // Da, i.e. g/mol
  m_K_ow = K_ow;
  m_pKa = pKa;
}
