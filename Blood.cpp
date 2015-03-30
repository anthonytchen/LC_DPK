#include "stdafx.h"
#include "Blood.h"

/*
 */
void Blood::Init(double frac_unbound, double body_mass, char gender)
{	
  double cardiac_output, frac_skin, frac_blood, blood_density;

  /* set up some constant values */

  switch (gender) {
  case 'M' :
    cardiac_output = 5.6* 1e-3/60; // 5.6 L/min, converted to m^3/s
    break;
  case 'F' :
    cardiac_output = 4.9* 1e-3/60; // 5.6 L/min, converted to m^3/s
    break;
  default:
    exit(-1);
    break;
  }
  frac_skin = 0.05; // fraction of skin blood flow as in total cardiac output
  m_flow_capil = frac_skin * cardiac_output;

  frac_blood = 0.07; // fraction of blood as in body mass
  blood_density = 1060; // kg / m^3
  m_vol_blood_body = frac_blood * body_mass / blood_density;

  m_concChem = .0;
  m_k_clear = 20e-6; // .0;
  m_f_unbound = frac_unbound;

}

void Blood::Release()
{
  return;
}


/* functions for computing the right-hand size of the odes */

// This function will be called after calculations in dermis
void Blood::updateMassInOutDermis(double massIn, double massOut, double factor)
{
  m_mass_into_dermis = massIn*factor;
  m_mass_outfrom_dermis = massOut*factor;
}

void Blood::compODE_dydt (double t, const double y[], double f[])
{
  f[0] = (-m_mass_into_dermis + m_mass_outfrom_dermis)/m_vol_blood_body;
  f[0] -= m_k_clear * y[0] / m_vol_blood_body;
}

/* ----------------- */

	
/*  +++  I/O functions +++++++++ */

void Blood::displayGrids()
{
  printf("# of grids: [x] 1, [y] 1 in blood\n");
  printf("B\n");
  fflush(stdout);
}

void Blood::saveConc(bool b_1st_time, const char fn[])
{
  FILE *file = NULL;
  int i, j, idx;

  if ( b_1st_time )
    file = fopen(fn, "w");
  else 
    file = fopen(fn, "a");
	
  fprintf(file, "%.5e\n", m_concChem);
  fclose(file);
}

/*  ------------ END <I/O functions> -------------------- */
