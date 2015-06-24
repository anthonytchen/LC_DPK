#include "stdafx.h"
#include "Blood.h"

/*
 */
void Blood::Init(double frac_unbound, double k_clear, double body_mass, char gender)
{	
  double cardiac_output, frac_skin, frac_blood, blood_density, vol_blood_per_kg;

  /* set up some constant values */

  // Gender dependent calculation of 
  //    -- cardiac output (flow rate)
  //    -- blood volume (Morgan, Mikhail, Murray. Clinical Anesthesiology. 3rd Edition.)
  //       Adult male: 75 mL/Kg, adult female: 65, infants: 80, neonates: 85, premature neonates: 96
  switch (gender) {
  case 'M' :
    cardiac_output = 5.6* 1e-3/60; // 5.6 L/min, converted to m^3/s
    vol_blood_per_kg = 75* 1e-6; // 75 mL/Kg, converted to m^3/Kg
    break;
  case 'F' :
    cardiac_output = 4.9* 1e-3/60; // 4.9 L/min, converted to m^3/s
    vol_blood_per_kg = 65* 1e-6; // 65 mL/Kg, converted to m^3/Kg
    break;
  default:
    exit(-1);
    break;
  }
  frac_skin = 0.05; // fraction of skin blood flow as in total cardiac output
  m_flow_capil = frac_skin * cardiac_output;

  /* This is one way of calculation of blood with no regard of gender
     Now we switch to a gender-dependent calculation as above
  frac_blood = 0.07; // fraction of blood as in body mass
  blood_density = 1060; // kg / m^3
  m_vol_blood_body = frac_blood * body_mass / blood_density;
  */
  m_vol_blood_body = vol_blood_per_kg * body_mass;

  m_concChem = .0;
  m_k_clear = k_clear;
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
