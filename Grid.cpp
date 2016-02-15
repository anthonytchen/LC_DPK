#include "stdafx.h"
#include "Grid.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* functions relating to points for creating grids */

void setPoint(struct Point& pt, double x_coord, double y_coord, double dx, double dy, const char x_type[], const char y_type[])
{
  pt.x_coord = x_coord;
  pt.y_coord = y_coord;
  pt.dx = dx;
  pt.dy = dy;
  strcpy(pt.x_type, x_type);
  strcpy(pt.y_type, y_type);
}
void cpyPoint(struct Point& dst, struct Point& src)
{
  dst.x_coord = src.x_coord;
  dst.y_coord = src.y_coord;
  dst.dx = src.dx;
  dst.dy = src.dy;
  strcpy(dst.x_type, src.x_type);
  strcpy(dst.y_type, src.y_type);
}
/* ------------ */

void Grid::operator=(const Grid &other)
{
  m_Kw = other.m_Kw;
  m_D = other.m_D;
  m_concChem = other.m_concChem;
  m_x_coord = other.m_x_coord;
  m_y_coord = other.m_y_coord;
  m_dx = other.m_dx;
  m_dy = other.m_dy;
  m_dz = other.m_dz;	
}


/*  +++++++++++++++++++++++++++++++
	Initiators
	+++++++++++++++++++++++++++++++ */

/* generic init */
void Grid::Init(const char name[], Chemical chem, double concChem, double x_coord, double y_coord,
		double dx, double dy, double dz)
{
  strcpy(m_name, name);
  m_chemical = chem;
  m_concChem = concChem;

  m_x_coord = x_coord; m_y_coord = y_coord;
  m_dx = dx;  m_dy = dy; m_dz = dz;
}

/* init vehicle */
void Grid::InitVH(const char name[], Chemical chem, double concChem, double x_coord, double y_coord,
		  double dx, double dy, double dz, double T, double eta, double K_vw, double diff_vh)
{
  Init(name, chem, concChem, x_coord, y_coord, dx, dy, dz);

  
  double K;

  if (diff_vh<0) { /* calculation of diffusivity in water */
    K = 1.3806488 * 1E-23; // Boltzmann constant, Kg m^2 s^{-2}
    m_D = K*T/6/M_PI/eta/chem.m_r_s; // diffusivity in water, Stoke-Eistein equation
  }
  else
    m_D = diff_vh;

  m_Kw = K_vw;
}

/* init sebum */
void Grid::InitSB(const char name[], Chemical chem, double concChem, double x_coord, double y_coord,
		  double dx, double dy, double dz, double K_sw, double D)
{
  Init(name, chem, concChem, x_coord, y_coord, dx, dy, dz);

  
  double K;

  if (D<0) { /* calculate diffusivity from QSPR model */
    SayBye("QSPR for calculating sebum diffusion coefficient not yet implemented");
  }
  else
    m_D = D;

  if (K_sw<0) { /* calculate partition coefficient from QSPR model */
    SayBye("QSPR for calculating sebum partition coefficient not yet implemented");
  }
  else
    m_Kw = K_sw;
}

/* init sink */
void Grid::InitSK(const char name[], Chemical chem, double concChem, double x_coord, double y_coord,
		  double dx, double dy, double dz)
{
  Init(name, chem, concChem, x_coord, y_coord, dx, dy, dz);

  m_D = 1; // very big value, essentially means anything in the sink can be quickly removed
  m_Kw = 1; // arbitrary value
}

void Grid::InitSK()
{
  strcpy(m_name, "SK");
  m_concChem = 0;

  m_x_coord = m_y_coord = 0;
  m_dx = m_dy = m_dz = 0;

  m_D = 1; // very big value, essentially means anything in the sink can be quickly removed
  m_Kw = 1; // arbitrary value, doesn't matter in flux calculation
}

/*  Init  stratum corneum
    T - temperature, default is 309 K (36 deg C)
    eta - water viscosity (default 0.0071 P at 36 deg C) */
void Grid::InitSC(const char name[], Chemical chem, double concChem, double mass_frac_water, double mass_frac_water_sat,
		double V_mortar_geometry, double V_brick_geometry, double V_all_geometry,
		double rou_lipid, double rou_keratin, double rou_water, double T, double eta, 
		double x_coord, double y_coord, double dx, double dy, double dz)
{
  /* 1. set up member variables */

  Init(name, chem, concChem, x_coord, y_coord, dx, dy, dz);

  /* 2. calulcate various mass/volume fractions needed */

  double f_l, f_k, mass_lipid, mass_keratin, theta_b, phi_b,
    V_all, V_lipid, V_keratin, V_water_mortar, V_water_brick;
	
  f_l = 0.125; // dry mass fraction of SC lipid and keratin
  f_k = 1 - f_l;

  // mass fraction of lipid and keratin
  mass_lipid = (1 - mass_frac_water) * f_l;
  mass_keratin = (1 - mass_frac_water) * f_k;

  V_all = mass_lipid/rou_lipid + mass_keratin/rou_keratin + mass_frac_water/rou_water;
  V_lipid = mass_lipid/rou_lipid / V_all;
  V_keratin = mass_keratin/rou_keratin / V_all;

  V_water_mortar = V_mortar_geometry/V_all_geometry - V_lipid;
  V_water_brick = V_brick_geometry/V_all_geometry - V_keratin;

  //	assert( fabs( V_water_mortar+V_water_brick+V_lipid+V_keratin - 1.0 ) < 1e-3 );
  theta_b = V_water_brick / V_brick_geometry * V_all_geometry;


  // do the same for saturated water
  //	TODO: move this part to a separate function
  // mass fraction of lipid and keratin
  mass_lipid = (1 - mass_frac_water_sat) * f_l;
  mass_keratin = (1 - mass_frac_water_sat) * f_k;

  V_all = mass_lipid/rou_lipid + mass_keratin/rou_keratin + mass_frac_water_sat/rou_water;
  V_lipid = mass_lipid/rou_lipid / V_all;
  V_keratin = mass_keratin/rou_keratin / V_all;

  V_water_mortar = V_mortar_geometry/V_all_geometry - V_lipid;
  V_water_brick = V_brick_geometry/V_all_geometry - V_keratin;

  //	assert( fabs( V_water_mortar+V_water_brick+V_lipid+V_keratin - 1.0 ) < 1e-3 );
  phi_b = V_water_brick / V_brick_geometry * V_all_geometry;	


  /* 3. calculate diffusivity and partition coefficient */

  double K, r_f, Dw, r_s_inA, alpha, beta, lambda, gamma, K_kw, r_f_inA, phi_f, k, S;

  K = 1.3806488 * 1E-23; // Boltzmann constant, Kg m^2 s^{-2}
  r_f = 3.5e-9; // keratin microfibril radius, 3.5 nm
  Dw = K*T/6/M_PI/eta/chem.m_r_s; // diffusivity in water, Stoke-Eistein equation

  if ( !strcmp(m_name, "LP") ) {	// lipid	

    //m_Kw = pow(m_K_ow, 0.7);
    m_Kw = rou_lipid / rou_water * pow(chem.m_K_ow,0.69);

    if (chem.m_mw <= 380){		
      r_s_inA = chem.m_r_s*1e10; // unit in Angstrom
      m_D = 2 * 1E-9 * exp(-0.46*r_s_inA*r_s_inA);
    } else {
      m_D = 3 * 1E-13;
    }

  } else if ( !strcmp(m_name, "CC") ) { // corneocyte

    /*
    if (m_K_ow>10)
      K_kw = 5.6 * pow(m_K_ow, 0.27);
    else 
      K_kw = 0.5* ( 1 + pow(m_K_ow, 0.7) );
    */
    K_kw = rou_keratin / rou_water * 4.2 * pow(chem.m_K_ow,0.31);
    m_Kw = (1-phi_b) * K_kw + theta_b;

    // empirically fitted parameters
    alpha = 9.47;
    beta = 9.32 * 1E-8;
    lambda = 1.09;
    gamma = -1.17;
		
    r_s_inA = chem.m_r_s*1e10; // unit in A
    r_f_inA = r_f*1e10; // unit in A
		
    phi_f = 1-theta_b;		
    k = beta*r_f_inA*r_f_inA* pow(phi_f, gamma);
    S = (r_s_inA+r_f_inA)/r_f_inA;
    S = phi_f * S*S;
	
    m_D = exp( -alpha*pow(S,lambda) ) / ( 1 + r_s_inA/sqrt(k) + r_s_inA*r_s_inA/3/k );
    m_D *= Dw;
  }
}


void Grid::InitVE_DE(const char name[], Chemical chem, double concChem, double x_coord, double y_coord, double dx, double dy, double dz)
{
  Init(name, chem, concChem, x_coord, y_coord, dx, dy, dz);

  double binding_factor, D_free;

  binding_factor = 0.68 + 0.32/chem.m_frac_unbound + 0.025*chem.m_frac_non_ion*pow(chem.m_K_ow, 0.7);

  // c.f. L. Chen's Phar. Res. paper; -8.15 used because of unit (m2/s)
  //  Kasting's original paper used -4.15 because of unit (cm2/s)
  D_free = pow(10, -8.15-0.655*log10(chem.m_mw));
  m_D = D_free / binding_factor;

  m_Kw = 0.7 * binding_factor;
}

/*  --------------------------------------------- */


/*  +++++++++++++++++++++++++++++++++++++++++++
	Main functions
	+++++++++++++++++++++++++++++++++++++++++++ */


// Compute the flux of solute from <other> to <this> grid
//  However, do not use concentration values in the Grid objects,
//  instead use <conc_this> and <conc_other>
double Grid::compFlux(Grid* other, double conc_this, double conc_other, 
		      double dist_this, double dist_other, double *deriv_this, double *deriv_other)
{
  double flux, K_other2this, tmp1, tmp2;
	
  K_other2this = other->m_Kw / this->m_Kw;
  flux = ( conc_other - K_other2this*conc_this );
  tmp1 = dist_other/other->m_D + K_other2this*dist_this/this->m_D;
  flux /= tmp1;
	
  if (deriv_this!=NULL)
    *deriv_this = -K_other2this / tmp1;
	
  if (deriv_other!=NULL)
    *deriv_other = 1 / tmp1;
		
  return flux;
}


/* -------------------------------------------- */

