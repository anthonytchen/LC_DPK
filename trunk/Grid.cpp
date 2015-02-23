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


/*  +++++++++++++++++++++++++++++++
	Initiators
	+++++++++++++++++++++++++++++++ */


//	args: mw - molecular weight
//		  T - temperature, default is 309 K (36 deg C)
//	      eta - water viscosity (default 0.0071 P at 36 deg C)
//  TODO: Read information from configuration files, instead of hard-coding here
void Grid::Init(const char name[], double mw, double mass_frac_water, double mass_frac_water_sat,
		double V_mortar_geometry, double V_brick_geometry, double V_all_geometry,
		double rou_lipid, double rou_keratin, double rou_water,
		double T, double eta, double K_ow, double x_coord, double y_coord,
		double dx, double dy, double dz)
{
  strcpy(m_name, name);
  m_concChem = 0;
  m_mass_diffused.bUp = m_mass_diffused.bLeft = m_mass_diffused.bRight = m_mass_diffused.bDown = false;
	
  double K;
	
  K = 1.3806488 * 1E-23; // Boltzmann constant, Kg m^2 s^{-2}
  m_r_f = 3.5e-9; // keratin microfibril radius, 3.5 nm

  m_mw = mw;
  m_mass_frac_water = mass_frac_water;
  m_r_s = pow( 0.9087 * mw * 3/4/M_PI, 1.0/3 )*1e-10; // from A to meter
  m_Dw = K*T/6/M_PI/eta/m_r_s;

  // calulcate various mass/volume fractions needed
  double f_l, f_k, mass_lipid, mass_keratin;
  double V_all, V_lipid, V_keratin, V_water_mortar, V_water_brick;
	
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
  m_theta_b = V_water_brick / V_brick_geometry * V_all_geometry;


  // do the same for saturated water
  //	TODO: move this part to a separate function
  // mass fraction of lipid and keratin
  mass_lipid = (1 - mass_frac_water_sat) * f_l;
  mass_keratin = (1 - mass_frac_water_sat) * f_k;

  V_all = mass_lipid/rou_lipid + mass_keratin/rou_keratin + mass_frac_water/rou_water;
  V_lipid = mass_lipid/rou_lipid / V_all;
  V_keratin = mass_keratin/rou_keratin / V_all;

  V_water_mortar = V_mortar_geometry/V_all_geometry - V_lipid;
  V_water_brick = V_brick_geometry/V_all_geometry - V_keratin;

  //	assert( fabs( V_water_mortar+V_water_brick+V_lipid+V_keratin - 1.0 ) < 1e-3 );
  m_phi_b = V_water_brick / V_brick_geometry * V_all_geometry;
	
  m_x_coord = x_coord; m_y_coord = y_coord;
  m_dx = dx;	m_dy = dy;	m_dz = dz;
	
  m_K_ow = K_ow;
  compDiffusivity();
  compKcoef();
}

void Grid::Init(const char name[], double concChem, double K_ow, 
		double dx, double dy, double dz, double D_vehicle)
{
  strcpy(m_name, name);
  m_concChem = concChem;

  m_dx = dx;
  m_dy = dy;
  m_dz = dz;
	
  m_K_ow = K_ow;
  compDiffusivity(D_vehicle);
  compKcoef();	
}

void Grid::InitVE(double mw, double Kow, double pKa, char acid_base,
		  double x_coord, double y_coord, double dx, double dy, double dz)
{
  int gsl_errno;

  strcpy(m_name, "VE");
  m_concChem = 0;
  m_mass_diffused.bUp = m_mass_diffused.bLeft = m_mass_diffused.bRight = m_mass_diffused.bDown = false;
	
  m_mw = mw;
  m_K_ow = Kow;
  m_pKa = pKa;

  m_x_coord = x_coord; m_y_coord = y_coord;
  m_dx = dx;	m_dy = dy;	m_dz = dz;

  // calculate the fraction of solute non-ionised at pH 7.4 (m_ve_fnon)
  //       and the fraction of unbound in a 2.7% albumin solution at pH 7.4 (m_ve_fu)
  // thus to calculate the binding factor in VE
  // Refs: 
  //     Florence AT, Attwood D (2006). Physicochemical Principles of Pharmacy, Pharmaceutical Press, London, p. 77.
  //     Yamazaki K, Kanaoka M (2004). Journal of Pharmaceutical Sciences, 93: 1480.
  switch (acid_base) {
  case 'A' : // weak acid
    m_ve_fnon = 1 / ( 1 + pow(10, 7.4-pKa) );
    m_ve_fu = ( 0.7936 * exp(log10(m_K_ow)) + 0.2239 ) / ( 0.7936 * exp(log10(m_K_ow)) + 1.2239 );
    break;
  case 'B' : // weak base
    m_ve_fnon = 1 / ( 1 + pow(10, pKa-7.4) );
    m_ve_fu = ( 0.5578 * exp(log10(m_K_ow)) + 0.0188 ) / ( 0.5578 * exp(log10(m_K_ow)) + 1.0188 );
    break;
  default :
    gsl_error ("Needs to provide whether it's acid or base", __FILE__, __LINE__, gsl_errno);
    exit(-1);
  }
  m_ve_binding_factor = 0.68 + 0.32/m_ve_fu + 0.025*m_ve_fnon*pow(m_K_ow, 0.7);

  compDiffusivity();
  compKcoef();
}

/*  END <Initiators>
	--------------------------------------------- */

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


void Grid::setConcFromDiffMass(void)
{
  double delConc;
  delConc = m_mass_diffused.up+m_mass_diffused.left+m_mass_diffused.right+m_mass_diffused.down;	
	
  m_concChem += delConc / (m_dx*m_dy*m_dz) ;
}

void Grid::set(Grid* other)
{
  m_Kw = other->m_Kw;
  m_D = other->m_D;
  m_Dw = other->m_Dw;
  m_K_ow = other->m_K_ow;
  m_mw = other->m_mw; 
  m_phi_b = other->m_phi_b;
  m_theta_b = other->m_theta_b;
  m_mass_frac_water = other->m_mass_frac_water;
  m_r_s = other->m_r_s;
  m_r_f = other->m_r_f;

  m_concChem = other->m_concChem;
  m_concWater = other->m_concWater;
  struct mass_diffused m_mass_diffused; // concentration of chemical
  m_x_coord = other->m_x_coord;
  m_y_coord = other->m_y_coord;
  m_dx = other->m_dx;
  m_dy = other->m_dy;
  m_dz = other->m_dz;	
}

/*  END <Main functions>
	--------------------------------------------- */

/*  ++++ Functions to calculate model parameters ++++ */

void Grid::compDiffusivity(double D_vehicle)
{	
  int gsl_errno;
  double alpha, beta, lambda, gamma, k, S, phi_f, r_s_inA, r_f_inA, 
    D_free, D_binding_factor;
	
  if ( !strcmp(m_name, "LP") ) {	// lipid	
    
    if (m_mw <= 380){		
      r_s_inA = m_r_s*1e10; // unit in Angstrom
      m_D = 2 * 1E-9 * exp(-0.46*r_s_inA*r_s_inA);
    } else {
      m_D = 3 * 1E-13;
    }
    //m_D = 9.6187e-12; // !!

  } else if ( !strcmp(m_name, "CC") ) { // corneocyte
		
    // empirically fitted parameters
    alpha = 9.47;
    beta = 9.32 * 1E-8;
    lambda = 1.09;
    gamma = -1.17;
		
    r_s_inA = m_r_s*1e10; // unit in A
    r_f_inA = m_r_f*1e10; // unit in A
		
    phi_f = 1-m_theta_b;		
    k = beta*r_f_inA*r_f_inA* pow(phi_f, gamma);
    S = (r_s_inA+r_f_inA)/r_f_inA;
    S = phi_f * S*S;
	
    m_D = exp( -alpha*pow(S,lambda) ) / ( 1 + r_s_inA/sqrt(k) + r_s_inA*r_s_inA/3/k );
    m_D *= m_Dw;
    //m_D = 1.3824e-015; // !!
			
  } else if ( !strcmp(m_name, "VE") ) { // viable epidermis

    // c.f. L. Chen's Phar. Res. paper; -8.15 should be -4.15 c.f. Kasting's original paper
    D_free = pow(10, -4.15-0.655*log10(m_mw));
    m_D = D_free / m_ve_binding_factor;

  } else if ( !strcmp(m_name, "SC") ) { // vehicle source
	
    assert( D_vehicle>0 ); // make sure some reasonable values are provided
    m_D = D_vehicle;
		
  } else if ( !strcmp(m_name, "SK") ) { // infinite sink
    
    m_D = 1; // very big value, essentially means anything in the sink can be quickly removed
    
  } else {
	
    gsl_error ("Grid name unknown", __FILE__, __LINE__, gsl_errno);
    exit(-1);		
    
  } // end if - else if - else
}


// Compute partition coefficient between this grid and water
void Grid::compKcoef()
{
  int gsl_errno;
  double K_kw;
	
  if ( !strcmp(m_name, "LP") ) {	// lipid	
    m_Kw = pow(m_K_ow, 0.7);
    //m_Kw = 0.9*pow(m_K_ow,0.69);
  } else if ( !strcmp(m_name, "CC") ) { // corneocyte
    if (m_K_ow>10)
      K_kw = 5.6 * pow(m_K_ow, 0.27);
    else 
      K_kw = 0.5* ( 1 + pow(m_K_ow, 0.7) );
    // K_kw = 1.37*4.2*pow(m_K_ow,0.31);
    m_Kw = (1-m_phi_b) * K_kw + m_theta_b;		
  } else if ( !strcmp(m_name, "VE") ) { // viable epidermis
    m_Kw = 0.7 * m_ve_binding_factor;
  } else if ( !strcmp(m_name, "SC") ) { // vehicle source
    m_Kw = 1.0;
  } else if ( !strcmp(m_name, "SK") ) { // sink
    m_Kw = 1.0;
  } else {
    gsl_error ("Grid name unknown", __FILE__, __LINE__, gsl_errno);
    exit(-1);
  }
}

/*  END <Functions to calculate model parameters>
	--------------------------------------------- */
