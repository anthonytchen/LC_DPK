#include "stdafx.h"
#include "Grid.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/*  +++++++++++++++++++++++++++++++
	Initiators
*/

//
//	args: mw - molecular weight
//		  T - temperature, default is 309 K (36 deg C)
//	      eta - water viscosity (default 0.0071 P at 36 deg C)
void Grid::Init(char name[], double mw, double mass_frac_water, double mass_frac_water_sat,
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

	assert( fabs( V_water_mortar+V_water_brick+V_lipid+V_keratin - 1.0 ) < 1e-3 );
	m_theta_b = V_water_brick / V_brick_geometry * V_all_geometry;
	// m_theta_b = 0.6302;
	// printf("m_theta_b = %lf\n", m_theta_b);

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

	assert( fabs( V_water_mortar+V_water_brick+V_lipid+V_keratin - 1.0 ) < 1e-3 );
	m_phi_b = V_water_brick / V_brick_geometry * V_all_geometry;
	// m_phi_b = 0.6302;
	// printf("m_phi_b = %lf\n", m_phi_b);
	
	m_x_coord = x_coord; m_y_coord = y_coord;
	m_dx = dx;	m_dy = dy;	m_dz = dz;
	
	m_K_ow = K_ow;
	compDiffusivity();
	compKcoef();
}

void Grid::Init(char name[], double concChem, double K_ow, 
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
// --------------------------------

/*
	Main functions
*/
void Grid::diffuse(Grid& gridUp, Grid& gridLeft, Grid& gridRight, Grid& gridDown, double dt)
{
	//	The four neighbouring grids are:
	//			 up 
	//	left		     right
	//			 down
	// x -- verticle axis (up to down), y -- lateral axis (left to right)

	m_mass_diffused.up = m_dy*m_dz*dt * compFlux(gridUp, gridUp.m_dx/2, m_dx/2);
	m_mass_diffused.bUp = true;
	if ( !strcmp( gridUp.m_name, "CC") || !strcmp( gridUp.m_name, "LP") ){
		gridUp.m_mass_diffused.bDown = true;
		gridUp.m_mass_diffused.down = -m_mass_diffused.up;
	}
	
	m_mass_diffused.left = m_dx*m_dz*dt * compFlux(gridLeft, gridLeft.m_dy/2, m_dy/2);
	m_mass_diffused.bLeft = true;
	if ( !strcmp( gridLeft.m_name, "CC") || !strcmp( gridLeft.m_name, "LP") ){
		gridLeft.m_mass_diffused.bRight = true;
		gridLeft.m_mass_diffused.right = -m_mass_diffused.left;
	}
		
	m_mass_diffused.right = m_dx*m_dz*dt * compFlux(gridRight, gridRight.m_dy/2, m_dy/2);
	m_mass_diffused.bRight = true;
	if ( !strcmp( gridRight.m_name, "CC") || !strcmp( gridRight.m_name, "LP") ){
		gridRight.m_mass_diffused.bLeft = true;
		gridRight.m_mass_diffused.left = -m_mass_diffused.right;
	}
		
	m_mass_diffused.down = m_dy*m_dz*dt * compFlux(gridDown, gridDown.m_dx/2, m_dx/2);
	m_mass_diffused.bDown = true;
	if ( !strcmp( gridDown.m_name, "CC") || !strcmp( gridDown.m_name, "LP") ){
		gridDown.m_mass_diffused.bUp = true;
		gridDown.m_mass_diffused.up = -m_mass_diffused.down;
	}		
}

void Grid::setConcFromDiffMass(void)
{
	double delConc;
	delConc = m_mass_diffused.up+m_mass_diffused.left+m_mass_diffused.right+m_mass_diffused.down;	
	
	m_concChem += delConc / (m_dx*m_dy*m_dz) ;
}

// Compute the flux of solute from <other> to <this> grid
double Grid::compFlux(Grid& other, double delA, double delB)
{
	double flux, K_AB;
	K_AB = other.m_Kw / m_Kw;
	
	flux = ( other.m_concChem - K_AB*m_concChem );
	//if ( flux / (fabs(other.m_concChem)+fabs(K_AB*m_concChem)) < 1e-3 )
	//	flux = 0;
	//else
		flux /= delA/other.m_D + K_AB*delB/m_D;
	
	return flux;
}

/*
	Rountines to calculate model parameters
*/

// computing diffusivity
void Grid::compDiffusivity(double D_vehicle)
{	
	double alpha, beta, lambda, gamma, k, S, phi_f, r_s_inA, r_f_inA;
	
	if ( !strcmp(m_name, "LP") ) {	// lipid	

		if (m_mw <= 380){		
			r_s_inA = m_r_s*1e10; // unit in A
			m_D = 2 * 1E-9 * exp(-0.46*r_s_inA*r_s_inA);
		} else {
			m_D = 3 * 1E-13;
		}

	} else if ( !strcmp(m_name, "CC") ) { // corneocyte
		
		// empirically fitted parameters
		alpha = 9.47;
		beta = 9.32 * 1E-8;
		lambda = 1.09;
		gamma = -1.17;
		
		r_s_inA = m_r_s*1e10; // unit in A
		r_f_inA = m_r_f*1e10; // unit in A
		
		phi_f = 1-m_theta_b;
		// phi_f = 1-0.6302;
		k = beta*r_f_inA*r_f_inA* pow(phi_f, gamma);
		S = (r_s_inA+r_f_inA)/r_f_inA;
		S = phi_f * S*S;
	
		m_D = exp( -alpha*pow(S,lambda) ) / ( 1 + r_s_inA/sqrt(k) + r_s_inA*r_s_inA/3/k );
		m_D *= m_Dw;
		// printf("m_D = %e\n", m_D);
	
	} else if ( !strcmp(m_name, "SC") ) { // vehicle source
	
		assert( D_vehicle>0 ); // make sure some reasonable values are provided
		m_D = D_vehicle;
		
	} else if ( !strcmp(m_name, "SK") ) { // infinite sink
	
		m_D = 1; // very big value, essentially means anything in the sink can be quickly removed
		
	} else {
	
		return;
		
	}
}


// computing partition coefficient between this grid and water
void Grid::compKcoef()
{
	double K_kw;
	
	if ( !strcmp(m_name, "LP") ) {	// lipid	
		m_Kw = pow(m_K_ow, 0.7);
	} else if ( !strcmp(m_name, "CC") ) { // corneocyte
		K_kw = 5.6 * pow(m_K_ow, 0.27);
		m_Kw = (1-m_phi_b) * K_kw + m_theta_b;
		// printf("m_Kw = %e\n", m_Kw);
	} else if ( !strcmp(m_name, "SC") ) { // vehicle source
		m_Kw = 1.0;
	} else if ( !strcmp(m_name, "SK") ) { // sink
		m_Kw = 1.0;
	} else {
		return;
	}
}

