#ifndef _H_GRID_
#define _H_GRID_

// #include "Parameters.h"

struct mass_diffused
{
	// indicates whether diffusion on this direction has been calculated
	bool bUp, bLeft, bRight, bDown;
	double up, left, right, down;
};

class Grid
{
public:
	char m_name[3]; // two-letter combination of the name of the grid
				 //		'CC' - corneocyte in stratum corneum
				 //		'LP' - lipid in stratum corneum
				 //		'SC' - source, usually vehicle
				 //		'SK' - sink
	
	double m_Kw,  // partition coefficient between this object and water
		m_D, // diffusivity of solute in this object
		m_Dw, // diffusivity of solute in water
		m_K_ow, // partition coefficient of the chemical between octanol and water
		m_mw, // solute molecular weight
		m_phi_b, // volume fraction of water in corneocyte at saturation
		m_theta_b, // volume fraction of water in corneocyte
		m_mass_frac_water, // mass fraction of water in this grid
		m_r_s, // solute radius
		m_r_f; // keratin microfibril radius, in nm

	double m_concChem, // concentration of chemical
		m_concWater; // concentration of water (hydration level)
	struct mass_diffused m_mass_diffused; // concentration of chemical
	double m_x_coord, m_y_coord; // the starting x and y coordinate for this grid
	double m_dx, // grid size in verticle direction
		m_dy, // grid size in lateral direction
		m_dz; // grid size in the 3rd dimension. Note this is 2D simulation, thus dz only used to calculate the diffusion area
		
	// Parameters paras; // A class to contain all parameters	

public:
	Grid(void) { };
	virtual ~Grid(void) {};
	
	void Init(char[], double, double, double, double, double, double D_vehicle=-1);
	void Init(char[], double, double, double, double, double, double, 
		double, double, double, double, double, double, double, double, double, double, double);
	void Release() {};
	

	double getConcChem() { return m_concChem; }
	char* getName() { return m_name; }
	
	void diffuse(Grid&, Grid&, Grid&, Grid&, double);
	void setConcFromDiffMass(void);
	double compFlux(Grid&, double, double);
	
	// Rountines to calculate model parameters
	void compD_lipid();
	void compD_corneocyte();
	void compDiffusivity(double D_vehicle=-1);
	void compKcoef(void);
};

#endif
