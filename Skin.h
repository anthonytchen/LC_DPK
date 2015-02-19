#ifndef _H_SKIN_
#define _H_SKIN_

#include "Chemical.h"
#include "StraCorn.h"
#include "ViaEpd.h"

class Skin
{
private:
  double m_rou_lipid, m_rou_keratin, m_rou_water, // the density of lipid, keratin and water
    m_mw, // molecular weight
    m_concSource, m_DSource, // source concentration and diffusivity
    m_K_ow,	//partition coefficient between octanol and water
    m_T, m_eta; // temperature, viscosity of water
  double m_dz, m_x_length, m_y_length, // the skin size in the z, x (verticle)
                                       //	and y (lateral) directions
    m_x_length_ve; // the depth of viable epidermis
  int m_nx, m_ny, // number of grids at x and y directions for SC
    m_nx_ve, // number of grids at x direction for viable epidermis (VE)
    m_nx_grids_lipid, m_nx_grids_cc, // number of grids for each lipid/corneocyte layer in the x direction
    m_ny_grids_lipid, m_ny_grids_cc_dh; // number of grids for each lipid/corneocyte layer in the y direction
  int m_boundary_cond;

  double *m_gsl_ode_Jacobian; // Jacobian matrix needed for GSL ODE solver

  Grid *m_grids, m_gridSource, m_gridSink, m_gridSinkLeft, m_gridSinkRight;
  Grid m_gridVehicle;

  // Chemical m_Chemical;
  StraCorn m_StraCorn;
  ViaEpd m_ViaEpd;

public:
  Skin(void) {};	
  ~Skin(void) {};
  void Init(Chemical, double, double, int, int, int, double);
  void Release();
  
  void diffuseMoL_cv(double t_start, double t_end); // method of lines using CVODE solver
  double compFlux_2sc();
  double compFlux_sc2ve();
  double compFlux_ve2sk();
	
  // Functions needed for ODE solver
  int compODE_dydt (double, const double[], double []);	
	
  // Functions needed for SUNDIALS' CVODE solver
  static int static_cvODE (double, N_Vector, N_Vector, void *);
  static int static_cvJacobian (long int, long int, long int, double, 
				N_Vector, N_Vector, DlsMat J, void *,
				N_Vector, N_Vector, N_Vector);
  static int static_cvJacobian (long int, double, N_Vector, N_Vector, DlsMat J, void *,
				N_Vector, N_Vector, N_Vector);
		
  int cvODE (double, const double[], double []);
	
  // I/O functions
  void displayGrids();
  void saveGrids(bool, const char []);
  void saveCoord(const char [], const char []);
};

#endif
