#ifndef _H_SKIN_
#define _H_SKIN_

#include "Grid.h"

class Skin
{
private:
  double m_rou_lipid, m_rou_keratin, m_rou_water, // the density of lipid, keratin and water
    m_mw, // molecular weight
    m_concSource, m_DSource, // source concentration and diffusivity
    m_K_ow,	//partition coefficient between octanol and water
    m_T, m_eta, // temperature, viscosity of water
    m_V_mortar, m_V_brick, m_V_all; // the volume of mortar, brick and sum (all) 
                                    //	in each element of stratum corneum dimensions
  double m_dz, m_x_length, m_y_length, // the skin size in the z, x (verticle)
                                       //	and y (lateral) directions
    m_x_length_ve; // the depth of viable epidermis
  int m_nx, m_ny, // number of grids at x and y directions for SC
    m_nx_ve, // number of grids at x direction for viable epidermis (VE)
    m_nx_grids_lipid, m_nx_grids_cc, // number of grids for each lipid/corneocyte layer in the x direction
    m_ny_grids_lipid, m_ny_grids_cc_dh; // number of grids for each lipid/corneocyte layer in the y direction
  int m_boundary_cond;
  double m_geom_g, m_geom_d, m_geom_s, m_geom_t, m_geom_dm, m_geom_dh, m_w, // geometry parameters
    m_offset_y; // offset at the left simulation boundary, relative to corneocyte
  double m_dt; // simulation time interval
  double *m_gsl_ode_Jacobian; // Jacobian matrix needed for GSL ODE solver

  Grid *m_grids, m_gridSource, m_gridSink, m_gridSinkLeft, m_gridSinkRight;

public:
  Skin(void) {};	
  ~Skin(void) {};
  void Init(double, double, double, double, double, double, double, double, int, int, double, double);
  void Release();
  
  void createGrids();
  void createGridsVE(double, double, double, double, double, int, int);
  void diffuseMoL_cv(double t_start, double t_end); // method of lines using CVODE solver
  /*
  double compFlux(Grid*, Grid*, double, double, double, double, 
		  double *deriv_this=NULL, double *deriv_other=NULL);
  */
  double compFlux_2sc();
  double compFlux_sc2ve();
  double compFlux_ve2sk();
	
  // Functions needed for ODE solver
  // 	The prefix "gsl" is for historical reasons
  //	i.e. the gsl ODE solvers were initially used, but changed to CVODE solvers
  static void* static_gslODE_threads(void *paras);
  int gslODE (double, const double[], double []);	
  void gslODE (double, const double[], double [], int, int, int, int);	
	
  // Functions needed for SUNDIALS' CVODE solver
  static int static_cvODE (double, N_Vector, N_Vector, void *);
  static int static_cvJacobian (long int, long int, long int, double, 
				N_Vector, N_Vector, DlsMat J, void *,
				N_Vector, N_Vector, N_Vector);
  static int static_cvJacobian (long int, double, N_Vector, N_Vector, DlsMat J, void *,
				N_Vector, N_Vector, N_Vector);
		
  int cvODE (double, const double[], double []);
	
  // Setting/copying functions
  void setPoint(struct Point&, double, double, double, double, const char [], const char []);
  void cpyPoint(struct Point&, struct Point&);
  
  // I/O functions
  void displayGrids();
  void saveGrids(bool, const char []);
  void saveCoord(const char [], const char []);
};

#endif
