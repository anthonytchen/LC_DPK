/* The header file for viable epidermis  */
#ifndef _H_VIAEPD_
#define _H_VIAEPD_

#include "Grid.h"

class ViaEpd
{
 public:
  double m_rou_lipid, m_rou_keratin, m_rou_water, // the density of lipid, keratin and water
    m_T, m_eta, // temperature, viscosity of water
    m_V_mortar, m_V_brick, m_V_all; // the volume of mortar, brick and sum (all) 
                                    //	in each element of stratum corneum dimensions
  double m_dz, m_x_length, m_y_length; // the skin size in the z, x (verticle)
                                       //	and y (lateral) directions
  int m_nx, m_ny; // number of grids at x and y directions for SC
  int m_boundary_cond;

  double m_mass_in, m_mass_out; // mass transferred in and out of VE
  bool m_bUseBdyUp;

  double *m_ode_Jacobian; // Jacobian matrix needed for GSL ODE solver

  Grid *m_grids, m_gridBdyUp, m_gridBdyDown, m_gridBdyLeft, m_gridBdyRight;

public:
  ViaEpd(void) {};	
  ~ViaEpd(void) {};
  void Init(int, int, double, int);
  void Release();
  
  void createGrids(double, double, double, double);
  void updateBoundary(Grid*, Grid*, Grid*, Grid*, double mass_in=0);
	
  // Functions needed for computing ODE's right hand side (i.e. dy/dt)
  void compODE_dydt (double, const double[], double []);
  static void* static_compODE_dydt_block_threads(void *);
  void compODE_dydt_block (double, const double[], double [], int, int, int, int);
  
  // I/O functions
  void displayGrids();
  void saveGrids(bool, const char []);
  void saveCoord(const char [], const char []);
};

#endif
