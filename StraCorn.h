/* The header file for stratum corneum  */
#ifndef _H_STRACORN_
#define _H_STRACORN_

#include "Grid.h"

class StraCorn
{
 public:
  double m_rou_lipid, m_rou_keratin, m_rou_water, // the density of lipid, keratin and water
    m_T, m_eta, // temperature, viscosity of water
    m_V_mortar, m_V_brick, m_V_all; // the volume of mortar, brick and sum (all) 
                                    //	in each element of stratum corneum dimensions
  double m_dz, m_x_length, m_y_length; // the skin size in the z, x (verticle)
                                       //	and y (lateral) directions
  int m_nx, m_ny, // number of grids at x and y directions for SC
    m_nx_grids_lipid, m_nx_grids_cc, // number of grids for each lipid/corneocyte layer in the x direction
    m_ny_grids_lipid, m_ny_grids_cc_dn, // number of grids for each lipid/corneocyte layer in the y direction
    m_n_layer_x; // number of layers in the x direction
  int m_boundary_cond;
  double m_geom_g, m_geom_d, m_geom_s, m_geom_t, m_geom_dm, m_geom_dn, m_w, // geometry parameters
    m_offset_y; // offset at the left simulation boundary, relative to corneocyte
  double m_mass_in, m_mass_out; // the mass transferred in and out of stratum corneum
  double *m_conc1D, *m_coord1D; // 1-d concentration and corresponding coordinates

  double *m_ode_Jacobian;

  Grid *m_grids, m_gridBdyUp, m_gridBdyDown, m_gridBdyLeft, m_gridBdyRight;

public:
  StraCorn(void) {};	
  ~StraCorn(void) {};
  void Init(double, double, double, double, double, int, int, double, int);
  void Release();
  
  void createGrids(double, double, double, double, double);
  void updateBoundary(Grid*, Grid*, Grid*, Grid*);
	
  // Functions needed for computing ODE's right hand side (i.e. dy/dt)
  void compODE_dydt (double, const double[], double []);
  static void* static_compODE_dydt_block_threads(void *);
  void compODE_dydt_block (double, const double[], double [], int, int, int, int);
  
  // I/O functions
  void displayGrids();
  void getGridsConc(double*, int);
  void getAmount(double*, double*, double*);
  void comp1DConc();
  void saveGrids(bool, const char []);
  void getXCoord(double*, int);
  void getYCoord(double*, int);
  void saveCoord(const char [], const char []);
};

#endif
