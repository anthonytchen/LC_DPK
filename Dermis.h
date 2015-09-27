/* The header file for dermis  */
#ifndef _H_DERMIS_
#define _H_DERMIS_

#include "Grid.h"

class Dermis
{
 public:
  double m_dz, m_x_length, m_y_length; // the skin size in the z, x (verticle)
                                       //	and y (lateral) directions
  int m_nx, m_ny; // number of grids at x and y directions
  int m_boundary_cond;

  double m_mass_in, m_mass_out; // mass transferred in and out of VE
  bool m_bUseBdyUp, m_bToBlood;

  double *m_ode_Jacobian; // Jacobian matrix needed for GSL ODE solver

  Grid *m_grids, m_gridBdyUp, m_gridBdyDown, m_gridBdyLeft, m_gridBdyRight;

  // blood involved: 
  double m_dermis_totalV, // total volume of dermis
    m_bld_skin_flow, // total blood flow rate in skin
    m_bld_concChem, // current concentration in blood
    m_bld_fu, // fraction of unbounded solute in blood
    m_par_de2blood, // partition coefficient from dermis to blood
    m_mass_into_dermis, m_mass_outof_dermis; // totoal mass (or molar) flow rates

public:
  Dermis(void) {};	
  ~Dermis(void) {};
  void Init(double, double, double, int, bool b_has_blood=1);
  void InitDermisBlood(double, double, double, double bld_concChem=0, double skin_area=1.8);
  void Release();
  
  void createGrids(double, double, double, double, double, char, double);
  void updateBoundary(Grid*, Grid*, Grid*, Grid*, double mass_in=0);
  void updateBlood(double);
	
  // Functions needed for computing ODE's right hand side (i.e. dy/dt)
  void compODE_dydt (double, const double[], double []);
  static void* static_compODE_dydt_block_threads(void *);
  void compODE_dydt_block (double, const double[], double [], int, int, int, int);
  
  // I/O functions
  void displayGrids();
  void getGridsConc(double*, int);
  double getAmount();
  void saveGrids(bool, const char []);
  void getXCoord(double*, int);
  void getYCoord(double*, int);
  void saveCoord(const char [], const char []);
};

#endif
