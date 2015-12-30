/* The header file for Comp -- Compartment.
   This will be the parent class that combines grids into a compartment.
   From this class the compartments of stratum corneum, viable epidermis,
   dermis, blood, hair follicle, sebum layer on top of skin, vehicle etc.
   will be derived.
*/
#ifndef _H_COMP_
#define _H_COMP_

#include "Grid.h"

/* The coordinate system (only 2D considered)
     -- Cartesian: 
     -- Cylindrical: axisymmetric                  */
enum CoordSys {Cartesian, Cylindrical}; 

/* The boundary conditions
    -- ZeroFlux: no flux across the boundary
    -- ZeroConc: boundary is a sink (N.A. for the up boundary)
    -- Periodic: left-right connected (N.A. for up/down boundary)
    -- FromOther: the boundary is in contact with other compartment
    -- ODE: to indicated boundary conditions are not needed if the compartment is ODE
 */
enum BdyCond {ZeroFlux, ZeroConc, Periodic, FromOther, ODE };

class Comp
{
 public:

  /* member variables to do with geometry */
  CoordSys m_coord_sys;
  double m_x_length, m_y_length, // compartment size in the x (verticle) and y (lateral) directions
                                 // if cylindrical coordinate (axisymmetric), y starts from centre to right
    m_dz_dtheta; // compartment size in the z (width) direction, either dz (Cartesian, metre) or dtheta (Cylindrical, degree)
  int m_nx, m_ny; // number of grids at x and y directions of the compartment; within each compartment the grids are regular

  /* member variables to do with boundaries */
  int m_n_gridsBdyRight, m_n_gridsBdyDown;
  BdyCond m_BdyCond_up, m_BdyCond_left, m_BdyCond_right, m_BdyCond_down;
  Grid *m_gridsBdyUp, *m_gridsBdyLeft, *m_gridsBdyRight, *m_gridsBdyDown; // grids of the 4 boundaries
  double *m_MassIn_up, *m_MassIn_left, // to avoid duplicate calculations, only record mass-in from up and left
    *m_MassOut_right, *m_MassOut_down; //  and mass-out to right and down

  /* the main grids */
  Grid *m_grids;

  /* member variables for I/O and MISC */
  double *m_conc1D, *m_coord1D; // 1-d concentration and corresponding coordinates, verticle direction

public:
  Comp(void) {};	
  ~Comp(void) {};
  void Init(CoordSys, double, BdyCond, BdyCond, BdyCond, BdyCond);
  void Release();
  
  void createGrids() {}; // not used; meshing should be done in specific compartments
  void createBoundary(int n_gridsBdyRight=0, int n_gridsBdyDown=0);
  void setBoundaryGrids(Grid *gridsBdyRight=NULL, Grid *gridsBdyDown=NULL);
  void setBoundaryConc(double *concBdyRight=NULL, double *concBdyDown=NULL);
  void setBdyMassInOutZero();
  void passBdyMassOut(Comp *bdyRight, Comp *bdydown);

  double compInterArea(Grid gridThiis, int direction);
  double compVolume(Grid gridThiis);
  double compTotalArea(int);
  double compTotalVolume();

  double compMassIrregGridsRight(Grid gridThiis, double conc_this);
  double compMassIrregGridsDown(Grid gridThiis, double conc_this);
	
  // Functions needed for computing ODE's right hand side (i.e. dy/dt)
  void compODE_dydt (double, const double[], double []);
  static void* static_compODE_dydt_block_threads(void *);
  void compODE_dydt_block (double, const double[], double [], int, int, int, int);
  
  // I/O functions
  void displayGrids();
  double getAmount();
  void getGridsConc(double*, int);
  void saveGrids(bool, const char []);
  void getXCoord(double*, int);
  void getYCoord(double*, int);
  void saveCoord(const char [], const char [], const char []);
};

/* Structure and definition for parallel computing */

struct pthread_struct {
	Comp *comp_obj;
	double t;
	const double *y;
	double *f;
	int x_start, x_end, y_start, y_end;
};
/* ------------------- */

#endif
