#ifndef _H_SKIN_
#define _H_SKIN_

#include "Chemical.h"
#include "Config.h"
#include "Vehicle.h"
#include "SurSebum.h"
//#include "Sebum.h"
#include "StraCorn.h"
#include "ViaEpd.h"
#include "Dermis.h"
#include "Blood.h"

struct Reaction
{
  int idx_substrate, idx_product; // the index pointing to the objects in Chemical array
  double Vmax, Km; // model parameters
};

struct CompIdx
{
  CompType type;
  Comp **pComp;
};

class Skin
{
public:
  double *m_concVehicleInit; // Initial concentration in the vehicle
  double m_dz_dtheta, m_x_length, m_y_length, // the skin size in the z, x (verticle)
                                       //	and y (lateral) directions
    m_x_length_ve; // the depth of viable epidermis

  int m_dim_vh, m_dim_sc, m_dim_ve, m_dim_de, m_dim_bd, m_dim_sb_sur, m_dim_sb_har, m_dim_all;
  int m_nChem; // No. of chemical species considered

  double m_Vehicle_area; // The dimensions in m_gridVehicle is for the microscopic grid used for simulation.
                         // The actual vehicle application area is contained here.
  bool m_bInfSrc, 
    m_b_has_blood; // whether has blood compartment

  CoordSys m_coord_sys;
  struct Reaction m_React;

  /* The amount of mass in individual compartments, row dominating arrangement */
  int m_n_amount;
  double *m_amount;
  
  /* The compartments */
  Vehicle *m_Vehicle;
  Sebum *m_Sebum;
  SurSebum *m_SurSebum;
  StraCorn *m_StraCorn;
  ViaEpd *m_ViaEpd;
  Dermis *m_Dermis;
  Blood *m_Blood;
  
  // the number of compartments in each type of compartment for each chemical species
  int m_nVehicle, m_nSurSebum, m_nSebum, m_nStraCorn, m_nViaEpd, m_nDermis, m_nBlood;
  int m_nxComp, m_nyComp;

  CompIdx **m_CompIdx; // 2D array to contain the compartment matrix

public:
  Skin(void) {};	
  ~Skin(void) {};
  void Init();
  void InitReaction(int, int, double, double); // initialisation for reaction parameters
  void Release();

  // functions to create individual compartments
  //     they return end-of-compartment x and y coordinates (x_end_coord, y_end_coord)
  //     that can be passed to subsequent compartments
  void createCompMatrix(int, int);
  void releaseCompMatrix();
  void createVH(const Chemical*, const double*, const double*, const double*, double, double, double, double, double, bool, BdyCondStr, double *x_end_coord, double *y_end_coord); // vehicle
  void createSurSB(const Chemical*, double, double, double, double, int, int, BdyCondStr, double *x_end_coord, double *y_end_coord, int, Crystal, double, double, double, double); // surface sebum
  void createSB(const Chemical*, double, double, double, double, int, int, BdyCondStr, double *x_end_coord, double *y_end_coord, int); // sebum
  void createSC(const Chemical*, double, double, int, int, double, BdyCondStr, double *x_end_coord, double *y_end_coord); // stratum corneum
  void createVE(const Chemical*, double, double, double, double, int, int, BdyCondStr, double *x_end_coord, double *y_end_coord); // viable epidermis
  void createDE(const Chemical*, double, double, double, double, int, int, bool, BdyCondStr, double *x_end_coord, double *y_end_coord); // dermis
  void createBD(const double*, const double*); // blood

  // functions to link compartments
  int getSizeBdyRight(int, int, int);
  Comp* getConcBdyRight(const double [], int, int, int, int, double*);
  int getSizeBdyDown(int, int, int);
  Comp* getConcBdyDown(const double [], int, int, int, int, double*);
  
  void diffuseMoL(double t_start, double t_end); // method of lines using CVODE solver
  void resetVehicle(double[], double[], double[]); // reset vehicle concentration, partition coefficient, diffusivity
  void removeVehicle(); 

  double compCompartAmount();
  void compFlux_2sc(double *flux);
  void compFlux_sc2down(double *flux);
  void compFlux_ve2down(double *flux);
  void compFlux_de2down(double *flux);
	
  // Functions needed for ODE solver
  int compODE_dydt (double, const double[], double []);
  void compReaction(); // add reaction with diffusion
	
  // Functions needed for SUNDIALS' CVODE solver
  static int static_cvODE (double, N_Vector, N_Vector, void *);
		
  // I/O functions

  // extracting information from stratum corneum
  double getSCYlen();
  int getNLayerXSc() { return m_StraCorn[0].m_n_layer_x; };
  int getNGridsXSc() { return m_StraCorn[0].m_nx; };
  int getNGridsYSc() { return m_StraCorn[0].m_ny; };
  void get1DConcSC(double *ret, int dim_ret, int idx_chem=0);
  void get1DCoordSC(double *ret, int dim_ret, int idx_chem=0);

  void getGridsConc(double *ret, int dim_ret, int idx_chem=0);
  void getLayersAmount(double *ret, int dim_ret, int idx_chem=0);

  void displayGrids();
  void saveGrids(bool, const char []);
  void saveAmount(bool, const char []);
  void getXCoord(double *coord_x, int dim);
  void getYCoord(double *coord_y, int dim);
  void saveCoord(const char [], const char []);
};

#endif
