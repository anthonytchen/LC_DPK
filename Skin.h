#ifndef _H_SKIN_
#define _H_SKIN_

#include "Chemical.h"
#include "Vehicle.h"
#include "Sebum.h"
#include "StraCorn.h"
#include "ViaEpd.h"
#include "Dermis.h"
#include "Blood.h"

struct Reaction
{
  int idx_substrate, idx_product; // the index pointing to the objects in Chemical array
  double Vmax, Km; // model parameters
};

class Skin
{
public:
  double *m_concVehicleInit; // Initial concentration in the vehicle
  double m_dz_dtheta, m_x_length, m_y_length, // the skin size in the z, x (verticle)
                                       //	and y (lateral) directions
    m_x_length_ve; // the depth of viable epidermis

  int m_dim_vh, m_dim_sc, m_dim_ve, m_dim_de, m_dim_bd, m_dim_all;
  int m_nChem; // No. of chemical species considered

  double m_Vehicle_area; // The dimensions in m_gridVehicle is for the microscopic grid used for simulation.
                         // The actual vehicle application area is contained here.
  bool m_bInfSrc, 
    m_b_has_SC, m_b_has_VE, m_b_has_DE, // whether has stratum corneum, viable epidermis, dermis
    m_b_has_blood; // whether has blood compartment

  struct Reaction m_React;

  // Chemical m_Chemical;
  Vehicle *m_Vehicle;
  Sebum *m_Sebum;
  StraCorn *m_StraCorn;
  ViaEpd *m_ViaEpd;
  Dermis *m_Dermis;
  Blood *m_Blood;

public:
  Skin(void) {};	
  ~Skin(void) {};
  void Init(Chemical*, int, const bool [], double*, double*, double*, double*, double*, double, double, double, double, int, int, int, int, double, bool);
  void InitReaction(int, int, double, double); // initialisation for reaction parameters
  void Release();

  // functions to create individual compartments
  //     they return end-of-compartment x and y coordinates (x_end_coord, y_end_coord)
  //     that can be passed to subsequent compartments
  void createVH(Chemical*, double*, double*, double*, double, double, double, double, double, bool, double *x_end_coord, double *y_end_coord); // vehicle
  void createSB(double *x_end_coord, double *y_end_coord); // sebum
  void createSC(double *x_end_coord, double *y_end_coord); // stratum corneum
  void createVE(double *x_end_coord, double *y_end_coord); // viable epidermis
  void createDE(double *x_end_coord, double *y_end_coord); // dermis
  void createBD(); // blood

  // functions to link compartments
  
  void diffuseMoL(double t_start, double t_end); // method of lines using CVODE solver
  void resetVehicle(double[], double[], double[]); // reset vehicle concentration, partition coefficient, diffusivity
  void removeVehicle(); 

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
  int getNLayerXSc() { return m_StraCorn[0].m_n_layer_x; };
  int getNGridsXSc() { return m_StraCorn[0].m_nx; };
  int getNGridsYSc() { return m_StraCorn[0].m_ny; };
  void get1DConcSC(double *ret, int dim_ret, int idx_chem=0);
  void get1DCoordSC(double *ret, int dim_ret, int idx_chem=0);

  void getGridsConc(double *ret, int dim_ret, int idx_chem=0);
  void getLayersAmount(double *ret, int dim_ret, int idx_chem=0);

  void displayGrids();
  void saveGrids(bool, const char []);
  void getXCoord(double *coord_x, int dim);
  void getYCoord(double *coord_y, int dim);
  void saveCoord(const char [], const char []);
};

#endif
