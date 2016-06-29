/* The header file for sebum on skin surface */
#ifndef _H_SURSEBUM_
#define _H_SURSEBUM_

#include "Sebum.h"

enum CryShape { Sphere, Cube, HyperRect };
struct Crystal {
  CryShape shape;
  double density; // kg / m^3
  double area; // area as in the 2D simulation
  int dim; // number of items in len[]
  double len[2]; // the length on each characteristic dimension
};

class SurSebum : public Sebum
{
public:
  double m_k_disv, m_k_rect, // rate constant of dissolution, reaction
    m_k_disv_per_area,
    m_Csat; // satuation concentration in sebum
  double m_rho_solid, // solid density
    m_radius_solid, m_V_solid; // radius and volume of solid
  Crystal m_crystal;
  bool m_b_has_react;

public:
  SurSebum(void) {};	
  ~SurSebum(void) {};
  void Init(double, double, double, int, int, CoordSys, BdyCond, BdyCond, BdyCond, BdyCond, Crystal, double init_mass_solid=-1, double k_disv_per_area = -1, double k_rect=-1, double Csat=-1);
  
  //  void createGrids(Chemical, double, double);
  void compODE_dydt (double, const double [], double []);

  void updateKdisv(CryShape);
  void updateKdisv(CryShape, double);

  void getGridsConc(double*, int);
  void setGridsConc(const double [], int);

  // I/O functions
  //  void saveCoord(const char [], const char []);
  void saveGrids(bool, const char []);
};

#endif
