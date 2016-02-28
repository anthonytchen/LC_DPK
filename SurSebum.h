/* The header file for sebum on skin surface */
#ifndef _H_SURSEBUM_
#define _H_SURSEBUM_

#include "Sebum.h"

enum CryShape { Sphere, Cube };

class SurSebum : public Sebum
{
public:
  double m_k_disv, m_k_rect, // rate constant of dissolution, reaction
    m_k_disv_per_area,
    m_Csat; // satuation concentration in sebum
  double m_mass_reacted, // accumulated mass due to reaction-induced breakdown
    m_mass_solid, // remaining mass in the solid
    m_rho_solid, // solid density
    m_radius_solid, m_V_solid; // radius and volume of solid
  bool m_b_has_solid, m_b_has_react;

public:
  SurSebum(void) {};	
  ~SurSebum(void) {};
  void Init(double, double, double, int, int, CoordSys, BdyCond, BdyCond, BdyCond, BdyCond, double init_mass_solid=-1, double k_disv_per_area = -1, double k_rect=-1, double Csat=-1);
  
  //  void createGrids(Chemical, double, double);
  void compODE_dydt (double, const double [], double []);

  void updateKdisv(CryShape shape);
  // I/O functions
  //  void saveCoord(const char [], const char []);
};

#endif
