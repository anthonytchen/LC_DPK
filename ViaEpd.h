/* The header file for viable epidermis  */
#ifndef _H_VIAEPD_
#define _H_VIAEPD_

#include "Comp.h"

class ViaEpd : public Comp
{
 public:

  //  double m_mass_in, m_mass_out; // mass transferred in and out of VE

public:
  ViaEpd(void) {};	
  ~ViaEpd(void) {};
  void Init(double, double, double, int, int, CoordSys, BdyCond, BdyCond, BdyCond, BdyCond);
  
  void createGrids(Chemical, double);
 
  // I/O functions
  void saveCoord(const char [], const char []);
};

#endif
