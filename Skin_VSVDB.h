/*
  Class definition for Skin_VSVDB
  This class is derived from Skin and meant to contain
     vehicle (V, or VH), stratum corneum (S, or SC), viable epidermi (V, or VE), 
     dermis (D, or DE) and blood (B, or BD) for predicting
     in vivo experiments
*/

#ifndef _H_SKIN_VSVDB_
#define _H_SKIN_VSVDB_

#include "Skin.h"

class Skin_VSVDB : public Skin
{
public:

public:
  void Init(Chemical*, int, double*, double*, double*, double, double, int, int, double, double, int, double, int, double*, double*, bool);
  void Release();

};

#endif
