/*
  Class definition for Skin_VS
  This class is derived from Skin and meant to contain
     vehicle (V, or VH) and stratum corneum (S, or SC) only for predicting
     in vitro / ex vivo experiments
*/

#ifndef _H_SKIN_VS_
#define _H_SKIN_VS_

#include "Skin.h"

class Skin_VS : public Skin
{
public:

public:
  void Init(Chemical*, int, double*, double*, double*, double, double, int, int, double, bool);
  void Release();

};

#endif
