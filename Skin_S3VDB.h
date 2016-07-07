/*!
  Class definition for Skin_S3VDB
  This class is derived from Skin and meant to contain
     surface sebum (S, or SB), stratum corneum (S, or SC), hair sebum (S, or SB), 
     viable epidermi (V, or VE), dermis (D, or DE) and blood (B, or BD) for predicting
     in vivo experiments
*/

#ifndef _H_SKIN_S3VDB_
#define _H_SKIN_S3VDB_

#include "Skin.h"

class Skin_S3VDB : public Skin
{
public:

public:
  void Init(Chemical*, int, double*, double*, double*, double, int, int, int, double, int, int, int, double, int, int, double, int, double*, double*);
  void InitConfig(Chemical *, Config&);
  void Release();

};

#endif
