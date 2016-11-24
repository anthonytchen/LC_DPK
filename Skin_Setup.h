/*!
  Class definition for Skin_Setup
  This class is derived from Skin and meant to set up the compartments
  in simulation as instructed by user

  Valid options:

  1. Vector compartments
     This means that all compartments are arranged one below the other (thus the term "vector")
     Bear in mind that within each compartment, we still consider 2D diffusion

     Possibilities include (where receptor is automatically added for in vitro set ups):
	1.1. VS: vehicle, stratum corneum, (then receptor)
	1.2. VSV: vehicle, stratum corneum, viable epidermis, (then receptor)
	1.3. VSVD: vehicle, stratum cornrum, viable epidermis, dermis, (then receptor)
	1.4. VSVDB: vehicle, stratum cornrum, viable epidermis, dermis, blood

  2. Matrix compartments
     This means that compartments are arranged into a 2D matrix

	2.1. VSH: vehicle, sc, hair
	2.2. VSVH: vehicle, stratum corneum, viable epidermis, hair
	2.3. VSVDH: vehicle, stratum cornrum, viable epidermis, dermis, hair
	2.4. VSVDBH: vehicle, stratum cornrum, viable epidermis, dermis, blood, hair

	(The following includes a sebum compartment on skin surface, thus no vehicle -- an ad hoc choice)
	2.5. SSH: sebum, stratum corneum, hair
	2.6. SSVH: sebum, stratum corneum, viable epidermis, hair
	2.7. SSVDH: sebum, stratum corneum, viable epidermis, dermis, hair
	2.8. SSVDBH: sebum, stratum corneum, viable epidermis, dermis, blood, hair
*/

#ifndef _H_SKIN_SETUP_
#define _H_SKIN_SETUP__

#include "Skin.h"

class Skin_Setup : public Skin
{
public:
  char m_sComps[1024];

public:
  void InitConfig(Chemical *, Config&);
  void Release();
  
  void InitVecCompart(Chemical*, int, double*, double*, double*, double, double, int, int, double, double, int, double, int, double*, double*, bool);
  void InitMtxCompart(Chemical*, int, double*, double*, double*, double, double, int, int, double, double, int, double, int, double*, double*, bool);

  bool saveFlux(bool, const char [], double*);

};

#endif
