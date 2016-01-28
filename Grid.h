#ifndef _H_GRID_
#define _H_GRID_

#include "Chemical.h"

/*
struct mass_diffused
{
  // indicates whether diffusion on this direction has been calculated
  bool bUp, bLeft, bRight, bDown;
  double up, left, right, down;
};
*/

/* structure relating to points for creating grids */
struct Point
{
  double x_coord, y_coord, dx, dy;
  char x_type[3], y_type[3];
};
void setPoint(struct Point&, double, double, double, double, const char[], const char[]);
void cpyPoint(struct Point&, struct Point&);
/* ------------ */


class Grid
{
public:
  char m_name[3]; // two-letter combination of the name of the grid
  //            'VH' - vehicle
  //		'CC' - corneocyte in stratum corneum
  //		'LP' - lipid in stratum corneum
  //		'VE' - viable epidermis
  //            'DE' - dermis
  //		'SK' - sink

  Chemical m_chemical;

  double m_Kw,  // partition coefficient between this object and water
    m_D, // diffusivity of solute in this object
    m_Dw, // diffisivity of solute in water
    m_concChem; // concentration of chemical

  double m_x_coord, m_y_coord; // the starting x and y coordinate for this grid
  double m_dx, m_dy, m_dz; // grid size in verticle direction, lateral direction, and the 3rd dimension. 
                           //    This is 2D simulation, thus dz only used to calculate the diffusion area
	

public:
  Grid(void) { };
  virtual ~Grid(void) {};
  void operator=(const Grid &other);

  void Init(const char[], Chemical, double, double, double, double, double, double);
  void InitVH(const char[], Chemical, double, double, double, double, double, double, double, double, double, double diff_vh=-1);
  void InitSK(const char[], Chemical, double, double, double, double, double, double);
  void InitSK();
  void InitSC(const char[], Chemical, double, double, double, double, double, double, 
	    double, double, double, double, double, double, double, double, double, double);
  void InitVE_DE(const char[], Chemical, double, double, double, double, double, double); // Init for viable epidermis or dermis grid
  void Release() {};
	

  double getConcChem() { return m_concChem; }
  char* getName() { return m_name; }
  void setName(const char name[]) { strcpy(m_name, name); };

  double compFlux(Grid*, double, double, double, double, double*, double*);

};

#endif
