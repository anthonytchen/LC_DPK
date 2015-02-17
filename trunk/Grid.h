#ifndef _H_GRID_
#define _H_GRID_

struct mass_diffused
{
  // indicates whether diffusion on this direction has been calculated
  bool bUp, bLeft, bRight, bDown;
  double up, left, right, down;
};

/* structure and functions relating to points for creating grids */

struct Point
{
  double x_coord, y_coord, dx, dy;
  char x_type[3], y_type[3];
};
void setPoint(struct Point& pt, double x_coord, double y_coord, double dx, double dy, const char x_type[], const char y_type[])
{
  pt.x_coord = x_coord;
  pt.y_coord = y_coord;
  pt.dx = dx;
  pt.dy = dy;
  strcpy(pt.x_type, x_type);
  strcpy(pt.y_type, y_type);
}
void cpyPoint(struct Point& dst, struct Point& src)
{
  dst.x_coord = src.x_coord;
  dst.y_coord = src.y_coord;
  dst.dx = src.dx;
  dst.dy = src.dy;
  strcpy(dst.x_type, src.x_type);
  strcpy(dst.y_type, src.y_type);
}
/* ------------ */


class Grid
{
public:
  char m_name[3]; // two-letter combination of the name of the grid
  //		'CC' - corneocyte in stratum corneum
  //		'LP' - lipid in stratum corneum
  //		'SC' - source, usually vehicle
  //		'SK' - sink
	
  double m_Kw,  // partition coefficient between this object and water
    m_D, // diffusivity of solute in this object
    m_Dw, // diffusivity of solute in water
    m_K_ow, // partition coefficient of the chemical between octanol and water
    m_mw, // solute molecular weight
    m_pKa, // solute ionisation
    m_phi_b, // volume fraction of water in corneocyte at saturation
    m_theta_b, // volume fraction of water in corneocyte
    m_mass_frac_water, // mass fraction of water in this grid
    m_r_s, // solute radius
    m_r_f; // keratin microfibril radius, in nm

  double m_concChem, // concentration of chemical
    m_concWater; // concentration of water (hydration level)
  struct mass_diffused m_mass_diffused; 
  double m_x_coord, m_y_coord; // the starting x and y coordinate for this grid
  double m_dx, m_dy, m_dz; // grid size in verticle direction, lateral direction, and the 3rd dimension. 
                           //    This is 2D simulation, thus dz only used to calculate the diffusion area
	

public:
  Grid(void) { };
  virtual ~Grid(void) {};
	
  void Init(const char[], double, double, double, double, double, double D_vehicle=-1);
  void Init(const char[], double, double, double, double, double, double, 
	    double, double, double, double, double, double, double, double, double, double, double);
  void InitVE(double, double, double, double, double, double, double, double); // Initialisation for a viable epidermis grid
  void Release() {};
	

  double getConcChem() { return m_concChem; }
  char* getName() { return m_name; }

  void setConcFromDiffMass(void);	
  void set(Grid*);
  double compFlux(Grid*, double, double, double, double, double*, double*);

  // Functions to calculate model parameters
  void compDiffusivity(double D_vehicle=-1);
  void compKcoef(void);
};

#endif
