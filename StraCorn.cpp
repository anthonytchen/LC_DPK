#include "stdafx.h"
#include "StraCorn.h"


/*
  Geometric dimensions:
    g: vertical size of lipid layer (typical: 0.075 micron)
    d: lateral size of corneocyte (typical: 40 micron)
    s: lateral space between two corneocytes (typical: 0.075 micron)
    t: vertical size of corneocyte (typical: 0.8 micron)
    dz: depth (perpendicular to x-y domain) of the simulation domain
 */
void StraCorn::Init(double g, double d, double s, double t, double dz_dtheta,
		    int n_layer_x, int n_layer_y, double offset_y, 
		    CoordSys coord_sys, BdyCond bdy_cond_up, BdyCond bdy_cond_left, BdyCond bdy_cond_right, BdyCond bdy_cond_down)
{	
  // call Init of the base class Comp
  Comp::Init (coord_sys, dz_dtheta, bdy_cond_up, bdy_cond_left, bdy_cond_right, bdy_cond_down);
  
  /* set up some constant values */

  // density of lipid, keratin and water
  m_rou_lipid = 0.9e3; // kg m^{-3}
  m_rou_keratin = 1.37e3; // kg m^{-3}
  m_rou_water = 1e3; // kg m^{-3}

  // dimension related; c.f. Readme.docx for more details
  m_w = 8.0; // offset ratio, 8.0
  m_nx_grids_lipid = 2; // # of x-grids for lipid layer, 2
  m_nx_grids_cc = 4; // # of x-grids for corneocyte layer, 4
  m_ny_grids_lipid = 2; // # of y-grids for lipid layer, 2
  m_ny_grids_cc_dn = 2; // # of y-grids for dn-part of the offset corneocyte layer, 2

  m_T = 309; // temperature (Kelvin)
  m_eta = 7.1E-4; // water viscosity at above temperature (Pa s),

  /* ---- */

  assert( s < d ); // inter-corneocyte gap must be less than the corneocyte width
  m_geom_g = g; 
  m_geom_d = d; 
  m_geom_s = s; 
  m_geom_t = t; 

  m_geom_dm = m_w*(m_geom_d-m_geom_s)/(1+m_w);
  m_geom_dn = m_geom_d-m_geom_s-m_geom_dm;
	
  // Vertical direction, lipid layer is at both top and bottom of the stratum corneum
  m_nx = (m_nx_grids_lipid+m_nx_grids_cc)*n_layer_x + m_nx_grids_lipid; 	
  // Lateral direction, [dh] [s] [dm] [s], here d=dh+dm+s, w=dm/dh
  m_ny = (int) ( m_ny_grids_lipid*2 + m_ny_grids_cc_dn + m_ny_grids_cc_dn*m_w ) * n_layer_y;
  m_dim = m_nx * m_ny;

  m_x_length = n_layer_x*(g+t)+g;
  m_y_length = n_layer_y*(d+s);

  // for the volume fraction calculation, we assume a Cartesian coordinate
  //   and the z-direction width m_dz_dtheta is directly used.
  //   This won't affect if other coordinates are used
  m_V_mortar = ( g*(d+s)+t*s ) * m_dz_dtheta;
  m_V_brick = d*t * m_dz_dtheta;
  m_V_all = m_V_mortar + m_V_brick;

  m_offset_y =  offset_y;

  // setup the array for 1D concentration and coordinates
  m_conc1D = new double [n_layer_x];
  m_coord1D = new double [n_layer_x];
  for ( int i = 0; i < n_layer_x-1; i ++ )
    m_coord1D[i] = 0.5*( i*(g+t) + (i+1)*(g+t) );
  m_coord1D[n_layer_x-1] = 0.5*( (n_layer_x-1)*(g+t) + n_layer_x*(g+t)+g );
  m_n_layer_x = n_layer_x;
}

void StraCorn::Release()
{
  delete [] m_conc1D;
  delete [] m_coord1D;

  Comp::Release();
}


void StraCorn::createGrids(Chemical chem, double water_frac_surface, double coord_x_start, double coord_y_start)
{
  bool bOffset = false;
  int i, j, idx, idx_x, idx_y, idx_y_offset, cc_subtype_offset;
  int cc_subtype = 0; // 0 = d_n; 1 = s; 2 = d_m, 3 = s;
  double dx_lipid, dx_cc, dy_lipid, dy_cc, dy_offset, dy_dm, coord_x, coord_y;

	
  dx_lipid = m_geom_g/m_nx_grids_lipid;
  dx_cc = m_geom_t/m_nx_grids_cc;
  dy_lipid = m_geom_s/m_ny_grids_lipid;
  dy_cc = m_geom_dn/m_ny_grids_cc_dn;
	
  m_grids = new Grid[m_nx*m_ny]; // organised in row dominant

  // work out the starting point from given offset
  
  double len = m_geom_d+m_geom_s, len_vec[4];
  len_vec[0]=m_geom_dn; len_vec[1]=m_geom_s; len_vec[2]=m_geom_dm; len_vec[3]=m_geom_s;

  while (m_offset_y > len) m_offset_y -= len;
  len = 0;
  for ( i = 0; i < 4; i ++ ){
    len += len_vec[i];
    if (m_offset_y < len) break;
  }
  cc_subtype_offset = i;
  len = m_offset_y - len + len_vec[i];

    
  switch (cc_subtype_offset) {
  case 0 :
  case 2 : // corneocyte width 
    idx_y_offset = (int) floor( len / dy_cc );
    dy_offset = dy_cc;
    break;
  case 1 :
  case 3 : // lipid width
    idx_y_offset = (int) floor( len / dy_lipid );
    dy_offset = dy_lipid;
    break;
  default :
    SayBye ("subtype name unknown");    
  }

  double water_frac, water_frac_sat = 0.55; // saturated water content (w/w)
  double water_increment_per_x = (water_frac_sat - water_frac_surface) / m_x_length;
    
  idx_x = 0; idx_y = idx_y_offset;
  cc_subtype = cc_subtype_offset;

  coord_x = coord_x_start; coord_y = coord_y_start;
  struct Point current_point;
  setPoint(current_point, coord_x, coord_y, dx_lipid, dy_offset, "LP", "LP"); // starting from lipid on the top layer

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    
    water_frac = water_frac_surface + (current_point.x_coord - coord_x_start) * water_increment_per_x;
    // printf("water fraction %lf at %.3e\n", water_frac, current_point.x_coord);

    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
			
      idx = i*m_ny + j;

      // assign type
      if ( !strcmp(current_point.x_type, "LP") || !strcmp(current_point.y_type, "LP") ) { 
	// entire lipid layer (1st strcmp) or lateral lipid between two coreneocytes
	m_grids[idx].InitSC("LP", chem, 0, water_frac, water_frac_sat, m_V_mortar, m_V_brick, m_V_all,
			  m_rou_lipid, m_rou_keratin, m_rou_water, m_T, m_eta, 
			  current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, m_dz_dtheta);
      } else {
	m_grids[idx].InitSC("CC", chem, 0, water_frac, water_frac_sat, m_V_mortar, m_V_brick, m_V_all,
			  m_rou_lipid, m_rou_keratin, m_rou_water, m_T, m_eta, 
			  current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, m_dz_dtheta);
      }

      // update current_point
      if (j==m_ny-1) { // last element in the lateral direction, move down
				
	idx_x ++;
	coord_y = coord_y_start;

	if ( !strcmp(current_point.x_type, "LP") ) {
	  coord_x += dx_lipid;
	  if ( idx_x == m_nx_grids_lipid ) {
	    switch (cc_subtype_offset) {
	    case 0 :
	    case 2 :
	      setPoint(current_point, coord_x, coord_y, dx_cc, dy_offset, "CC", "CC");
	      break;
	    case 1 :
	      if (!bOffset)   setPoint(current_point, coord_x, coord_y, dx_cc, dy_offset, "CC", "CC");
	      else            setPoint(current_point, coord_x, coord_y, dx_cc, dy_offset, "CC", "LP");
	      break;
	    case 3 :
	      if (!bOffset)   setPoint(current_point, coord_x, coord_y, dx_cc, dy_offset, "CC", "LP");
	      else            setPoint(current_point, coord_x, coord_y, dx_cc, dy_offset, "CC", "CC");
	      break;
	    default :
	      break;
	    }	      
	    idx_x = 0;
	  } else {
	    setPoint(current_point, coord_x, coord_y, dx_lipid, dy_offset, "LP", "LP");
	  }
	} else if ( !strcmp(current_point.x_type, "CC") ) {
	  coord_x += dx_cc;
	  if ( idx_x == m_nx_grids_cc ) {
	    setPoint(current_point, coord_x, coord_y, dx_lipid, dy_offset, "LP", "LP");
	    idx_x = 0;
	    bOffset = !bOffset;
	  } else {
	    switch (cc_subtype_offset) {
	    case 0 :
	    case 2 :
	      setPoint(current_point, coord_x, coord_y, dx_cc, dy_offset, "CC", "CC");
	      break;
	    case 1 :
	      if (!bOffset)   setPoint(current_point, coord_x, coord_y, dx_cc, dy_offset, "CC", "CC");
	      else            setPoint(current_point, coord_x, coord_y, dx_cc, dy_offset, "CC", "LP");
	      break;
	    case 3 :
	      if (!bOffset)   setPoint(current_point, coord_x, coord_y, dx_cc, dy_offset, "CC", "LP");
	      else            setPoint(current_point, coord_x, coord_y, dx_cc, dy_offset, "CC", "CC");
	      break;
	    default :
	      break;
	    }
	  }
	}
	// coord_y = 0;
	idx_y = idx_y_offset;
	cc_subtype = cc_subtype_offset;

      } else { // not the last element in the lateral direction, thus move to the right
			
	idx_y ++;
				
	if ( !strcmp(current_point.x_type, "LP") ) { // current row is lipid					
	  switch (cc_subtype) {
	  case 0 : // now within dh
	    coord_y += dy_cc;
	    if (idx_y==m_ny_grids_cc_dn){
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_lipid, "LP", "LP");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_cc, "LP", "LP");
	    }
	    break;
	  case 1 : // now within s
	    coord_y += dy_lipid;
	    if (idx_y==m_ny_grids_lipid){
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_cc, "LP", "LP");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_lipid, "LP", "LP");
	    }
	    break;
	  case 2 : // now wtihin dm
	    coord_y += dy_cc;
	    if (idx_y==m_ny_grids_cc_dn*m_w){
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_lipid, "LP", "LP");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_cc, "LP", "LP");
	    }
	    break;
	  case 3 : // now within the 2nd s
	    coord_y += dy_lipid;
	    if (idx_y==m_ny_grids_lipid){
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_cc, "LP", "LP");
	      idx_y = 0;
	      cc_subtype = 0;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_lipid, dy_lipid, "LP", "LP");
	    }
	    break;
	  default :
	    SayBye("cc_subtype not implemented");
	    exit(-1);
	  } // switch-case
	} else if ( !strcmp(current_point.x_type, "CC") ) { // current row is corneocyte
	  switch (cc_subtype) {
	  case 0 : // now within dh
	    coord_y += dy_cc;
	    if (idx_y==m_ny_grids_cc_dn){
	      if (bOffset)
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "LP");
	      else
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "CC");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_cc, dy_cc, "CC", "CC");
	    }
	    break;
	  case 1 : // now within s
	    coord_y += dy_lipid;
	    if (idx_y==m_ny_grids_lipid){
	      setPoint(current_point, coord_x, coord_y, dx_cc, dy_cc, "CC", "CC");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      if (bOffset)
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "LP");
	      else
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "CC");
	    }
	    break;
	  case 2 : // now wtihin dm
	    coord_y += dy_cc;
	    if (idx_y==m_ny_grids_cc_dn*m_w){
	      if (bOffset)
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "CC");
	      else
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "LP");
	      idx_y = 0;
	      cc_subtype ++;
	    } else {
	      setPoint(current_point, coord_x, coord_y, dx_cc, dy_cc, "CC", "CC");
	    }
	    break;
	  case 3 : // now within the 2nd s
	    coord_y += dy_lipid;
	    if (idx_y==m_ny_grids_lipid){
	      setPoint(current_point, coord_x, coord_y, dx_cc, dy_cc, "CC", "CC");
	      idx_y = 0;
	      cc_subtype = 0;
	    } else {
	      if (bOffset)
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "CC");
	      else
		setPoint(current_point, coord_x, coord_y, dx_cc, dy_lipid, "CC", "LP");
	    }
	    break;
	  default :
	    SayBye("cc_subtype not implemented\n");
	  } // switch-case
	} // if-else-if
      } // if (j==m_ny-1)
      
    } // for j
  } // for i
  //  exit(0);
}


/*  ++++++++++++++++++++++++++++++++++
	I/O functions
	++++++++++++++++++++++++++++++++++ */

void StraCorn::getAmount(double *amount_total, double *amount_lipid, double *amount_corneocyte)
{
  assert( m_grids );

  *amount_total = *amount_lipid = *amount_corneocyte = .0;

  int i, j, idx, gsl_errno;
  double volume, amount;
	
  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
       idx = i*m_ny + j;

       volume = Comp::compVolume(m_grids[idx]);
       amount = m_grids[idx].getConcChem() * volume;

       if ( !strcmp(m_grids[idx].m_name, "LP") )
	 *amount_lipid += amount;
       else if (!strcmp(m_grids[idx].m_name, "CC") )
	 *amount_corneocyte += amount;
       else
	 SayBye ("subtype name unknown");       
     }
  }

  *amount_total = *amount_lipid + *amount_corneocyte;
}


void StraCorn::comp1DConc()
{
  int i, j, idx, nx_grids_lipid, nx_grids_cc, idx_conc1D;
  double conc_sum, area, area_sum;

  assert( m_n_layer_x > 1 );

  idx_conc1D = 0;

  conc_sum = 0; area_sum = 0;
  nx_grids_lipid = m_nx_grids_lipid;
  nx_grids_cc = m_nx_grids_cc;

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right

       idx = i*m_ny + j;

       area = m_grids[idx].m_dx * m_grids[idx].m_dy;
       conc_sum += m_grids[idx].getConcChem() * area;
       area_sum += area;

    } // for j

    if ( nx_grids_lipid > 0 )
      nx_grids_lipid --;
    else
      nx_grids_cc --;

    if ( !nx_grids_lipid && !nx_grids_cc ) { // finished calculation for 1 lipid and 1 corneocyte layer
      m_conc1D[idx_conc1D] = conc_sum / area_sum;
      idx_conc1D ++;
      conc_sum = 0; area_sum = 0;
      if ( idx_conc1D == m_n_layer_x-1 )
	nx_grids_lipid = m_nx_grids_lipid * 2; // also count the botton lipid layer
      else 
	nx_grids_lipid = m_nx_grids_lipid;
      nx_grids_cc = m_nx_grids_cc;
    }

  } // for i

}

void StraCorn::saveCoord(const char fn_x[], const char fn_y[])
{
  Comp::saveCoord(fn_x, fn_y, ".sc");
}
/*  ---- END <I/O functions> ---- */
