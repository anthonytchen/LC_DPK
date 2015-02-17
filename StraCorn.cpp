#include "stdafx.h"
#include "StraCorn.h"

/* Structure and definition for parallel computing */
#define NTHREADS 8 // number of threads for parallel computing
struct pthread_struct {
	StraCorn *sc_obj;
	double t;
	const double *y;
	double *f;
	int x_start, x_end, y_start, y_end;
};
/* ------------------- */



/*
  Geometric dimensions:
    g: vertical size of lipid layer (typical: 0.075 micron)
    d: lateral size of corneocyte (typical: 40 micron)
    s: lateral space between two corneocytes (typical: 0.075 micron)
    t: vertical size of corneocyte (typical: 0.8 micron)
    dz: depth (perpendicular to x-y domain) of the simulation domain
 */
void StraCorn::Init(double g, double d, double s, double t, double dz,
		    int n_layer_x, int n_layer_y, double offset_y)
{	

  /* set up some constant values */

  // density of lipid, keratin and water
  m_rou_lipid = 1e3; // kg m^{-3}
  m_rou_keratin = 1.2e3; // kg m^{-3}
  m_rou_water = 1e3; // kg m^{-3}

  // dimension related; c.f. Readme.docx for more details
  m_w = 8.0; // offset ratio, 8.0
  m_nx_grids_lipid = 2; // # of x-grids for lipid layer
  m_nx_grids_cc = 4; // # of x-grids for corneocyte layer
  m_ny_grids_lipid = 2; // # of y-grids for lipid layer
  m_ny_grids_cc_dn = 2; // # of y-grids for dn-part of the offset corneocyte layer

  m_T = 309; // temperature (Kelvin)
  m_eta = 7.1E-4; // water viscosity at above temperature (Pa s),

  /* ---- */

  m_grids = NULL;
  m_ode_Jacobian = NULL;

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

  m_x_length = n_layer_x*(g+t)+g;
  m_y_length = n_layer_y*(d+s);

  m_V_mortar = ( g*(d+s)+t*s ) * dz;	
  m_V_brick = d*t * dz;
  m_V_all = m_V_mortar + m_V_brick;

  m_offset_y =  offset_y;
}

void StraCorn::Release(void)
{
  if (!m_grids)
    delete [] m_grids;
}


void StraCorn::createGrids(double MW, double Kow, double water_frac)
{
  bool bOffset = false;
  int i, j, idx, idx_x, idx_y, idx_y_offset, cc_subtype_offset, gsl_errno;
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
    gsl_error ("subtype name unknown", __FILE__, __LINE__, gsl_errno);    
  }

  struct Point current_point;
  setPoint(current_point, 0, 0, dx_lipid, dy_offset, "LP", "LP"); // starting from lipid on the top layer
    
  idx_x = 0; idx_y = idx_y_offset;
  coord_x = coord_y = 0;
  cc_subtype = cc_subtype_offset;

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
			
      idx = i*m_ny + j;

      // assign type
      if ( !strcmp(current_point.x_type, "LP") || !strcmp(current_point.y_type, "LP") ) { 
	// entire lipid layer (1st strcmp) or lateral lipid between two coreneocytes
	m_grids[idx].Init("LP", MW, water_frac, water_frac, m_V_mortar, m_V_brick, m_V_all,
			  m_rou_lipid, m_rou_keratin, m_rou_water, m_T, m_eta, Kow, 
			  current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, m_dz);
      } else {
	m_grids[idx].Init("CC", MW, water_frac, water_frac, m_V_mortar, m_V_brick, m_V_all,
			  m_rou_lipid, m_rou_keratin, m_rou_water, m_T, m_eta, Kow, 
			  current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, m_dz);
      }

      // update current_point
      if (j==m_ny-1) { // last element in the lateral direction, move down
				
	idx_x ++;

	if ( !strcmp(current_point.x_type, "LP") ) {
	  coord_x += dx_lipid;
	  if ( idx_x == m_nx_grids_lipid ) {
	    switch (cc_subtype_offset) {
	    case 0 :
	    case 2 :
	      setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      break;
	    case 1 :
	      if (!bOffset)   setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      else            setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "LP");
	      break;
	    case 3 :
	      if (!bOffset)   setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "LP");
	      else            setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      break;
	    default :
	      break;
	    }	      
	    idx_x = 0;
	  } else {
	    setPoint(current_point, coord_x, 0, dx_lipid, dy_offset, "LP", "LP");
	  }
	} else if ( !strcmp(current_point.x_type, "CC") ) {
	  coord_x += dx_cc;
	  if ( idx_x == m_nx_grids_cc ) {
	    setPoint(current_point, coord_x, 0, dx_lipid, dy_offset, "LP", "LP");
	    idx_x = 0;
	    bOffset = !bOffset;
	  } else {
	    switch (cc_subtype_offset) {
	    case 0 :
	    case 2 :
	      setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      break;
	    case 1 :
	      if (!bOffset)   setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      else            setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "LP");
	      break;
	    case 3 :
	      if (!bOffset)   setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "LP");
	      else            setPoint(current_point, coord_x, 0, dx_cc, dy_offset, "CC", "CC");
	      break;
	    default :
	      break;
	    }
	  }
	}
	coord_y = 0;
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
	    gsl_error("cc_subtype not implemented", __FILE__, __LINE__, gsl_errno);
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
	    printf("error: cc_subtype not implemented\n");
	    exit(0);
	  } // switch-case
	} // if-else-if
      } // if (j==m_ny-1)
      
    } // for j
  } // for i

}

/* functions for computing the right-hand size of the odes */

// todo: remove repeated calculation of mass transfer
void StraCorn::compODE_dydt (double t, const double y[], double f[])
{
  int i, rc;
	
  if (NTHREADS==1) {
    compODE_dydt_block (t, y, f, 0, m_nx, 0, m_ny);
  } else {		
    struct pthread_struct p[NTHREADS];
    pthread_t threads[NTHREADS];
			
    for ( i=0; i < NTHREADS; i++ ) {
      p[i].sc_obj = this;
      p[i].t=t; p[i].y=y; p[i].f=f;			
      p[i].x_start=0; p[i].x_end=m_nx;
			
      p[i].y_start=i*m_ny/NTHREADS; p[i].y_end=(i+1)*m_ny/NTHREADS;
      pthread_create(&threads[i], NULL, static_compODE_dydt_block_threads, (void *) &p[i]);
    }
		
    for ( i=0; i<NTHREADS; i++ ) {
      rc = pthread_join(threads[i], NULL);
      assert(rc==0);
    }		
  }
}

// the static container function needed for using multipe threads
void* StraCorn::static_compODE_dydt_block_threads(void *paras)
{
  struct pthread_struct p = *((struct pthread_struct *) paras);
  p.sc_obj->compODE_dydt_block (p.t, p.y, p.f, p.x_start, p.x_end, p.y_start, p.y_end);
}

// the actual funtion to calculate dy/dy
void StraCorn::compODE_dydt_block (double t, const double y[], double f[], 
				   int idx_x_start, int idx_x_end, int idx_y_start, int idx_y_end)
{
  int i, j, idx_this, idx_other, dim;
  double flux, mass_transfer_rate, conc_this, conc_other, deriv_this, deriv_other,
    deriv_this_sum, volume_this;
	
  dim = m_nx*m_ny;
	
  assert(idx_x_start>=0 && idx_x_start<=m_nx);
  assert(idx_x_end>idx_x_start && idx_x_end<=m_nx);
  assert(idx_y_start>=0 && idx_y_start<=m_ny);
  assert(idx_y_end>idx_y_start && idx_y_end<=m_ny);
	
  /* Assumptions:
     -- The vehicle is a infinite source.
     -- Left/right/down boundaries are sink (concentration always zero)
  */
  Grid *gridThiis, *gridUp, *gridLeft, *gridRight, *gridDown;
	
  gridThiis = gridUp = gridLeft = gridRight = gridDown = NULL;
	
  // Calculate diffused mass
  for ( i=idx_x_start; i<idx_x_end; i++ ) { // x direction up to down
    for ( j=idx_y_start; j<idx_y_end; j++ ) { // y direction left to right
				
      mass_transfer_rate = 0;
      idx_this = i*m_ny+j;
			
      gridThiis = &m_grids[idx_this];
      conc_this = y[idx_this];
      volume_this = gridThiis->m_dx * gridThiis->m_dy * gridThiis->m_dz;
			
      // Setup the neighbouring grids
      // 	and calculate the mass transfer rate.
      // If Jacobian required, calculate the Jacobian in the following order:
      //	up, left, self, right, down.
			
      if ( i==0 ) { // top layer, its top is source
	gridUp = &m_gridSource;
	conc_other = m_gridSource.m_concChem;
      } else {
	idx_other = (i-1)*m_ny+j;
	gridUp = &m_grids[idx_other];
	conc_other = y[idx_other];
      }
      flux = gridThiis->compFlux( gridUp, conc_this, conc_other,
				  gridThiis->m_dx/2, gridUp->m_dx/2, &deriv_this, &deriv_other);
      /*
      if ( i==m_nx )
	printf("i=%d, j=%d, mass=%e\t", 
	       i, j, gridThiis->m_dy*gridThiis->m_dz * flux);
      */			
      mass_transfer_rate += gridThiis->m_dy*gridThiis->m_dz * flux;			
      if (m_ode_Jacobian!=NULL) {
	deriv_this_sum = deriv_this / gridThiis->m_dx;
	if ( i!=0 ) 
	  m_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dx;
      }

      if ( j==0 ) { // leftmost layer, its left is sink
	
	if ( m_boundary_cond == 0 ) {
	  gridLeft = &m_gridSinkLeft;
	  conc_other = m_gridSinkLeft.m_concChem;
	  flux = 0; // left impermeable
	} else if ( m_boundary_cond == 1 ) {
	  idx_other = i*m_ny+m_ny-1;
	  gridLeft = &m_grids[idx_other];
	  conc_other = y[idx_other];
	  flux = gridThiis->compFlux( gridLeft, conc_this, conc_other, 
				      gridThiis->m_dy/2, gridLeft->m_dy/2, &deriv_this, &deriv_other);	  
	}
	
      } else {
	idx_other = i*m_ny+j-1;
	gridLeft = &m_grids[idx_other];
	conc_other = y[idx_other];
	flux = gridThiis->compFlux( gridLeft, conc_this, conc_other, 
				    gridThiis->m_dy/2, gridLeft->m_dy/2, &deriv_this, &deriv_other);
      }	
      /*
      if ( i==m_nx )
	printf("mass=%e\t", gridThiis->m_dx*gridThiis->m_dz * flux);
      */
	
      mass_transfer_rate += gridThiis->m_dx*gridThiis->m_dz * flux;
      if (m_ode_Jacobian!=NULL) {
	deriv_this_sum += deriv_this / gridThiis->m_dy;
	if ( j!=0 ) 
	  m_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dy;
      }

      if ( j==m_ny-1 ) { // right layer, its right is sink

	if ( m_boundary_cond == 0 ) {
	  gridRight = &m_gridSinkRight;
	  conc_other = m_gridSinkRight.m_concChem;
	  flux = 0; // left impermeable
	} else if ( m_boundary_cond == 1 ) {
	  idx_other = i*m_ny;
	  gridRight = &m_grids[idx_other];
	  conc_other = y[idx_other];
	  flux = gridThiis->compFlux( gridRight, conc_this, conc_other, 
				      gridThiis->m_dy/2, gridRight->m_dy/2, &deriv_this, &deriv_other);
	}
	
      } else {
	idx_other = i*m_ny+j+1;
	gridRight = &m_grids[idx_other];
	conc_other = y[idx_other];
	flux = gridThiis->compFlux( gridRight, conc_this, conc_other, 
				    gridThiis->m_dy/2, gridRight->m_dy/2, &deriv_this, &deriv_other);
      }	

      /*
      if ( i==m_nx )
	printf("mass=%e\t", gridThiis->m_dx*gridThiis->m_dz * flux);
      */

      mass_transfer_rate += gridThiis->m_dx*gridThiis->m_dz * flux;
      if (m_ode_Jacobian!=NULL) {
	deriv_this_sum += deriv_this / gridThiis->m_dy;
	if ( j!=m_ny-1 ) 
	  m_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dy;
      }

			
      if ( i==m_nx-1 ) { // bottom layer, its down is sink
	gridDown = &m_gridSink;
	conc_other = m_gridSink.m_concChem;				
      } else {
	idx_other = (i+1)*m_ny+j;
	gridDown = &m_grids[idx_other];
	conc_other = y[idx_other];
      }
      flux = gridThiis->compFlux( gridDown, conc_this, conc_other, 
				  gridThiis->m_dx/2, gridDown->m_dx/2, &deriv_this, &deriv_other);
      /*
      if ( i==m_nx )
	printf("mass=%e\n", gridThiis->m_dy*gridThiis->m_dz * flux);
      */

      mass_transfer_rate += gridThiis->m_dy*gridThiis->m_dz * flux;
      if (m_ode_Jacobian!=NULL) {
	deriv_this_sum += deriv_this / gridThiis->m_dx;
	if ( i!=m_nx-1 ) 
	  m_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dx;
      }

			
      f[idx_this] = mass_transfer_rate / volume_this;
      if (m_ode_Jacobian!=NULL) 
	m_ode_Jacobian[ idx_this*dim + idx_this ] = deriv_this_sum;

    } // for j
    //printf("\n");
  } // for i
	
}
/* ----------------- */

	
/*  ++++++++++++++++++++++++++++++++++
	I/O functions
	++++++++++++++++++++++++++++++++++ */

void StraCorn::displayGrids()
{
  assert( m_grids );

  int i, j, idx, gsl_errno;;
  printf("# of grids: [x] %d, [y] %d\n", m_nx, m_ny);

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right	

      idx = i*m_ny + j;			
      if ( !strcmp(m_grids[idx].m_name, "LP") )
	printf("L ");
      else if ( !strcmp(m_grids[idx].m_name, "CC") )
	printf("C ");
      else
	gsl_error ("subtype name unknown", __FILE__, __LINE__, gsl_errno); 
				
    } // for j
    printf("\n");
  } // for i
  fflush(stdout);
}

void StraCorn::saveGrids(bool b_1st_time, const char fn[])
{
	assert( m_grids );

	FILE *file = NULL;
	int i, j, idx;

	// save grids
	if ( b_1st_time )
		file = fopen(fn, "w");
	else 
		file = fopen(fn, "a");
	
	for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
		for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right		

			idx = i*m_ny + j;
			fprintf(file, "%.5e\t", m_grids[idx].getConcChem());
			
		} // for j
		fprintf(file, "\n");
	} // for i

	fclose(file);
}

void StraCorn::saveCoord(const char fn_x[], const char fn_y[])
{
  assert( m_grids );

  FILE *file_x, *file_y;
  int i, j, idx;

  // save grids
  file_x = fopen(fn_x, "w");
  file_y = fopen(fn_y, "w");

  for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
    for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right		

      idx = i*m_ny + j;
      fprintf(file_x, "%.5e\t", m_grids[idx].m_x_coord);
      fprintf(file_y, "%.5e\t", m_grids[idx].m_y_coord);
			
    } // for j
    fprintf(file_x, "\n");
    fprintf(file_y, "\n");
  } // for i

  fclose(file_x);
  fclose(file_y);
}
/*  END <I/O functions>
	------------------------------ */