#include "stdafx.h"
#include "Skin.h"


/* Initialisation for possibly multiple chemicals */
void Skin::Init(Chemical *chemSolute, int nChem, const bool b_has_compartments[],
		double *conc_vehicle, double *partition_vehicle, double *diffu_vehicle, 
		double *par_dermis2blood, double *blood_k_clear,
		double dx_vehicle, double area_vehicle, double x_len_ve, double x_len_de,
		int n_layer_x_sc, int n_layer_y_sc, int n_grids_x_ve, int n_grids_x_de, double offset_y_sc,
		bool bInfSrc)
{

  int i;

  m_dz_dtheta = 0.01; // fixing dz, the dimension perpendicular to x-y domain
  BdyCond bdy_left_right = Periodic;
  double x_len_sc, y_len_sc, y_len_ve, y_len_de;

  m_nChem = nChem;
  m_b_has_SC = b_has_compartments[0]; 
  m_b_has_VE = b_has_compartments[1]; 
  m_b_has_DE = b_has_compartments[2];
  m_b_has_blood = b_has_compartments[3];

  /* now set up each compartment */

  // 1. set up stratum corneum using fixed geometry

  assert(m_b_has_SC); // currently we have to have SC to make the code work

  if (m_b_has_SC) {
    m_StraCorn = new StraCorn[nChem];

    double g, d, s, t, water_frac_surface;
    g=.075e-6; d=40e-6; s=0.075e-6; t=0.8e-6;
    water_frac_surface = 0.55; // mass fraction of water in stratum corneum

    for (i=0; i<m_nChem; i++) {
      if (m_b_has_VE) {
	m_StraCorn[i].Init(g, d, s, t, m_dz_dtheta, n_layer_x_sc, n_layer_y_sc, offset_y_sc, 
			   Cartesian, FromOther, bdy_left_right, bdy_left_right, FromOther); // bdy conditions: u/l/r/d
	m_StraCorn[i].createBoundary(0, m_ViaEpd[i].m_ny);    
	m_StraCorn[i].setBoundaryGrids(NULL, m_ViaEpd[i].m_grids);
      }
      else {
	m_StraCorn[i].Init(g, d, s, t, m_dz_dtheta, n_layer_x_sc, n_layer_y_sc, offset_y_sc, 
			   Cartesian, FromOther, bdy_left_right, bdy_left_right, ZeroConc); // bdy conditions: u/l/r/d
	m_StraCorn[i].createBoundary(0, 0);    
	m_StraCorn[i].setBoundaryGrids(NULL, NULL);
      }
      m_StraCorn[i].createGrids(chemSolute[i], water_frac_surface);     
    }
    m_dim_sc = m_StraCorn[0].m_nx * m_StraCorn[0].m_ny;
    x_len_sc = n_layer_x_sc*(g+t) + g;
    y_len_sc = n_layer_y_sc*(d+s);
  }
  else {
    m_dim_sc = 0;
    x_len_sc = 0;
  }

  // 2. set up viable epidermis 

  if (m_b_has_VE) {
    assert(m_b_has_SC);
    m_ViaEpd = new ViaEpd[nChem];

    y_len_ve = y_len_sc; // depends on the setup for stratum corneum

    for (i=0; i<m_nChem; i++) {
      if (m_b_has_DE) {
	m_ViaEpd[i].Init(x_len_ve, y_len_ve, m_dz_dtheta, n_grids_x_ve, 1, Cartesian,
			 FromOther, bdy_left_right, bdy_left_right, FromOther);
	m_ViaEpd[i].createBoundary(0, m_Dermis[i].m_ny);    
	m_ViaEpd[i].setBoundaryGrids(NULL, m_Dermis[i].m_grids);
      }
      else {
	m_ViaEpd[i].Init(x_len_ve, y_len_ve, m_dz_dtheta, n_grids_x_ve, 1, Cartesian,
			 FromOther, bdy_left_right, bdy_left_right, ZeroConc);
	m_ViaEpd[i].createBoundary(0, 0);
	m_ViaEpd[i].setBoundaryGrids(NULL, NULL);
      }
      m_ViaEpd[i].createGrids(chemSolute[i], x_len_sc);
    }
    m_dim_ve = m_ViaEpd[0].m_nx * m_ViaEpd[0].m_ny;
  }
  else
    m_dim_ve = 0;

  // 3. set up dermis 

  if (m_b_has_DE) {
    assert(m_b_has_SC && m_b_has_VE);
    m_Dermis = new Dermis[nChem];

    y_len_de = y_len_ve; // depends on the setup for ve

    for (i=0; i<m_nChem; i++) {
      m_Dermis[i].Init(x_len_de, y_len_de, m_dz_dtheta, n_grids_x_de, 1, m_b_has_blood, Cartesian,
		       FromOther, bdy_left_right, bdy_left_right, ZeroFlux);
      m_Dermis[i].createGrids(chemSolute[i], x_len_sc+x_len_ve);

      // setup boundaries
      m_Dermis[i].createBoundary(0, 0);    m_Dermis[i].setBoundaryGrids(NULL, NULL);
    }
    m_dim_de = m_Dermis[0].m_nx * m_Dermis[0].m_ny;
  }
  else
    m_dim_de = 0;

  // 4. set up blood
  if (m_b_has_blood) {
    assert(m_b_has_SC && m_b_has_VE && m_b_has_DE);
    m_Blood = new Blood[nChem];

    for (i=0; i<m_nChem; i++) {
      m_Blood[i].Init(m_Dermis[i].m_grids->m_chemical.m_frac_unbound, blood_k_clear[i], 70, 'M');
      m_Dermis[i].InitDermisBlood(m_Blood[i].m_flow_capil, m_Blood[i].m_f_unbound, par_dermis2blood[i]);
    }
    m_dim_bd = 2;
  }
  else
    m_dim_bd = 0;

  // 5. set up vehicle 

  m_concVehicleInit = new double[nChem];
  m_Vehicle = new Vehicle[nChem];

  m_Vehicle_area = area_vehicle;
  m_bInfSrc = bInfSrc; // whether the vehicle is a infinite source

  for (i=0; i<m_nChem; i++) {
    m_concVehicleInit[i] = conc_vehicle[i];
    if (m_b_has_SC) {
      m_Vehicle[i].Init(dx_vehicle, y_len_sc, m_dz_dtheta, 1, 1, conc_vehicle[i], partition_vehicle[i], diffu_vehicle[i],
			Cartesian, ZeroFlux, bdy_left_right, bdy_left_right, FromOther);
      m_Vehicle[i].createBoundary(0, m_StraCorn[i].m_ny);    
      m_Vehicle[i].setBoundaryGrids(NULL, m_StraCorn[i].m_grids);
    }
    else {
      m_Vehicle[i].Init(dx_vehicle, y_len_sc, m_dz_dtheta, 1, 1, conc_vehicle[i], partition_vehicle[i], diffu_vehicle[i],
			Cartesian, ZeroFlux, bdy_left_right, bdy_left_right, ZeroFlux);
      m_Vehicle[i].createBoundary(0, 0);    
      m_Vehicle[i].setBoundaryGrids(NULL, NULL);
    }
    m_Vehicle[i].createGrids(chemSolute[i], -dx_vehicle); // coordinate 0 starts from stratum corneum
  }
  m_dim_vh = m_Vehicle[0].m_nx * m_Vehicle[0].m_ny;

  // overall dimension
  m_dim_all =  m_dim_vh + m_dim_sc + m_dim_ve + m_dim_de +  m_dim_bd;

  // If InitReaction() is not called, set m_React.idx_substrate to -1 to indicate no reaction
  m_React.idx_substrate = -1;
}

void Skin::InitReaction(int idx_substrate, int idx_product, double Vmax, double Km)
{
  /* Now only one reaction is implemented. 
     m_React needs to be expanded into an array to consider multiple reactions */

  m_React.idx_substrate = idx_substrate;
  m_React.idx_product = idx_product;
  m_React.Vmax = Vmax;
  m_React.Km = Km;
}

void Skin::Release(void)
{
  int i;

  delete [] m_concVehicleInit;
  for (i<0; i<m_nChem; i++)
    m_Vehicle[i].Release();
  delete [] m_Vehicle;

  if(m_b_has_SC) {
    for (i<0; i<m_nChem; i++) 
      m_StraCorn[i].Release();
    delete [] m_StraCorn;
  }

  if(m_b_has_VE) {
    for (i<0; i<m_nChem; i++)
      m_ViaEpd[i].Release();
    delete [] m_ViaEpd;
  }

  if(m_b_has_DE) {
    for (i<0; i<m_nChem; i++)
      m_Dermis[i].Release();
    delete [] m_Dermis;
  }

  if(m_b_has_blood) {
    for (i<0; i<m_nChem; i++)
      m_Blood[i].Release();
    delete [] m_Blood;
  }

}


int Skin::static_cvODE (double t, N_Vector y, N_Vector dydt, void *paras)
{
  double *p_y, *p_dydt;
	
  p_y = NV_DATA_S(y);
  p_dydt = NV_DATA_S(dydt);
  return ((Skin*)paras)->compODE_dydt(t, p_y, p_dydt);
}


int Skin::compODE_dydt (double t, const double y[], double f[])
{
  int i, j, idx;
  double *concBdy = NULL;

  /* y and f are organised as (for each chemical):
     vehicle (1 grid), sc (row dominant), ve (row dominant), de (row dominant), bd */

  for (i=0; i<m_nChem; i++) { // for each chemical

    idx = i*m_dim_all;

    /* re-set all boundary in-out mass transfer to 0 */
    m_Vehicle[i].setBdyMassInOutZero();
    if (m_b_has_SC)
      m_StraCorn[i].setBdyMassInOutZero();
    if (m_b_has_VE)
      m_ViaEpd[i].setBdyMassInOutZero();
    if (m_b_has_DE)
      m_Dermis[i].setBdyMassInOutZero();


    /* compute for vehicle */

    // update the concentration in boundary grids according to y[]
    if (m_b_has_SC) {
      concBdy = new double[m_StraCorn[i].m_ny];
      memcpy(concBdy, y+idx+m_dim_vh, sizeof(double)*m_StraCorn[i].m_ny);
      m_Vehicle[i].setBoundaryConc(NULL, concBdy);
      delete [] concBdy;
    }
    m_Vehicle[i].compODE_dydt(t, y+idx, f+idx); // calculate dydt
    if (m_bInfSrc)                              // infinite source, concentration doesn't change, but above dydt calculation is still needed
      memset(f, 0, sizeof(double)*m_dim_vh);    //    since calling compODE_dydt will calculate the flux across boundaries properly

    if (m_b_has_SC)
      m_Vehicle[i].passBdyMassOut(NULL, &m_StraCorn[i]); // pass the right/down boundary mass transfer to appropriate compartment

    /* compute for stratum corneum */

    if (m_b_has_SC) {
      idx += m_dim_vh;

      // update the concentration in boundary grids according to y[]
      if (m_b_has_VE) {
	concBdy = new double[m_ViaEpd[i].m_ny];
	memcpy(concBdy, y+idx+m_dim_sc, sizeof(double)*m_ViaEpd[i].m_ny);
	m_StraCorn[i].setBoundaryConc(NULL, concBdy);
	delete [] concBdy;
      }
      m_StraCorn[i].compODE_dydt(t, y+idx, f+idx);
      if (m_b_has_VE)
	m_StraCorn[i].passBdyMassOut(NULL, &m_ViaEpd[i]); // pass the right/down boundary mass transfer to appropriate compartment
    }

    /* compute for viable epidermis */

    if (m_b_has_VE) {
      idx += m_dim_sc;

      // update the concentration in boundary grids according to y[]
      if (m_b_has_DE) {
	concBdy = new double[m_Dermis[i].m_ny];
	memcpy(concBdy, y+idx+m_dim_ve, sizeof(double)*m_ViaEpd[i].m_ny);
	m_ViaEpd[i].setBoundaryConc(NULL, concBdy);
	delete [] concBdy;
      }
      m_ViaEpd[i].compODE_dydt(t, y+idx, f+idx);
      if (m_b_has_DE)
	m_ViaEpd[i].passBdyMassOut(NULL, &m_Dermis[i]); // pass the right/down boundary mass transfer to appropriate compartment
    }

    /* compute for dermis */

    if (m_b_has_DE) {
      idx += m_dim_ve;
      if (m_b_has_blood)
	m_Dermis[i].updateBlood(y[idx+m_dim_de]);
      m_Dermis[i].compODE_dydt(t, y+idx, f+idx);
      m_Dermis[i].passBdyMassOut(NULL, NULL); // pass the right/down boundary mass transfer to appropriate compartment
    }

    /* compute for blood */

    if (m_b_has_blood){
      // simulation is for a small skin area, but needs to multiple
      //  the mass transport due to blood flow by the actual topical application area
      double factor = m_Vehicle_area / m_Dermis[i].compTotalArea(0);
      idx += m_dim_de;

      m_Blood[i].updateMassInOutDermis(m_Dermis[i].m_mass_into_dermis, m_Dermis[i].m_mass_outof_dermis, factor);
      m_Blood[i].compODE_dydt(t, y+idx, f+idx);
    }

  } // for i, each chemical species


  /* Compute for reaction */

  int idx_s, idx_p;
  double C_s, C_p, mass_reacted;

  if (m_React.idx_substrate != -1) { // if idx_substrate not -1, reaction needs to be considered

    idx_s = m_React.idx_substrate * m_dim_all + m_dim_vh + m_dim_sc;
    idx_p = m_React.idx_product * m_dim_all + m_dim_vh + m_dim_sc;

    // For each grid in viable epidermis
    for ( i = 0; i < m_ViaEpd[0].m_nx; i++ ){ // verticle direction up to down
      for ( j = 0; j < m_ViaEpd[0].m_ny; j++ ){ // lateral direction left to right

	idx = i*m_ViaEpd[0].m_ny + j;
	C_s = y[idx_s+idx];
	C_p = y[idx_p+idx];

	mass_reacted = m_React.Vmax * C_s / (m_React.Km + C_s);
	f[idx_s+idx] -= mass_reacted;
	f[idx_p+idx] += mass_reacted;
      }
    }


    idx_s += m_dim_ve;
    idx_p += m_dim_ve;

    // For each grid in dermis
    for ( i = 0; i < m_Dermis[0].m_nx; i++ ){ // verticle direction up to down
      for ( j = 0; j < m_Dermis[0].m_ny; j++ ){ // lateral direction left to right

	idx = i*m_Dermis[0].m_ny + j;
	C_s = y[idx_s+idx];
	C_p = y[idx_p+idx];

	mass_reacted = m_React.Vmax * C_s / (m_React.Km + C_s);
	f[idx_s+idx] -= mass_reacted;
	f[idx_p+idx] += mass_reacted;
      }
    }

  }

  return 0;
}

/* add reaction with diffusion */
void Skin::compReaction()
{
  SayBye("not implemented yet");
}


/* Diffusion using method of lines (MoL) */
void Skin::diffuseMoL(double t_start, double t_end)
{		
  int i, j, k, idx, gsl_status, flag;
  double reltol, abstol, t;
  double *y = NULL;
  N_Vector y0;
  void *cvode_mem = NULL;

  /* get current concentration, and set as initial conditions */

  y = new double[m_dim_all*m_nChem];

  for (k=0; k<m_nChem; k++) {
  
    idx = k*m_dim_all;

    // from vehicle
    for ( i=0; i<m_Vehicle[k].m_nx; i++ ) // x direction up to down
      for ( j=0; j<m_Vehicle[k].m_ny; j++ ) // y direction left to right	
	y[idx + i*m_Vehicle[k].m_ny + j] = m_Vehicle[k].m_grids[i*m_Vehicle[k].m_ny + j].m_concChem;

    // from stratum corneum
    if (m_b_has_SC) {
      idx += m_dim_vh;
      for ( i=0; i<m_StraCorn[k].m_nx; i++ ) // x direction up to down
	for ( j=0; j<m_StraCorn[k].m_ny; j++ ) // y direction left to right	
	  y[idx + i*m_StraCorn[k].m_ny + j] = m_StraCorn[k].m_grids[i*m_StraCorn[k].m_ny + j].m_concChem;
    }

    // from viable epidermis
    if (m_b_has_VE) {
      idx += m_dim_sc;
      for ( i=0; i<m_ViaEpd[k].m_nx; i++ ) // x direction up to down
	for ( j=0; j<m_ViaEpd[k].m_ny; j++ ) // y direction left to right	
	  y[idx + i*m_ViaEpd[k].m_ny + j] = m_ViaEpd[k].m_grids[i*m_ViaEpd[k].m_ny + j].m_concChem;
    }

    // from dermis
    if (m_b_has_DE) {
      idx += m_dim_ve;
      for ( i=0; i<m_Dermis[k].m_nx; i++ ) // x direction up to down
	for ( j=0; j<m_Dermis[k].m_ny; j++ ) // y direction left to right	
	  y[idx + i*m_Dermis[k].m_ny + j] = m_Dermis[k].m_grids[i*m_Dermis[k].m_ny + j].m_concChem;
    }

    // from blood
    if (m_b_has_blood) {
      idx += m_dim_de;
      y[idx] = m_Blood[k].m_concChem;
      y[idx+1] = m_Blood[k].m_concCleared;
    }

  } // for k, each chemical

  /* ------------- */

  y0 = N_VMake_Serial(m_dim_all*m_nChem, y);
	
  // Call CVodeCreate to create the solver memory and specify the 
  //	 Backward Differentiation Formula and the use of a Newton iteration
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
		
  // Call CVodeInit to initialize the integrator memory and specify the
  //	right hand side function in y'=f(t,y), the inital time t_start, and the initial condition y0.
  CVodeInit(cvode_mem, static_cvODE, t_start, y0);

  // Call CVodeSStolerances to specify the scalar relative tolerance and scalar absolute tolerance.
  reltol=1e-6; abstol=1e-10;
  CVodeSStolerances(cvode_mem, reltol, abstol);
	
  // Set the pointer to user-defined data
  CVodeSetUserData(cvode_mem, this);
  CVodeSetMaxNumSteps(cvode_mem, 10000*5);
	
  // Call CVBand to specify the CVBAND band linear solver
  //	m_ny is the bandwidth of the banded Jacobian matrix.
  //	Then setup the Jacobian function
  //CVBand(cvode_mem, dim, m_ny, m_ny);
  //CVDlsSetBandJacFn(cvode_mem, static_cvJacobian);
  
  CVDense(cvode_mem, m_dim_all*m_nChem);
  // CVDlsSetDenseJacFn(cvode_mem, static_cvJacobian);
  
  // prepare for the memory
  // m_gsl_ode_Jacobian = new double [dim*dim];
  // memset(m_gsl_ode_Jacobian, 0, sizeof(double)*dim*dim);
  
  CVode(cvode_mem, t_end, y0, &t, CV_NORMAL);
	
	
  /* Extract calculated values */
		  
  for (k=0; k<m_nChem; k++) {

    idx = k*m_dim_all;

    // for vehicle
    for ( i=0; i<m_Vehicle[k].m_nx; i++ ) // x direction up to down
      for ( j=0; j<m_Vehicle[k].m_ny; j++ ) // y direction left to right	
	m_Vehicle[k].m_grids[i*m_Vehicle[k].m_ny + j].m_concChem = NV_Ith_S(y0, idx + i*m_Vehicle[k].m_ny + j); 


    // for stratum corneum
    if (m_b_has_SC) {
      idx += m_dim_vh;
      for ( i=0; i<m_StraCorn[k].m_nx; i++ ) // x direction up to down
	for ( j=0; j<m_StraCorn[k].m_ny; j++ ) // y direction left to right	
	  m_StraCorn[k].m_grids[i*m_StraCorn[k].m_ny + j].m_concChem = NV_Ith_S(y0, idx + i*m_StraCorn[k].m_ny + j);
    }
    
    // for viable epidermis
    if (m_b_has_VE) {
      idx += m_dim_sc;
      for ( i=0; i<m_ViaEpd[k].m_nx; i++ ) // x direction up to down
	for ( j=0; j<m_ViaEpd[k].m_ny; j++ ) // y direction left to right	
	  m_ViaEpd[k].m_grids[i*m_ViaEpd[k].m_ny + j].m_concChem = NV_Ith_S(y0, idx + i*m_ViaEpd[k].m_ny + j);
    }

    // for dermis
    if (m_b_has_DE) {
      idx += m_dim_ve;
      for ( i=0; i<m_Dermis[k].m_nx; i++ ) // x direction up to down
	for ( j=0; j<m_Dermis[k].m_ny; j++ ) // y direction left to right	
	  m_Dermis[k].m_grids[i*m_Dermis[k].m_ny + j].m_concChem = NV_Ith_S(y0, idx + i*m_Dermis[k].m_ny + j);
    }

    // for blood
    if (m_b_has_blood){
      idx += m_dim_de;
      m_Blood[k].m_concChem = NV_Ith_S(y0, idx);
      m_Blood[k].m_concCleared = NV_Ith_S(y0, idx+1);
    }

  } // for k, each chemical

  /* ------------------------- */
  

  N_VDestroy_Serial(y0);  
  CVodeFree(&cvode_mem);
  delete [] y;
}

 // reset vehicle concentration, partition coefficient, diffusivity
void Skin::resetVehicle(double concChem[], double partition_coef[], double diffu_coef[])
{
  int i, j, k, idx;
  for (k=0; k<m_nChem; k++) {
    
    for (i=0; i<m_Vehicle[k].m_nx; i++) {
      for (j=0; j<m_Vehicle[k].m_ny; j++) {
	
	idx = i*m_Vehicle[k].m_ny + j;
	m_Vehicle[k].m_grids[idx].m_concChem = concChem[k];
	m_Vehicle[k].m_grids[idx].m_Kw = partition_coef[k];
	m_Vehicle[k].m_grids[idx].m_D = diffu_coef[k];
      }
    }

  } // for k
}

void Skin::removeVehicle()
{
  double *concChem = new double[m_nChem];
  double *partition_coef = new double[m_nChem];
  double *diffu_coef = new double[m_nChem];

  for (int i=0; i<m_nChem; i++) {
    concChem[i] = 0;
    partition_coef[i] = 1e10; // very large partition coefficient
    diffu_coef[i] = 1e-100;   //   and very small diffusivity effectively remove the vehicle
  }
  resetVehicle(concChem, partition_coef, diffu_coef);

  delete [] concChem;
  delete [] partition_coef;
  delete [] diffu_coef;
}

// Compute flux into stratum corneum
void Skin::compFlux_2sc(double *flux)
{
  int i, j, idx;
  double f, mass_transfer_rate, total_area;
  double *concBdy = NULL;

  for ( i=0; i<m_nChem; i++ ) {

    m_Vehicle[i].setBdyMassInOutZero();

    concBdy = new double[m_StraCorn[i].m_ny];
    for (j=0; j<m_StraCorn[i].m_ny; j++)
      concBdy[j] = m_StraCorn[i].m_grids[j].m_concChem;
    m_Vehicle[i].setBoundaryConc(NULL, concBdy);
    delete [] concBdy;

    idx = (m_Vehicle[i].m_nx-1)*m_Vehicle[i].m_ny;
    mass_transfer_rate = 0;

    for ( j=0; j<m_Vehicle[i].m_ny; j++ ) { // y direction left to right    
      mass_transfer_rate += - m_Vehicle[i].compMassIrregGridsDown( m_Vehicle[i].m_grids[idx+j], m_Vehicle[i].m_grids[idx+j].m_concChem );
    }

    total_area = m_Vehicle[i].compTotalArea(3);
    flux[i] = mass_transfer_rate / total_area;
  } // for i, each chemical

}

// Compute flux from stratum corneum to the layer down (e.g. viable epidermis)
void Skin::compFlux_sc2down(double *flux)
{
  int i, j, idx;
  double f, mass_transfer_rate, area, total_area, conc_this, deriv_this, deriv_other;
  Grid *gridThiis;
  double *concBdy = NULL;

  for ( i=0; i<m_nChem; i++ ) {

    if (m_b_has_VE) {
      m_StraCorn[i].setBdyMassInOutZero();

      concBdy = new double[m_ViaEpd[i].m_ny];
      for (j=0; j<m_ViaEpd[i].m_ny; j++)
	concBdy[j] = m_ViaEpd[i].m_grids[j].m_concChem;
      m_StraCorn[i].setBoundaryConc(NULL, concBdy);
      delete [] concBdy;
    }

    idx = (m_StraCorn[i].m_nx-1)*m_StraCorn[i].m_ny;
    mass_transfer_rate = 0;

    for ( j=0; j<m_StraCorn[i].m_ny; j++ ) { // y direction left to right

      gridThiis = &m_StraCorn[i].m_grids[idx+j];
      conc_this = gridThiis->m_concChem;      

      if (m_b_has_VE)
	mass_transfer_rate += - m_StraCorn[i].compMassIrregGridsDown( *gridThiis, conc_this );
      else { // sc down is sink
	f = gridThiis->compFlux( &m_StraCorn[i].m_gridSink, conc_this, 0, gridThiis->m_dx/2, 0, &deriv_this, &deriv_other);
	area = m_StraCorn[i].compInterArea(*gridThiis, 3);
	mass_transfer_rate += - f*area;
      }
    }

    total_area = m_StraCorn[i].compTotalArea(3);
    flux[i] = mass_transfer_rate / total_area;

  } // for i, each chemical

}


/* Compute the mass (or mol, depending on concentration unit used) of solute in 
   each layer
 */
void Skin::getLayersAmount(double *fLayersAmount, int dim, int idx_chem)
{
  // order of values (9 items):
  //    [initial amount in vehicle] [current amount in vehicle]
  //    [current amount in stratum corneum] [current amount in lipid of stratum corneum] [current amount in corneocyte of stratum corneum]
  //    [current amount in viable epidermis]
  //    [current amount in dermis]
  //    [current amount in blood]
  //    [current amount in skin, i.e. cumulative amount cleared from blood]
  assert( fLayersAmount && dim==9 );

  int i = idx_chem;

  // Compute the amount in vehicle
  fLayersAmount[0] = m_concVehicleInit[i] * m_Vehicle[i].compTotalVolume();
  fLayersAmount[1] = m_Vehicle[i].getAmount();

  m_StraCorn[i].getAmount(fLayersAmount+2, fLayersAmount+3, fLayersAmount+4);
  fLayersAmount[5] = m_ViaEpd[i].getAmount();
  fLayersAmount[6] = m_Dermis[i].getAmount();
  fLayersAmount[7] = m_Blood[i].getAmount();
  fLayersAmount[8] = m_Blood[i].getClearedAmount();
}

	
/*  ++++++++++++++++++++++++++++++++++
	I/O functions
	++++++++++++++++++++++++++++++++++ */

void Skin::getGridsConc(double *fGridsConc, int dim, int idx_chem)
{
  assert( fGridsConc && dim==m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de+m_dim_bd-1 );

  int i = idx_chem;

  // extracting concentration values from the various compartments
  m_Vehicle[i].getGridsConc( fGridsConc, m_dim_vh );
  m_StraCorn[i].getGridsConc( fGridsConc+m_dim_vh, m_dim_sc );
  m_ViaEpd[i].getGridsConc( fGridsConc+m_dim_vh+m_dim_sc, m_dim_ve );
  m_Dermis[i].getGridsConc( fGridsConc+m_dim_vh+m_dim_sc+m_dim_ve, m_dim_de );
  fGridsConc[m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de] = m_Blood[i].getConcChem();

}

void Skin::get1DConcSC(double *ret, int dim_ret, int idx_chem)
{
  assert( dim_ret == m_StraCorn[idx_chem].m_n_layer_x );
  m_StraCorn[idx_chem].comp1DConc();
  for ( int i = 0; i < dim_ret; i ++ )
    ret[i] = m_StraCorn[idx_chem].m_conc1D[i];
}


void Skin::get1DCoordSC(double *ret, int dim_ret, int idx_chem)
{
  assert( dim_ret == m_StraCorn[idx_chem].m_n_layer_x );
  for ( int i = 0; i < dim_ret; i ++ )
    ret[i] = m_StraCorn[idx_chem].m_coord1D[i];
}


void Skin::displayGrids()
{
  m_Vehicle[0].displayGrids();
  if (m_b_has_SC)
    m_StraCorn[0].displayGrids();
  if (m_b_has_VE)
    m_ViaEpd[0].displayGrids();
  if (m_b_has_DE)
    m_Dermis[0].displayGrids();
  if (m_b_has_blood)
    m_Blood[0].displayGrids();

  fflush(stdout);
}


void Skin::saveGrids(bool b_1st_time, const char fn[])
{
  int i;
  char fn_tmp[1024];

  for (i=0; i<m_nChem; i++) {

    sprintf(fn_tmp, "%s_vh_chem%d.txt", fn, i);
    m_Vehicle[i].saveGrids(b_1st_time, fn_tmp);

    if (m_b_has_SC) {
      sprintf(fn_tmp, "%s_sc_chem%d.txt", fn, i);
      m_StraCorn[i].saveGrids(b_1st_time, fn_tmp);
    }

    if (m_b_has_VE) {
      sprintf(fn_tmp, "%s_ve_chem%d.txt", fn, i);
      m_ViaEpd[i].saveGrids(b_1st_time, fn_tmp);
    }

    if (m_b_has_DE) {
      sprintf(fn_tmp, "%s_de_chem%d.txt", fn, i);
      m_Dermis[i].saveGrids(b_1st_time, fn_tmp);
    }

    if (m_b_has_blood) {
      sprintf(fn_tmp, "%s_bd_chem%d.txt", fn, i);
      m_Blood[i].saveConc(b_1st_time, fn_tmp);
    }
  }
}


void Skin::getXCoord(double *coord_x, int dim)
{
  assert( coord_x && dim==m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de+m_dim_bd-1 );

  // extracting X coordinates from the various compartments

  coord_x[0] = -1e10; // vehicle
  m_StraCorn[0].getXCoord( coord_x+m_dim_vh, m_dim_sc );
  m_ViaEpd[0].getXCoord( coord_x+m_dim_vh+m_dim_sc, m_dim_ve );
  m_Dermis[0].getXCoord( coord_x+m_dim_vh+m_dim_sc+m_dim_ve, m_dim_de );
  coord_x[m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de] = 1e10; // blood
}

void Skin::getYCoord(double *coord_y, int dim)
{
  assert( coord_y && dim==m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de+m_dim_bd-1 );

  // extracting Y coordinates from the various compartments

  coord_y[0] = 0; // vehicle
  m_StraCorn[0].getYCoord( coord_y+m_dim_vh, m_dim_sc );
  m_ViaEpd[0].getYCoord( coord_y+m_dim_vh+m_dim_sc, m_dim_ve );
  m_Dermis[0].getYCoord( coord_y+m_dim_vh+m_dim_sc+m_dim_ve, m_dim_de );
  coord_y[m_dim_vh+m_dim_sc+m_dim_ve+m_dim_de] = 0; // blood
}

void Skin::saveCoord(const char fn_x[], const char fn_y[])
{
  if (m_b_has_SC)
    m_StraCorn[0].saveCoord(fn_x, fn_y);
  if (m_b_has_VE)
    m_ViaEpd[0].saveCoord(fn_x, fn_y); 
  if (m_b_has_DE)
    m_Dermis[0].saveCoord(fn_x, fn_y); 
}
/*  END <I/O functions>
	------------------------------ */
