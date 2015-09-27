#include "stdafx.h"
#include "Skin.h"

#define NTHREADS 8

struct pthread_struct {
	Skin *skin_obj;
	double t;
	const double *y;
	double *f;
	int x_start, x_end, y_start, y_end;
};

/* Initialisation for possibly multiple chemicals */
void Skin::Init(Chemical *chemSolute, int nChem, bool b_has_blood,
		double *conc_vehicle, double *diffu_vehicle, double *partition_vehicle, 
		double *par_dermis2blood, double *blood_k_clear,
		double dx_vehicle, double area_vehicle, double x_len_ve, double x_len_de,
		int n_layer_x_sc, int n_layer_y_sc, int n_grids_x_ve, int n_grids_x_de, double offset_y_sc,
		bool bInfSrc)
{

  int i;

  m_nChem = nChem;
  m_b_has_blood = b_has_blood;

  m_concVehicleInit = new double[nChem];
  m_gridVehicle = new Grid[nChem];
  m_StraCorn = new StraCorn[nChem];
  m_ViaEpd = new ViaEpd[nChem];
  m_Dermis = new Dermis[nChem];
  if (m_b_has_blood) 
    m_Blood = new Blood[nChem];
  m_gridSink = new Grid[nChem];

  m_dz = 0.01; // fixing dz, the dimension perpendicular to x-y domain

  m_boundary_cond = 1; // Boundary condition, [0] - zero flux at left/right sides
	               //    [1] - mirror flux at left/right sides (i.e. periodic boundary condition)

  /* set up stratum corneum using fixed geometry */

  double g, d, s, t, water_frac;
  g=.075e-6; d=40e-6; s=0.075e-6*5.21; t=0.8e-6*5.21;
  water_frac = 0.55; // mass fraction of water in stratum corneum

  for (i=0; i<m_nChem; i++) {
    m_StraCorn[i].Init(g, d, s, t, m_dz, n_layer_x_sc, n_layer_y_sc, offset_y_sc, m_boundary_cond);
    m_StraCorn[i].createGrids(chemSolute[i].m_mw, chemSolute[i].m_K_ow, water_frac, conc_vehicle[i], diffu_vehicle[i]);
  }


  /* set up viable epidermis using fixed geometry */

  double x_len_sc, y_len_ve;
  x_len_sc = n_layer_x_sc*(g+t) + g;
  y_len_ve = n_layer_y_sc*(d+s); // depends on the setup for stratum corneum

  for (i=0; i<m_nChem; i++) {
    m_ViaEpd[i].Init(x_len_ve, y_len_ve, m_dz, n_grids_x_ve);
    m_ViaEpd[i].createGrids(chemSolute[i].m_mw, chemSolute[i].m_K_ow, chemSolute[i].m_pKa, chemSolute[i].m_frac_non_ion, chemSolute[i].m_frac_unbound, chemSolute[i].m_acid_base, x_len_sc);
  }

  /* set up dermis using fixed geometry */

  double y_len_de;
  y_len_de = y_len_ve; // depends on the setup for ve (thus also on stratum corneum)

  for (i=0; i<m_nChem; i++) {
    m_Dermis[i].Init(x_len_de, y_len_de, m_dz, n_grids_x_de, m_b_has_blood);
    m_Dermis[i].createGrids(chemSolute[i].m_mw, chemSolute[i].m_K_ow, chemSolute[i].m_pKa, chemSolute[i].m_frac_non_ion, chemSolute[i].m_frac_unbound, chemSolute[i].m_acid_base, x_len_sc+x_len_ve);
  }

  /* set up blood compartment, then set up the blood properties in dermis */
  if (m_b_has_blood)  {
    for (i=0; i<m_nChem; i++) {
      m_Blood[i].Init(m_Dermis[i].m_grids->m_ve_fu, blood_k_clear[i], 70, 'M');
      m_Dermis[i].InitDermisBlood(m_Blood[i].m_flow_capil, m_Blood[i].m_f_unbound, par_dermis2blood[i]);
    }
  }

  /* set up vehicle using fixed geometry */

  m_Vehicle_area = area_vehicle;

  for (i=0; i<m_nChem; i++) {
    m_concVehicleInit[i] = conc_vehicle[i];
    m_gridVehicle[i].Init("SC", conc_vehicle[i], chemSolute[i].m_K_ow, dx_vehicle, y_len_ve, m_dz, diffu_vehicle[i], partition_vehicle[i]);
    m_gridSink[i].Init("SK", 0, chemSolute[i].m_K_ow, 0, 0, m_dz);
  }
  m_bInfSrc = bInfSrc; // whether the vehicle is a infinite source

  // setup the dimensions
  m_dim_vh = 1;
  m_dim_sc = m_StraCorn[0].m_nx * m_StraCorn[0].m_ny;
  m_dim_ve = m_ViaEpd[0].m_nx * m_ViaEpd[0].m_ny;
  m_dim_de = m_Dermis[0].m_nx * m_Dermis[0].m_ny;
  if (m_b_has_blood) 
    m_dim_bd = 2;
  else
    m_dim_bd = 0;
  m_dim_all =  m_dim_vh + m_dim_sc + m_dim_ve + m_dim_de +  m_dim_bd;

  // If InitReaction is not called, set m_React.idx_substrate to -1 to indicate no reaction
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
  for (i<0; i<m_nChem; i++) {
    m_gridVehicle[i].Release();
    m_gridSink[i].Release();
    m_StraCorn[i].Release();
    m_ViaEpd[i].Release();
    m_Dermis[i].Release();
    if (m_b_has_blood)
      m_Blood[i].Release();
  }
  
  delete [] m_concVehicleInit;
  delete [] m_gridVehicle;
  delete [] m_gridSink;
  delete [] m_StraCorn;
  delete [] m_ViaEpd;
  delete [] m_Dermis;
  if (m_b_has_blood)
    delete [] m_Blood;

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
  Grid gridUp, gridDown;

  /* y and f are organised as (for each chemical):
     vehicle (1 grid), sc (row dominant), ve (row dominant), de (row dominant), bd */

  for (i=0; i<m_nChem; i++) { // for each chemical

    idx = i*m_dim_all;

    /* compute for stratum corneum */
    gridUp.set(&m_gridVehicle[i]);  gridUp.m_concChem = y[idx];
    gridDown.set(m_ViaEpd[i].m_grids); gridDown.m_concChem = y[idx+m_dim_vh+m_dim_sc];
    m_StraCorn[i].updateBoundary(&gridUp, &gridDown, NULL, NULL);  // update top (vehicle) and bottom (ve) boundary
    m_StraCorn[i].compODE_dydt(t, y+idx+m_dim_vh, f+idx+m_dim_vh);

    /* compute for vehicle */
    //     this is after the computation for SC because mass transfer from vehicle to SC is calculated in the SC
    if (m_bInfSrc)
      f[idx] = 0; // infinite source, concentration does not change
    else
      f[idx] = -m_StraCorn[i].m_mass_in/(m_gridVehicle[i].m_dx*m_gridVehicle[i].m_dy*m_gridVehicle[i].m_dz);

    /* compute for viable epidermis */
    idx += m_dim_vh+m_dim_sc;
    gridDown.set(m_Dermis[i].m_grids); gridDown.m_concChem = y[idx+m_dim_ve];
    m_ViaEpd[i].updateBoundary(NULL, &gridDown, NULL, NULL, m_StraCorn[i].m_mass_out); // update top (sc) and bottom (de) boundary
    m_ViaEpd[i].compODE_dydt(t, y+idx, f+idx);

    /* compute for dermis */
    idx += m_dim_ve;
    gridDown.set(&m_gridSink[i]);
    m_Dermis[i].updateBoundary(NULL, &gridDown, NULL, NULL, m_ViaEpd[i].m_mass_out); // update top (ve) and bottom (sink) boundary
    if (m_b_has_blood)
      m_Dermis[i].updateBlood(y[idx+m_dim_de]);
    m_Dermis[i].compODE_dydt(t, y+idx, f+idx);

    /* compute for blood */

    if (m_b_has_blood){
      // simulation is for a small skin area, but needs to multiple
      //  the mass transport due to blood flow by the actual topical application area
      double factor = m_Vehicle_area / (m_Dermis[i].m_y_length*m_Dermis[i].m_dz);
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

  return GSL_SUCCESS;	
}

/* add reaction with diffusion */
void Skin::compReaction()
{
}


/* Diffusion using method of lines (MoL) */
void Skin::diffuseMoL(double t_start, double t_end)
{		
  bool bJacobianRequired = false;
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
    y[idx] = m_gridVehicle[k].m_concChem; 

    // from stratum corneum
    idx += m_dim_vh;
    for ( i=0; i<m_StraCorn[k].m_nx; i++ ) // x direction up to down
      for ( j=0; j<m_StraCorn[k].m_ny; j++ ) // y direction left to right	
	y[idx + i*m_StraCorn[k].m_ny + j] = m_StraCorn[k].m_grids[i*m_StraCorn[k].m_ny + j].m_concChem;


    // from viable epidermis
    idx += m_dim_sc;
    for ( i=0; i<m_ViaEpd[k].m_nx; i++ ) // x direction up to down
      for ( j=0; j<m_ViaEpd[k].m_ny; j++ ) // y direction left to right	
	y[idx + i*m_ViaEpd[k].m_ny + j] = m_ViaEpd[k].m_grids[i*m_ViaEpd[k].m_ny + j].m_concChem;

    // from dermis
    idx += m_dim_ve;
    for ( i=0; i<m_Dermis[k].m_nx; i++ ) // x direction up to down
      for ( j=0; j<m_Dermis[k].m_ny; j++ ) // y direction left to right	
	y[idx + i*m_Dermis[k].m_ny + j] = m_Dermis[k].m_grids[i*m_Dermis[k].m_ny + j].m_concChem;

    if (m_b_has_blood) {
      // from blood
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
    m_gridVehicle[k].m_concChem = NV_Ith_S(y0, idx); 

    // for stratum corneum
    idx += m_dim_vh;
    for ( i=0; i<m_StraCorn[k].m_nx; i++ ) // x direction up to down
      for ( j=0; j<m_StraCorn[k].m_ny; j++ ) // y direction left to right	
	m_StraCorn[k].m_grids[i*m_StraCorn[k].m_ny + j].m_concChem = NV_Ith_S(y0, idx + i*m_StraCorn[k].m_ny + j);

    
    // for viable epidermis
    idx += m_dim_sc;
    for ( i=0; i<m_ViaEpd[k].m_nx; i++ ) // x direction up to down
      for ( j=0; j<m_ViaEpd[k].m_ny; j++ ) // y direction left to right	
	m_ViaEpd[k].m_grids[i*m_ViaEpd[k].m_ny + j].m_concChem = NV_Ith_S(y0, idx + i*m_ViaEpd[k].m_ny + j);

    // for dermis
    idx += m_dim_ve;
    for ( i=0; i<m_Dermis[k].m_nx; i++ ) // x direction up to down
      for ( j=0; j<m_Dermis[k].m_ny; j++ ) // y direction left to right	
	m_Dermis[k].m_grids[i*m_Dermis[k].m_ny + j].m_concChem = NV_Ith_S(y0, idx + i*m_Dermis[k].m_ny + j);

    if (m_b_has_blood){
      // for blood
      idx += m_dim_de;
      m_Blood[k].m_concChem = NV_Ith_S(y0, idx);
      m_Blood[k].m_concCleared = NV_Ith_S(y0, idx+1);
    }

  } // for k, each chemical

  /* ------------------------- */
  

  N_VDestroy_Serial(y0);  
  CVodeFree(&cvode_mem);
  //  delete [] m_gsl_ode_Jacobian;
  delete [] y;
}

 // reset vehicle concentration, partition coefficient, diffusivity
void Skin::resetVehicle(double concChem[], double partition_coef[], double diffu_coef[])
{
  for (int i=0; i<m_nChem; i++) {
    m_gridVehicle[i].Init("SC", concChem[i], m_gridVehicle[i].m_K_ow, 
			  m_gridVehicle[i].m_dx, m_gridVehicle[i].m_dy, m_gridVehicle[i].m_dz, 
			  diffu_coef[i], partition_coef[i]);
    m_StraCorn[i].updateBoundary(&m_gridVehicle[i], NULL, NULL, NULL);  // update top (vehicle) boundary for stratum corneum
  }
}

void Skin::removeVehicle()
{
  double *concChem = new double[m_nChem];
  double *partition_coef = new double[m_nChem];
  double *diffu_coef = new double[m_nChem];

  for (int i=0; i<m_nChem; i++) {
    concChem[i] = m_gridVehicle[i].m_concChem;
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
  int i, j;
  double f, conc_this, conc_other, deriv_this, deriv_other, mass_transfer_rate, area, total_area;
  Grid *gridThiis = NULL;


  for ( i=0; i<m_nChem; i++ ) {
    
    mass_transfer_rate = total_area = 0;
  
    for ( j=0; j<m_StraCorn[i].m_ny; j++ ) { // y direction left to right
				
      gridThiis = &m_StraCorn[i].m_grids[j]; 
      conc_this = gridThiis->m_concChem;
      conc_other = m_gridVehicle[i].m_concChem;
      
      f = gridThiis->compFlux( &m_gridVehicle[i], conc_this, conc_other, 
			       gridThiis->m_dx/2, m_gridVehicle[i].m_dx/2, &deriv_this, &deriv_other);
      area = gridThiis->m_dy*gridThiis->m_dz;
    
      total_area += area;
      mass_transfer_rate += f*area;
      //printf("j %d, area %e, total area %e\n", j, area, total_area);    
    }

    flux[i] = mass_transfer_rate / total_area;
  } // for i, each chemical
}

// Compute flux from stratum corneum to the layer down (e.g. viable epidermis)
void Skin::compFlux_sc2down(double *flux)
{
  int i, j, idx;
  double f, conc_this, conc_other, deriv_this, deriv_other, mass_transfer_rate, area, total_area;
  Grid *gridThiis, *gridOther;

  for ( i=0; i<m_nChem; i++ ) {
    
    mass_transfer_rate = total_area = 0;
    idx = (m_StraCorn[i].m_nx-1)*m_StraCorn[i].m_ny;

    //  printf("\n compflux_sc2ve\n");
    assert(m_ViaEpd[i].m_ny==1);
    gridOther = m_ViaEpd[i].m_grids;
    conc_other = gridOther->m_concChem;

    for ( j=0; j<m_StraCorn[i].m_ny; j++ ) { // y direction left to right

      gridThiis = &m_StraCorn[i].m_grids[idx+j]; 
      conc_this = gridThiis->m_concChem;   

      f = gridThiis->compFlux( gridOther, conc_this, conc_other, 
			       gridThiis->m_dx/2, gridOther->m_dx/2, &deriv_this, &deriv_other);
      area = gridThiis->m_dy*gridThiis->m_dz;

      total_area += area;
      mass_transfer_rate += f*area;
    }

    //  printf("\t total area sc2down %e\n", total_area);

    flux[i] = -mass_transfer_rate / total_area;
  } // for i, each chemical

}

// Compute flux from viable epidermis to the layer down (e.g. dermis)
void Skin::compFlux_ve2down(double *flux)
{
  int i, j, idx;
  double f, conc_this, conc_other, deriv_this, deriv_other, mass_transfer_rate, area, total_area;
  Grid *gridThiis, *gridOther;

  for ( i=0; i<m_nChem; i++ ) {

    idx = (m_ViaEpd[i].m_nx-1)*m_ViaEpd[i].m_ny;
    mass_transfer_rate = total_area = 0;

    //  printf("\n compflux_ve2sk\n");
    assert(m_Dermis[i].m_ny==1);
    gridOther = m_Dermis[i].m_grids;
    conc_other = gridOther->m_concChem;

    for ( j=0; j<m_ViaEpd[i].m_ny; j++ ) { // y direction left to right				

      gridThiis = &m_ViaEpd[i].m_grids[idx+j];
      conc_this = gridThiis->m_concChem;

      f = gridThiis->compFlux( gridOther, conc_this, conc_other, 
			       gridThiis->m_dx/2, gridOther->m_dx/2, &deriv_this, &deriv_other);
      area = gridThiis->m_dy*gridThiis->m_dz;

      total_area += area;
      mass_transfer_rate += f*area;
    }
    //  printf("\t total area ve2down %e\n", total_area);
    flux[i] = -mass_transfer_rate / total_area;
  } // for i, each chemical

}


// Compute flux from dermis to the layer down (e.g. deep tissue)
void Skin::compFlux_de2down(double *flux)
{
  int i, j, idx;
  double f, conc_this, conc_other, deriv_this, deriv_other, mass_transfer_rate, area, total_area;
  Grid *gridThiis, *gridOther;

  for ( i=0; i<m_nChem; i++ ) {

    idx = (m_Dermis[i].m_nx-1)*m_Dermis[i].m_ny;
    mass_transfer_rate = total_area = 0;

    assert(m_Dermis[i].m_ny==1);
    gridOther = &m_gridSink[i];
    conc_other = gridOther->m_concChem;

    for ( j=0; j<m_Dermis[i].m_ny; j++ ) { // y direction left to right				

      gridThiis = &m_Dermis[i].m_grids[idx+j];
      conc_this = gridThiis->m_concChem;

      f = gridThiis->compFlux( gridOther, conc_this, conc_other, 
			       gridThiis->m_dx/2, gridOther->m_dx/2, &deriv_this, &deriv_other);
      area = gridThiis->m_dy*gridThiis->m_dz;

      total_area += area;
      mass_transfer_rate += f*area;  
    }

    // printf("\t total area de2down %e\n", total_area);

    flux[i] = -mass_transfer_rate / total_area;
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
  fLayersAmount[0] = m_concVehicleInit[i] * (m_gridVehicle[i].m_dx*m_gridVehicle[i].m_dy*m_gridVehicle[i].m_dz);
  fLayersAmount[1] = m_gridVehicle[i].getConcChem() * (m_gridVehicle[i].m_dx*m_gridVehicle[i].m_dy*m_gridVehicle[i].m_dz);

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
  fGridsConc[0] = m_gridVehicle[i].getConcChem();
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
  m_StraCorn[0].displayGrids();
  m_ViaEpd[0].displayGrids();
  m_Dermis[0].displayGrids();
  m_Blood[0].displayGrids();

  fflush(stdout);
}


void Skin::saveGrids(bool b_1st_time, const char fn[])
{
  int i;
  char fn_tmp[1024];

  for (i=0; i<m_nChem; i++) {

    sprintf(fn_tmp, "%s_vh_chem%d.txt", fn, i);
    saveVehicle(b_1st_time, fn_tmp, i);

    sprintf(fn_tmp, "%s_sc_chem%d.txt", fn, i);
    m_StraCorn[i].saveGrids(b_1st_time, fn_tmp);

    sprintf(fn_tmp, "%s_ve_chem%d.txt", fn, i);
    m_ViaEpd[i].saveGrids(b_1st_time, fn_tmp);

    sprintf(fn_tmp, "%s_de_chem%d.txt", fn, i);
    m_Dermis[i].saveGrids(b_1st_time, fn_tmp);

    sprintf(fn_tmp, "%s_bd_chem%d.txt", fn, i);
    m_Blood[i].saveConc(b_1st_time, fn_tmp);

  }
}

void Skin::saveVehicle(bool b_1st_time, const char fn[], int idx_chem)
{
  FILE *file = NULL;
  if ( b_1st_time )
    file = fopen(fn, "w");
  else 
    file = fopen(fn, "a");
	
  fprintf(file, "%.5e\n", m_gridVehicle[idx_chem].getConcChem());
  fclose(file);
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
  m_StraCorn[0].saveCoord(fn_x, fn_y);
  m_ViaEpd[0].saveCoord(fn_x, fn_y); 
  m_Dermis[0].saveCoord(fn_x, fn_y); 
}
/*  END <I/O functions>
	------------------------------ */
