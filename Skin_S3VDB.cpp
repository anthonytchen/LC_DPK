#include "stdafx.h"
#include "Skin_S3VDB.h"

/*!
 */
void Skin_S3VDB::Init(Chemical *chemSolute, int nChem, 
		      double *conc_vehicle, double *partition_vehicle, double *diffu_vehicle, 
		      //double x_len_sb_sur, int n_grids_x_sb_sur, double y_len_sb_har, int n_grids_y_sb_har,
		      double x_len_sb_sur, int n_grids_x_sb_sur, int n_grids_y_sb_sur,
		      int n_grids_x_sb_har, double y_len_sb_har, int n_grids_y_sb_har,		      
		      int n_layer_x_sc, int n_layer_y_sc, double offset_y_sc,
		      int x_len_ve, int n_grids_x_ve, double x_len_de, int n_grids_x_de,
		      double *par_dermis2blood=NULL, double *blood_k_clear=NULL)
{

  Skin::Init();
  
  int i, nxComp, nyComp, dim;
  double coord_x_start, coord_y_start, coord_x_end, coord_y_end;
  double x_len_sc, y_len_sc;

  m_dz_dtheta = 0.01; // fixing dz, the dimension perpendicular to x-y domain

  m_coord_sys = Cylindrical;
  

  // setup compartment matrix
  //  [ SurSB, SurSB ]
  //  [ SC   , HarSB ]
  nxComp = 2; nyComp = 2;
  m_nStraCorn = m_nSebum = 1;
  m_nSurSebum = 2;
  m_nVehicle = m_nViaEpd = m_nDermis = m_nBlood = 0;
  createCompMatrix(nxComp, nyComp);

  // n_grids_y_sb_sur = 5;

  // boundary condition: up, left, right, down
  BdyCondStr bdys_sb_sur1 = {ZeroFlux,  ZeroFlux, FromOther, FromOther};
  BdyCondStr bdys_sb_sur2 = {ZeroFlux,  FromOther, ZeroFlux, FromOther};
#ifdef _DEBUG_ // in debug mode, switch off the sink below sb_har & sc
  BdyCondStr bdys_sb_har = {FromOther,  FromOther, ZeroFlux, ZeroFlux};
  BdyCondStr bdys_sc = {FromOther, ZeroFlux, FromOther, ZeroFlux};
#else
  BdyCondStr bdys_sb_har = {FromOther,  FromOther, ZeroFlux, ZeroConc};
  BdyCondStr bdys_sc = {FromOther, ZeroFlux, FromOther, ZeroConc};
#endif
  
  //BdyCondStr bdys_ve = bdys_sc;
  //BdyCondStr bdys_de = {FromOther, Periodic, Periodic, ZeroFlux};

  m_nChem = nChem;

  /*  set up the compartments */

  y_len_sc = getSCYlen() * n_layer_y_sc;
  dim = 0;

  m_SurSebum = new SurSebum[m_nChem*2];
  m_Sebum = new Sebum[m_nChem*1];

  // SB, surface sebum

  double sursb_init_mass_solid, sursb_k_disv_per_area, sursb_k_rect, sursb_Csat;
  double d_particle, S_particle;
  /*
  d_particle = 1.4e-6; // calculated from BET surface area measured in Davies 1985
  S_particle = 6*d_particle*d_particle; // cuboid  
  
  sursb_k_disv_per_area = 34.3e-6 / S_particle;
  */
  Crystal crystal;
  //crystal.shape = HyperRect;
  crystal.shape = BottomOnly;
  crystal.density = 1.782e3; // kg/m^3, ZnPT, from EC Opinion on ZnPT
  crystal.dim = 2;
  crystal.len[0] = 0.5e-6;// 10e-6; // 0.5e-6;
  crystal.len[1] = 0.5e-6; //0.01e-6;   // 0.5e-6;
  //crystal.len[0] = 10e-6; 
  //crystal.len[1] = 0.01e-6;
  crystal.area = crystal.len[0]*crystal.len[1];

  if (m_coord_sys == Cartesian)
    sursb_init_mass_solid = crystal.area*m_dz_dtheta * crystal.density;
  else if (m_coord_sys == Cylindrical ){
    sursb_init_mass_solid = M_PI * (crystal.len[0]*crystal.len[0]) * m_dz_dtheta / 360; // top view area
    sursb_init_mass_solid *= crystal.len[1] * crystal.density;
  }
  else
    SayBye("Coordinate system not implemented");
  
  sursb_k_disv_per_area = 5e-7;
  sursb_k_rect = 4.235e-6; // from Unilever, coverted to 1/s, note large error because the reaction is not first order
  // sursb_k_rect = 1.736e-5; // from Unilever, only use the first 2 data points (i.e. up to day 1) for estimation
  //sursb_k_rect *= 0; //0.1;
  // sursb_Csat = 40 * 1e-6/(1e-3*0.9105); // 40 ppm, sebum specific gravity is 0.9105
  sursb_Csat = 40 * 1e-3; // 40 ppm, converted to kg/m3

  coord_x_start = 0; coord_y_start = 0;
  createSurSB(chemSolute, coord_x_start, coord_y_start, x_len_sb_sur, y_len_sc, n_grids_x_sb_sur, n_grids_y_sb_sur,
	      bdys_sb_sur1, &coord_x_end, &coord_y_end, 0,
	      crystal, sursb_init_mass_solid, sursb_k_disv_per_area, sursb_k_rect, sursb_Csat);
  // m_SurSebum[0].setGridConc(1.0, 0, 0); // give the first grid some initial concentration
  m_CompIdx[0][0].type = emSurSB;
  m_CompIdx[0][0].pComp = new Comp*[1];
  m_CompIdx[0][0].pComp[0] = &m_SurSebum[0];
  dim += m_SurSebum[0].m_dim;
  // m_SurSebum[0].m_grids[19].m_concChem = 1;

  coord_x_start = 0; coord_y_start = y_len_sc;
  createSurSB(chemSolute, coord_x_start, coord_y_start, x_len_sb_sur, y_len_sb_har, n_grids_x_sb_sur, 1,
	      bdys_sb_sur2, &coord_x_end, &coord_y_end, 1,
	      crystal, -1, -1, sursb_k_rect, -1);
  m_CompIdx[0][1].type = emSurSB;
  m_CompIdx[0][1].pComp = new Comp*[1];
  m_CompIdx[0][1].pComp[0] = &m_SurSebum[1];
  dim += m_SurSebum[1].m_dim;

  // SC
  coord_x_start = coord_x_end; coord_y_start = 0;
  createSC(chemSolute, coord_x_start, coord_y_start, n_layer_x_sc, n_layer_y_sc, offset_y_sc, 
	   bdys_sc, &coord_x_end, &coord_y_end);
  m_CompIdx[1][0].type = emSC;
  m_CompIdx[1][0].pComp = new Comp*[1];
  m_CompIdx[1][0].pComp[0] = &m_StraCorn[0]; 
  dim += m_StraCorn[0].m_dim;

  // SB, hair sebum
  // double n_grids_x_sb_har = 10;
  coord_x_start = coord_x_start; coord_y_start = y_len_sc;
  createSB(chemSolute, coord_x_start, coord_y_start, x_len_sb_sur, y_len_sb_har, n_grids_x_sb_har, n_grids_y_sb_har,
	   bdys_sb_har, &coord_x_end, &coord_y_end, 0);
  m_CompIdx[1][1].type = emSB;
  m_CompIdx[1][1].pComp = new Comp*[1];
  m_CompIdx[1][1].pComp[0] = &m_Sebum[0];
  dim += m_Sebum[0].m_dim;

  /*
  // VE
  coord_x_start = coord_x_end; coord_y_start = 0;
  createVE(chemSolute, coord_x_start, coord_y_start, x_len_ve, y_len_sc, n_grids_x_ve, 1,
	   bdys_ve, &coord_x_end, &coord_y_end);

  // DE
  coord_x_start = coord_x_end; coord_y_start = 0;
  createDE(chemSolute, coord_x_start, coord_y_start, x_len_de, y_len_sc, n_grids_x_de, 1, true,
	   bdys_de, &coord_x_end, &coord_y_end);

  // BD
  createBD(par_dermis2blood, blood_k_clear);
  */

  /* link the compartments through boundary setting
     bdy conditions: up/left/right/down */

  for (i=0; i<m_nChem; i++) {

    m_SurSebum[i*2+0].createBoundary(m_SurSebum[i*2+1].m_nx, m_StraCorn[i].m_ny);
    m_SurSebum[i*2+0].setBoundaryGrids(m_SurSebum[i*2+1].m_grids, m_StraCorn[i].m_grids);

    m_SurSebum[i*2+1].createBoundary(0, m_Sebum[i].m_ny);
    m_SurSebum[i*2+1].setBoundaryGrids(NULL, m_Sebum[i].m_grids);

    m_StraCorn[i].createBoundary(m_Sebum[i].m_nx, 0);
    m_StraCorn[i].setBoundaryGrids(m_Sebum[i].m_grids, NULL);

    m_Sebum[i].createBoundary(0, 0);
    m_Sebum[i].setBoundaryGrids(NULL, NULL);

  }

  // overall dimension
  m_dim_all =  dim;

  // If InitReaction() is not called, set m_React.idx_substrate to -1 to indicate no reaction
  m_React.idx_substrate = -1;
}

void Skin_S3VDB::InitConfig(Chemical *chemSolute, Config &conf)
			    /*int nChem, 
		      double *conc_vehicle, double *partition_vehicle, double *diffu_vehicle, 
		      double x_len_sb_sur, int n_grids_x_sb_sur, double y_len_sb_har, int n_grids_y_sb_har,
		      int n_layer_x_sc, int n_layer_y_sc, double offset_y_sc,
		      int x_len_ve, int n_grids_x_ve, double x_len_de, int n_grids_x_de,
		      double *par_dermis2blood=NULL, double *blood_k_clear=NULL)*/
{
  Init(chemSolute, conf.m_nChem, &conf.m_conc_vehicle, &conf.m_partition_vehicle, &conf.m_diffu_vehicle,       
       conf.m_x_len_sb_sur, conf.m_n_grids_x_sb_sur, conf.m_n_grids_y_sb_sur,
       conf.m_n_grids_x_sb_har, conf.m_y_len_sb_har, conf.m_n_grids_y_sb_har,
       conf.m_n_layer_x_sc, conf.m_n_layer_y_sc, conf.m_offset_y_sc,
       conf.m_x_len_ve, conf.m_n_grids_x_ve, conf.m_x_len_de, conf.m_n_grids_x_de);
}

void Skin_S3VDB::Release(void)
{
  Skin::Release();
}

