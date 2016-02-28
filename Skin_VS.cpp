#include "stdafx.h"
#include "Skin_VS.h"

void Skin_VS::Init(Chemical *chemSolute, int nChem, 
		      double *conc_vehicle, double *partition_vehicle, double *diffu_vehicle, 
		      double dx_vehicle, double area_vehicle, 
		      int n_layer_x_sc, int n_layer_y_sc, double offset_y_sc,
		      bool bInfSrc)
{
  Skin::Init();

  int i, nxComp, nyComp;
  double coord_x_start, coord_y_start, coord_x_end, coord_y_end;
  double y_len_sc, y_len_ve, y_len_de;

  m_dz_dtheta = 0.01; // fixing dz, the dimension perpendicular to x-y domain

  m_nVehicle = m_nStraCorn = 1;

  // setup compartment matrix
  nxComp = 2; nyComp = 1;
  m_nVehicle = m_nStraCorn = 1;
  createCompMatrix(nxComp, nyComp);

  // boundary condition: up, left, right, down
  BdyCondStr bdys_vh = {ZeroFlux,  Periodic, Periodic, FromOther};
  BdyCondStr bdys_sc = {FromOther, Periodic, Periodic, ZeroConc};

  m_nChem = nChem;

  /* set up the compartments */

  y_len_sc = getSCYlen() * n_layer_y_sc;

  coord_x_start = 0; coord_y_start = 0;
  createVH(chemSolute, conc_vehicle, partition_vehicle, diffu_vehicle,
	   coord_x_start, coord_y_start, dx_vehicle, y_len_sc, area_vehicle,
	   bInfSrc, bdys_vh, &coord_x_end, &coord_y_end);
  m_CompIdx[0][0].type = emVH;
  m_CompIdx[0][0].pComp = new Comp*[1];
  m_CompIdx[0][0].pComp[0] = &m_Vehicle[0];

  coord_x_start = dx_vehicle; coord_y_start = 0;
  createSC(chemSolute, coord_x_start, coord_y_start, n_layer_x_sc, n_layer_y_sc, offset_y_sc, 
	   bdys_sc, &coord_x_end, &coord_y_end);
  m_CompIdx[1][0].type = emSC;
  m_CompIdx[1][0].pComp = new Comp*[1];
  m_CompIdx[1][0].pComp[0] = &m_StraCorn[0];


  /* link the compartments through boundary setting
     bdy conditions: up/left/right/down */

  for (i=0; i<m_nChem; i++) {

    m_Vehicle[i].createBoundary(0, m_StraCorn[i].m_ny);    
    m_Vehicle[i].setBoundaryGrids(NULL, m_StraCorn[i].m_grids);

    m_StraCorn[i].createBoundary(0, 0);
    m_StraCorn[i].setBoundaryGrids(NULL, NULL);

  }

  // set up dimensions 
  m_dim_all =  m_dim_vh + m_dim_sc;
  
  // If InitReaction() is not called, set m_React.idx_substrate to -1 to indicate no reaction
  m_React.idx_substrate = -1;
 
}

void Skin_VS::Release(void)
{
  Skin::Release();
}

