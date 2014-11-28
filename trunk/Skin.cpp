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

void Skin::Init(double g, double d, double s, double t, 
				double K_ow, double concSource, double DSource,
				int n_layer_x, int n_layer_y, double tinv)
{	
	double w;
	
	m_rou_lipid = 1e3; // kg m^{-3}
	m_rou_keratin = 1.2e3; // kg m^{-3}
	m_rou_water = 1e3; // kg m^{-3}
	m_dz = 0.01; // in metre
	m_T = 309; // in Kelvin
	m_eta = 7.05E-4; // Pa s,

	m_mw = 119.12; // Da, i.e. g/mol

	m_grids = NULL;
	m_gsl_ode_Jacobian = NULL;

	assert( s < d ); // inter-corneocyte gap must be less than the corneocyte width
	m_geom_g = g; // * 1e6;
	m_geom_d = d; // * 1e6;
	m_geom_s = s; // * 1e6;
	m_geom_t = t; // * 1e6;
	m_w = 8.0; // offset ratio
	m_geom_dm = m_w*(m_geom_d-m_geom_s)/(1+m_w);
	m_geom_dh = m_geom_d-m_geom_s-m_geom_dm;
	
	m_nx_grids_lipid = 2;
	m_nx_grids_cc = 4;
	m_ny_grids_lipid = 2;
	m_ny_grids_cc_dh = 2; 
	
	int n1, n2;
	// Vertical direction, lipid layer is at both top and bottom of the stratum corneum
	m_nx = (m_nx_grids_lipid+m_nx_grids_cc)*n_layer_x + m_nx_grids_lipid; 	
	// Lateral direction, [dh] [s] [dm] [s], here d=dh+dm+s, w=dm/dh
	m_ny = (int) ( m_ny_grids_lipid*2 + m_ny_grids_cc_dh + m_ny_grids_cc_dh*m_w ) * n_layer_y;

	printf("# of grids: [x] %d, [y] %d\n", m_nx, m_ny);
	fflush(stdout);
	
	m_x_length = n_layer_x*(g+t)+g;
	m_y_length = n_layer_y*(d+s);
	m_dt = tinv;

	m_V_mortar = ( g*(d+s)+t*s ) * m_dz;	
	m_V_brick = d*t * m_dz;
	m_V_all = m_V_mortar + m_V_brick;

	m_K_ow = K_ow;
	m_concSource = concSource * m_mw; // converting from mol/L to kg/m^3 (same as g/l)
	m_DSource = DSource;
}

void Skin::Release(void)
{
	if (!m_grids)
		delete [] m_grids;
}


void Skin::createGrids()
{
	bool bOffset = false;
	int i, j, idx, idx_x, idx_y;
	int cc_subtype = 0; // 0 = d_h; 1 = s; 2 = d_m, 3 = s; this order matters
	double dx_lipid, dx_cc, dy_lipid, dy_cc, dy_dm, coord_x, coord_y;
	
	dx_lipid = m_geom_g/m_nx_grids_lipid;
	dx_cc = m_geom_t/m_nx_grids_cc;
	dy_lipid = m_geom_s/m_ny_grids_lipid;
	dy_cc = m_geom_dh/m_ny_grids_cc_dh;
	
	m_grids = new Grid[m_nx*m_ny]; // organised in row dominant

	// Initialise source and sink grids; 
	//	the distance perpendicular to diffusion is zero
	//	todo: the setting of dimensions may need change
	m_gridSource.Init("SC", m_concSource, m_K_ow, 0, 0, m_dz, m_DSource);
	m_gridSink.Init("SK", 0, m_K_ow, 0, 0, m_dz);
	m_gridSinkLeft.Init("SK", 0, m_K_ow, 0, 0, m_dz);
	m_gridSinkRight.Init("SK", 0, m_K_ow, 0, 0, m_dz);

	struct Point current_point;
	setPoint(current_point, 0, 0, dx_lipid, dy_cc, "LP", "LP"); // starting from lipid on the top layer

	idx_x = idx_y = 0;
	coord_x = coord_y = 0;
	double *y = new double[m_nx*m_ny];
	for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
		for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right
			
			idx = i*m_ny + j;
#ifdef _DEBUG_3_			
			if ( i < 2 ){
				if ( j >= m_ny-2 )
					y[idx] = m_concSource/10;
				else
					y[idx] = m_concSource;
			}
			else
				y[idx] = 0;
#endif
			// assign type
			if ( !strcmp(current_point.x_type, "LP") || !strcmp(current_point.y_type, "LP") ) { 
				// entire lipid layer (1st strcmp) or lateral lipid between two coreneocytes
				m_grids[idx].Init("LP", m_mw, 0.55, 0.55,
					m_V_mortar, m_V_brick, m_V_all,
					m_rou_lipid, m_rou_keratin, m_rou_water, m_T, m_eta, m_K_ow, 
					current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, m_dz);
			} else {
				m_grids[idx].Init("CC", m_mw, 0.55, 0.55,
					m_V_mortar, m_V_brick, m_V_all,
					m_rou_lipid, m_rou_keratin, m_rou_water, m_T, m_eta, m_K_ow, 
					current_point.x_coord, current_point.y_coord, current_point.dx, current_point.dy, m_dz);
			}

			// update current_point
			if (j==m_ny-1) { // last element in the lateral direction, move down
				
				idx_x ++;

				if ( !strcmp(current_point.x_type, "LP") ) {
					coord_x += dx_lipid;
					if ( idx_x == m_nx_grids_lipid ) {
						setPoint(current_point, coord_x, 0, dx_cc, dy_cc, "CC", "CC");
						idx_x = 0;
					} else {
						setPoint(current_point, coord_x, 0, dx_lipid, dy_cc, "LP", "LP");
					}
				} else if ( !strcmp(current_point.x_type, "CC") ) {
					coord_x += dx_cc;
					if ( idx_x == m_nx_grids_cc ) {
						setPoint(current_point, coord_x, 0, dx_lipid, dy_cc, "LP", "LP");
						idx_x = 0;
						bOffset = !bOffset;
					} else {
						setPoint(current_point, coord_x, 0, dx_cc, dy_cc, "CC", "CC");
					}
				}
				coord_y = 0;
				idx_y = 0;
				cc_subtype = 0;

			} else { // not the last element in the lateral direction, thus move to the right
			
				idx_y ++;
				
				if ( !strcmp(current_point.x_type, "LP") ) { // current row is lipid					
					switch (cc_subtype) {
						case 0 : // now within dh
							coord_y += dy_cc;
							if (idx_y==m_ny_grids_cc_dh){
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
							if (idx_y==m_ny_grids_cc_dh*m_w){
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
							printf("error: cc_subtype not implemented\n");
							exit(0);
					} // switch-case
				} else if ( !strcmp(current_point.x_type, "CC") ) { // current row is corneocyte
					switch (cc_subtype) {
						case 0 : // now within dh
							coord_y += dy_cc;
							if (idx_y==m_ny_grids_cc_dh){
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
							if (idx_y==m_ny_grids_cc_dh*m_w){
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
	
#ifdef _DEBUG_3_
	double *f = new double[m_nx*m_ny];
	memset(f, 0, sizeof(double)*m_nx*m_ny);	
	
	gslODE (0, y, f, 0, m_nx, 0, m_ny);
	delete [] y;
	delete [] f;
	exit(0);
#endif
}

// Cellular automaton implementation of diffusion
void Skin::diffuseCA()
{
	int i, j;

	/* Assumptions:
		-- The vehicle is a infinite source.
		-- Left/right/down boundaries are sink (concentration always zero)
	*/
	Grid *gridUp, *gridLeft, *gridRight, *gridDown;
	
	gridUp = gridLeft = gridRight = gridDown = NULL;

	// Calculate diffused mass
	for ( i=0; i<m_nx; i++ ) { // x direction up to down
		for ( j=0; j<m_ny; j++ ) { // y direction left to right			
			
			// setup the neighbouring grids
			if ( i==0 ) { // top layer, its top is source
				gridUp = &m_gridSource;
			} else {
				gridUp = &m_grids[(i-1)*m_ny+j];
			}			

			if ( j==0 ) { // left layer, its left is sink
				gridLeft = &m_gridSinkLeft;
				// gridLeft->m_Kw = m_Kw; // so that the sink has the same partition coefficient as this object
			} else {
				gridLeft = &m_grids[i*m_ny+j-1];
			}

			if ( j==m_ny-1 ) { // right layer, its right is sink
				gridRight = &m_gridSinkRight;
				// gridRight->m_Kw = m_Kw;
			} else {
				gridRight = &m_grids[i*m_ny+j+1];
			}

			if ( i==m_nx-1 ) { // bottom layer, its down is sink
				gridDown = &m_gridSink;
				// gridDown->m_Kw = m_Kw;
			} else {
				gridDown = &m_grids[(i+1)*m_ny+j];
			}

			m_grids[i*m_ny + j].diffuse(*gridUp, *gridLeft, *gridRight, *gridDown, m_dt);
			
		} // for j
	} // for i

	// Update concentrations
	for ( i=0; i<m_nx; i++ ) { // x direction up to down
		for ( j=0; j<m_ny; j++ ) { // y direction left to right
			m_grids[i*m_ny + j].setConcFromDiffMass();
		} // for j
	} // for i

}

// Method of lines implementation of diffusion
void Skin::diffuseMoL(double t_start, double t_end)
{
	bool bJacobianRequired = false;
	int i, j, gsl_status, dim;
	double *y = NULL;
	gsl_odeiv2_driver *driver = NULL;
	
	dim = m_nx*m_ny;
	
	// setup GSL ODE solver system
	if (bJacobianRequired){ // most implicit methods need Jacobian to solve the stiff problem
		m_gsl_ode_Jacobian = new double [dim*dim];
		memset(m_gsl_ode_Jacobian, 0, sizeof(double)*dim*dim);
		gsl_odeiv2_system sys = {static_gslODE, static_gslJacobian, m_nx*m_ny, this};
		driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk1imp, 1e-5, 1e-4, 0.0);
		// gsl_odeiv2_step_msbdf, gsl_odeiv2_step_bsimp
	} else { // rk8pd (and other explicit methods) does not need Jacobian
		gsl_odeiv2_system sys = {static_gslODE, NULL, m_nx*m_ny, this};
		driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-5, 1e-4, 0.0);
	} // gsl_odeiv2_step_rk8pd
	
	// get current concentration
	y = new double[dim];
	for ( i=0; i<m_nx; i++ ){ // x direction up to down
		for ( j=0; j<m_ny; j++ ){ // y direction left to right	
			y[i*m_ny + j] = m_grids[i*m_ny + j].m_concChem;
		}
	}
	
	gsl_status = gsl_odeiv2_driver_apply(driver, &t_start, t_end, y);
	assert(gsl_status==GSL_SUCCESS);
	
	for ( i=0; i<m_nx; i++ ){ // x direction up to down
		for ( j=0; j<m_ny; j++ ){ // y direction left to right	
			m_grids[i*m_ny + j].m_concChem = y[i*m_ny + j];
		}
	}
	
	gsl_odeiv2_driver_free(driver);
	delete [] y;
	if (bJacobianRequired)
		delete [] m_gsl_ode_Jacobian;
}

int Skin::static_gslJacobian(double t, const double y[], double *dfdy, double dfdt[], void *paras)
{
	return ((Skin*)paras)->gslJacobian(t, y, dfdy, dfdt);
}

int Skin::gslJacobian(double t, const double y[], double *dfdy, double dfdt[])
{
	int i, dim;
	dim = m_nx*m_ny;
	
	// memset(dfdy, 0, sizeof(double)*dim*dim);
	assert(m_gsl_ode_Jacobian!=NULL);
	gslODE(t, y, dfdt); // dfdt here is a dummy place holder
	memcpy(dfdy, m_gsl_ode_Jacobian, sizeof(double)*dim*dim);
	
	// df/dt = 0
	memset(dfdt, 0, sizeof(double)*dim);
	// for ( i=0; i< m_nx*m_ny; i++ )
	// 	dfdt[i] = 0;
	//printf("Jacobian called\n");
		
	return GSL_SUCCESS;
}


int Skin::static_cvJacobian (long int N, double t, 
		N_Vector y, N_Vector dydt, DlsMat Jac, void *paras,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	int i, j, k, left, right, nx, ny;
	double *p_y, *p_dydt, *p_fy, *kthCol;
	
	p_fy = new double[N];
	
	p_y = NV_DATA_S(y);
	p_dydt = NV_DATA_S(dydt);
	
	// todo: a lot of repeated calculations; to be improved
	((Skin*)paras)->gslODE(t, p_y, p_fy); // p_fy here is a dummy place holder
	/*	
	for ( i = 0; i < N; i ++ ) {
		for ( j = 0; j < N; j ++ ) {
			printf("%e\t", BAND_ELEM(Jac, i, j));
		}
		printf("\n");
	}
	printf("============ \n ===========");
	*/
	for ( i = 0; i < N; i ++ ) {
/*
		for ( j=left; j<=right; j++ ) {
			BAND_ELEM(Jac, i, j) = ((Skin*)paras)->m_gsl_ode_Jacobian[i*N + j];
			if (i == 0)
				printf("i=%d, j=%d, %e\n", i, j, BAND_ELEM(Jac, i, j));
		}
*/		
		for ( j = 0; j < N; j ++ ) {
			DENSE_ELEM(Jac, i, j) = ((Skin*)paras)->m_gsl_ode_Jacobian[i*N + j];
			//printf("%e\t", BAND_ELEM(Jac, i, j));
		}
		//printf("\n");
	}
	
	printf("------ ---- \n");
	for ( i = 0; i < N; i ++ ) {
		for ( j = 0; j < N; j ++ ) {
			printf("%e\t", DENSE_ELEM(Jac, i, j));
		}
		printf("\n");
	}
	
	printf("------ ---- \n");
	
	for ( i = 0; i < N; i ++ ) {
		for ( j = 0; j < N; j ++ ) {
			printf("%e\t", ((Skin*)paras)->m_gsl_ode_Jacobian[i*N + j]);
		}
		printf("\n");
	}

	
	delete [] p_fy;
	// exit(0);
}

int Skin::static_cvJacobian (long int N, long int mupper, long int mlower, double t, 
		N_Vector y, N_Vector dydt, DlsMat Jac, void *paras,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	int i, j, k, left, right, nx, ny;
	double *p_y, *p_dydt, *p_fy, *kthCol;
	
	p_fy = new double[N];
	
	p_y = NV_DATA_S(y);
	p_dydt = NV_DATA_S(dydt);
	
	// todo: a lot of repeated calculations; to be improved
	((Skin*)paras)->gslODE(t, p_y, p_fy); // p_fy here is a dummy place holder
	nx = ((Skin*)paras)->m_nx;
	ny = ((Skin*)paras)->m_nx;
	/*	
	for ( i = 0; i < N; i ++ ) {
		for ( j = 0; j < N; j ++ ) {
			printf("%e\t", BAND_ELEM(Jac, i, j));
		}
		printf("\n");
	}
	printf("============ \n ===========");
	*/
	for ( i = 0; i < N; i ++ ) {
		left = i-mlower;
		if ( left < 0 ) left = 0;
		right = i+mupper;
		if ( right > N-1 ) right = N-1;
/*
		for ( j=left; j<=right; j++ ) {
			BAND_ELEM(Jac, i, j) = ((Skin*)paras)->m_gsl_ode_Jacobian[i*N + j];
			if (i == 0)
				printf("i=%d, j=%d, %e\n", i, j, BAND_ELEM(Jac, i, j));
		}
*/		
		for ( j = 0; j < N; j ++ ) {
			BAND_ELEM(Jac, i, j) = ((Skin*)paras)->m_gsl_ode_Jacobian[i*N + j];
			//printf("%e\t", BAND_ELEM(Jac, i, j));
		}
		//printf("\n");
	}
/*	
	printf("------ ---- \n");
	for ( i = 0; i < N; i ++ ) {
		for ( j = 0; j < N; j ++ ) {
			printf("%e\t", BAND_ELEM(Jac, i, j));
		}
		printf("\n");
	}
	
	printf("------ ---- \n");
	
	for ( i = 0; i < N; i ++ ) {
		for ( j = 0; j < N; j ++ ) {
			printf("%e\t", ((Skin*)paras)->m_gsl_ode_Jacobian[i*N + j]);
		}
		printf("\n");
	}
*/
	
	delete [] p_fy;
	// exit(0);
}

int Skin::static_gslODE(double t, const double y[], double f[], void *paras)
{
	return ((Skin*)paras)->gslODE(t, y, f);
}

int Skin::static_cvODE (double t, N_Vector y, N_Vector dydt, void *paras)
{
	double *p_y, *p_dydt;
	
	p_y = NV_DATA_S(y);
	p_dydt = NV_DATA_S(dydt);
	return ((Skin*)paras)->gslODE(t, p_y, p_dydt);
}
	
void* Skin::static_gslODE_threads(void *paras)
{
	struct pthread_struct p = *((struct pthread_struct *) paras);
	p.skin_obj->gslODE (p.t, p.y, p.f, p.x_start, p.x_end, p.y_start, p.y_end);
}

// todo: remove repeated calculation of mass transfer
int Skin::gslODE (double t, const double y[], double f[])
{
	int i, rc;
	
	if (NTHREADS==1) {
		gslODE (t, y, f, 0, m_nx, 0, m_ny);
	} else {		
		struct pthread_struct p[NTHREADS];
		pthread_t threads[NTHREADS];
			
		for ( i=0; i < NTHREADS; i++ ) {
			p[i].skin_obj = this;
			p[i].t=t; p[i].y=y; p[i].f=f;			
			p[i].x_start=0; p[i].x_end=m_nx;
			
			p[i].y_start=i*m_ny/NTHREADS; p[i].y_end=(i+1)*m_ny/NTHREADS;
			pthread_create(&threads[i], NULL, static_gslODE_threads, (void *) &p[i]);
		}
		
		for ( i=0; i<NTHREADS; i++ ) {
			rc = pthread_join(threads[i], NULL);
			assert(rc==0);
		}		
	}
	
	return GSL_SUCCESS;
}

void Skin::gslODE (double t, const double y[], double f[], 
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
#ifdef _DEBUG_3_
	printf("\n =============== t = %lf ================\n", t);
	printf("\n === y ====\n");
	for ( i=idx_x_start; i<idx_x_end; i++ ) { // x direction up to down
		for ( j=idx_y_start; j<idx_y_end; j++ ) { // y direction left to right
			printf("%e\t", y[i*m_ny+j]);
		}
		printf("\n");
	}
	printf("\n === f ====\n");
#endif	
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
			flux = compFlux( gridThiis, gridUp, conc_this, conc_other, 
				gridThiis->m_dx/2, gridUp->m_dx/2, &deriv_this, &deriv_other);
			//if ( i==2 && j < 3 )
			//	printf("flux=%e\t", flux);
				
			mass_transfer_rate += gridThiis->m_dy*gridThiis->m_dz * flux;			
			if (m_gsl_ode_Jacobian!=NULL) {
				deriv_this_sum = deriv_this / gridThiis->m_dx;
				if ( i!=0 ) 
					m_gsl_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dx;
			}

			if ( j==0 ) { // left layer, its left is sink
				gridLeft = &m_gridSinkLeft;
				conc_other = m_gridSinkLeft.m_concChem;
				flux = 0; // left impermeable
			} else {
				idx_other = i*m_ny+j-1;
				gridLeft = &m_grids[idx_other];
				conc_other = y[idx_other];
				flux = compFlux( gridThiis, gridLeft, conc_this, conc_other, 
					gridThiis->m_dy/2, gridLeft->m_dy/2, &deriv_this, &deriv_other);
			}			
			//if ( i==2 && j < 3 )
			//	printf("flux=%e\t", flux);
			mass_transfer_rate += gridThiis->m_dx*gridThiis->m_dz * flux;
			if (m_gsl_ode_Jacobian!=NULL) {
				deriv_this_sum += deriv_this / gridThiis->m_dy;
				if ( j!=0 ) 
					m_gsl_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dy;
			}

			if ( j==m_ny-1 ) { // right layer, its right is sink
				gridRight = &m_gridSinkRight;
				conc_other = m_gridSinkRight.m_concChem;
				flux = 0; // left impermeable
			} else {
				idx_other = i*m_ny+j+1;
				gridRight = &m_grids[idx_other];
				conc_other = y[idx_other];
				flux = compFlux( gridThiis, gridRight, conc_this, conc_other, 
					gridThiis->m_dy/2, gridRight->m_dy/2, &deriv_this, &deriv_other);
			}		
			//if ( i==2 && j < 3 )
			//	printf("flux=%e\t", flux);			
			mass_transfer_rate += gridThiis->m_dx*gridThiis->m_dz * flux;
			if (m_gsl_ode_Jacobian!=NULL) {
				deriv_this_sum += deriv_this / gridThiis->m_dy;
				if ( j!=m_ny-1 ) 
					m_gsl_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dy;
			}

			
			if ( i==m_nx-1 ) { // bottom layer, its down is sink
				gridDown = &m_gridSink;
				conc_other = m_gridSink.m_concChem;				
			} else {
				idx_other = (i+1)*m_ny+j;
				gridDown = &m_grids[idx_other];
				conc_other = y[idx_other];
			}
			flux = compFlux( gridThiis, gridDown, conc_this, conc_other, 
				gridThiis->m_dx/2, gridDown->m_dx/2, &deriv_this, &deriv_other);
			//if ( i==2 && j < 3 )
			//	printf("flux=%e\t", flux);
			mass_transfer_rate += gridThiis->m_dy*gridThiis->m_dz * flux;
			if (m_gsl_ode_Jacobian!=NULL) {
				deriv_this_sum += deriv_this / gridThiis->m_dx;
				if ( i!=m_nx-1 ) 
					m_gsl_ode_Jacobian[ idx_this*dim + idx_other ] = deriv_other / gridThiis->m_dx;
			}

			
			f[idx_this] = mass_transfer_rate / volume_this;
			if (m_gsl_ode_Jacobian!=NULL) 
				m_gsl_ode_Jacobian[ idx_this*dim + idx_this ] = deriv_this_sum;
#ifdef _DEBUG_3_
			printf("%e\t", y[idx_this]);
			
			if ( i==2 && j < 3 )
				printf("t=%lf, i=%d, j=%d, mass=%e, f=%e\n", t, i, j, mass_transfer_rate, f[idx_this]);
#endif
		} // for j
		//printf("\n");
	} // for i
	
}

void Skin::diffuseMoL_cv(double t_start, double t_end)
{		
	bool bJacobianRequired = false;
	int i, j, gsl_status, dim, flag;
	double reltol, abstol, t;
	double *y = NULL;
	gsl_odeiv2_driver *driver = NULL;
	N_Vector y0;
	void *cvode_mem = NULL;
	
	dim = m_nx*m_ny;	

	// get current concentration, and set as initial conditions
	y = new double[dim];
	for ( i=0; i<m_nx; i++ ){ // x direction up to down
		for ( j=0; j<m_ny; j++ ){ // y direction left to right	
			y[i*m_ny + j] = m_grids[i*m_ny + j].m_concChem;
		}
	}
	y0 = N_VMake_Serial(dim, y);
	
	// Call CVodeCreate to create the solver memory and specify the 
    //	 Backward Differentiation Formula and the use of a Newton iteration
	// cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
		
	// Call CVodeInit to initialize the integrator memory and specify the
	//	user's right hand side function in y'=f(t,y), the inital time t_start, and
	//	the initial condition y0.
	CVodeInit(cvode_mem, static_cvODE, t_start, y0);

	// Call CVodeSStolerances to specify the scalar relative tolerance
	//	 and scalar absolute tolerance.
	reltol=1e-4; abstol=1e-5;
	// reltol=1; abstol=1;
	CVodeSStolerances(cvode_mem, reltol, abstol);
	
	// Set the pointer to user-defined data
	CVodeSetUserData(cvode_mem, this);
	CVodeSetMaxNumSteps(cvode_mem, 10000*5);
	
	// Call CVBand to specify the CVBAND band linear solver
	//	m_ny is the bandwidth of the banded Jacobian matrix.
	//	Then setup the Jacobian function
	//CVBand(cvode_mem, dim, m_ny, m_ny);
	//CVDlsSetBandJacFn(cvode_mem, static_cvJacobian);
  
	CVDense(cvode_mem, dim);
	// CVDlsSetDenseJacFn(cvode_mem, static_cvJacobian);
  
	// prepare for the memory
	m_gsl_ode_Jacobian = new double [dim*dim];
	memset(m_gsl_ode_Jacobian, 0, sizeof(double)*dim*dim);
	
	CVode(cvode_mem, t_end, y0, &t, CV_NORMAL);
	
	
	//y = NV_DATA_S(y0);
		
	for ( i=0; i<m_nx; i++ ){ // x direction up to down
		for ( j=0; j<m_ny; j++ ){ // y direction left to right	
			m_grids[i*m_ny + j].m_concChem = NV_Ith_S(y0, i*m_ny+j); //y[i*m_ny + j];
		}
	}
	
	N_VDestroy_Serial(y0);  
	CVodeFree(&cvode_mem);
  	delete [] m_gsl_ode_Jacobian;
	delete [] y;
}


// Compute the flux of solute from <other> to <this> grid
//	However, do not use concentration values in the Grid objects,
//	instead use <conc_this> and <conc_other>
double Skin::compFlux(Grid* thiis, Grid* other, double conc_this, double conc_other, 
	double dist_this, double dist_other, double *deriv_this, double *deriv_other)
{
	double flux, K_other2this, tmp1, tmp2;
	
	K_other2this = other->m_Kw / thiis->m_Kw;	
	flux = ( conc_other - K_other2this*conc_this );
	tmp1 = dist_other/other->m_D + K_other2this*dist_this/thiis->m_D;
	flux /= tmp1;
	
	if (deriv_this!=NULL)
		*deriv_this = -K_other2this / tmp1;
	
	if (deriv_other!=NULL)
		*deriv_other = 1 / tmp1;
		
	return flux;
}

/*
	Setting/copying functions
*/

void Skin::setPoint(struct Point& pt, double x_coord, double y_coord, double dx, double dy, char x_type[], char y_type[])
{
	pt.x_coord = x_coord;
	pt.y_coord = y_coord;
	pt.dx = dx;
	pt.dy = dy;
	strcpy(pt.x_type, x_type);
	strcpy(pt.y_type, y_type);
}

void Skin::cpyPoint(struct Point& dst, struct Point& src)
{
	dst.x_coord = src.x_coord;
	dst.y_coord = src.y_coord;
	dst.dx = src.dx;
	dst.dy = src.dy;
	strcpy(dst.x_type, src.x_type);
	strcpy(dst.y_type, src.y_type);
}


/*
	I/O functions
*/

void Skin::displayGrids()
{
	assert( m_grids );

	int i, j, idx;

	for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
		for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right		

			idx = i*m_ny + j;
			// check type
			if ( !strcmp(m_grids[idx].m_name, "LP") )
				printf("L ");
			else
				printf("C ");
			
		} // for j
		printf("\n");
	} // for i

}

void Skin::saveGrids(bool b_1st_time, char fn[])
{
	assert( m_grids );

	FILE *file = NULL;
	int i, j, idx;

	// save grids
	if ( b_1st_time ) {
		file = fopen(fn, "w");
	} else {
		file = fopen(fn, "a");
	}

	for ( i = 0; i < m_nx; i++ ){ // verticle direction up to down
		for ( j = 0; j < m_ny; j++ ){ // lateral direction left to right		

			idx = i*m_ny + j;
			fprintf(file, "%.5e\t", m_grids[idx].getConcChem());
			
		} // for j
		fprintf(file, "\n");
	} // for i

	fclose(file);
}

void Skin::saveCoord(char fn_x[], char fn_y[])
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
