// LCP.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "skin.h"

/*
extern void __stdcall MEBDFSO(N,T0,HO,Y0,TOUT,TEND,MF,IDID,LOUT,
     LWORK, WORK,LIWORK,IWORK,MAXDER,ITOL,RTOL,ATOL,F,PDERV,
     NSP,NONZ,IPAR,RPAR,IERR);
*/
// extern void __stdcall PYTHAGORAS (float a, float b, float *c);

int main (int argc, char* argv[])
{
	double g, d, s, t, concSource, K_ow, DSource,
		t_simu, t_end, t_inv;
	bool b_1st_save = true;
	int i, n_layer_x, n_layer_y;
	char fn[] = "conc.txt";
	
	clock_t start, end;
	double cpu_time_used;

	g=0.075E-6; d=40E-6; s=0.075E-6; t=0.8E-6;
	// g=0.075E-6; d=40E-6; s=0.075E-6; t=0.8E-6; // these are the real values	
	n_layer_x = 16; //16
	//x_length = 16*(g+t);
	n_layer_y = 2; // 2
	// y_length = 1.2*(d+s); // 
	
	t_end = 900; t_inv = 10; // simulation time and interval ()in seconds

	K_ow = pow(10,1.6); // partition coefficient between octanol and water
	concSource = 0.11; // in mol/L; converted to mass concnetration when initialising Skin
	DSource = 9.12e-10; // diffusivity of solute in vehicle

	Skin _skin;
	_skin.Init( g, d, s, t, K_ow, concSource, DSource, n_layer_x, n_layer_y, t_inv );
	_skin.createGrids();
	_skin.saveCoord("coord_x.txt", "coord_y.txt");
	_skin.displayGrids();
	
	
/*	start = clock();
//	_skin.diffuseMoL(0, 1e-3);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("cpu time for GSL = %e s \n", cpu_time_used);
	

	start = clock();
	_skin.diffuseMoL_cv(0, 1e-7);
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("cpu time for CVODE = %e s \n", cpu_time_used);
	_skin.saveGrids(b_1st_save, fn);
	
	return 0;
*/		
	fflush(stdout);
	


	for ( t_simu=.0; t_simu<t_end; t_simu+=t_inv ){
		start = clock();
		_skin.diffuseMoL_cv(t_simu, t_simu+t_inv);
		//_skin.diffuseMoL(t_simu, t_simu+t_inv);
		
		_skin.saveGrids(b_1st_save, fn);
		if ( b_1st_save )
			b_1st_save = !b_1st_save;
		
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		
		printf("Simulation time is %e, cpu time = %e s \n", t_simu+t_inv, cpu_time_used);
		fflush(stdout);
	}


	_skin.Release();

	return 0;
}

