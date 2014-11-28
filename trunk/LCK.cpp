// LCP.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "arg.h"
#include "skin.h"


int main (int argc, char* argv[])
{
	double g, d, s, t, concSource, K_ow, DSource,
		t_simu, t_end, t_inv;
	bool b_1st_save = true;
	int i, n_layer_x, n_layer_y;
	char fn_coord_x[20], fn_coord_y[20];
	
	static int nDis = 1;
	static char *fn_conc = "conc.txt";
	static char *pre_coord = "coord";
	
	t_end = 900; t_inv = 10; // simulation time and interval ()in seconds
	g=0.075E-6; d=40E-6; s=0.075E-6; t=0.8E-6;
	n_layer_x = 16; //16
	n_layer_y = 2; // 2

	K_ow = pow(10,1.6); // partition coefficient between octanol and water
	concSource = 0.11; // in mol/L; converted to mass concnetration when initialising Skin
	DSource = 9.12e-10; // diffusivity of solute in vehicle
	
	// Provide a command line user interface
    static Config_t params[] = {

      "tinv", "Simulation time interval (s)",
      "-tinv", DOUBLE, (caddr_t)&t_inv,

      "tend", "Simulation end time (s)",
      "-tend", DOUBLE, (caddr_t)&t_end,  
	  
      "fn_conc", "File name to store concentration",
      "-fn_conc", STRING, (caddr_t)&fn_conc,
	  
	  "pre_coord", "Prefix (max 10 chars) to store coordinates to <fn_coord>_x.txt and <fn_coord>_y.txt",
      "-pre_coord", STRING, (caddr_t)&pre_coord,
	  
	  "dis", "Display options; 0 - most parsimonious; 3 - most verbose",
      "-dis", INT, (caddr_t)&nDis,
	  
      0, 0, 0, NOTYPE, 0

    };
    if ( argc < 2) {
      pusage ( argv[0], params );
      exit (0);
    }
    if ( ppconf (argc, argv, params, true) ) {
      pusage (argv[0], params);
      exit (1);
    }
		
	clock_t start, end;
	double cpu_time_used;

	Skin _skin;
	
	_skin.Init( g, d, s, t, K_ow, concSource, DSource, n_layer_x, n_layer_y, t_inv );
	_skin.createGrids();
	
	strcpy(fn_coord_x, pre_coord); strcat(fn_coord_x, "_x.txt");
	strcpy(fn_coord_y, pre_coord); strcat(fn_coord_y, "_y.txt");
	_skin.saveCoord( fn_coord_x, fn_coord_y );
	
	if ( nDis > 1 ) {
		_skin.displayGrids();	
		fflush(stdout);
	}

	for ( t_simu=.0; t_simu<t_end; t_simu+=t_inv ){
		start = clock();
		_skin.diffuseMoL_cv(t_simu, t_simu+t_inv);	
		_skin.saveGrids(b_1st_save, fn_conc);
		if ( b_1st_save )
			b_1st_save = !b_1st_save;
		
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		
		if ( nDis > 0 ) {
			printf("Simulation time is %e, cpu time = %e s \n", t_simu+t_inv, cpu_time_used);
			fflush(stdout);
		}
	}

	_skin.Release();
	return 0;
}

