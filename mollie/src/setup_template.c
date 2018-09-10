/* Setup file for radiative transfer simulation. This
   file is part of a C program and C syntax is required. */

/* Pathname to the directory where most of the data will be 
 * stored. This directory should be accessible 
 * to the machine where the master is running. 
 * Other files are written by the slaves in the
 * directory where the executable resides. This directory
 * must be cross-mounted to all the machines that are
 * running the slaves. At least "./" is required. */

/*        strcpy(pathname,"/data2/rest_3/keto/L1544/R24/data"); */

        strcpy(pathname,"./");


/* Select how many iterations of lambda iteration. For 
molecular cloud cores of moderate optical depth 10 - 100
try 10 to 30 iterations. Check the convergence with
the print out and the plotting program.*/               

/* NH3 and CH3CN are full LTE only. No point in more than
   one iteration for these 2 molecules. */

	iterations = 1;

/* Set the initial radiation field to either CMB or ZERO or LTE.  
   ZERO means initial Jbar = 0 and I0 = 0
   CMB  means initial Jbar = CMB and I0 = CMB
   LTE  means initial Jbar = Planck function and I0 = CMB */

	initial_radiation = CMB;

/* ----------------------------------------------------*/

/* Select the molecule for the simulation.
Choices:

	H2Oortho
	CO		checked OK with accel 1-8-09
	C13O		checked OK with accel 1-7-09
	C17O		checked OK with accel 1-8-09
	C18O		checked OK with accel 1-7-09
	CS		checked OK with accel 1-8-09
	N2Hplus		checked OK with accel 1-7-09
	N2Dplus		checked OK with accel 1-8-09
	HCOplus		checked OK with accel 1-19-09
	H13COplus	checked OK with accel 1-19-09
	SiO		checked OK with accel 1-19-09
	HCN		checked OK with accel 1-19-09

All the molecules above have either no hyperfine structure or 
use the approximation of hyperfine statistical equilibrium (LTE hyperfines)

N2Hhyp below is for N2H+ with non-LTE hyperfines. N2Hplus above
is N2H+ with the LTE hyperfine approximation
	N2Hhyp		checked OK with accel 1-7-09
	HCNvib

NH3 and CH3CN are LTE for both the rotational and hyperfine transitions
	NH3		checked OK LTE 1-9-09
	CH3CN  		checked OK LTE 1-8-09

*/

	molecule = CH3CN;

/* Lambda acceleration works most of the time. 
   Tests with HCO+ and H13CO+ failed in that the 
   computed source function was much higher than 
   the gas temperature. Other molecules generally OK. 
   Not relevant for LTE molecules NH3 and CH3CN. */

	acceleration = NO;

	dust_continuum = NO;

/* Photodissociation and depletion are applicable to 
   clouds such as starless cores */

/* Normally we photodissociate all the molecules at the edge
   of the cloud */
	photodissociation = NO ;

/* Normally we deplete all the C molecules, but not the N molecules */
	depletion = NO;
	if (molecule == N2Hplus ||
	    molecule == N2Hhyp  ||
	    molecule == N2Dplus ||
	    molecule == NH3 ) depletion = NO;

/* Some standard abundances for some molecules */
/* These abundances are not automatically used. Depends
   on how define_model.c is written. Also set_model_parameters
   may overwrite these abundances depending on how it is 
   written */
	if (molecule == CO)	std_abundance = 5.625e-5;
	if (molecule == C13O)	std_abundance = 5.625e-5 / 77.;
	if (molecule == C18O)	std_abundance = 5.625e-5 / 77. / 7.3;
        if (molecule == C18O)   std_abundance = 5.625e-5 / 77. / 7.3 * 1.8 ;
        if (molecule == C17O)   std_abundance = 5.625e-5 / 77. / 7.3 / 4.0 * 1.8;
/*	if (molecule == CS)	std_abundance = 3.e-9; */
	if (molecule == CS)	std_abundance = 5.e-10;  /* use this for C34S */
	if (molecule == N2Hplus)std_abundance = 3.e-10;
	if (molecule == N2Hhyp )std_abundance = 2.e-10;
	if (molecule == N2Dplus)std_abundance = 3.e-10;
	if (molecule == NH3)    std_abundance = 0.5*2.4e-8;
        if (molecule == HCN)    std_abundance = 1.5e-8;
        if (molecule == HCNvib) std_abundance = 1.5e-8;
        if (molecule == CH3CN)  std_abundance = 6.e-8;

/* HCN abundance,
        Irvine, Goldsmith, & Hjalmarson 1987
        Paglione et al 1988
        1.5 x 10-8

   CH3CN abundance from Zhang's paper on AFGL5192
*/


/* Set the number of lines and number of levels. For linear rotors
   this is straight forward, and the number of lines should be 
   less than the number of levels. The lines and levels also have
   to be set in the 2 files nlines_C. and nlines_F77.h */

/* There are some special cases that require a specific combination of
   states and lines.

	NH3: requires		n_state=12	n_lines=6

	N2Hhyp can be run with one of the following pairs

	n_state = 10;
	n_lines = 1;

	n_state = 19;
	n_lines = 2;

	n_state = 28;
	n_lines = 3;

	n_state = 37;
	n_lines = 4;

	n_state = 64;
	n_lines = 7;

	CH3CN: requires 	nstate = 120	nlines = 105
*/

/* This is for Rainer's HCN
	n_state = 62;
	n_lines = 125;
*/

/* This is for NH3
	n_state = 12;
	n_lines = 6;
*/

/* This is for CH3CN  */

	n_state = 120;
	n_lines = 105;


/* A nice pair for rotors */
/*
	n_state = 8;
	n_lines = 7;
*/

/* H2Oortho can be run with 8 states and 11 lines 
   There are only collision rates for 8 states.

These are the combinations that can be used.
   8 states and 11 lines
   7 states and 9 lines
   6 states and 7 lines
   5 states and 5 lines
   4 states and 4 lines
   3 states and 2 lines
   2 states and 1 lines

The upper states of H2O have high excitation energies and
are probably not populated in cold clouds. 
If you get negative opacities in the upper levels, try
running with fewer states.
*/

/* Set the channel width in cm/s */     
/* Set the velocity range in cm/s */    
/* Choose Hanning smoothing YES or NO */          

/* 
 * velocity resolution and range 
*/
	chanwd = ${chwidth}; //cm/s
	//velrange = {velran}; //cm/s
	velrange = ${nchan}*chanwd; //cm/s

        hanning = YES ;

/* Set the minimum linewidth to decide which hyperfine lines to compute. */
/* Set the maximum linewidth for computing a range of line profiles. */

/* These lwmin,lwmax numbers are the gaussian widths (sigma) of the lines */

/* For CS 
        lwmin =  3000.;
        lwmax =  45000.;

*/

/* For CH3CN */
		// for precalculations of exponentials for gaussians
        lwmin =  ${minlwidth}; //cm/s
        lwmax =  ${maxlwidth}; //cm/s


/* ----------------------------------------------------*/

/* Number of nested boxes. Numbers greater than 1 have not
   been tested. */

	nbox = 3;
	//nbox = 2;

/* For each box set the dimensions. There is only one box
   and it is number zero. */

/* You cannot run with nx < 3. ny and nz can be as small as 1 
   If you want to use the beam convolution feature, then
   each dimension must be a power of 2 because the convolution
   is done by FFT. If beamx = beamy = 0 below, then no
   convolution is done and the grid does not have to be a power
   of 2.
*/

	ngrid = 150; //;

	ibox = 0;
	nx[ibox] = 80;  //38;//80;  //ngrid;
	ny[ibox] = 80;  //38;//80;  //ngrid;
	nz[ibox] = 80;  //38;//80;  //ngrid;
                                  
    ibox = 1;                     
    nx[ibox] = 100; //40;//100; //ngrid;
    ny[ibox] = 100; //40;//100; //ngrid;
    nz[ibox] = 100; //40;//100; //ngrid-120;//ngrid;
                                  
    ibox = 2;                     
    nx[ibox] = 100; //64;//100; //ngrid;
    ny[ibox] = 100; //64;//100; //ngrid;
    nz[ibox] = 50;  //14;//60;  //ngrid;
   


/* Radius of model sphere in pc. The calculations are done inside
  a sphere which can be smaller or larger than the input grid.

  There is only one cellsize for each box, 
  in other words, only square voxels.  Normally the 
  diameter of the model sphere (2*radius) would 
  be the same as the grid width. 

  This will put the sphere entirely inside the model grid.
  In this case, even on angles aligned parallel to the one
  axis of the model grid, not all the pathlengths through the 
  model will be the same. Also the pixels in the corners of the
  model grid will have no rays through them
	radius = 0.100 ;
        cellsize[ibox] = 2.*radius/nx[ibox];

  This case will make the sphere just large enough to contain
  the entire model grid. In this case, all the pathlengths
  on angles aligned parallel to the grid will be the same.
  All pixels, including those in the corners will have rays
  through them.
	radius = 0.200 ;
        cellsize[ibox] = radius/nx[ibox];

  In both these examples, the model grid is the same size,
  0.100 pc^3.

*/      

//	radius = 0.05 ;
    radius = ${radius} ;

	ibox = 0;
        cellsize[ibox] = (2.*radius)/nx[ibox];

    ibox = 1;
            //cellsize[ibox] = 2.*50000./206265./nx[ibox];//radius/nx[ibox];
            cellsize[ibox] = 2.*25000./206265./nx[ibox];//radius/nx[ibox];
   
	ibox = 2;
            //cellsize[ibox] = 2.*500./206265./nx[ibox];//(radius/2.)/nx[ibox];
	        cellsize[ibox] = 2.*1000./206265./nx[ibox];//(radius/2.)/nx[ibox];
	        //cellsize[ibox] = 2.*2500./206265./nx[ibox];//(radius/2.)/nx[ibox];
	
	
	ibox = 0;
/* Set the number of angles through the cloud. These are the
   angles along which radiative transfer will be done. 
   Use nlng = nlat = 1 for the LTE molecules, NH3 and CH3CN
*/
/*
   These are rotation angles and the rotation by
   latitude is applied first. The order makes a
   difference. For example, 90 lat followed 0 long
   is a view down the north pole along the z axis
   with x,y as the plane. */

/* This is the commonly used "6 ray approximation" */
	// for integral to find term for stimulated emission Jbar
       nlng = ${nphi};
       nlat = ${nth};

// For non-lte
//       nlng = 3;
//       nlat = 3;

/* ----------------------------------------------------*/
/* ----------------------------------------------------*/

/* Output parameters */

/* Set the number of lines to view in number_lines.
   number_lines must be 1 or more.
   Then select the lines. First line to view is
   linelist[0], next line (if number_lines > 1) is
   linelist[1]. For dipole molecules the (1-0)
   transition is number 0, the (2-1) transition is
   number 1, etc. */

/*        number_lines = 9; /* output lines to view */

/*
	for (i=0;i<number_lines-1;i++){
        linelist[i] = i;        
	}
*/

	// linelist[0] = 0;
	// linelist[1] = 2;

 /* For the CH3CN(12,0 -> 11,0) you need the following
    specific lines to be written. */


        number_lines = ${nlines_rt};      
        for (i=0;i<number_lines;i++){
                linelist[i] = i+${line_id}-1;
        }



/* If the model is cylindrically symmetric, alng can
   only be 0, but set it here as 0. These parameters
   follow the style of number_lines and linelist above. 
   But the numbers to be entered here are the degrees
   of the viewing angles. There is a section in cmain.c
   that picks the closest grid angle to each angle 
   specified here.*/


        outviews = 1;
          outlng[0] =  ${phi} ;
          outlat[0] =  90.-${th} ; // 0. is edge-on!
/*          alng[1] =  90. ;
          alat[1] =  0. ; */




/* Set the parameters for the output grid */


/* set output grid to write out every other grid point (after smoothing) */
/* Output axes are image like (x and y are axes on the plane of the sky), 
   but model axes are physics like (x and y are the equatorial plane.)
   These are output axes */

//        noutx = nouty = ngrid/2 + 1;
/* or easier and safer to write every point */
        noutx = nouty = ${npix};//4.*ngrid;

        //cellx = cellsize[0];
        //celly = cellx;

        cellx = celly = ${pixsize}; //2.*radius/noutx;
        //celly = 2.*radius/nouty;


/* beamx and beamy are the gaussian widths (sigma) used to weight 
   the output data. The relation between beam and FWHM is
   beam = FWHM/(2.sqrt(ln(2)) The beam size is in pc. */

/* ; N2H+ beam is from 37m Haystack, 27" FWHM
   ; At 150 pc,
   ; sigma = 150.*tan(27./3600*PI/180.)/2.35 pc
   ; sigma = 0.00835530
*/
//	beamx = 0.0; //
//    beamx = 0.25*4200./206265./2.35;
//	beamx = 0.8*1700./206265./2.35;
	//beamx = 0.00133
//        beamy = beamx ;


/* ----------------------------------------------------*/
/* ----------------------------------------------------*/

/* Comparison can be either  NO,  ANNEAL, or  SIMPLEX, 
 * Simplex function is not working. YES is equivalent to
 * ANNEAL */

	comparison = NO;

/* Set the number of models to compute. This should be 1 unless
   the code is set up to read data and compare model spectra
   against the data. */

        number_models = 1;

/* The number of adjustable parameters for each model. This is
   irrelevant unless the code is set up for comparisons. Has to be
   set up with the function set_model(). */

	nParameters = 6;
 
