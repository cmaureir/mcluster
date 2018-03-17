/***************************************************************************
 *   Copyright (C) 2009 by Andreas H.W. Kuepper                            *
 *   akuepper@astro.uni-bonn.de                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/***************************************************************************
 *   Compile using the command: cc -lm -o mcluster mcluster.c              *
 *   or use the Makefile, i.e. type: make mcluster or make mcluster_sse    *
 ***************************************************************************/
/***************************************************************************
 * Main array with stellar parameters:                                     *
 *      star[][0] = mass(t = epoch)                                        *
 *      star[][1-3] = {x, y, z}                                            *
 *      star[][4-6] = {vx, vy, vz}                                         *
 *      star[][7] = mass(t = 0)                                            *
 *      star[][8] = kstar, stellar type in (SSE/BSE)                       *
 *      star[][9] = epoch1, age within evolutionary phase (SSE/BSE)        *
 *      star[][10] = ospin, spin of star (SSE/BSE)                         *
 *      star[][11] = rstar, radius of star (SSE/BSE)                       *
 *      star[][12] = lstar, luminosity (SSE/BSE)                           *
 *      star[][13] = epochstar, age of star (SSE/BSE)                      *
 *      star[][14] = Zstar, metallicity of star                            *
 ***************************************************************************/

#include "include/functions.h"



int main (int argv, char **argc) {

    /*******************
     * Input variables *
     *******************/

    // Basic physical parameters

    // Number of stars, Mcl will be set to 0 if specified!
    int N = 0;
    // Total mass of the cluster, only used when N is set to 0, necessary for
    // usage of maximum stellar mass relation of Weidner & Kroupa (2006)
    double Mcl = 1000.0;
    // Density profile;
    // =0 Plummer model,
    // =1 King model (based on king0.f by D.C. Heggie),
    // =2 mass segregated Subr model (Subr et al. 2007),
    // =3 2-dimensional EFF template (Elson, Fall & Freeman 1987) or
    // Nuker template (Lauer et al. 1995),
    // =-1 no density gradient
    int profile = 0;
    // King's W0 paramter [0.3-12.0]
    double W0 = 5.0;
    // Fraction of mass segregation for profile; =0.0 unsegregated,
    // =1.0 completely segregated (take maximally S=0.99 for profile=2)
    double S = 0.0;
    // Fractal dimension; =3.0 unfractal, =2.6 2/8 fractal, =2.0 4/8 fractal,
    // =1.6 6/8 fractal, (2^D children per parent, following Goodwin &
    // Whitworth 2004)
    double D = 3.0;
    // Initial virial ratio; =0.5 virial equilibrium, <0.5 collapsing,
    // >0.5 expanding
    double Q = 0.5;
    // Half-mass radius [pc], ignored if profile = 3, set =-1 for using Marks
    // & Kroupa (2012) Mcl-Rh relation
    double Rh = 0.8;
    // Power-law slopes of EFF/Nuker templates (outer slope, inner slope,
    // transition); set gamma[1] = 0.0 and gamma[2] = 2.0 for EFF (profile = 3)
    double gamma[3] = {2.0, 0.0, 2.0};
    // Scale radius of EFF/Nuker template (profile = 3) [pc]
    double a = 1.0;
    // Cut-off radius for EFF/Nuker template (profile = 3) [pc]
    double Rmax = 100.0;
    // Simulation time [N-body units (Myr in Nbody6 custom)]
    double tcrit = 100.0;
    // Tidal field: =0 no tidal field,
    // =1 Near-field approximation,
    // =2 point-mass galaxy,
    // =3 Allen & Santillan (1991) MW potential
    // (or Sverre's version of it)
    int tf = 3;
    // Initial Galactic coordinates of the cluster [pc]
    double RG[3] = {8500.0,0.0,0.0};
    // Initial velocity of the cluster [km/s]
    double VG[3] = {0.0,220.0,0.0};

    //Mass function parameters

    // 0 = single mass stars;
    // 1 = use Kroupa (2001) mass function;
    // 2 = use multi power law (based on mufu.c by L.Subr)
    int mfunc = 1;
    // Stellar mass in case of single-mass cluster
    double single_mass = 1.0;
    // Lower mass limit for mfunc = 1 & mfunc = 4
    double mlow = 0.08;
    // Upper mass limit for mfunc = 1 & mfunc = 4
    double mup = 100.0;
    // alpha slopes for mfunc = 2
    double alpha[MAX_AN] = {-1.35, -2.35, -2.7, 0.0, 0.0};
    // mass limits for mfunc = 2
    double mlim[MAX_MN] = {0.08, 0.5, 4.0, 100.0, 0.0, 0.0};
    // alpha slope for mfunc = 4 (L3 mass function, Maschberger 2012)
    double alpha_L3 = 2.3;
    // beta slope for mfunc = 4
    double beta_L3 = 1.4;
    // mu parameter for mfunc = 4
    double mu_L3 = 0.2;
    // Usage of Weidner & Kroupa (2006) relation for most massive star;
    // =0 off, =1 on
    int weidner = 0;
    //Stellar evolution; 0 = off, 3 = Eggleton, Tout & Hurley [KZ19]
    int mloss = 3;
    // Use random kick velocity and present-day escape velocity to determine
    // retention of compact remnants & evolved binary components
    // (only for SSE/BSE version); =0 off, =1 on
    int remnant = 1;
    // Age of the cluster, i.e. star burst has been ... Myr before
    // [e.g. 1000.0, default = 0.0] [needs special compiling and SSE by
    // Hurley, Pols & Tout]
    double epoch = 0.0;
    // Metallicity [0.0001-0.03, 0.02 = solar]
    double Z = 0.02;
    // Metallicity [Fe/H], only used when Z is set to 0
    double FeH = -1.41;
    // Usage of Prantzos (2007) relation for the life-times of stars.
    // Set upper mass limit to Lifetime(mup) >= epoch
    int prantzos = 0;

    // Binary parameters

    // Number of primordial binary systems
    int nbin = 0;
    // Primordial binary fraction, number of binary systems = 0.5*N*fbin,
    // only used when nbin is set to 0
    double fbin = 0.0;
    // Pairing of binary components; 0= random pairing,
    // 1= ordered pairing for components with masses M>msort,
    // 2= random but separate pairing for components with masses m>Msort;
    // 3= Uniform distribution of mass ratio (0.1<q<1.0) for m>Msort and
    // random pairing for remaining (Kiminki & Kobulnicky 2012;
    // Sana et al., 2012; Kobulnicky et al. 2014; implemented by Long Wang)
    int pairing = 3;
    // Stars with masses > msort will be sorted and preferentially paired into
    // binaries if pairing = 1
    double msort = 5.0;
    // Semi-major axis distribution; 0= flat ranging from amin to amax,
    // 1= based on Kroupa (1995) period distribution,
    // 2= based on Duquennoy & Mayor (1991) period distribution,
    // 3= based on Kroupa (1995) period distribution for M<Msort;
    // based on Sana et al. (2012); Oh, S., Kroupa, P., & Pflamm-Altenburg,
    // J. (2015) period distribution for M>Msort (implemented by Long Wang)
    int adis = 3;
    // Use period distribution for massive binaries with M_primary > msort
    // from Sana & Evans (2011) if OBperiods = 1
    int OBperiods = 1;
    // Minimum semi-major axis for adis = 0 [pc]
    double amin = 0.0001;
    // Maximum semi-major axis for adis = 0 [pc]
    double amax = 0.01;
#ifdef SSE
    // Use Kroupa (1995) eigenevolution for pre-main sequence short-period
    // binaries;
    // =0 off,
    // =1 on [use either eigenevolution or BSE;
    // BSE recommended when using SSE]
    int eigen = 0;
    // Apply binary star evolution using BSE (Hurley, Tout & Pols 2002)
    // =0 off,
    // =1 on [use either eigenevolution or BSE;
    // BSE recommended when using SSE]
    int BSE = 1;
#else
    // Use Kroupa (1995) eigenevolution for pre-main sequence short-period
    // binaries;
    // =0 off,
    // =1 on [use either eigenevolution or BSE;
    // BSE recommended when using SSE]
    int eigen = 1;
    // Apply binary star evolution using BSE (Hurley, Tout & Pols 2002)
    // [needs special compiling and BSE];
    // =0 off,
    // =1 on [use either eigenevolution or BSE;
    // BSE recommended when using SSE]
    int BSE = 0;
#endif

    // Gas parameters (only used for Nbody6 input)

    // external Plummer (gas) sphere mass [Msun]
    double extmass = 0.0;
    // external Plummer (gas) sphere scale factor [pc]
    double extrad = 0.0;
    // decay time for gas expulsion (set 0 for no decay) [Myr]
    double extdecay = 0.0;
    // delay time for start of gas expulsion [Myr]
    double extstart = 0.0;

    // Code parameters

    // Nbody version:
    // =0 Nbody6,
    // =1 Nbody4,
    // =2 Nbody6 custom,
    // =3 only create output list of stars,
    // =4 Nbody7 (not yet fully functional),
    // =5 Nbody6++GPU
    int code = 3;
    // Number seed for random number generator;
    // =0 for randomization by local time
    unsigned int seed = 0;
    // Name of output files
    char *output = "test";
    // DTADJ [N-body units (Myr in Nbody6 custom)], energy-check time step
    double dtadj = 1.0;
    // DELTAT [N-body units (Myr in Nbody6 custom)],
    // output interval, must be multiple of DTADJ
    double dtout = 1.0;
    // DTPLOT [N-body units (Myr in Nbody6 custom)],
    // output of HRdiagnostics, should be multiple of DTOUT,
    // set to zero if output not desired
    double dtplot = 1.0;
    // Use of GPU,
    // 0= off,
    // 1= on
    int gpu = 0;
    // Update of regularization parameters during computation;
    // 0 = off,
    // 0 > on
    int regupdate = 0;
    // Update of ETAI & ETAR during computation;
    // 0 = off,
    // 0 > on
    int etaupdate = 0;
    // Removal of escapers;
    // 0 = no removal,
    // 1 = regular removal at 2*R_tide;
    // 2 = removal and output in ESC
    int esc = 2;
    // Units of McLuster output;
    // 0= Nbody-Units,
    // 1= astrophysical units
    int units = 1;

    // McLuster internal parameters

    // Make cluster half-mass radius exactly match the expected/desired value;
    // =0 off,
    // =1 on (recommended)
    int match = 1;
    // Force spherical symmetry for fractal clusters;
    // =0 off,
    // =1 on (recommended)
    int symmetry = 1;
    // Make energy check at end of McLuster;
    // =0 off,
    // =1 on
    int check = 0;
    // Creates a radial density profile and prints it to the screen;
    // =0 off,
    // =1 on
    int create_radial_profile = 1;
    // Creates a radial cumulative profile and prints it to the screen;
    // =0 off,
    // =1 on
    int create_cumulative_profile = 1;
    // Distance of cluster from sun for artificial CMD with
    // observational errors [pc]
    double Rgal = 10000.0;
    // Solar metallicity
    double Zsun = 0.02;
    // Maximum number of stars & orbits allowed in McLuster
    int NMAX = 1500000;
    // Maximum number of neighbours allowed in NBODY6
    int NNBMAX_NBODY6 = 500;
    // Maximum stellar mass allowed in McLuster [Msun]
    double upper_IMF_limit = 150.0;
    // Counter for number of alpha slopes for mfunc = 2
    int an = 0;
    // Counter for number of mass limits for mfunc = 1, 2 & 4
    int mn = 0;
    // Counter for components of galactocentric radius vector
    int xn = 0;
    // Counter for components of cluster velocity vector
    int vn = 0;
    // Counter for external potential input parameters
    int xx = 0;
    // Counter for EFF/Nuker profile parameters
    int gn = 0;
    // Input array for external potential parameters
    double extgas[4];


    // SSE internal parameters (see Hurley, Pols & Tout 2000)

    // Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally)
    value1_.neta = 0.5;
    // Binary enhanced mass loss parameter (inactive for single)
    value1_.bwind = 0.0;
    // Helium star mass loss factor (1.0 normally)
    value1_.hewind = 1.0;
    // Maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1)
    value1_.mxns = 3.0;
    // Time-step parameter in evolution phase: MS (0.05)
    points_.pts1 = 0.05;
    // Time-step parameter in evolution phase: GB, CHeB, AGB, HeGB (0.01)
    points_.pts2 = 0.01;
    // Time-step parameter in evolution phase: HG, HeMS (0.02)
    points_.pts3 = 0.02;
    // Kick velocities
    value4_.sigma = 190.0;
    // bhflag > 0 allows velocity kick at BH formation
    value4_.bhflag = 1;

    // BSE internal parameters (see Hurley, Pols & Tout 2002)

    // ceflag > 0 activates spin-energy correction in common-envelope (0)
    // #defunct#
    flags_.ceflag = 0;
    // tflag > 0 activates tidal circularisation (1)
    flags_.tflag = 1;
    // ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0)
    flags_.ifflag = 0;
    // nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ,572,407(1)
    flags_.nsflag = 1;
    // wdflag > 0 uses modified-Mestel cooling for WDs (0)
    flags_.wdflag = 1;

    // beta is wind velocity factor: proportional to vwind**2 (1/8)
    value5_.beta = 1.0/8.0;
    // xi is the wind accretion efficiency factor (1.0)
    value5_.xi = 1.0;
    // acc2 is the Bondi-Hoyle wind accretion factor (3/2)
    value5_.acc2 = 3.0/2.0;
    // epsnov is the fraction of accreted matter retained in nova
    // eruption (0.001)
    value5_.epsnov = 0.001;
    // eddfac is Eddington limit factor for mass transfer (1.0)
    value5_.eddfac = 1.0;
    // gamma is the angular momentum factor for mass lost during Roche (-1.0)
    value5_.gamma = -1.0;
    // alpha1 is the common-envelope efficiency parameter (1.0)
    value2_.alpha1 = 1.0;
    // lambda is the binding energy factor for common envelope evolution (0.5)
    value2_.lambda = 0.5;


    // Command line input
    int option;
    while ((option = getopt(argv, argc, "N:M:P:W:R:r:c:g:S:D:T:Q:C:A:O:G:o:f:a:m:B:b:p:s:t:e:Z:X:x:V:u:h:?")) != -1) switch (option)
    {
        case 'N': N = atoi(optarg); Mcl = 0.0; break;
        case 'M': Mcl = atof(optarg); N = 0; break;
        case 'P': profile = atoi(optarg); break;
        case 'W': W0 = atof(optarg); break;
        case 'R': Rh = atof(optarg); break;
        case 'r': a = atof(optarg); break;
        case 'c': Rmax = atof(optarg); break;
        case 'g':
            if (gn < 3) { gamma[gn] = atof(optarg); gn++; break; }
        case 'S': S = atof(optarg); break;
        case 'D': D = atof(optarg); break;
        case 'T': tcrit = atof(optarg); break;
        case 'Q': Q = atof(optarg); break;
        case 'C': code = atoi(optarg); break;
        case 'A': dtadj = atof(optarg); break;
        case 'O': dtout = atof(optarg); break;
        case 'G': gpu = atoi(optarg); break;
        case 'o': output = optarg; break;
        case 'f': mfunc = atoi(optarg); break;
        case 'a' :
            if (an < MAX_AN) {
                alpha[an] = atof(optarg);
                if (an == 0) alpha_L3 = atof(optarg);
                if (an == 1) beta_L3 = atof(optarg);
                if (an == 2) mu_L3 = atof(optarg);
                an++;
                break;
            } else { printf("\nError: Number of alphas exceeded maximum limit of %d\n", MAX_AN); return 1; }
        case 'm' :
            if (mn < MAX_MN) {
                mlim[mn] = atof(optarg);
                if (mn == 0) mlow = atof(optarg);
                if (mn == MAX_MN-1) mup = atof(optarg);
                mn++;
                break;
            } else { printf("\nError: Number of mass params exceded maximum limit of %d\n", MAX_MN); return 1; }
        case 'B': nbin = atoi(optarg); break;
        case 'b': fbin = atof(optarg); break;
        case 'p': pairing = atoi(optarg); break;
        case 's': seed = atoi(optarg); break;
        case 't': tf = atoi(optarg); break;
        case 'e': epoch = atof(optarg); break;
        case 'Z': Z = atof(optarg); break;
        case 'X' :
            if (xn < 3) { RG[xn] = atof(optarg); xn++; break; }
        case 'V' :
            if (vn < 3) { VG[vn] = atof(optarg); vn++; break; }
        case 'x' :
            if (xx < 4) { extgas[xx] = atof(optarg); xx++; break; }
        case 'u': units = atoi(optarg); break;
        case ':':
        case 'h':    help(msort); return 1;
        case '?':    help(msort); return 1;
    };

    if (mn-1 > 0) mup = mlim[mn-1];

    // print summary of input parameters to .info file
    info(output, N, Mcl, profile, W0, S, D, Q, Rh, gamma, a, Rmax, tcrit, tf,
        RG, VG, mfunc, single_mass, mlow, mup, alpha, mlim, alpha_L3, beta_L3,
        mu_L3, weidner, mloss, remnant, epoch, Z, prantzos, nbin, fbin,
        pairing, msort, adis, amin, amax, eigen, BSE, extmass, extrad,
        extdecay, extstart, code, seed, dtadj, dtout, dtplot, gpu, regupdate,
        etaupdate, esc, units, match, symmetry, OBperiods);

    /*********
     * Start *
     *********/

    printf("\n-----START----         \n");

#ifdef NOOMP
    clock_t t1, t2;
    // start stop-watch
    t1 = clock();
#else
    double t1, t2;
#pragma omp parallel
    {
        if (omp_get_thread_num() == 0)
            if (omp_get_num_threads() > 1)
                printf("\n\nUsing OpenMP with %d threads\n",
                    omp_get_num_threads());

        t1 = omp_get_wtime(); //start stop-watch
    }
#endif

    if (seed) {
        // initialize random number generator by seed
        srand48(seed);
    }
    else {
        // initialize random number generator by local time
        seed = (unsigned) time(NULL);
        seed %= 100000;
        srand48(seed);
    }

    printf ("\n\nRandom seed = %i\n\n\n", seed);
    if (seed) {
        // idum is the random number seed used in the kick routine.
        value3_.idum = seed;
    }
    else {
        value3_.idum = 10000000.0*drand48();
    }

    int i,j;
    // Total mass [M_sun]
    double M;
    // Mean stellar mass [M_sun]
    double mmean;
    // Maximum neighbour number (Nbody6 only)
    int NNBMAX;
    // Initial radius of neighbour sphere [pc], Nbody6 only
    double RS0;
    // Tidal radius [pc]
    double rtide;
    // Angular velocity of cluster around the galaxy
    double omega;
    // Virial radius [pc]
    double rvir=0;
    // For CoM correction
    double cmr[7];
    // King-, Plummer radius
    double rking=0, rplummer=0;
    // most massive star
    double MMAX;
    // time-scale factor
    double tscale;
    // kinetic energy
    double ekin = 0.0;
    // potential energy
    double epot = 0.0;
    // velocity dispersion
    double sigma = 0.0;
    // KZ(22) parameter (binaries yes/no)
    int bin;
    // (evolved stellar population yes/no)
    int sse;
    // mass function parameters for mfunc = 2
    double submass[MAX_AN], subcount[MAX_AN], norm[MAX_AN], N_tmp, M_tmp;
    // actual 2D/3D half-mass radius of the model
    double Rh2D, Rh3D;

    // set half-mass radius temporarily to scale radius for computation of
    // escape velocity
    if (profile == 3)
        Rh = a;

    if ((Mcl) && (N)) {
        printf("\nWARNING:"
            "specify either Mcl (-M) or N (-N)!\n\n");
        exit (1);
    } else if ((!Mcl) && (mfunc == 3)) {
        printf("\nWARNING:"
            "specify Mcl (-M) when using optimal sampling (-f 3)!\n\n");
        exit (1);
    }

    if (xx > 0) {
        if ((xx == 4) && (extgas[2])) {
            extmass = extgas[0];
            extrad = extgas[1];
            extdecay = extgas[2];
            extstart = extgas[3];
        } else {
            printf("\nWARNING: Insufficient or incorrect parameters specified"
                " for external Plummer potential!\n\n");
            exit (1);
        }
    }


    /***********************
     * Generate star array *
     ***********************/

    int columns = 15;
    double **star;
    star = (double **)calloc(NMAX,sizeof(double *));
    for (j=0;j<NMAX;j++){
        star[j] = (double *)calloc(columns,sizeof(double));
        if (star[j] == NULL) {
            printf("\nMemory allocation failed!\n");
            return 0;
        }
    }



    /*******************************************
     * Evaluate Z from [Fe/H] if Z is set to 0 *
     *******************************************/

     if (!Z) {
        // Bertelli, Bressan, Chiosi, Fagotto, Nasi, 1994, A&AS, 106, 275
        Z = pow(10.0, 0.977*FeH)*Zsun;
        printf("\nUsing Bertelli et al. (1994) relation to convert FeH = %.3f"
            " into Z = %.3f\n", FeH, Z);
    }



    /**********************************
     * Calculate maximum stellar mass *
     **********************************/

    //always use Weidner relation when using optimal sampling!
     if (mfunc == 3)
        weidner = 1;

    if (!N && weidner && mfunc) {
        mup = upper_IMF_limit;

        printf("\nUsing maximum stellar mass-cluster mass relation for upper"
            " stellar mass limit\n");
        printf("\n(Weidner & Kroupa 2006, Pflamm-Altenburg & Kroupa 2007)\n");

        // analytic fit to the observational data from
        // Pflamm-Altenburg & Kroupa (2007), implementation by M. Kruckow
        MMAX = pow(10,2.56*log10(Mcl) * pow(pow(3.82,9.17) + pow(log10(Mcl),
            9.17),-1/9.17)-0.38);

        /*
        //three-part power-law fit to the observational data
        if (Mcl < 1000.0) {
            MMAX = (log10(Mcl)*0.540563175) - 0.14120167;
            MMAX = pow(10.0,MMAX);
        } else if (Mcl < 3300.0) {
            MMAX = (log10(Mcl)*0.19186051) + 0.9058611;
            MMAX = pow(10.0,MMAX);
        } else {
            MMAX = (log10(Mcl)*0.360268003) + 0.313342031;
            MMAX = pow(10.0,MMAX);
        }*/
    } else {
        MMAX = mup;
    }

    if (mfunc && epoch && prantzos) {
        printf("\nUsing Prantzos (2007) relation to reduce upper mass"
            " limit to Lifetime(mup) > epoch\n");
        while (Lifetime(MMAX) < sqrt(pow(epoch,2))) {
            MMAX -= 0.01;
        }
    }


    /*******************
     * Generate masses *
     *******************/

    printf("\n\n-----GENERATE MASSES-----     \n");

    // This loop has to be called multiple times with different N (or Mcl),
    // Z and epoch in order to create multiple stellar populations  make sure
    // that epoch and Z are stored together with the stars, then concatenate
    // all arrays for now, all Z and epochs are the same,
    // i.e. a single stellar population
    if (mfunc == 1) {
        printf("\nMaximum stellar mass set to: %.2f\n", MMAX);
        generate_m1(&N, star, mlow, mup, &M, &mmean, MMAX, Mcl, epoch, Z, Rh,
            remnant);
    } else if (mfunc == 2) {

        if (mn) {
            for (i = mn+1; i < MAX_MN; i++) mlim[i] = 0.0;
        } else {
            for (i=0; i<MAX_MN; i++) {
                if (mlim[i]) mn++;
            }
        }

        if (an) {
            for (i = an+1; i < MAX_AN; i++) alpha[i] = 0.0;
        } else {
            for (i=0; i<MAX_AN; i++) {
                if (alpha[i]) an++;
            }
        }

        if (an >= mn)
            an = mn - 1;

        mn = an + 1;

        if (!mn){
            printf("\nError: at least one mass limit has to be specified\n");
            return 1;
        } else if (mn == 1) {
            single_mass = mlim[0];
            printf("\nSetting stellar masses to %g solar mass\n",single_mass);

            if (!N)
                N = Mcl/single_mass;

            for (j=0;j<N;j++)
                star[j][0] = 1.0/N;

            mmean = single_mass;
            M = N*mmean;
            printf("\nM = %g\n", M);
        } else {
            printf("\nMaximum stellar mass set to: %.2f\n", MMAX);
            // normalization factor of integral
            norm[an-1] = 1.;
            // integrated number of stars in interval [mlim[an-1]:mlim[an]]
            N_tmp = subcount[an-1] = subint(mlim[an-1], mlim[an],
                alpha[an-1] + 1.);
            // integrated mass of stars in interval [mlim[an-1]:mlim[an]]
            M_tmp = submass[an-1] = subint(mlim[an-1], mlim[an],
                alpha[an-1] + 2.);

            for (i = an - 2; i >= 0; i--) {
                norm[i] = norm[i+1] * pow(mlim[i+1], alpha[i+1] - alpha[i]);
                subcount[i] = norm[i] * subint(mlim[i], mlim[i+1],
                    alpha[i] + 1.);
                N_tmp += subcount[i];
                submass[i] = norm[i] * subint(mlim[i], mlim[i+1],
                    alpha[i] + 2.);
                M_tmp += submass[i];
            }

            generate_m2(an, mlim, alpha, Mcl, M_tmp, subcount, &N, &mmean, &M,
                star, MMAX, epoch, Z, Rh, remnant);
        }
    } else if (mfunc == 3) {
        printf("\nMaximum stellar mass set to: %.2f\n",MMAX);
        generate_m3(&N, star, mlow, mup, &M, &mmean, MMAX, Mcl);
        randomize(star, N);
    } else if (mfunc == 4) {
        printf("\nMaximum stellar mass set to: %.2f\n", MMAX);
        printf("\nUsing L3 IMF (Maschberger 2012)\n");
        generate_m4(&N, star, alpha_L3, beta_L3, mu_L3, mlow, mup, &M, &mmean,
            MMAX, Mcl, epoch, Z, Rh, remnant);
    } else {
        printf("\nSetting stellar masses to %.1f solar mass\n", single_mass);
        if (!N) N = Mcl/single_mass;
        for (j=0;j<N;j++) {
            star[j][0] = single_mass;
            star[j][7] = single_mass;
            star[j][8] = 0;
            star[j][9] = 0.0;
            star[j][10] = 0.0;
            star[j][11] = 0.0;
            star[j][12] = 0.0;
            star[j][13] = 0.0;
            star[j][14] = 0.0;
        }
        mmean = single_mass;
        M = N*mmean;
        printf("\nM = %g\n", M);
        mloss = 0;
    }

    // set all stars to the same metallicity and age for now
    double epochstar, zstar;
    // age compared to the oldest stars in the cluster [Myr]
    epochstar = 0.0;
    zstar = Z;

    for (j=0;j<N;j++) {
        star[j][13] = epochstar;
        star[j][14] = zstar;
    }

    // Pair binary masses and convert to centre-of-mass particles
    int Nstars;
    if (!nbin)
        nbin = 0.5*N*fbin;

    // component mass & stellar evol parameter array
    double **mbin;
    mbin = (double **)calloc(nbin,sizeof(double *));
    for (j=0;j<nbin;j++) {
        mbin[j] = (double *)calloc(20,sizeof(double));

        if (mbin[j] == NULL) {
            printf("\nMemory allocation failed!\n");
            return 0;
        }
    }

    Nstars = N;

    if (nbin) {
        printf("\nPreparing binary components.\n");
        // Specify component pairing
        if (pairing) {

            if (pairing == 1) {
                printf("\nApplying ordered pairing for stars with masses"
                    " > %.1f Msun.\n", msort);
            }
            else if (pairing == 2) {
                printf("\nApplying random pairing for stars with masses"
                    " > %.1f Msun.\n", msort);
            }
            else if (pairing == 3) {
                printf("\nApplying uniform mass ratio distribution for stars"
                    " with masses > %.1f Msun.\n", msort);
            }

            order(star, N, M, msort, pairing);
        } else {
            randomize(star, N);
            printf("\nApplying random pairing.\n");
        }
        N -= nbin;
        for (j=0;j<nbin;j++) {
            // system mass
            mbin[j][0] = star[2*j][0]+star[2*j+1][0];
            // primary mass
            mbin[j][1] = star[2*j][0];
            // secondary mass
            mbin[j][2] = star[2*j+1][0];
            // primary m0
            mbin[j][3] = star[2*j][7];
            // primary kw
            mbin[j][4] = star[2*j][8];
            // primary epoch1
            mbin[j][5] = star[2*j][9];
            // primary spin
            mbin[j][6] = star[2*j][10];
            // primary r
            mbin[j][7] = star[2*j][11];
            // primary lum
            mbin[j][8] = star[2*j][12];
            // secondary m0
            mbin[j][9] = star[2*j+1][7];
            // secondary kw
            mbin[j][10] = star[2*j+1][8];
            // secondary epoch1
            mbin[j][11] = star[2*j+1][9];
            // secondary spin
            mbin[j][12] = star[2*j+1][10];
            //secondary r
            mbin[j][13] = star[2*j+1][11];
            // secondary lum
            mbin[j][14] = star[2*j+1][12];
            // identifier
            mbin[j][15] = 1000+j;
            // primary epochstar
            mbin[j][16] = star[2*j][13];
            // secondary epochstar
            mbin[j][17] = star[2*j+1][13];
            // primary zstar
            mbin[j][18] = star[2*j][14];
            // secondary zstar
            mbin[j][19] = star[2*j+1][14];

            //s ystem mass
            star[2*j][0] += star[2*j+1][0];
            star[2*j+1][0] = 0.0;
            // identifier
            star[2*j][7] = 1000+j;
            // identifier
            star[2*j+1][7] = 0.0;
            // system luminosity
            star[2*j][12] += star[2*j+1][12];
            star[2*j+1][12] = 0.0;
        }

        order(star, Nstars, M, 0.0, 0);
        randomize(star, N);
    }


    // prepare mass segregation

    // search lowest mass star
    double mlowest = MMAX;
    //search highest mass star
    double mhighest = 0;
    double mmeancom = 0.0;

    for (i=0;i<N;i++) {
        // scale masses to Nbody units
        star[i][0] /= M;
        mmeancom += star[i][0];

        if (star[i][0] < mlowest)
            mlowest = star[i][0];

        if (star[i][0] > mhighest)
            mhighest = star[i][0];
    }
    mmeancom /= N;
    // number of necessary pos & vel pairs for Baumgardt et al. (2008)
    // mass segregation routine
    int Nseg = ceil(N*mmeancom/mlowest);
    int Nunseg = N;

    double *Mcum;
    Mcum = (double *)calloc(N,sizeof(double));

    // sort masses when mass segregation parameter > 0
    if ((S) && !(profile == 2)) {
        printf("\nApplying mass segregation with S = %f\n",S);
        segregate(star, N, S);

        //calculate cumulative mass function Mcum
        for (i=0;i<N;i++) {
            Mcum[i] = 0.0;
            for (j=0;j<=i;j++) Mcum[i] = Mcum[i] + star[j][0];
        }
        N = Nseg;
    }


    /*************************************
     * Generate positions and velocities *
     *************************************/

    printf("\n\n-----GENERATE POSITIONS & VELOCITIES-----   \n");

    // calculate half-mass radius according to Marks & Kroupa 2012
    // if Rh is set to -1
    if ((Rh == -1) && (Mcl)) {
        // Marks & Kroupa (2012), implementation by M. Kruckow
        Rh = 0.1*pow(Mcl,0.13);
        printf("\nUsing Marks & Kroupa (2012) relation to derive half-mass"
            " radius from cluster mass: %g (pc)\n", Rh);
    }

    // evaluate approximate tidal radius assuming circular orbit
    if (tf == 3) {
        // in the case of Allen & Santillan potential,
        // assume kappa = 1.4omega (eq. 9 in Kuepper et al. 2010)
        omega = sqrt(VG[0]*VG[0]+VG[1]*VG[1]+VG[2]*VG[2])/sqrt(RG[0]*RG[0]+
            RG[1]*RG[1]+RG[2]*RG[2]);
        rtide = pow(G*M/(2.0*omega*omega),1.0/3.0);
    } else if (!tf) {
        rtide = 1.0E5;
    } else if ((tf == 1) && (code == 0 || code == 4 || code == 5)) {
        // in case of Sverre's Nbody6 standard tidal field
        rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*8000.0;
    } else {
        // in the case of a point mass potential or near field approximation
        rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*sqrt(RG[0]*RG[0]+RG[1]*
            RG[1]+RG[2]*RG[2]);
    }
    printf("\nApproximate tidal radius: %g (pc)\n", rtide);


    // generate scaled pos & vel, postpone scaling for Plummer and King in
    // case of mass segregation
    double rhtemp, rvirtemp;

    if (profile == 1) {
        printf("\nGenerating King model with parameters: "
            "N = %i\t W0 = %g\t Rh = %.3f\t D = %.2f\n",N, W0, Rh, D);
        generate_king(N, W0, star, &rvirtemp, &rhtemp, &rking, D, symmetry);
    } else if (profile == 2) {
        N = Nunseg;
        printf("\nGenerating segregated Subr model with parameters: "
            "N = %i\t S = %g\t Rh = %.3f\n",N, S, Rh);
        // value provided by L. Subr
        rvir = Rh/0.76857063065978;
        generate_subr(N, S, star, rtide, rvir);
        printf ("\nrvir = %.5f\t rh = %.5f\t rtide = %.5f (pc)\n",
            rvir, Rh, rtide);
    } else if (profile == 3) {
        if (gamma[1] == 0.0 && gamma[2] == 2.0) {
            printf("\nGenerating EFF model with parameters: "
                "N = %i\t a = %.1f\t gamma = %.3f\t Rmax = %.1f\t D = %.2f\n",
                N, a, gamma[0], Rmax, D);
        }
        else {
            printf("\nGenerating Nuker model with parameters: "
                "N = %i\t a = %.1f\t gamma (outer) = %.3f\t gamma (inner) "
                "= %.3f\t transition = %.3f\t Rmax = %.1f\t D = %.2f\n",
                N, a, gamma[0],gamma[1],gamma[2], Rmax, D);
        }

        double p[6];

        // rho0, will be scaled according to Rmax and Mtot
        p[1] = 1.0;
        //scale radius
        p[2] = a;
        //outer power-law slope (>0.5)
        p[3] = gamma[0];
        //inner power-law slope (0.0 for EFF template)
        p[4] = gamma[1];
        //transition parameter (2.0 for EFF template)
        p[5] = gamma[2];

        generate_profile(N, star, Rmax, M, p, &Rh, D, symmetry);
        printf("\nRh = %.1f pc\n", Rh);
    } else if (profile == -1) {
        printf("\nGenerating fractal distribution with parameters: "
            "N = %i\t Rh = %.3f\t D = %.2f\n", N, Rh, D);
        fractalize(D, N, star, 0, symmetry);
        rvir = Rh;
    } else {
        printf("\nGenerating Plummer model with parameters: "
            "N = %i\t Rh = %.3f\t D = %.2f\n", N, Rh, D);
        rvir = Rh/0.772764;
        rplummer = Rh/1.305;
        generate_plummer(N, star, rtide, rvir, D, symmetry);
    }

    // Apply Baumgardt et al. (2008) mass segregation
    if (!(profile == 2) && (S)) {
        double **m_temp;
        m_temp = (double **)calloc(Nunseg,sizeof(double *));

        for (j=0;j<Nunseg;j++) {
            m_temp[j] = (double *)calloc(9,sizeof(double));

            if (m_temp[j] == NULL) {
                printf("\nMemory allocation failed!\n");
                return 0;
            }
        }

        // temporarily store the masses & stellar evol parameters
        for (i=0;i<Nunseg;i++) {
            m_temp[i][0] = star[i][0];
            m_temp[i][1] = star[i][7];
            m_temp[i][2] = star[i][8];
            m_temp[i][3] = star[i][9];
            m_temp[i][4] = star[i][10];
            m_temp[i][5] = star[i][11];
            m_temp[i][6] = star[i][12];
            m_temp[i][7] = star[i][13];
            m_temp[i][8] = star[i][14];
        }

        printf("\nOrdering orbits by energy.\n");
        energy_order(star, N, Nstars);

        int nlow, nhigh, nrandom;
        for (i=0;i<Nunseg;i++) {
            nhigh = Nseg*Mcum[i];
            if (i) {
                nlow = Nseg*Mcum[i-1];
            } else {
                nlow = 0;
            }
            nrandom = (nhigh-nlow)*drand48()+nlow;
            star[i][0] = m_temp[i][0];
            star[i][1] = star[nrandom][1];
            star[i][2] = star[nrandom][2];
            star[i][3] = star[nrandom][3];
            star[i][4] = star[nrandom][4];
            star[i][5] = star[nrandom][5];
            star[i][6] = star[nrandom][6];
            star[i][7] = m_temp[i][1];
            star[i][8] = m_temp[i][2];
            star[i][9] = m_temp[i][3];
            star[i][10] = m_temp[i][4];
            star[i][11] = m_temp[i][5];
            star[i][12] = m_temp[i][6];
            star[i][13] = m_temp[i][7];
            star[i][14] = m_temp[i][8];
        }


        for (j=0;j<Nunseg;j++)
            free (m_temp[j]);

        free(m_temp);
        N = Nunseg;
    }

    // CoM correction
    printf("\nApplying centre-of-mass correction.\n");
    for (j=0; j<7; j++)
        cmr[j] = 0.0;

    for (j=0; j<N; j++) {
        for (i=1;i<7;i++)
            cmr[i] += star[j][0]*star[j][i];
    }

    for (j=0; j<N; j++) {
        for (i=1;i<7;i++)
            star[j][i] -= cmr[i];
    }


    // Apply scaling to Nbody-units
    if (profile == 0) {
        printf("\nRe-scaling of orbits (dt ~ N^2!)\n");
        double ke = 0.0;
        double pe = 0.0;
        double sx, sv, r2;

        #ifdef GPU
        gpupot(N,star,&pe);
        #ifndef NOOMP
        #pragma omp parallel shared(N, star)  private(i)
        {
            #pragma omp for reduction(+: ke) schedule(dynamic)
        #endif
            for (i=0;i<N;i++) {
                ke += star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+
                    pow(star[i][6],2));
            }
        #ifndef NOOMP
        }
        #endif
        //        printf("PE_GPU %lg ke %lg\n",pe,ke);
        #else
        #ifndef NOOMP
        #pragma omp parallel shared(N, star)  private(i, j, r2)
        {
            #pragma omp for reduction(+: pe, ke) schedule(dynamic)
        #endif

        for (i=0;i<N;i++) {
            if (i) {
                for (j=0;j<i-1;j++) {
                    r2=(star[i][1] - star[j][1]) * (star[i][1] - star[j][1]) +
                       (star[i][2] - star[j][2]) * (star[i][2] - star[j][2]) +
                       (star[i][3] - star[j][3]) * (star[i][3] - star[j][3]) ;
                    pe -= star[i][0]*star[j][0]/sqrt(r2);
                }
            }
            ke += star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+
                pow(star[i][6],2));
        }
        #ifndef NOOMP
        }
        #endif
        //        printf("PE_HOST %lg ke %lg\n",pe,ke);
        #endif

        rvir = -GNBODY*pow(MNBODY,2)/(2.0*pe);
        sx = 1.0/rvir;
        ke *= 0.5;
        sv = sqrt(4.0*ke);

        for (i=0;i<N;i++) {
            star[i][1] *= sx;
            star[i][2] *= sx;
            star[i][3] *= sx;
            star[i][4] /= sv;
            star[i][5] /= sv;
            star[i][6] /= sv;
        }

        ke = 0;
        for (i=0;i<N;i++) {
            ke += M*star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+
                pow(star[i][6],2));
        }
        ke /= 0.5*N;
        printf("Dynamical temperature of centre-of-mass particles kT = "
            "%lf\n\n", ke);

        // make half-mass radius of the system match the desired one
        radial_profile(star, N, rvir, M, 0, 0, code, &NNBMAX, &RS0, &Rh2D,
            &Rh3D, NNBMAX_NBODY6);

        if (match) {
            printf("\nmeasuring half-mass radius: %.7f \t %.7f (should/is)\n"
                "and correcting for this factor\n", Rh, Rh3D);
            rvir = rvir *Rh/Rh3D;
        }

        printf ("\nrvir = %.5f\t rh = %.5f\t rplummer = %.5f\t rtide = "
            "%.5f (pc)\n", rvir, Rh, rplummer, rtide);
    } else if ((profile == 1) || (profile == 3) || (profile == -1)) {
        printf("\nRe-scaling of orbits (dt ~ N^2!)\n");

        double pe = 0.0;
        double ke = 0.0;
        double r2, vscale;

        #ifdef GPU
        gpupot(N,star,&pe);
        #ifndef NOOMP
        #pragma omp parallel shared(N, star)  private(i)
        {
        #pragma omp for reduction(+: ke) schedule(dynamic)
        #endif
          for (i=0;i<N;i++) {
            ke += star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+
                pow(star[i][6],2));
          }
        #ifndef NOOMP
        }
        #endif
        //        printf("PE_GPU %lg ke %lg\n",pe,ke);
        #else
        #ifndef NOOMP
        #pragma omp parallel shared(N, star)  private(i, j, r2)
        {
        #pragma omp for reduction(+: pe, ke) schedule(dynamic)
        #endif
        for (i=0;i<N;i++) {
            for (j=0;j<i-1;j++) {
                r2 = (star[i][1] - star[j][1]) * (star[i][1] - star[j][1]) +
                     (star[i][2] - star[j][2]) * (star[i][2] - star[j][2]) +
                     (star[i][3] - star[j][3]) * (star[i][3] - star[j][3]) ;

                pe -=  star[i][0]*star[j][0]/sqrt(r2);
            }
            ke += star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+pow(star[i][6],2));
        }
        #ifndef NOOMP
        }
        #endif
        //        printf("PE_HOST %lg ke %lg\n",pe,ke);
        #endif
        ke *= 0.5;
        rvir = -GNBODY*pow(MNBODY,2)/(2.0*pe);
        vscale = sqrt(4.0*ke);

        ke = 0;
        for (i=0;i<N;i++) {

            for (j=0;j<3;j++) {
                star[i][j+1] /= rvir;
                star[i][j+4] /= vscale;
            }
            ke += M*star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+
                pow(star[i][6],2));
        }
        ke /= 0.5*N;
        printf("Dynamical temperature of centre-of-mass particles kT ="
            "%lf\n\n",ke);

        // make half-mass radius of the system match the desired one
        radial_profile(star, N, rvir, M, 0, 0, code, &NNBMAX, &RS0, &Rh2D,
            &Rh3D, NNBMAX_NBODY6);

        if (match) {
            printf("\nmeasuring half-mass radius: %.7f \t %.7f (should/is)\n"
                "and correcting for this factor\n", Rh, Rh3D);
            rvir = rvir *Rh/Rh3D;
        }
        //printf("\nrvir = %.5f\t rh = %.5f\t rtide = %.5f (pc)\n", rvir, Rh,
        //  rtide);
        //printf("Edge radius (King units) = %g\t(Nbody units) = %g\n", rking,
        //  rking/rvirtemp);
        //printf("Core radius (King units) = %g\t(Nbody units) = %g\n\n", 1.0,
        //  1.0/rvirtemp);
        //printf("Concentration = %g\n", log10(rking));
    } else if (profile == 2) {
        double ke = 0;
        #ifndef NOOMP
        #pragma omp parallel shared(N, star)  private(i)
        {
        #pragma omp for reduction(+: ke) schedule(dynamic)
        #endif
        for (i=0;i<N;i++) {
            ke += M*star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+
                pow(star[i][6],2));
        }
        #ifndef NOOMP
        }
        #endif

        ke /= 0.5*N;
        printf("Dynamical temperature of centre-of-mass particles kT = "
            "%lf\n\n",ke);
        // make half-mass radius of the system match the desired one
        radial_profile(star, N, rvir, M, 0, 0, code, &NNBMAX, &RS0, &Rh2D,
            &Rh3D, NNBMAX_NBODY6);
        printf("\nmeasuring half-mass radius: %.7f \t %.7f (should/is)\n"
            "and correcting for this factor\n", Rh, Rh3D);
        rvir = rvir *Rh/Rh3D;
    }

    // Calculate radial density profile, estimate NNBMAX and RS0
    // (important for Nbody6 only)
    radial_profile(star, N, rvir, M, create_radial_profile,
        create_cumulative_profile, code, &NNBMAX, &RS0, &Rh2D, &Rh3D,
        NNBMAX_NBODY6);
    printf("\nActual half-mass radius of the cluster ="
        "(%.4f / %.4f) pc (3D / 2D)\n", Rh3D, Rh2D);

    // Scale RS0 to nbody units for Nbody6
    RS0 /= 1.0*rvir;

    /*********************
     * Generate Binaries *
     *********************/

    printf("\n\n-----GENERATE BINARIES-----   \n");

    if ((!fbin) && (!nbin)) {
        printf("\nNo primordial binaries!\n");
    } else {
        // re-create original array with Nstars (original N) entries
        columns = 15;
        double **star_temp;
        star_temp = (double **)calloc(Nstars,sizeof(double *));

        for (j=0;j<Nstars;j++) {
            star_temp[j] = (double *)calloc(columns,sizeof(double));

            if (star_temp[j] == NULL) {
                printf("\nMemory allocation failed!\n");
                return 0;
            }
        }

        double **mbin_index; //sort mbin by identifier
        mbin_index = (double **)calloc(nbin,sizeof(double *));

        for (j=0;j<nbin;j++) {
            mbin_index[j] = (double *)calloc(2,sizeof(double));

            if (mbin_index[j] == NULL) {
                printf("\nMemory allocation failed!\n");
                return 0;
            }
        }

        for (j=0;j<nbin;j++) {
            mbin_index[j][0] = mbin[j][15];
            mbin_index[j][1] = j;
        }

        shellsort(mbin_index,nbin,2);

        // sort star by identifier
        double **star_index;
        star_index = (double **)calloc(N,sizeof(double *));

        for (j=0;j<N;j++){
            star_index[j] = (double *)calloc(2,sizeof(double));

            if (star_index[j] == NULL) {
                printf("\nMemory allocation failed!\n");
                return 0;
            }
        }

        for (j=0;j<N;j++) {
            star_index[j][0] = star[j][7];
            star_index[j][1] = j;
        }

        shellsort(star_index,N,2);

        int p;
        for (j=0;j<nbin;j++) {
            int mbin_itmp = mbin_index[j][1];
            int star_itmp = star_index[j][1];

            if (mbin[(int)mbin_itmp][15] == star[(int)star_itmp][7]) {

                // primary mass
                star_temp[2*j][0] = mbin[mbin_itmp][1]/(1.0*M);
                // secondary mass
                star_temp[2*j+1][0] = mbin[mbin_itmp][2]/(1.0*M);
                // primary m0
                star_temp[2*j][7] = mbin[mbin_itmp][3];
                // secondary m0
                star_temp[2*j+1][7] = mbin[mbin_itmp][9];
                // primary kw
                star_temp[2*j][8] = mbin[mbin_itmp][4];
                // secondary kw
                star_temp[2*j+1][8] = mbin[mbin_itmp][10];
                // primary epoch
                star_temp[2*j][9] = mbin[mbin_itmp][5];
                //secondary epoch
                star_temp[2*j+1][9] = mbin[mbin_itmp][11];
                // primary spin
                star_temp[2*j][10] = mbin[mbin_itmp][6];
                // secondary spin
                star_temp[2*j+1][10] = mbin[mbin_itmp][12];
                // primary r
                star_temp[2*j][11] = mbin[mbin_itmp][7];
                // secondary r
                star_temp[2*j+1][11] = mbin[mbin_itmp][13];
                // primary lum
                star_temp[2*j][12] = mbin[mbin_itmp][8];
                // secondary lum
                star_temp[2*j+1][12] = mbin[mbin_itmp][14];
                // primary epochstar
                star_temp[2*j][13] = mbin[mbin_itmp][16];
                // secondary epochstar
                star_temp[2*j+1][13] = mbin[mbin_itmp][17];
                // primary Zstar
                star_temp[2*j][14] = mbin[mbin_itmp][18];
                // secondary Zstar
                star_temp[2*j+1][14] = mbin[mbin_itmp][19];

                for (p = 1; p < 7; p++) {
                    star_temp[2*j][p] = star[star_itmp][p];
                    star_temp[2*j+1][p] = star[star_itmp][p];
                }
            }
        }

        for (j=nbin;j<N;j++) {
            for (p=0;p<columns;p++) {
                star_temp[j+nbin][p] = star[(int)star_index[j][1]][p];
            }
        }

        N += nbin;

        for (j=0;j<N;j++) {
            for (p=0;p<columns;p++) {
                //copy temporary array back to original
                star[j][p] = star_temp[j][p];
            }
        }

        for (j=0;j<Nstars;j++)
            free (star_temp[j]);

        free(star_temp);

        printf("\nCreating %i primordial binary systems, fraction: "
            "%6.2f percent.\n", nbin, 2.0*nbin/N*100.0);

        if (seed)
            srand48(seed);

        get_binaries(nbin, star, M, rvir, pairing, &N, adis, amin, amax, Rh,
            eigen, BSE, epoch, Z, remnant, OBperiods, msort);
    }

    // Specify KZ(22) & the sse parameter
    #ifdef SSE
    if (epoch) {
        // If feeding an evolved stellar population to Nbody6, KZ(12) has to
        // be =2 in order to read-in fort.12
        sse = 1;
    }
    else
        sse = 0;
    #else
    sse = 0;
    #endif

    // KZ(22)
    bin = 4;


    /***********
     * Scaling *
     ***********/

    printf("\n\n-----SCALING-----      \n");

    // scale masses, pos & vel to astrophysical units or Nbody units
    tscale = sqrt(rvir*rvir*rvir/(G*M));

    if (units) {
        printf("\nScaling to astrophysical units.\n");

        for (j=0; j<N; j++)
            star[j][0] *= M;

        for (j=0; j<N; j++) {
            for (i=1;i<4;i++)
                star[j][i] *= rvir;
        }

        for (j=0; j<N; j++) {
            for (i=4;i<7;i++)
                star[j][i] *= rvir/tscale;
        }

        bin = -1; //KZ(22)
    } else {
        printf("\nScaling to Nbody units.\n");
    }

    // Scale mass, radius and decay time of external (gas) potential to
    // Nbody units
    if (!(code == 3 && units)) {

        if (extmass)
            extmass /= M;

        if (extrad)
            extrad /= rvir;

        if (extdecay)
            extdecay = 1.0/(extdecay/tscale);

        if (extstart)
            extstart = extstart/tscale;
    }


    /**********
     * Output *
     **********/

    printf("\n\n-----OUTPUT-----      \n");

    if (code == 0) {
        output0(output, N, NNBMAX, RS0, dtadj, dtout, tcrit, rvir, mmean, tf,
            regupdate, etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch,
            dtplot, Z, nbin, Q, RG, VG, rtide, gpu, star, sse, seed, extmass,
            extrad, extdecay, extstart);
    }
    else if (code == 1) {
        output1(output, N, dtadj, dtout, tcrit, rvir, mmean, tf, regupdate,
            etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch, Z, nbin, Q,
            RG, VG, rtide, gpu, star);
    }
    else if (code == 2) {
        output2(output, N, NNBMAX, RS0, dtadj, dtout, tcrit, rvir, mmean, tf,
            regupdate, etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch,
            dtplot, Z, nbin, Q, RG, VG, rtide, gpu, star, sse, seed, extmass,
            extrad, extdecay, extstart);
    }
    else if (code == 3) {
        output3(output, N, rvir, Rh, mmean, M, epoch, Z, RG, VG, rtide, star,
            Rgal, extmass, extrad);
    }
    else if (code == 4) {
        output4(output, N, NNBMAX, RS0, dtadj, dtout, tcrit, rvir, mmean, tf,
            regupdate, etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch,
            dtplot, Z, nbin, Q, RG, VG, rtide, gpu, star, sse, seed, extmass,
            extrad, extdecay, extstart);
    }
    else if (code == 5) {
        output5(output, N, NNBMAX, RS0, dtadj, dtout, tcrit, rvir, mmean, tf,
            regupdate, etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch,
            dtplot, Z, nbin, Q, RG, VG, rtide, gpu, star, sse, seed, extmass,
            extrad, extdecay, extstart);
    }

    /**********************
     * Final energy check *
     **********************/

    printf("\n\n-----FINISH-----  \n");

    if (check) {
        printf("\nMaking final energy check... "
            "(may take a while but can be aborted by pressing CTRL+c)\n");

        #ifndef NOOMP
        #pragma omp parallel shared(N, star)  private(i, j)
        {
        #pragma omp for reduction(+: ekin, epot, sigma) schedule(dynamic)
        #endif
            for (j=0; j<N; j++) {
                ekin += star[j][0]*((star[j][4]*star[j][4])+
                    (star[j][5]*star[j][5])+(star[j][6]*star[j][6]));
                if (j) {
                    for (i=0;i<j-1;i++)
                        epot -= star[i][0]*star[j][0]/sqrt(
                            (star[i][1]-star[j][1])*(star[i][1]-star[j][1])+
                            (star[i][2]-star[j][2])*(star[i][2]-star[j][2])+
                            (star[i][3]-star[j][3])*(star[i][3]-star[j][3]));
                }
                sigma += star[j][4]*star[j][4]+star[j][5]*star[j][5]+
                    star[j][6]*star[j][6];
            }
        #ifndef NOOMP
        }
        #endif

        if (units)
            epot *= G;

        ekin *= 0.5;
        sigma = sqrt(sigma/N);
        tscale = sqrt(rvir*rvir*rvir/(G*M));

        printf("\nEkin = %g\t Epot = %g\t Etot = %g \t kT = %g",
            ekin, epot, ekin+epot, ekin/(N-nbin));
        printf("\nVel.Disp. = %g\tCross.Time = %g \n", sigma, 2.0/sigma);

        if (units) {
            printf("Vel.Disp. = %g\tCross.Time = %g (Nbody units)\n",
                sigma/rvir*tscale, 2.0/sigma/tscale);
        }
        else {
            printf("Vel.Disp. = %g\tCross.Time = "
                "%g (physical units, km/s, Myr)\n",
                sigma*rvir/tscale, 2.0/sigma*tscale);
        }
    }

    free(Mcum);

    for (j=0;j<nbin;j++)
        free (mbin[j]);

    free(mbin);

    for (j=0;j<NMAX;j++)
        free (star[j]);

    free(star);

    #ifdef NOOMP
    // stop stop-watch
    t2 = clock();
    // print stopped time
    printf("\nElapsed time: %g sec\n",
        (double)(t2-t1)/CLOCKS_PER_SEC);
    #else
    #pragma omp parallel
    {
        // stop stop-watch
        t2 = omp_get_wtime();
    }
    //print stopped time
    printf("\nElapsed time: %g sec\n",t2-t1);
    #endif

    return 0;
}
