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

#include "include/functions.h"


int main (int argv, char **argc) {

    struct options opt;

    // SSE internal parameters (see Hurley, Pols & Tout 2000)
    struct value1 value1_;
    struct points points_;

    // BSE internal parameters (see Hurley, Pols & Tout 2002)
    struct flags flags_;
    struct value5 value5_;
    struct value2 value2_;
    struct value3 value3_;


    // Command line input
    int option;
    while ((option = getopt(argv, argc, "N:M:P:W:R:r:c:g:S:D:T:Q:C:A:O:G:o:f:a:m:B:b:p:s:t:e:Z:X:x:V:u:h:?")) != -1) switch (option)
    {
        case 'N': opt.N = atoi(optarg); opt.Mcl = 0.0; break;
        case 'M': opt.Mcl = atof(optarg); opt.N = 0; break;
        case 'P': opt.profile = atoi(optarg); break;
        case 'W': opt.W0 = atof(optarg); break;
        case 'R': opt.Rh = atof(optarg); break;
        case 'r': opt.a = atof(optarg); break;
        case 'c': opt.Rmax = atof(optarg); break;
        case 'g':
            if (opt.gn < 3) { opt.gamma[opt.gn] = atof(optarg); opt.gn++; break; }
        case 'S': opt.S = atof(optarg); break;
        case 'D': opt.D = atof(optarg); break;
        case 'T': opt.tcrit = atof(optarg); break;
        case 'Q': opt.Q = atof(optarg); break;
        case 'C': opt.code = atoi(optarg); break;
        case 'A': opt.dtadj = atof(optarg); break;
        case 'O': opt.dtout = atof(optarg); break;
        case 'G': opt.gpu = atoi(optarg); break;
        case 'o': opt.output = optarg; break;
        case 'f': opt.mfunc = atoi(optarg); break;
        case 'a' :
            if (opt.an < MAX_AN) {
                opt.alpha[opt.an] = atof(optarg);
                if (opt.an == 0) opt.alpha_L3 = atof(optarg);
                if (opt.an == 1) opt.beta_L3 = atof(optarg);
                if (opt.an == 2) opt.mu_L3 = atof(optarg);
                opt.an++;
                break;
            } else { printf("\nError: Number of alphas exceeded maximum limit of %d\n", MAX_AN); return 1; }
        case 'm' :
            if (opt.mn < MAX_MN) {
                opt.mlim[opt.mn] = atof(optarg);
                if (opt.mn == 0) opt.mlow = atof(optarg);
                if (opt.mn == MAX_MN-1) opt.mup = atof(optarg);
                opt.mn++;
                break;
            } else { printf("\nError: Number of mass params exceded maximum limit of %d\n", MAX_MN); return 1; }
        case 'B': opt.nbin = atoi(optarg); break;

        case 'b': opt.fbin = atof(optarg); break;
        case 'p': opt.pairing = atoi(optarg); break;
        case 's': opt.seed = atoi(optarg); break;
        case 't': opt.tf = atoi(optarg); break;
        case 'e': opt.epoch = atof(optarg); break;
        case 'Z': opt.Z = atof(optarg); break;
        case 'X' :
            if (opt.xn < 3) { opt.RG[opt.xn] = atof(optarg); opt.xn++; break; }
        case 'V' :
            if (opt.vn < 3) { opt.VG[opt.vn] = atof(optarg); opt.vn++; break; }
        case 'x' :
            if (opt.xx < 4) { opt.extgas[opt.xx] = atof(optarg); opt.xx++; break; }
        case 'u': opt.units = atoi(optarg); break;
        case ':':
        case 'h':    help(opt.msort); return 1;
        case '?':    help(opt.msort); return 1;
    };

    if (opt.mn-1 > 0) opt.mup = opt.mlim[opt.mn-1];

    // print summary of input parameters to .info file
    info(opt);

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

    if (opt.seed) {
        // initialize random number generator by seed
        srand48(opt.seed);
    }
    else {
        // initialize random number generator by local time
        opt.seed = (unsigned) time(NULL);
        opt.seed %= 100000;
        srand48(opt.seed);
    }

    printf ("\n\nRandom seed = %i\n\n\n", opt.seed);
    if (opt.seed) {
        // idum is the random number seed used in the kick routine.
        value3_.idum = opt.seed;
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
    if (opt.profile == 3)
        opt.Rh = opt.a;

    if ((opt.Mcl) && (opt.N)) {
        printf("\nWARNING:"
            "specify either Mcl (-M) or N (-N)!\n\n");
        exit (1);
    } else if ((!opt.Mcl) && (opt.mfunc == 3)) {
        printf("\nWARNING:"
            "specify Mcl (-M) when using optimal sampling (-f 3)!\n\n");
        exit (1);
    }

    if (opt.xx > 0) {
        if ((opt.xx == 4) && (opt.extgas[2])) {
            opt.extmass = opt.extgas[0];
            opt.extrad = opt.extgas[1];
            opt.extdecay = opt.extgas[2];
            opt.extstart = opt.extgas[3];
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

    struct star_data star[opt.NMAX];
    for (j = 0; j < opt.NMAX; j++) {
        star[j] = {0};
    }



    /*******************************************
     * Evaluate Z from [Fe/H] if Z is set to 0 *
     *******************************************/

     if (!opt.Z) {
        // Bertelli, Bressan, Chiosi, Fagotto, Nasi, 1994, A&AS, 106, 275
        opt.Z = pow(10.0, 0.977*opt.FeH)*opt.Zsun;
        printf("\nUsing Bertelli et al. (1994) relation to convert FeH = %.3f"
            " into Z = %.3f\n", opt.FeH, opt.Z);
    }



    /**********************************
     * Calculate maximum stellar mass *
     **********************************/

    //always use Weidner relation when using optimal sampling!
     if (opt.mfunc == 3)
        opt.weidner = 1;

    if (!opt.N && opt.weidner && opt.mfunc) {
        opt.mup = opt.upper_IMF_limit;

        printf("\nUsing maximum stellar mass-cluster mass relation for upper"
            " stellar mass limit\n");
        printf("\n(Weidner & Kroupa 2006, Pflamm-Altenburg & Kroupa 2007)\n");

        // analytic fit to the observational data from
        // Pflamm-Altenburg & Kroupa (2007), implementation by M. Kruckow
        MMAX = pow(10,2.56*log10(opt.Mcl) * pow(pow(3.82,9.17) + pow(log10(opt.Mcl),
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
        MMAX = opt.mup;
    }

    if (opt.mfunc && opt.epoch && opt.prantzos) {
        printf("\nUsing Prantzos (2007) relation to reduce upper mass"
            " limit to Lifetime(mup) > epoch\n");
        while (Lifetime(MMAX) < sqrt(pow(opt.epoch,2))) {
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
    if (opt.mfunc == 1) {
        printf("\nMaximum stellar mass set to: %.2f\n", MMAX);
        generate_m1(&opt, star, &M, &mmean, MMAX);
    } else if (opt.mfunc == 2) {

        if (opt.mn) {
            for (i = opt.mn+1; i < MAX_MN; i++) opt.mlim[i] = 0.0;
        } else {
            for (i=0; i<MAX_MN; i++) {
                if (opt.mlim[i]) opt.mn++;
            }
        }

        if (opt.an) {
            for (i = opt.an+1; i < MAX_AN; i++) opt.alpha[i] = 0.0;
        } else {
            for (i=0; i<MAX_AN; i++) {
                if (opt.alpha[i]) opt.an++;
            }
        }

        if (opt.an >= opt.mn)
            opt.an = opt.mn - 1;

        opt.mn = opt.an + 1;

        if (!opt.mn){
            printf("\nError: at least one mass limit has to be specified\n");
            return 1;
        } else if (opt.mn == 1) {
            opt.single_mass = opt.mlim[0];
            printf("\nSetting stellar masses to %g solar mass\n",opt.single_mass);

            if (!opt.N)
                opt.N = opt.Mcl/opt.single_mass;

            for (j=0;j<opt.N;j++)
                star[j].mass_epoch = 1.0/opt.N;

            mmean = opt.single_mass;
            M = opt.N*mmean;
            printf("\nM = %g\n", M);
        } else {
            printf("\nMaximum stellar mass set to: %.2f\n", MMAX);
            // normalization factor of integral
            norm[opt.an-1] = 1.;
            // integrated number of stars in interval [mlim[an-1]:mlim[an]]
            N_tmp = subcount[opt.an-1] = subint(opt.mlim[opt.an-1], opt.mlim[opt.an],
                opt.alpha[opt.an-1] + 1.);
            // integrated mass of stars in interval [mlim[an-1]:mlim[an]]
            M_tmp = submass[opt.an-1] = subint(opt.mlim[opt.an-1], opt.mlim[opt.an],
                opt.alpha[opt.an-1] + 2.);

            for (i = opt.an - 2; i >= 0; i--) {
                norm[i] = norm[i+1] * pow(opt.mlim[i+1], opt.alpha[i+1] - opt.alpha[i]);
                subcount[i] = norm[i] * subint(opt.mlim[i], opt.mlim[i+1],
                    opt.alpha[i] + 1.);
                N_tmp += subcount[i];
                submass[i] = norm[i] * subint(opt.mlim[i], opt.mlim[i+1],
                    opt.alpha[i] + 2.);
                M_tmp += submass[i];
            }

            generate_m2(&opt, star, M_tmp, subcount, &mmean, &M, MMAX);
        }
    } else if (opt.mfunc == 3) {
        printf("\nMaximum stellar mass set to: %.2f\n",MMAX);
        generate_m3(&opt, star, &M, &mmean, MMAX);
        randomize(star, opt.N);
    } else if (opt.mfunc == 4) {
        printf("\nMaximum stellar mass set to: %.2f\n", MMAX);
        printf("\nUsing L3 IMF (Maschberger 2012)\n");
        generate_m4(&opt, star, &mmean, MMAX, &M);
    } else {
        printf("\nSetting stellar masses to %.1f solar mass\n", opt.single_mass);
        if (!opt.N) opt.N = opt.Mcl/opt.single_mass;
        for (j=0;j<opt.N;j++) {
            star[j].mass_epoch = opt.single_mass;
            star[j].mass = opt.single_mass;
            star[j].kstar = 0;
            star[j].epoch1 = 0.0;
            star[j].ospin = 0.0;
            star[j].rstar = 0.0;
            star[j].lstar = 0.0;
            star[j].epochstar = 0.0;
            star[j].zstar = 0.0;
        }
        mmean = opt.single_mass;
        M = opt.N*mmean;
        printf("\nM = %g\n", M);
        opt.mloss = 0;
    }

    // set all stars to the same metallicity and age for now
    double epochstar, zstar;
    // age compared to the oldest stars in the cluster [Myr]
    epochstar = 0.0;
    zstar = opt.Z;

    for (j=0;j<opt.N;j++) {
        star[j].epochstar = epochstar;
        star[j].zstar = zstar;
    }

    // Pair binary masses and convert to centre-of-mass particles
    int Nstars;
    if (!opt.nbin)
        opt.nbin = 0.5*opt.N*opt.fbin;

    // component mass & stellar evol parameter array
    double **mbin;
    mbin = (double **)calloc(opt.nbin,sizeof(double *));
    for (j=0;j<opt.nbin;j++) {
        mbin[j] = (double *)calloc(20,sizeof(double));

        if (mbin[j] == NULL) {
            printf("\nMemory allocation failed!\n");
            return 0;
        }
    }

    Nstars = opt.N;

    if (opt.nbin) {
        printf("\nPreparing binary components.\n");
        // Specify component opt.pairing
        if (opt.pairing) {

            if (opt.pairing == 1) {
                printf("\nApplying ordered pairing for stars with masses"
                    " > %.1f Msun.\n", opt.msort);
            }
            else if (opt.pairing == 2) {
                printf("\nApplying random pairing for stars with masses"
                    " > %.1f Msun.\n", opt.msort);
            }
            else if (opt.pairing == 3) {
                printf("\nApplying uniform mass ratio distribution for stars"
                    " with masses > %.1f Msun.\n", opt.msort);
            }

//            order(star, opt.N, M, opt.msort, opt.pairing);
        } else {
//            randomize(star, opt.N);
            printf("\nApplying random pairing.\n");
        }
        opt.N -= opt.nbin;
        for (j=0;j<opt.nbin;j++) {
            // system mass
            mbin[j][0] = star[2*j].mass_epoch+star[2*j+1].mass_epoch;
            // primary mass
            mbin[j][1] = star[2*j].mass_epoch;
            // secondary mass
            mbin[j][2] = star[2*j+1].mass_epoch;
            // primary m0
            mbin[j][3] = star[2*j].mass;
            // primary kw
            mbin[j][4] = star[2*j].kstar;
            // primary epoch1
            mbin[j][5] = star[2*j].epoch1;
            // primary spin
            mbin[j][6] = star[2*j].ospin;
            // primary r
            mbin[j][7] = star[2*j].rstar;
            // primary lum
            mbin[j][8] = star[2*j].lstar;
            // secondary m0
            mbin[j][9] = star[2*j+1].mass;
            // secondary kw
            mbin[j][10] = star[2*j+1].kstar;
            // secondary epoch1
            mbin[j][11] = star[2*j+1].epoch1;
            // secondary spin
            mbin[j][12] = star[2*j+1].ospin;
            //secondary r
            mbin[j][13] = star[2*j+1].rstar;
            // secondary lum
            mbin[j][14] = star[2*j+1].lstar;
            // identifier
            mbin[j][15] = 1000+j;
            // primary epochstar
            mbin[j][16] = star[2*j].epochstar;
            // secondary epochstar
            mbin[j][17] = star[2*j+1].epochstar;
            // primary zstar
            mbin[j][18] = star[2*j].zstar;
            // secondary zstar
            mbin[j][19] = star[2*j+1].zstar;

            //s ystem mass
            star[2*j].mass_epoch += star[2*j+1].mass_epoch;
            star[2*j+1].mass_epoch = 0.0;
            // identifier
            star[2*j].mass = 1000+j;
            // identifier
            star[2*j+1].mass = 0.0;
            // system luminosity
            star[2*j].lstar += star[2*j+1].lstar;
            star[2*j+1].lstar = 0.0;
        }

//        order(star, Nstars, M, 0.0, 0);
//        randomize(star, opt.N);
    }


    // prepare mass segregation

    // search lowest mass star
    double mlowest = MMAX;
    //search highest mass star
    double mhighest = 0;
    double mmeancom = 0.0;

    for (i=0;i<opt.N;i++) {
        // scale masses to Nbody units
        star[i].mass_epoch /= M;
        mmeancom += star[i].mass_epoch;

        if (star[i].mass_epoch < mlowest)
            mlowest = star[i].mass_epoch;

        if (star[i].mass_epoch > mhighest)
            mhighest = star[i].mass_epoch;
    }
    mmeancom /= opt.N;
    // number of necessary pos & vel pairs for Baumgardt et al. (2008)
    // mass segregation routine
    int Nseg = ceil(opt.N*mmeancom/mlowest);
    int Nunseg = opt.N;

    double *Mcum;
    Mcum = (double *)calloc(opt.N,sizeof(double));

    // sort masses when mass segregation parameter > 0
    if ((opt.S) && !(opt.profile == 2)) {
        printf("\nApplying mass segregation with S = %f\n",opt.S);
//        segregate(star, opt.N, opt.S);

        //calculate cumulative mass function Mcum
        for (i=0;i<opt.N;i++) {
            Mcum[i] = 0.0;
            for (j=0;j<=i;j++) Mcum[i] = Mcum[i] + star[j].mass_epoch;
        }
        opt.N = Nseg;
    }


    /*************************************
     * Generate positions and velocities *
     *************************************/

    printf("\n\n-----GENERATE POSITIONS & VELOCITIES-----   \n");

    // calculate half-mass radius according to Marks & Kroupa 2012
    // if Rh is set to -1
    if ((opt.Rh == -1) && (opt.Mcl)) {
        // Marks & Kroupa (2012), implementation by M. Kruckow
        opt.Rh = 0.1*pow(opt.Mcl,0.13);
        printf("\nUsing Marks & Kroupa (2012) relation to derive half-mass"
            " radius from cluster mass: %g (pc)\n", opt.Rh);
    }

    // evaluate approximate tidal radius assuming circular orbit
    if (opt.tf == 3) {
        // in the case of Allen & Santillan potential,
        // assume kappa = 1.4omega (eq. 9 in Kuepper et al. 2010)
        omega = sqrt(opt.VG[0]*opt.VG[0]+opt.VG[1]*opt.VG[1]+opt.VG[2]*opt.VG[2])/sqrt(opt.RG[0]*opt.RG[0]+
            opt.RG[1]*opt.RG[1]+opt.RG[2]*opt.RG[2]);
        rtide = pow(G*M/(2.0*omega*omega),1.0/3.0);
    } else if (!opt.tf) {
        rtide = 1.0E5;
    } else if ((opt.tf == 1) && (opt.code == 0 || opt.code == 4 || opt.code == 5)) {
        // in case of Sverre's Nbody6 standard tidal field
        rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*8000.0;
    } else {
        // in the case of a point mass potential or near field approximation
        rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*sqrt(opt.RG[0]*opt.RG[0]+opt.RG[1]*
            opt.RG[1]+opt.RG[2]*opt.RG[2]);
    }
    printf("\nApproximate tidal radius: %g (pc)\n", rtide);


    // generate scaled pos & vel, postpone scaling for Plummer and King in
    // case of mass segregation
    double rhtemp, rvirtemp;

    if (opt.profile == 1) {
        printf("\nGenerating King model with parameters: "
            "N = %i\t W0 = %g\t Rh = %.3f\t D = %.2f\n",opt.N, opt.W0, opt.Rh, opt.D);
//        generate_king(opt.N, opt.W0, star, &rvirtemp, &rhtemp, &rking, opt.D, opt.symmetry);
    } else if (opt.profile == 2) {
        opt.N = Nunseg;
        printf("\nGenerating segregated Subr model with parameters: "
            "N = %i\t S = %g\t Rh = %.3f\n",opt.N, opt.S, opt.Rh);
        // value provided by L. Subr
        rvir = opt.Rh/0.76857063065978;
//        generate_subr(opt.N, opt.S, star, rtide, rvir);
        printf ("\nrvir = %.5f\t rh = %.5f\t rtide = %.5f (pc)\n",
            rvir, opt.Rh, rtide);
    } else if (opt.profile == 3) {
        if (opt.gamma[1] == 0.0 && opt.gamma[2] == 2.0) {
            printf("\nGenerating EFF model with parameters: "
                "N = %i\t a = %.1f\t gamma = %.3f\t Rmax = %.1f\t D = %.2f\n",
                opt.N, opt.a, opt.gamma[0], opt.Rmax, opt.D);
        }
        else {
            printf("\nGenerating Nuker model with parameters: "
                "N = %i\t a = %.1f\t gamma (outer) = %.3f\t gamma (inner) "
                "= %.3f\t transition = %.3f\t Rmax = %.1f\t D = %.2f\n",
                opt.N, opt.a, opt.gamma[0],opt.gamma[1],opt.gamma[2], opt.Rmax, opt.D);
        }

        double p[6];

        // rho0, will be scaled according to Rmax and Mtot
        p[1] = 1.0;
        //scale radius
        p[2] = opt.a;
        //outer power-law slope (>0.5)
        p[3] = opt.gamma[0];
        //inner power-law slope (0.0 for EFF template)
        p[4] = opt.gamma[1];
        //transition parameter (2.0 for EFF template)
        p[5] = opt.gamma[2];

//        generate_profile(opt.N, star, opt.Rmax, M, p, &opt.Rh, opt.D, opt.symmetry);
        printf("\nRh = %.1f pc\n", opt.Rh);
    } else if (opt.profile == -1) {
        printf("\nGenerating fractal distribution with parameters: "
            "N = %i\t Rh = %.3f\t D = %.2f\n", opt.N, opt.Rh, opt.D);
//        fractalize(opt.D, opt.N, star, 0, opt.symmetry);
        rvir = opt.Rh;
    } else {
        printf("\nGenerating Plummer model with parameters: "
            "N = %i\t Rh = %.3f\t D = %.2f\n", opt.N, opt.Rh, opt.D);
        rvir = opt.Rh/0.772764;
        rplummer = opt.Rh/1.305;
//        generate_plummer(opt.N, star, rtide, rvir, opt.D, opt.symmetry);
    }

    // Apply Baumgardt et al. (2008) mass segregation
    if (!(opt.profile == 2) && (opt.S)) {
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
            m_temp[i][0] = star[i].mass_epoch;
            m_temp[i][1] = star[i].mass;
            m_temp[i][2] = star[i].kstar;
            m_temp[i][3] = star[i].epoch1;
            m_temp[i][4] = star[i].ospin;
            m_temp[i][5] = star[i].rstar;
            m_temp[i][6] = star[i].lstar;
            m_temp[i][7] = star[i].epochstar;
            m_temp[i][8] = star[i].zstar;
        }

        printf("\nOrdering orbits by energy.\n");
//        energy_order(star, opt.N, Nstars);

        int nlow, nhigh, nrandom;
        for (i=0;i<Nunseg;i++) {
            nhigh = Nseg*Mcum[i];
            if (i) {
                nlow = Nseg*Mcum[i-1];
            } else {
                nlow = 0;
            }
            nrandom = (nhigh-nlow)*drand48()+nlow;
            star[i].mass_epoch = m_temp[i][0];
            star[i].rx = star[nrandom].rx;
            star[i].ry = star[nrandom].ry;
            star[i].rz = star[nrandom].rz;
            star[i].vx = star[nrandom].vx;
            star[i].vy = star[nrandom].vy;
            star[i].vz = star[nrandom].vz;
            star[i].mass = m_temp[i][1];
            star[i].kstar = m_temp[i][2];
            star[i].epoch1 = m_temp[i][3];
            star[i].ospin = m_temp[i][4];
            star[i].rstar = m_temp[i][5];
            star[i].lstar = m_temp[i][6];
            star[i].epochstar = m_temp[i][7];
            star[i].zstar = m_temp[i][8];
        }


        for (j=0;j<Nunseg;j++)
            free (m_temp[j]);

        free(m_temp);
        opt.N = Nunseg;
    }

    // CoM correction
    printf("\nApplying centre-of-mass correction.\n");
    for (j=0; j<7; j++)
        cmr[j] = 0.0;

    for (j=0; j<opt.N; j++) {
        cmr[1] += star[j].mass_epoch*star[j].rx;
        cmr[2] += star[j].mass_epoch*star[j].ry;
        cmr[3] += star[j].mass_epoch*star[j].rz;
        cmr[4] += star[j].mass_epoch*star[j].vx;
        cmr[5] += star[j].mass_epoch*star[j].vy;
        cmr[6] += star[j].mass_epoch*star[j].vz;
    }

    for (j=0; j<opt.N; j++) {
            star[j].rx -= cmr[i];
            star[j].ry -= cmr[i];
            star[j].rz -= cmr[i];
            star[j].vx -= cmr[i];
            star[j].vy -= cmr[i];
            star[j].vz -= cmr[i];
    }


    // Apply scaling to Nbody-units
    if (opt.profile == 0) {
        printf("\nRe-scaling of orbits (dt ~ N^2!)\n");
        double ke = 0.0;
        double pe = 0.0;
        double sx, sv, r2;

        #ifdef GPU
        gpupot(opt.N,star,&pe);
        #ifndef NOOMP
        int Ntemp = opt.N;
        #pragma omp parallel shared(Ntemp, star)  private(i)
        {
            #pragma omp for reduction(+: ke) schedule(dynamic)
        #endif
            for (i=0;i<opt.N;i++) {
                ke += star[i].mass_epoch*(pow(star[i].vx,2)+pow(star[i].vy,2)+
                    pow(star[i].vz,2));
            }
        #ifndef NOOMP
        }
        #endif
        //        printf("PE_GPU %lg ke %lg\n",pe,ke);
        #else
        #ifndef NOOMP
        int Ntemp = opt.N;
        #pragma omp parallel shared(Ntemp, star)  private(i, j, r2)
        {
            #pragma omp for reduction(+: pe, ke) schedule(dynamic)
        #endif

        for (i=0;i<opt.N;i++) {
            if (i) {
                for (j=0;j<i-1;j++) {
                    r2=(star[i].rx - star[j].rx) * (star[i].rx - star[j].rx) +
                       (star[i].ry - star[j].ry) * (star[i].ry - star[j].ry) +
                       (star[i].rz - star[j].rz) * (star[i].rz - star[j].rz) ;
                    pe -= star[i].mass_epoch*star[j].mass_epoch/sqrt(r2);
                }
            }
            ke += star[i].mass_epoch*(pow(star[i].vx,2)+pow(star[i].vy,2)+
                pow(star[i].vz,2));
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

        for (i=0;i<opt.N;i++) {
            star[i].rx *= sx;
            star[i].ry *= sx;
            star[i].rz *= sx;
            star[i].vx /= sv;
            star[i].vy /= sv;
            star[i].vz /= sv;
        }

        ke = 0;
        for (i=0;i<opt.N;i++) {
            ke += M*star[i].mass_epoch*(pow(star[i].vx,2)+pow(star[i].vy,2)+
                pow(star[i].vz,2));
        }
        ke /= 0.5*opt.N;
        printf("Dynamical temperature of centre-of-mass particles kT = "
            "%lf\n\n", ke);

        // make half-mass radius of the system match the desired one
//        radial_profile(star, opt, rvir, M, 0, 0, &NNBMAX, &RS0, &Rh2D,
//            &Rh3D);

        if (opt.match) {
            printf("\nmeasuring half-mass radius: %.7f \t %.7f (should/is)\n"
                "and correcting for this factor\n", opt.Rh, Rh3D);
            rvir = rvir *opt.Rh/Rh3D;
        }

        printf ("\nrvir = %.5f\t rh = %.5f\t rplummer = %.5f\t rtide = "
            "%.5f (pc)\n", rvir, opt.Rh, rplummer, rtide);
    } else if ((opt.profile == 1) || (opt.profile == 3) || (opt.profile == -1)) {
        printf("\nRe-scaling of orbits (dt ~ N^2!)\n");

        double pe = 0.0;
        double ke = 0.0;
        double r2, vscale;

        #ifdef GPU
        gpupot(opt.N,star,&pe);
        #ifndef NOOMP
        #pragma omp parallel shared(opt.N, star)  private(i)
        {
        #pragma omp for reduction(+: ke) schedule(dynamic)
        #endif
          for (i=0;i<opt.N;i++) {
            ke += star[i].mass_epoch*(pow(star[i].vx,2)+pow(star[i].vy,2)+
                pow(star[i].vz,2));
          }
        #ifndef NOOMP
        }
        #endif
        //        printf("PE_GPU %lg ke %lg\n",pe,ke);
        #else
        #ifndef NOOMP
        int Ntemp = opt.N;
        #pragma omp parallel shared(Ntemp, star)  private(i, j, r2)
        {
        #pragma omp for reduction(+: pe, ke) schedule(dynamic)
        #endif
        for (i=0;i<opt.N;i++) {
            for (j=0;j<i-1;j++) {
                r2 = (star[i].rx - star[j].rx) * (star[i].rx - star[j].rx) +
                     (star[i].ry - star[j].ry) * (star[i].ry - star[j].ry) +
                     (star[i].rz - star[j].rz) * (star[i].rz - star[j].rz) ;

                pe -=  star[i].mass_epoch*star[j].mass_epoch/sqrt(r2);
            }
            ke += star[i].mass_epoch*(pow(star[i].vx,2)+pow(star[i].vy,2)+pow(star[i].vz,2));
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
        for (i=0;i<opt.N;i++) {

            star[i].rx /= rvir;
            star[i].ry /= rvir;
            star[i].rz /= rvir;
            star[i].vx /= vscale;
            star[i].vy /= vscale;
            star[i].vz /= vscale;

            ke += M*star[i].mass_epoch*(pow(star[i].vx,2)+pow(star[i].vy,2)+
                pow(star[i].vz,2));
        }
        ke /= 0.5*opt.N;
        printf("Dynamical temperature of centre-of-mass particles kT ="
            "%lf\n\n",ke);

        // make half-mass radius of the system match the desired one
//        radial_profile(star, opt.N, rvir, M, 0, 0, opt.code, &NNBMAX, &RS0, &Rh2D,
//            &Rh3D, opt.NNBMAX_NBODY6);

        if (opt.match) {
            printf("\nmeasuring half-mass radius: %.7f \t %.7f (should/is)\n"
                "and correcting for this factor\n", opt.Rh, Rh3D);
            rvir = rvir *opt.Rh/Rh3D;
        }
        //printf("\nrvir = %.5f\t rh = %.5f\t rtide = %.5f (pc)\n", rvir, Rh,
        //  rtide);
        //printf("Edge radius (King units) = %g\t(Nbody units) = %g\n", rking,
        //  rking/rvirtemp);
        //printf("Core radius (King units) = %g\t(Nbody units) = %g\n\n", 1.0,
        //  1.0/rvirtemp);
        //printf("Concentration = %g\n", log10(rking));
    } else if (opt.profile == 2) {
        double ke = 0;
        #ifndef NOOMP
        int Ntemp = opt.N;
        #pragma omp parallel shared(Ntemp, star)  private(i)
        {
        #pragma omp for reduction(+: ke) schedule(dynamic)
        #endif
        for (i=0;i<opt.N;i++) {
            ke += M*star[i].mass_epoch*(pow(star[i].vx,2)+pow(star[i].vy,2)+
                pow(star[i].vz,2));
        }
        #ifndef NOOMP
        }
        #endif

        ke /= 0.5*opt.N;
        printf("Dynamical temperature of centre-of-mass particles kT = "
            "%lf\n\n",ke);
        // make half-mass radius of the system match the desired one
//        radial_profile(star, opt.N, rvir, M, 0, 0, opt.code, &NNBMAX, &RS0, &Rh2D,
//            &Rh3D, opt.NNBMAX_NBODY6);
        printf("\nmeasuring half-mass radius: %.7f \t %.7f (should/is)\n"
            "and correcting for this factor\n", opt.Rh, Rh3D);
        rvir = rvir *opt.Rh/Rh3D;
    }

    // Calculate radial density profile, estimate NNBMAX and RS0
    // (important for Nbody6 only)
//    radial_profile(star, opt.N, rvir, M, opt.create_radial_profile,
//        opt.create_cumulative_profile, opt.code, &NNBMAX, &RS0, &Rh2D, &Rh3D,
//        opt.NNBMAX_NBODY6);
    printf("\nActual half-mass radius of the cluster ="
        "(%.4f / %.4f) pc (3D / 2D)\n", Rh3D, Rh2D);

    // Scale RS0 to nbody units for Nbody6
    RS0 /= 1.0*rvir;

    /*********************
     * Generate Binaries *
     *********************/

    printf("\n\n-----GENERATE BINARIES-----   \n");

    if ((!opt.fbin) && (!opt.nbin)) {
        printf("\nNo primordial binaries!\n");
    } else {
        // re-create original array with Nstars (original N) entries
        columns = 15;
        struct star_data star_temp[Nstars];

        double **mbin_index; //sort mbin by identifier
        mbin_index = (double **)calloc(opt.nbin,sizeof(double *));

        for (j=0;j<opt.nbin;j++) {
            mbin_index[j] = (double *)calloc(2,sizeof(double));

            if (mbin_index[j] == NULL) {
                printf("\nMemory allocation failed!\n");
                return 0;
            }
        }

        for (j=0;j<opt.nbin;j++) {
            mbin_index[j][0] = mbin[j][15];
            mbin_index[j][1] = j;
        }

//        shellsort(mbin_index,opt.nbin,2);

        // sort star by identifier
        double **star_index;
        star_index = (double **)calloc(opt.N,sizeof(double *));

        for (j=0;j<opt.N;j++){
            star_index[j] = (double *)calloc(2,sizeof(double));

            if (star_index[j] == NULL) {
                printf("\nMemory allocation failed!\n");
                return 0;
            }
        }

        for (j=0;j<opt.N;j++) {
            star_index[j][0] = star[j].mass;
            star_index[j][1] = j;
        }

//        shellsort(star_index,opt.N,2);

        for (j=0;j<opt.nbin;j++) {
            int mbin_itmp = mbin_index[j][1];
            int star_itmp = star_index[j][1];

            if (mbin[mbin_itmp][15] == star[star_itmp].mass) {

                // primary mass
                star_temp[2*j].mass_epoch = mbin[mbin_itmp][1]/(1.0*M);
                // secondary mass
                star_temp[2*j+1].mass_epoch = mbin[mbin_itmp][2]/(1.0*M);
                // primary m0
                star_temp[2*j].mass = mbin[mbin_itmp][3];
                // secondary m0
                star_temp[2*j+1].mass = mbin[mbin_itmp][9];
                // primary kw
                star_temp[2*j].kstar = mbin[mbin_itmp][4];
                // secondary kw
                star_temp[2*j+1].kstar = mbin[mbin_itmp][10];
                // primary epoch
                star_temp[2*j].epoch1 = mbin[mbin_itmp][5];
                //secondary epoch
                star_temp[2*j+1].epoch1 = mbin[mbin_itmp][11];
                // primary spin
                star_temp[2*j].ospin = mbin[mbin_itmp][6];
                // secondary spin
                star_temp[2*j+1].ospin = mbin[mbin_itmp][12];
                // primary r
                star_temp[2*j].rstar = mbin[mbin_itmp][7];
                // secondary r
                star_temp[2*j+1].rstar = mbin[mbin_itmp][13];
                // primary lum
                star_temp[2*j].lstar = mbin[mbin_itmp][8];
                // secondary lum
                star_temp[2*j+1].lstar = mbin[mbin_itmp][14];
                // primary epochstar
                star_temp[2*j].epochstar = mbin[mbin_itmp][16];
                // secondary epochstar
                star_temp[2*j+1].epochstar = mbin[mbin_itmp][17];
                // primary Zstar
                star_temp[2*j].zstar = mbin[mbin_itmp][18];
                // secondary Zstar
                star_temp[2*j+1].zstar = mbin[mbin_itmp][19];


                star_temp[2*j].rx = star[star_itmp].rx;
                star_temp[2*j].ry = star[star_itmp].ry;
                star_temp[2*j].rz = star[star_itmp].rz;
                star_temp[2*j].vx = star[star_itmp].vx;
                star_temp[2*j].vy = star[star_itmp].vy;
                star_temp[2*j].vz = star[star_itmp].vz;
                star_temp[2*j].mass = star[star_itmp].mass;

                star_temp[2*j+1].rx = star[star_itmp].rx;
                star_temp[2*j+1].ry = star[star_itmp].ry;
                star_temp[2*j+1].rz = star[star_itmp].rz;
                star_temp[2*j+1].vx = star[star_itmp].vx;
                star_temp[2*j+1].vy = star[star_itmp].vy;
                star_temp[2*j+1].vz = star[star_itmp].vz;
                star_temp[2*j+1].mass = star[star_itmp].mass;
            }
        }

        for (j=opt.nbin;j<opt.N;j++) {
            int istar = star_index[j][1];

            star_temp[j+opt.nbin] = star[istar];
        }

        opt.N += opt.nbin;

        for (j=0;j<opt.N;j++) {
            star[j] = star_temp[j];
        }

        printf("\nCreating %i primordial binary systems, fraction: "
            "%6.2f percent.\n", opt.nbin, 2.0*opt.nbin/opt.N*100.0);

        if (opt.seed)
            srand48(opt.seed);

//        get_binaries(opt, star, M, rvir, &opt.N);
    }

    // Specify KZ(22) & the sse parameter
    #ifdef SSE
    if (opt.epoch) {
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

    if (opt.units) {
        printf("\nScaling to astrophysical units.\n");

        for (j=0; j<opt.N; j++)
            star[j].mass_epoch *= M;

        for (j=0; j<opt.N; j++) {
            star[j].rx *= rvir;
            star[j].ry *= rvir;
            star[j].rz *= rvir;
        }

        for (j=0; j<opt.N; j++) {
            star[j].vx *= rvir/tscale;
            star[j].vy *= rvir/tscale;
            star[j].vz *= rvir/tscale;
        }

        bin = -1; //KZ(22)
    } else {
        printf("\nScaling to Nbody units.\n");
    }

    // Scale mass, radius and decay time of external (gas) potential to
    // Nbody units
    if (!(opt.code == 3 && opt.units)) {

        if (opt.extmass)
            opt.extmass /= M;

        if (opt.extrad)
            opt.extrad /= rvir;

        if (opt.extdecay)
            opt.extdecay = 1.0/(opt.extdecay/tscale);

        if (opt.extstart)
            opt.extstart = opt.extstart/tscale;
    }


    /**********
     * Output *
     **********/

    printf("\n\n-----OUTPUT-----      \n");

    if (opt.code == 0) {
//        output0(opt, NNBMAX, RS0, rvir, mmean, bin, M, MMAX, rtide, star, sse);
    }
    else if (opt.code == 1) {
//        output1(opt, rvir, mmean, bin, M, MMAX, rtide, star);
    }
    else if (opt.code == 2) {
//        output2(opt, NNBMAX, RS0, rvir, mmean, bin, M, MMAX, rtide, star, sse);
    }
    else if (opt.code == 3) {
//        output3(opt, rvir, mmean, M, rtide, star);
    }
    else if (opt.code == 4) {
//        output4(opt, NNBMAX, RS0, rvir, mmean, bin, M, MMAX, rtide, star, sse);
    }
    else if (opt.code == 5) {
//        output5(opt, NNBMAX, RS0, rvir, mmean, bin, M, MMAX,  rtide, star, sse);
    }

    /**********************
     * Final energy check *
     **********************/

    printf("\n\n-----FINISH-----  \n");

    if (opt.check) {
        printf("\nMaking final energy check... "
            "(may take a while but can be aborted by pressing CTRL+c)\n");

        #ifndef NOOMP
        int Ntemp = opt.N;
        #pragma omp parallel shared(Ntemp, star)  private(i, j)
        {
        #pragma omp for reduction(+: ekin, epot, sigma) schedule(dynamic)
        #endif
            for (j=0; j<opt.N; j++) {
                ekin += star[j].mass_epoch*((star[j].vx*star[j].vx)+
                    (star[j].vy*star[j].vy)+(star[j].vz*star[j].vz));
                if (j) {
                    for (i=0;i<j-1;i++)
                        epot -= star[i].mass_epoch*star[j].mass_epoch/sqrt(
                            (star[i].rx-star[j].rx)*(star[i].rx-star[j].rx)+
                            (star[i].ry-star[j].ry)*(star[i].ry-star[j].ry)+
                            (star[i].rz-star[j].rz)*(star[i].rz-star[j].rz));
                }
                sigma += star[j].vx*star[j].vx+star[j].vy*star[j].vy+
                    star[j].vz*star[j].vz;
            }
        #ifndef NOOMP
        }
        #endif

        if (opt.units)
            epot *= G;

        ekin *= 0.5;
        sigma = sqrt(sigma/opt.N);
        tscale = sqrt(rvir*rvir*rvir/(G*M));

        printf("\nEkin = %g\t Epot = %g\t Etot = %g \t kT = %g",
            ekin, epot, ekin+epot, ekin/(opt.N-opt.nbin));
        printf("\nVel.Disp. = %g\tCross.Time = %g \n", sigma, 2.0/sigma);

        if (opt.units) {
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

    for (j=0;j<opt.nbin;j++)
        free (mbin[j]);

    free(mbin);

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
