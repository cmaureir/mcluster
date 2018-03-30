#ifndef MAIN_H
#define MAIN_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<sys/stat.h>
#include<getopt.h>

#include <iostream>
#include <vector>

#ifdef GPU
#include <cuda.h>
#include <cuda_runtime.h>
#endif

// Use OpenMP if not specified otherwise
#ifndef NOOMP
#include<omp.h>
#endif
//Constants
#define G 0.0043009211 //in pc*km^2*s^-2*Msun
#define Pi 3.14159265
#define PI 3.14159265
#define TWOPI 6.2831853   /* 2PI */
#define GNBODY 1.0
#define MNBODY 1.0
#define PARSEC  3.08568E13      /* KM pro PC */
#define GBIN    1.327126E11     /* G in (km/sec)^3/Msun */
#define RSUN    6.96265E5       /* Solar radius in km */

//Functions
#define mmax(a,b)         (a < b) ?  (b) : (a)
#define mmin(a,b)         (a < b) ?  (a) : (b)
#define Lifetime(Mstar)  1.13E4*pow(Mstar,-3)+0.6E2*pow(Mstar,-0.75)+1.2 //Myr    [Prantzos 2007]


//MW potential consisting of bulge, disk & NFW halo - constants:
//based on Law, Majewski & Johnston (2009)
#define b1_LMJ 0.7000;       //[kpc]
#define M1_LMJ 3.4e10;       //[solar masses]
#define a2_LMJ 6.5000;       //[kpc]
#define b2_LMJ 0.2600;        //[kpc]
#define M2_LMJ 1.0e11;       //[solar masses]
#define R_NFW 42.9647;       //[kpc]
#define M_NFW 9.832538e+11;  //[solar masses]
#define q_halo 1.0;          //[dimensionless]

//Allen & Santillan (1991) MW potential - constants:
#define b1allen 0.3873            //kpc
#define M1allen 606.0*2.32e07    //solar masses
#define a2allen 5.3178
#define b2allen 0.2500
#define M2allen 3690.0*2.32e07
#define a3allen 12.0000
#define M3allen 4615.0*2.32e07

//Additional parameters for Sverre's MW potential
#define VCIRC 220.0            //circular velocity at RCIRC
#define RCIRC 8.500            //kpc

//new input parameters for bulge potential
#define GMB 0.0                //Central bulge mass (Msun).
#define AR 1.0                //Scale radius in gamma/eta model (kpc, Dehnen 1993).
#define GAM 1.0                //gamma (in the range 0 =< gamma < 3).

//Point-mass potential - constants:
#define M1pointmass 9.565439E+10    //solar masses

//Constants for user-defined profile
#define EPS 1.E-05  //precision of integrator
#define JMAX 5      //maximum number of integration steps pow(3, JMAX)
#define BIGNUMBER 1.0e30

//Mass function variables
#define MAX_AN  10
#define MAX_MN  11

//Plumix
#define SWAP(a,b) { temp = star1[a].mass; star1[a].mass = star1[b].mass; star1[b].mass = temp; temp = star1[a].m0; star1[a].m0 = star1[b].m0; star1[b].m0 = temp; temp = star1[a].kw; star1[a].kw = star1[b].kw; star1[b].kw = temp; temp = star1[a].epoch; star1[a].epoch = star1[b].epoch; star1[b].epoch = temp; temp = star1[a].spin; star1[a].spin = star1[b].spin; star1[b].spin = temp; temp = star1[a].rstar; star1[a].rstar = star1[b].rstar; star1[b].rstar = temp; temp = star1[a].lum; star1[a].lum = star1[b].lum; star1[b].lum = temp; temp = star1[a].epochstar; star1[a].epochstar = star1[b].epochstar; star1[b].epochstar = temp; temp = star1[a].zstar; star1[a].zstar = star1[b].zstar; star1[b].zstar = temp;};
#define _G      6.673e-8
#define _pc     3.0856e18
#define _M_sun  1.989e33
#define _R_sun  6.96265e10
#define _AU     1.49597870e13
#define _GMsun  1.3272597e26
#define RSCALE  0.58904862254809  // = 3pi / 16 = 3pi / (16 * |U_tot|)
#define BUFF_STEP 1024

struct t_star1
{
    double mass, U_tmp, Ui;
    double r[3];
    double m0, kw, epoch, spin, rstar, lum, epochstar, zstar;
};
struct t_star2
{
    double v[3];
    double U, U_sub, E;
    double UT_add, UT;
    double rad, M_sub;
    long   ntry;
};

//functions and COMMON blocks for SSE (Hurley, Pols & Tout, 2000, MNRAS, 315, 543) and BSE (Hurley, Tout & Pols, 2002, MNRAS, 329, 897)
#ifdef SSE
extern void zcnsts_(double *z, double *zpars);
extern void evolv1_(int *kw, double *mass, double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin, double *epoch, double *tms, double *tphys, double *tphysf, double *dtp, double *z, double *zpars, double *vkick);
extern void evolv2_(int *kw, double *mass, double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin, double *epoch, double *tms, double *tphys, double *tphysf, double *dtp, double *z, double *zpars, double *tb, double *ecc, double *vkick);
#endif

struct value1 {
    double neta;
    double bwind;
    double hewind;
    double mxns;

    value1() {
        // Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally)
        neta = 0.5;
        // Binary enhanced mass loss parameter (inactive for single)
        bwind = 0.0;
        // Helium star mass loss factor (1.0 normally)
        hewind = 1.0;
        // Maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1)
        mxns = 3.0;
    }
};

struct points {
    double pts1;
    double pts2;
    double pts3;
    points() {
        // Time-step parameter in evolution phase: MS (0.05)
        pts1 = 0.05;
        // Time-step parameter in evolution phase: GB, CHeB, AGB, HeGB (0.01)
        pts2 = 0.01;
        // Time-step parameter in evolution phase: HG, HeMS (0.02)
        pts3 = 0.02;
    }
};

struct value4 {
    double sigma;
    int bhflag;
    value4() {
        // Kick velocities
        sigma = 190.0;
        // bhflag > 0 allows velocity kick at BH formation
        bhflag = 1;
    }
};

struct value3 {
    int idum;
};

struct flags {
    int ceflag;
    int tflag;
    int ifflag;
    int nsflag;
    int wdflag;
    flags () {
        // ceflag > 0 activates spin-energy correction in common-envelope (0)
        // #defunct#
        ceflag = 0;
        // tflag > 0 activates tidal circularisation (1)
        tflag = 1;
        // ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0)
        ifflag = 0;
        // nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ,572,407(1)
        nsflag = 1;
        // wdflag > 0 uses modified-Mestel cooling for WDs (0)
        wdflag = 1;
    }

};

struct value5 {
    double beta;
    double xi;
    double acc2;
    double epsnov;
    double eddfac;
    double gamma;
    value5() {
        // beta is wind velocity factor: proportional to vwind**2 (1/8)
        beta = 1.0/8.0;
        // xi is the wind accretion efficiency factor (1.0)
        xi = 1.0;
        // acc2 is the Bondi-Hoyle wind accretion factor (3/2)
        acc2 = 3.0/2.0;
        // epsnov is the fraction of accreted matter retained in nova
        // eruption (0.001)
        epsnov = 0.001;
        // eddfac is Eddington limit factor for mass transfer (1.0)
        eddfac = 1.0;
        // gamma is the angular momentum factor for mass lost during Roche (-1.0)
        gamma = -1.0;
    }
};

struct value2 {
    double alpha1;
    double lambda;
    value2() {
        // alpha1 is the common-envelope efficiency parameter (1.0)
        alpha1 = 1.0;
        // lambda is the binding energy factor for common envelope evolution (0.5)
        lambda = 0.5;
    }
};

struct coord {
    double rx, ry, rz;
    double vx, vy, vz;
    coord() {
        rx = 0;
        ry = 0;
        rz = 0;
        vx = 0;
        vy = 0;
        vz = 0;
    }
};

//      star[][0] = mass(t = epoch)                                        *
//      star[][1-3] = {x, y, z}                                            *
//      star[][4-6] = {vx, vy, vz}                                         *
//      star[][7] = mass(t = 0)                                            *
//      star[][8] = kstar, stellar type in (SSE/BSE)                       *
//      star[][9] = epoch1, age within evolutionary phase (SSE/BSE)        *
//      star[][10] = ospin, spin of star (SSE/BSE)                         *
//      star[][11] = rstar, radius of star (SSE/BSE)                       *
//      star[][12] = lstar, luminosity (SSE/BSE)                           *
//      star[][13] = epochstar, age of star (SSE/BSE)                      *
//      star[][14] = Zstar, metallicity of star                            *
struct star_data {
    double mass_epoch;
    double rx;
    double ry;
    double rz;
    double vx;
    double vy;
    double vz;
    double mass;
    double kstar;
    double epoch1;
    double ospin;
    double rstar;
    double lstar;
    double epochstar;
    double zstar;
    star_data () {
        mass_epoch = 0.0;
        rx = 0.0;
        ry = 0.0;
        rz = 0.0;
        vx = 0.0;
        vy = 0.0;
        vz = 0.0;
        mass = 0.0;
        kstar = 0.0;
        epoch1 = 0.0;
        ospin = 0.0;
        rstar = 0.0;
        lstar = 0.0;
        epochstar = 0.0;
        zstar = 0.0;
    }
    struct star_data& operator=(const star_data& other) {
        mass_epoch = other.mass_epoch;
        rx = other.rx;
        ry = other.ry;
        rz = other.rz;
        vx = other.vx;
        vy = other.vy;
        vz = other.vz;
        mass = other.mass;
        kstar = other.kstar;
        epoch1 = other.epoch1;
        ospin = other.ospin;
        rstar = other.rstar;
        lstar = other.lstar;
        epochstar = other.epochstar;
        zstar = other.zstar;
        return *this;
    }
};

struct options {
    // Basic physical parameters
    // Number of stars, Mcl will be set to 0 if specified!
    int N;
    // Total mass of the cluster, only used when N is set to 0, necessary for
    // usage of maximum stellar mass relation of Weidner & Kroupa (2006)
    double Mcl;
    // Density profile;
    // =0 Plummer model,
    // =1 King model (based on king0.f by D.C. Heggie),
    // =2 mass segregated Subr model (Subr et al. 2007),
    // =3 2-dimensional EFF template (Elson, Fall & Freeman 1987) or
    // Nuker template (Lauer et al. 1995),
    // =-1 no density gradient
    int profile;
    // King's W0 paramter [0.3-12.0]
    double W0;
    // Fraction of mass segregation for profile; =0.0 unsegregated,
    // =1.0 completely segregated (take maximally S=0.99 for profile=2)
    double S;
    // Fractal dimension; =3.0 unfractal, =2.6 2/8 fractal, =2.0 4/8 fractal,
    // =1.6 6/8 fractal, (2^D children per parent, following Goodwin &
    // Whitworth 2004)
    double D;
    // Initial virial ratio; =0.5 virial equilibrium, <0.5 collapsing,
    // >0.5 expanding
    double Q;
    // Half-mass radius [pc], ignored if profile = 3, set =-1 for using Marks
    // & Kroupa (2012) Mcl-Rh relation
    double Rh;
    // Power-law slopes of EFF/Nuker templates (outer slope, inner slope,
    // transition); set gamma[1] = 0.0 and gamma[2] = 2.0 for EFF (profile = 3)
    double gamma[3];
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
    char const *output = "test";
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
    // Maximum number of stars & orbits allowed in McLuster
    int NMAX = 1500000;
    // Maximum number of neighbours allowed in NBODY6
    int NNBMAX_NBODY6 = 500;
    // Distance of cluster from sun for artificial CMD with
    // observational errors [pc]
    double Rgal = 10000.0;
    // Solar metallicity
    double Zsun = 0.02;
    // Maximum stellar mass allowed in McLuster [Msun]
    double upper_IMF_limit = 150.0;
    // Input array for external potential parameters
    double extgas[4];
    // DTADJ [N-body units (Myr in Nbody6 custom)], energy-check time step
    double dtadj = 1.0;
    // DELTAT [N-body units (Myr in Nbody6 custom)],
    // output interval, must be multiple of DTADJ
    double dtout = 1.0;
    // DTPLOT [N-body units (Myr in Nbody6 custom)],
    // output of HRdiagnostics, should be multiple of DTOUT,
    // set to zero if output not desired
    double dtplot = 1.0;
    // Gas parameters (only used for Nbody6 input)

    // external Plummer (gas) sphere mass [Msun]
    double extmass = 0.0;
    // external Plummer (gas) sphere scale factor [pc]
    double extrad = 0.0;
    // decay time for gas expulsion (set 0 for no decay) [Myr]
    double extdecay = 0.0;
    // delay time for start of gas expulsion [Myr]
    double extstart = 0.0;
};

#endif

