#ifndef FUNCTIONS_H
#define FUNCTIONS_H
/*************
 * Functions *
 *************/
#include "main.h"

int generate_m1(int *N, double **star, double mlow, double mup, double *M,
    double *mmean, double MMAX, double Mcl, double epoch, double Z, double Rh,
    int remnant);

int generate_m2(int an, double *mlim, double *alpha, double Mcl, double M_tmp,
    double *subcount, int *N, double *mmean, double *M, double **star,
    double MMAX, double epoch, double Z, double Rh, int remnant);

int generate_m3(int *N, double **star, double mlow, double mup, double *M,
    double *mmean, double MMAX, double Mcl);

double subint(double min, double max, double alpha);

double mlow(double mhigh, double alpha, double norma, double delta);

int generate_m4(int *N, double **star, double alpha, double beta, double mu,
    double mlow, double mup, double *M, double *mmean, double MMAX, double Mcl,
    double epoch, double Z, double Rh, int remnant);

double alogam(double x, int *ifault);

double betain(double x, double p, double q, double beta, int *ifault);

double r8_abs(double x);

int generate_plummer(int N, double **star, double rtide, double rvir, double D,
    int symmetry);

int generate_king(int N, double W0, double **star, double *rvir, double *rh,
    double *rking, double D, int symmetry);

double densty(double z);

int odeint(double ystart0, double ystart1, double x1, double x2, double den,
    int *kount, double *xp, double **yp, int M, int KMAX);

int derivs(double x, double *y, double *dydx, double den);

int rkqc(double *y,double *dydx, double *x, double *h, double den,
    double *yscal, double *hdid, double *hnext, double TOL);

int rk4(double x, double *y, double *dydx, double h, double *yout,
    double den);

int generate_subr(int N, double S, double **star, double rtide, double rvir);

void quick(int start, int stop);

void position(int id, double U_sub, double UT_sub, double S);

void find_beta(double *avg, double *beta);

void isorand(double *r);

int cmpmy(double *x1, double *x2);

int cmpmy_reverse(double *x1, double *x2);

double generate_profile (int N, double **star, double Rmax, double Mtot,
    double *p, double *Rh, double D, int symmetry);

double dfridr(double (*func)(double, double*), double x, double h, double *err,
    double *p);

double midexp(double (*func)(double, double*), double aa, double bb, int n,
    double *p);

double midsql(double (*func)(double, double*), double aa, double bb, int n,
    double *p);

double midsqu(double (*func)(double, double*), double aa, double bb, int n,
    double *p);

double midinf(double (*func)(double, double*), double aa, double bb, int n,
    double *p);

double midpnt(double (*func)(double, double*), double a, double b, int n,
    double *p);

double kernel (double x, double *p);

double rhoR (double x, double *p);

double rho(double r, double *p);

double rho_kernel (double x, double *p);

double M(double r, double *p);

double sigma_kernel (double x, double *p);

double sigma(double r, double *p);

double get_gauss(void);

double fractalize(double D, int N, double **star, int radial, int symmetry);

int get_binaries(int nbin, double **star, double M, double rvir, int pairing,
    int *N, int adis, double amin, double amax, double Rh, int eigen, int BSE,
    double epoch, double Z, int remnant, int OBperiods, double msort);

void shellsort(double **array, int N, int k);

void shellsort_reverse(double **array, int N, int k);

void shellsort_1d(double *array, int N);

void shellsort_reverse_1d(double *array, int N);

int order(double **star, int N, double M, double msort, int pairing);

int segregate(double **star, int N, double S);

int energy_order(double **star, int N, int Nstars);

int randomize(double **star, int N);

double rtnewt (double ecc, double ma);

int eigenevolution(double *m1, double *m2, double *ecc, double *abin);

int radial_profile(double **star, int N, double rvir, double M,
    int create_radial_profile, int create_cumulative_profile, int code,
    int *NNBMAX, double *RS0, double *Rh2D, double *Rh3D, int NNBMAX_NBODY6);

int cmd(double **star, int l, double Rgal, double *abvmag, double *vmag,
    double *BV, double *Teff, double *dvmag, double *dBV);

int output0(char *output, int N, int NNBMAX, double RS0, double dtadj,
    double dtout, double tcrit, double rvir, double mmean, int tf,
    int regupdate, int etaupdate, int mloss, int bin, int esc, double M,
    double mlow, double mup, double MMAX, double epoch, double dtplot,
    double Z, int nbin, double Q, double *RG, double *VG, double rtide,
    int gpu, double **star, int sse, int seed, double extmass, double extrad,
    double extdecay, double extstart);

int output1(char *output, int N, double dtadj, double dtout, double tcrit,
    double rvir, double mmean, int tf, int regupdate, int etaupdate, int mloss,
    int bin, int esc, double M, double mlow, double mup, double MMAX,
    double epoch, double Z, int nbin, double Q, double *RG, double *VG,
    double rtide, int gpu, double **star);

int output2(char *output, int N, int NNBMAX, double RS0, double dtadj,
    double dtout, double tcrit, double rvir, double mmean, int tf,
    int regupdate, int etaupdate, int mloss, int bin, int esc, double M,
    double mlow, double mup, double MMAX, double epoch, double dtplot,
    double Z, int nbin, double Q, double *RG, double *VG, double rtide,
    int gpu, double **star, int sse, int seed, double extmass, double extrad,
    double extdecay, double extstart);

int output3(char *output, int N, double rvir, double rh, double mmean,
    double M, double epoch, double Z, double *RG, double *VG, double rtide,
    double **star, double Rgal, double extmass, double extrad);

int output4(char *output, int N, int NNBMAX, double RS0, double dtadj,
    double dtout, double tcrit, double rvir, double mmean, int tf,
    int regupdate, int etaupdate, int mloss, int bin, int esc, double M,
    double mlow, double mup, double MMAX, double epoch, double dtplot,
    double Z, int nbin, double Q, double *RG, double *VG, double rtide,
    int gpu, double **star, int sse, int seed, double extmass, double extrad,
    double extdecay, double extstart);

int output5(char *output, int N, int NNBMAX, double RS0, double dtadj,
    double dtout, double tcrit, double rvir, double mmean, int tf,
    int regupdate, int etaupdate, int mloss, int bin, int esc, double M,
    double mlow, double mup, double MMAX, double epoch, double dtplot,
    double Z, int nbin, double Q, double *RG, double *VG, double rtide,
    int gpu, double **star, int sse, int seed, double extmass, double extrad,
    double extdecay, double extstart);

void info(char *output, int N, double Mcl, int profile, double W0, double S,
    double D, double Q, double Rh, double gamma[], double a, double Rmax,
    double tcrit, int tf, double RG[], double VG[], int mfunc,
    double single_mass, double mlow, double mup, double alpha[], double mlim[],
    double alpha_L3, double beta_L3, double mu_L3, int weidner, int mloss,
    int remnant, double epoch, double Z, int prantzos, int nbin, double fbin,
    int pairing, double msort, int adis, double amin, double amax, int eigen,
    int BSE, double extmass, double extrad, double extdecay, double extstart,
    int code, int seed, double dtadj, double dtout, double dtplot, int gpu,
    int regupdate, int etaupdate, int esc, int units, int match, int symmetry,
    int OBperiods);

void help(double msort);
#endif
