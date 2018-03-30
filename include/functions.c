/*************
 * Functions *
 *************/
#include "functions.h"

struct t_star1 *star1;
struct t_star2 *star2;

#ifndef SSE
void zcnsts_(double *z, double *zpars) {/*DUMMY*/};
void evolv1_(int *kw, double *mass, double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin, double *epoch, double *tms, double *tphys, double *tphysf, double *dtp, double *z, double *zpars, double *vkick) {/*DUMMY*/};
void evolv2_(int *kw, double *mass, double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin, double *epoch, double *tms, double *tphys, double *tphysf, double *dtp, double *z, double *zpars, double *tb, double *ecc, double *vkick) {/*DUMMY*/};
#endif

int generate_m1(struct options *opt, vector<struct star_data> star, double *M,
    double *mmean, double MMAX) {

    //int ty;
    int i;
    double alpha1, alpha2, c1, c2, k1, k2, xx, mth;

    // set up parameters and variables for SSE (Hurley, Pols & Tout 2002)

    // stellar type
    int kw;
    // initial mass
    double mass;
    // actual mass
    double mt;
    // radius
    double r = 0.0;
    // luminosity
    double lum = 0.0;
    // core mass
    double mc = 0.0;
    // core radius
    double rc = 0.0;
    // envelope mass
    double menv = 0.0;
    // envelope radius
    double renv = 0.0;
    // spin
    double ospin;
    // time spent in current evolutionary state
    double epoch1;
    // main-sequence lifetime
    double tms = 0.0;
    // initial age
    double tphys;
    // final age
    double tphysf = opt->epoch;
    // data store value, if dtp>tphys no data will be stored
    double dtp = opt->epoch+1;
    // metallicity
    double z = opt->Z;
    // metallicity parameters
    double zpars[20];
    // kick velocity for compact remnants
    double vkick;
    // opt.escape velocity of cluster
    double vesc;
    // number of ejected compact remnants
    int lostremnants = 0;
    // mass of ejected compact remnants
    double lostremnantsmass = 0.0;

    if (opt->Mcl && opt->remnant && (opt->Rh != -1)) {
        vesc = sqrt(2.0*G*opt->Mcl/opt->Rh);
        printf("Escape velocity of cluster = %.4f km/s\n", vesc);
    } else if ((!opt->remnant) || (opt->Rh == -1)) {
        vesc = 1.0E10;
        printf("Keeping all compact remnants\n");
    } else {
        vesc = sqrt(2.0*0.4*opt->N/opt->Rh);
        printf("Estimated opt.escape velocity of cluster assuming mean stellar "
            "mass of 0.4 Msun = %.4f km/s\n", vesc);
    }

    for (i=0; i<20; i++)
        zpars[i] = 0;

    // get metallicity parameters
    zcnsts_(&z,zpars);

    printf("\nSetting up stellar population with opt.Z = %.4f.\n",opt->Z);

    if (opt->epoch)
        printf("\nEvolving stellar population for %.1f Myr.\n",opt->epoch);

    // set up mass function parameters
    //ty = 2;
    alpha1 = 1.3;
    alpha2 = 2.3;

    c1 = 1.0-alpha1;
    c2 = 1.0-alpha2;

    k1 = 2.0/c1*(pow(0.5,c1)-pow(opt->mlow,c1));
    if (opt->mlow>0.5) {
        k1 = 0;
        k2 = 1.0/c2*(pow(opt->mup,c2)-pow(opt->mlow,c2));
    } else
        k2 = k1 + 1.0/c2*(pow(opt->mup,c2)-pow(0.5,c2));
    if (opt->mup<0.5) {
        k1 = 2.0/c1*(pow(opt->mup,c1)-pow(opt->mlow,c1));
        k2 = k1;
    }

    // determine theoretical mean mass from mass function
    c1 = 2.0-alpha1;
    c2 = 2.0-alpha2;

    if (opt->mlow != opt->mup) {
        int mlow_c1 = pow(opt->mlow, c1);
        int mlow_c2 = pow(opt->mlow, c2);
        int mup_c1  = pow(opt->mup, c1);
        int mup_c2  = pow(opt->mup, c2);

        if (opt->mlow>0.5) {
            mth = (1.0/c2*(mup_c2-mlow_c2))/k2;
        } else if (opt->mup<0.5) {
            mth = (2.0/c1*(mup_c1-mlow_c1))/k2;
        } else
            mth = (2.0/c1*(pow(0.5,c1)-mlow_c1)+1.0/c2*(mup_c2-pow(0.5,c2)))/k2;
    } else {
        mth = opt->mlow;
    }

    if (!opt->N) {
        opt->N = mmax(floor((opt->Mcl-MMAX)/mth), 1);
        if (!opt->epoch)
            printf("Estimated number of necessary stars: %i\n", opt->N);
        opt->N = 1;
    }

    c1 = 1.0-alpha1;
    c2 = 1.0-alpha2;
    *mmean = 0.0;
    *M = 0.0;
    double mostmassive = 0.0;

    for (i=0; i<opt->N; i++) {
        do{
            do {
                xx = drand48();
                if (xx<k1/k2)
                    star[i].mass_epoch = pow(0.5*c1*xx*k2+pow(opt->mlow,c1),1.0/c1);
                else
                    star[i].mass_epoch = pow(c2*(xx*k2-k1)+pow(mmax(0.5,opt->mlow),c2),1.0/c2);
            } while (star[i].mass_epoch > MMAX);

            // evolve star for deltat = opt.epoch with SSE
            // (Hurley, Pols & Tout 2002)
            tphys = 0.0;
            kw = 1;
            // initial mass
            mass = star[i].mass_epoch;
            // actual mass
            mt = mass;
            ospin = 0.0;
            tphysf = opt->epoch;
            epoch1 = 0.0;
            vkick = 0.0;
            //printf("MASS %.2f", mass);
            star[i].mass = mass;

            evolv1_(&kw, &mass, &mt, &r, &lum, &mc, &rc, &menv, &renv, &ospin,
                &epoch1, &tms, &tphys, &tphysf, &dtp, &z, zpars, &vkick);
            // printf("-> %.2f\n", mass);
            // if (vkick) printf("KICK: %.5f\n", vkick);
            lostremnants++;
            lostremnantsmass += mt;
        } while (vkick > vesc);

        lostremnants--;
        lostremnantsmass -= mt;

        star[i].mass_epoch = mt;
        star[i].kstar = kw;
        star[i].epoch1 = epoch1;
        star[i].ospin = ospin;
        star[i].rstar = r;
        star[i].lstar = lum;

        if (star[i].mass_epoch > mostmassive)
            mostmassive = star[i].mass_epoch;

        *M += star[i].mass_epoch;

        if ((i==opt->N-1) && (*M<opt->Mcl))
            opt->N += 1;
    }
    if (lostremnants) {
        printf("Number of ejected compact remnants: %i (%.1f Msun)\n",
            lostremnants, lostremnantsmass);
    }

    printf("Total mass: %g\t(%i stars)\n",*M,opt->N);
    *mmean = *M/ opt->N;
    printf("Most massive star: %g\n",mostmassive);
    if (!opt->epoch) {
        printf("Mean masses theoretical/data: %lf %f\n",mth,*mmean);
    }

    return 0;
}

int generate_m2(struct options *opt, vector<struct star_data> star, double M_tmp,
    double *subcount, double *mmean, double *M, double MMAX) {

    int i, j;

    // set up parameters and variables for SSE (Hurley, Pols & Tout 2002)

    // stellar type
    int kw;
    // initial mass
    double mass;
    // actual mass
    double mt;
    // radius
    double r = 0.0;
    // luminosity
    double lum = 0.0;
    // core mass
    double mc = 0.0;
    // core radius
    double rc = 0.0;
    // envelope mass
    double menv = 0.0;
    // envelope raius
    double renv = 0.0;
    // spin
    double ospin;
    // time of birth?
    double epoch1;
    // main-sequence lifetime
    double tms = 0.0;
    // initial age
    double tphys;
    // final age
    double tphysf = opt->epoch;
    // data store value, if dtp>tphys no data will be stored
    double dtp = opt->epoch+1;
    // metallicity
    double z = opt->Z;
    // metallicity parameters
    double zpars[20];
    // kick velocity for compact remnants
    double vkick;
    // opt.escape velocity of cluster
    double vesc;
    // number of ejected compact remnants
    int lostremnants = 0;
    // mass of ejected compact remnants
    double lostremnantsmass = 0.0;

    if (opt->Mcl && opt->remnant && (opt->Rh != -1)) {
        vesc = sqrt(2.0*G*opt->Mcl/opt->Rh);
        printf("Escape velocity of cluster = %.4f km/s\n", vesc);
    } else if ((!opt->remnant) || (opt->Rh == -1)) {
        vesc = 1.0E10;
        printf("Keeping all compact remnants\n");
    } else {
        vesc = sqrt(2.0*0.4* opt->N/opt->Rh);
        printf("Estimated opt.escape velocity of cluster assuming mean stellar "
            "mass of 0.4 Msun = %.4f km/s\n", vesc);
    }

    for (i=0; i<20; i++)
        zpars[i] = 0;

    // get metallicity parameters
    zcnsts_(&z,zpars);

    printf("\nSetting up stellar population with opt.Z = %.4f.\n",opt->Z);

    if (opt->epoch)
        printf("\nEvolving stellar population for %.1f Myr.\n", opt->epoch);

    double tmp, ml, mup;
    double mostmassive = 0.0;
    *mmean = 0.0;
    *M = 0.0;

    if (!opt->N)
        opt->N = 1;

    for (i = 0; i < opt->an; i++) {
        printf("# <%.2f , %.2f> .. %.2f\n", opt->mlim[i], opt->mlim[i+1],
            opt->alpha[i]);
    }

    for (i = 1; i < opt->an; i++)
        subcount[i] += subcount[i-1];

    for (i = 0; i < opt->N; i++) {
        do {
            do {
                tmp = drand48() * subcount[opt->an-1];
                for (j = 0; (j < opt->an) && (subcount[j] < tmp); j++);
                if (opt->alpha[j] != -1.) {
                    ml = pow(opt->mlim[j], 1. + opt->alpha[j]);
                    mup = pow(opt->mlim[j+1], 1. + opt->alpha[j]);
                } else {
                    ml = log(opt->mlim[j]);
                    mup = log(opt->mlim[j+1]);
                }
                tmp = ml + drand48() * (mup - ml);

                if (opt->alpha[j] != -1.)
                    star[i].mass_epoch = pow(tmp, 1. / (1. + opt->alpha[j]));
                else
                    star[i].mass_epoch = exp(tmp);

            } while (star[i].mass_epoch > MMAX);
            //printf("%8.4f\n", star[i].mass_epoch);

            // evolve star for deltat = opt.epoch with SSE
            // (Hurley, Pols & Tout 2002)
            tphys = 0.0;
            kw = 1;
            // initial mass
            mass = star[i].mass_epoch;
            // actual mass
            mt = mass;
            ospin = 0.0;
            vkick = 0.0;
            tphysf = opt->epoch;
            epoch1 = 0.0;
            //printf("MASS %.2f", mass);
            star[i].mass = mass;
            if (opt->epoch) {
                evolv1_(&kw, &mass, &mt, &r, &lum, &mc, &rc, &menv, &renv,
                    &ospin, &epoch1, &tms, &tphys, &tphysf, &dtp, &z, zpars,
                    &vkick);
            }
            //printf("-> %.2f\n", mass);
            //printf("KICK: %.5f\n", vkick);
            lostremnants++;
            lostremnantsmass += mt;
        } while (vkick > vesc);

        lostremnants--;
        lostremnantsmass -= mt;

        star[i].mass_epoch = mt;
        star[i].kstar = kw;
        star[i].epoch1 = epoch1;
        star[i].ospin = ospin;
        star[i].rstar = r;
        star[i].lstar = lum;

        if (star[i].mass_epoch > mostmassive)
            mostmassive = star[i].mass_epoch;

        *M += star[i].mass_epoch;

        if ((i==opt->N-1) && (*M<opt->Mcl))
            opt->N += 1;
    }

    if (lostremnants) {
        printf("Number of ejected compact remnants: %i (%.1f Msun)\n",
            lostremnants, lostremnantsmass);
    }

    printf("Total mass: %g\t(%i stars)\n",*M,opt->N);

    *mmean = *M/ opt->N;
    printf("Most massive star: %g\n",mostmassive);
    printf("Mean mass: %f\n",*mmean);

    return 0;
}

int generate_m3(struct options *opt, vector<struct star_data> star, double *M,
    double *mmean,  double MMAX) {

    int i;
    double a1, a2, k1, k2, mb, m;

    // set up mass function parameter
    a1 = 1.3;
    a2 = 2.3;
    mb = 0.5;

    k2 = opt->Mcl/(pow(mb, a1-a2)/(2.0-a1)*(pow(mb,2.0-a1) - pow(opt->mlow,2.0-a1)) +
        1.0/(2.0-a2)*(pow(MMAX,2.0-a2)-pow(mb,2.0-a2)));
    k1 = k2*pow(mb, a1-a2);

    printf( "Mcl = %f\tMMAX = %f\tk1 = %f\tk2 = %f\n\n",opt->Mcl,MMAX, k1, k2);

    *mmean = 0.0;
    *M = 0.0;
    opt->N = 0;
    i = 0;

    // set first star to MMAX
    star[i].mass_epoch = MMAX;
    star[i].mass = MMAX;
    star[i].kstar = 0;
    star[i].epoch1 = 0.0;
    star[i].ospin = 0.0;
    star[i].rstar = 0.0;
    star[i].lstar = 0.0;
    star[i].epochstar = 0.0;
    star[i].zstar = 0.0;

    *M += star[i].mass_epoch;

    do{
        i++;

        if (star[i-1].mass_epoch > mb) {
            m = pow(pow(star[i-1].mass_epoch, 2.0-a2) - star[i-1].mass_epoch/k2*(2.0-a2),
                1.0/(2.0-a2));

            if (m < mb) {
                m = pow(pow(mb, 2.0-a1) - (2.0-a1)/pow(mb, a1-a2)*
                    (star[i-1].mass_epoch/k2 + pow(mb, 2.0-a2)/(2.0-a2) - 1.0/(2.0-a2)*
                    (pow(star[i-1].mass_epoch, 2.0-a2))), 1.0/(2.0-a1));
            }
        } else {
            m = pow((a1-2.0)/k1*star[i-1].mass_epoch +
                pow(star[i-1].mass_epoch, 2.0-a1), 1.0/(2.0-a1));
        }

        if (m< opt->mlow) {
            i--;
            break;
        }

        star[i].mass_epoch = m;
        star[i].mass = m;
        star[i].kstar = 0;
        star[i].epoch1 = 0.0;
        star[i].ospin = 0.0;
        star[i].rstar = 0.0;
        star[i].lstar = 0.0;
        star[i].epochstar = 0.0;
        star[i].zstar = 0.0;

        *M += star[i].mass_epoch;
        printf("-------%f\n",*M);

    } while ( *M < opt->Mcl);

    opt->N = i+1;

    printf("Total mass: %g\t(%i stars)\n",*M,opt->N);
    *mmean = *M/ opt->N;
    printf("Most massive star: %g\n",MMAX);

    return 0;
}

double subint(double min, double max, double alpha) {
    if (alpha == 0.)
        return log(max / min);
    else
        return (pow(max, alpha) - pow(min, alpha)) / alpha;
}

double mlow(double mhigh, double alpha, double norma, double delta) {
    if (alpha == 0.)
        return mhigh / exp(delta / norma);
    else
        return pow(pow(mhigh, alpha) - delta * alpha / norma, 1. / alpha);
}

int generate_m4(struct options *opt, vector<struct star_data> star, double *mmean,
    double MMAX, double *M) {

    // L_3 IMF
    // needs functions
    //     double alogam ( double x, int *ifault );
    //     double betain ( double x, double p, double q, double beta,
    //          int *ifault );
    //     double r8_abs ( double x );
    // taken from the Applied Statistics Algorithm 63; ASA063 is a C library
    // which evaluates the incomplete Beta function, by KL Majumder and
    // G Bhattacharjee.

    // set up parameters and variables for SSE (Hurley, Pols & Tout 2002)
    double mu = opt->mup;

    //stellar type
    int kw;
    // initial mass
    double mass;
    // actual mass
    double mt;
    // radius
    double r = 0.0;
    // luminosity
    double lum = 0.0;
    // core mass
    double mc = 0.0;
    // core radius
    double rc = 0.0;
    // envelope mass
    double menv = 0.0;
    // envelope radius
    double renv = 0.0;
    // spin
    double ospin;
    // time spent in current evolutionary state
    double epoch1;
    // main-sequence lifetime
    double tms = 0.0;
    // initial age
    double tphys;
    // final age
    double tphysf = opt->epoch;
    // data store value, if dtp>tphys no data will be stored
    double dtp = opt->epoch+1;
    // metallicity
    double z = opt->Z;
    // metallicity parameters
    double zpars[20];
    // kick velocity for compact remnants
    double vkick;
    // escape velocity of cluster
    double vesc;
    // number of ejected compact remnants
    int lostremnants = 0;
    // mass of ejected compact remnants
    double lostremnantsmass = 0.0;
    double mostmassive = 0.0;

    double G_l;
    double G_u;
    double u;
    double x;
    int i;
    double mth;

    double gamma;
    double mpeak;
    // variables for the mean
    double B_l;
    double B_u;
    double a;
    double b;
    double beta_log;
    int ifault;

    printf("\nL3 Parameters:\n");
    printf("alpha = %1.2f beta = %1.2f mu = %1.2f m_low = %3.3f m_up = "
        "%3.3f\n",opt->alpha_L3,opt->beta_L3,opt->mu_L3,opt->mlow,opt->mup);

    mpeak = opt->mu_L3*pow((opt->beta_L3-1.0),1.0/(opt->alpha_L3-1.0));
    gamma = opt->alpha_L3 + opt->beta_L3*(1.0-opt->alpha_L3);

    printf("'Peak' mass = %3.3f   Effective low mass exponent gamma = "
        "%1.2f\n\n",mpeak,gamma);

    // normalisation
    G_l = pow(1.0 + pow(opt->mlow/opt->mu_L3, 1.0- opt->alpha_L3) , 1.0- opt->beta_L3);
    G_u = pow(1.0 + pow(opt->mup/opt->mu_L3, 1.0- opt->alpha_L3) , 1.0- opt->beta_L3);


    // Calculate the mean mass
    a = (2.0-opt->alpha_L3)/(1.0- opt->alpha_L3);
    b = opt->beta_L3 - a;

    // logarithm of the beta function necessary to calculate
    // the incomplete beta fn
    beta_log = alogam( a, &ifault ) + alogam( b, &ifault ) -
        alogam( a + b, &ifault );

    x = pow(opt->mlow/opt->mu_L3,1.0-opt->alpha_L3);
    x = x/(1.0+x);
    B_l = betain(x, a, b, beta_log, &ifault)*exp(beta_log);

    x = pow(opt->mup/opt->mu_L3,1.0-opt->alpha_L3);
    x = x/(1.0+x);
    B_u = betain(x, a, b, beta_log, &ifault)*exp(beta_log);
    mth = opt->mu_L3*(1.0-opt->beta_L3)*(B_u - B_l)/(G_u - G_l) ;

    if (opt->Mcl && opt->remnant && (opt->Rh != -1)) {
        vesc = sqrt(2.0*G*opt->Mcl/opt->Rh);
        printf("Escape velocity of cluster = %.4f km/s\n", vesc);
    } else if ((!opt->remnant) || (opt->Rh == -1)) {
        vesc = 1.0E10;
        printf("Keeping all compact remnants\n");
    } else {
        vesc = sqrt(2.0*0.4*opt->N/opt->Rh);
        printf("Estimated escape velocity of cluster assuming mean stellar "
            "mass of 0.4 Msun = %.4f km/s\n", vesc);
    }
    for (i=0; i<20; i++)
        zpars[i] = 0;

    zcnsts_(&z,zpars);  //get metallicity parameters

    printf("\nSetting up stellar population with Z = %.4f.\n",opt->Z);

    if (opt->epoch)
        printf("\nEvolving stellar population for %.1f Myr.\n",opt->epoch);

    if (!opt->N) {
        opt->N = mmax(floor((opt->Mcl-opt->mup)/mth), 1);
        if (!opt->epoch)
            printf("Estimated number of necessary stars: %i\n", opt->N);
        opt->N = 1;
    }

    // generate random masses
    for (i=0; i< opt->N; i++) {
        do {
            u = drand48();
            x = u*(G_u-G_l)+G_l;
            x = pow(x,1/(1-opt->beta_L3)) - 1.0;
            star[i].mass_epoch = opt->mu_L3*pow(x,1.0/(1.0-opt->alpha_L3));

            tphys = 0.0;
            kw = 1;
            // initial mass
            mass = star[i].mass_epoch;
            // actual mass
            mt = mass;
            ospin = 0.0;
            tphysf = opt->epoch;
            epoch1 = 0.0;
            vkick = 0.0;
            //printf("MASS %.2f", mass);
            star[i].mass = mass;
            evolv1_(&kw, &mass, &mt, &r, &lum, &mc, &rc, &menv, &renv, &ospin,
                &epoch1, &tms, &tphys, &tphysf, &dtp, &z, zpars, &vkick);
            //printf("-> %.2f\n", mass);
            //if (vkick) printf("KICK: %.5f\n", vkick);
            lostremnants++;
            lostremnantsmass += mt;
        } while (vkick > vesc);

        lostremnants--;
        lostremnantsmass -= mt;

        star[i].mass_epoch = mt;
        star[i].kstar = kw;
        star[i].epoch1 = epoch1;
        star[i].ospin = ospin;
        star[i].rstar = r;
        star[i].lstar = lum;

        if (star[i].mass_epoch > mostmassive)
            mostmassive = star[i].mass_epoch;

        *M += star[i].mass_epoch;

        if ((i==opt->N-1) && (*M<opt->Mcl))
            opt->N += 1;
    }
    if (lostremnants) {
        printf("Number of ejected compact remnants: %i (%.1f Msun)\n",
            lostremnants, lostremnantsmass);
    }

    *mmean = *M/(1.0* opt->N);

    printf("Total mass: %g\t(%i stars)\n",*M,opt->N);
    printf("Most massive star: %g\n",mostmassive);

    if (!opt->epoch)
        printf("Mean masses theoretical/data: %lf %f\n",mth,*mmean);

    return 0;
}

double alogam(double x, int *ifault) {
    /*
    Purpose:
    ALOGAM computes the logarithm of the Gamma function.

    Licensing:
    This code is distributed under the GNU LGPL license.

    Modified:
    20 October 2010

    Author:
    Original FORTRAN77 version by Malcolm Pike, David Hill.
    C version by John Burkardt.

    Reference:
    Malcolm Pike, David Hill,
    Algorithm 291:
    Logarithm of Gamma Function,
    Communications of the ACM,
    Volume 9, Number 9, September 1966, page 684.

    Parameters:
    Input, double X, the argument of the Gamma function.
    X should be greater than 0.

    Output, int *IFAULT, error flag.
    0, no error.
    1, X <= 0.

    Output, double ALOGAM, the logarithm of the Gamma
    function of X.
    */

    double f;
    double value;
    double y;
    double z;

    if (x <= 0.0)
    {
        *ifault = 1;
        value = 0.0;

        return value;
    }

    *ifault = 0;
    y = x;

    if ( x < 7.0 )
    {
        f = 1.0;
        z = y;

        while ( z < 7.0 )
        {
            f = f * z;
            z = z + 1.0;
        }
        y = z;
        f = - log ( f );
    }
    else
    {
        f = 0.0;
    }

    z = 1.0 / y / y;

    value = f + ( y - 0.5 ) * log ( y ) - y + 0.918938533204673 +
        (((- 0.000595238095238 * z + 0.000793650793651 ) * z -
        0.002777777777778 ) * z + 0.083333333333333 ) / y;

    return value;
}

double betain(double x, double p, double q, double beta, int *ifault) {
    /*
    Purpose:
    BETAIN computes the incomplete Beta function ratio.

    Licensing:
    This code is distributed under the GNU LGPL license.

    Modified:
    31 October 2010

    Author:
    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
    C version by John Burkardt.

    Reference:
    KL Majumder, GP Bhattacharjee,
    Algorithm AS 63:
    The incomplete Beta Integral,
    Applied Statistics,
    Volume 22, Number 3, 1973, pages 409-411.

    Parameters:
    Input, double X, the argument, between 0 and 1.

    Input, double P, Q, the parameters, which
    must be positive.

    Input, double BETA, the logarithm of the complete
    beta function.

    Output, int *IFAULT, error flag.
    0, no error.
    nonzero, an error occurred.

    Output, double BETAIN, the value of the incomplete
    Beta function ratio.
    */

    double acu = 0.1E-14;
    double ai;
    double cx;
    int indx;
    int ns;
    double pp;
    double psq;
    double qq;
    double rx;
    double temp;
    double term;
    double value;
    double xx;

    value = x;
    *ifault = 0;
    /*
     Check the input arguments.
     */
    if ( p <= 0.0 || q <= 0.0 )
    {
        *ifault = 1;
        return value;
    }

    if ( x < 0.0 || 1.0 < x )
    {
        *ifault = 2;
        return value;
    }
    /*
     Special cases.
     */
    if ( x == 0.0 || x == 1.0 )
    {
        return value;
    }
    /*
     Change tail if necessary and determine S.
     */
    psq = p + q;
    cx = 1.0 - x;

    if ( p < psq * x )
    {
        xx = cx;
        cx = x;
        pp = q;
        qq = p;
        indx = 1;
    }
    else
    {
        xx = x;
        pp = p;
        qq = q;
        indx = 0;
    }

    term = 1.0;
    ai = 1.0;
    value = 1.0;
    ns = ( int ) ( qq + cx * psq );
    /*
     Use the Soper reduction formula.
     */
    rx = xx / cx;
    temp = qq - ai;
    if ( ns == 0 )
    {
        rx = xx;
    }

    for ( ; ; )
    {
        term = term * temp * rx / ( pp + ai );
        value = value + term;;
        temp = r8_abs ( term );

        if ( temp <= acu && temp <= acu * value )
        {
            value = value * exp ( pp * log ( xx )
                                 + ( qq - 1.0 ) * log ( cx ) - beta ) / pp;

            if ( indx )
            {
                value = 1.0 - value;
            }
            break;
        }

        ai = ai + 1.0;
        ns = ns - 1;

        if ( 0 <= ns )
        {
            temp = qq - ai;
            if ( ns == 0 )
            {
                rx = xx;
            }
        }
        else
        {
            temp = psq;
            psq = psq + 1.0;
        }
    }

    return value;
}

double r8_abs(double x) {
    /*
     Purpose:

     R8_ABS returns the absolute value of an R8.

     Licensing:

     This code is distributed under the GNU LGPL license.

     Modified:

     07 May 2006

     Author:

     John Burkardt

     Parameters:

     Input, double X, the quantity whose absolute value is desired.

     Output, double R8_ABS, the absolute value of X.
     */

    double value;

    if ( 0.0 <= x )
    {
        value = + x;
    }
    else
    {
        value = - x;
    }
    return value;
}

int generate_plummer(int N, vector<struct star_data> star, double rtide, double rvir, double D,
    int symmetry){

    int i, h;
    double a[9], ri, sx, rcut;
    double r_norm, v_norm;

    // Scale length
    sx = 1.5*TWOPI/16.0;

    printf("Setting cut-off radius of Plummer sphere to approximate tidal "
        "radius\n");

    // cut-off radius for Plummer sphere = tidal radius in scaled length
    rcut = rtide/(sx*rvir);

    printf("\nGenerating Orbits:\n");

    if (D>=3.0) {
        for (i=0;i<N;i++) {

            if ((i/1000)*1000 == i)
                printf("Generating orbit #%i\n", i);

            // Positions
            do {
                do {
                    a[1] = drand48();
                } while (a[1]<1.0E-10);

                ri = 1.0/sqrt(pow(a[1],-2.0/3.0) - 1.0);

                a[2] = drand48();
                a[3] = drand48();

                star[i].rz = (1.0 - 2.0*a[2])*ri;
                star[i].rx = sqrt(ri*ri-pow(star[i].rz,2))*
                    cos(TWOPI*a[3]);
                star[i].ry = sqrt(ri*ri-pow(star[i].rz,2))*
                    sin(TWOPI*a[3]);

            // reject particles beyond tidal radius
            } while (sqrt(pow(star[i].rx,2)+pow(star[i].ry,2)+
                pow(star[i].rz,2))>rcut);

        // velocities
        do {
            a[4] = drand48();
            a[5] = drand48();
            a[6] = pow(a[4],2)*pow(1.0 - pow(a[4],2),3.5);
        } while (0.1*a[5]>a[6]);

        a[8] = a[4]*sqrt(2.0)/pow(1.0 + ri*ri,0.25);
        a[6] = drand48();
        a[7] = drand48();

        star[i].vz = (1.0 - 2.0*a[6])*a[8];
        star[i].vx = sqrt(a[8]*a[8] - pow(star[i].vz,2))*cos(TWOPI*a[7]);
        star[i].vy = sqrt(a[8]*a[8] - pow(star[i].vz,2))*sin(TWOPI*a[7]);
        }
    } else {
        double xcut;
        xcut = 0.000;
        while (1.0/sqrt(pow(xcut,-2.0/3.0) - 1.0)<=rcut) xcut+=0.00001;

        fractalize(D, N, star, 1, symmetry);
        for (i=0;i<N;i++) {

            if ((i/1000)*1000 == i)
                printf("Generating orbit #%i\n", i);

            ri = sqrt(pow(star[i].rx,2)+pow(star[i].ry,2)+
                pow(star[i].rz,2))*xcut;
            r_norm = 1.0/sqrt(pow(ri,-2.0/3.0) - 1.0);
            v_norm = sqrt(2.0)/pow(1.0 + ri*ri,0.25);

            star[i].rx *= r_norm/ri;
            star[i].ry *= r_norm/ri;
            star[i].rz *= r_norm/ri;

            star[i].vx *= v_norm;
            star[i].vy *= v_norm;
            star[i].vz *= v_norm;
        }
    }

    return 0;
}

int generate_king(int N, double W0, vector<struct star_data> star, double *rvir, double *rh,
    double *rking, double D, int symmetry) {

    // ODE variables

    // Number of interpolation points
    int M = 10001;
    // Maximum number of output steps of the integrator
    int KMAX = 10000;

    int i,j,k;
    double h;
    double den;
    double xstart, ystart0, ystart1;
    double x1, x2;
    double xp[KMAX], x[M];
    int kount = 0;
    double yking[M][2], mass[M];
    double rmin = 0.2;

    double **yp;
    yp = (double **)calloc(KMAX,sizeof(double *));

    for (i=0;i<KMAX;i++){
        yp[i] = (double *)calloc(2,sizeof(double));

        if (yp[i] == NULL) {
            printf("\nMemory allocation failed!\n");
            return 0;
        }
    }

    double pot = 0.0;
    double totmas;
    double zh;

    double ve, cg1;
    double fmass, r2;
    double costh, sinth, phi, r, xstar, ystar, zstar;
    double w1, w=0, wj1, wj;
    double vmax, speed, fstar, ustar, vstar, wstar;
    double ri;

    struct coord coord[N];

    if (W0>12.0) {
        printf("W0 too large\n");
        return 0;
    } else if (W0 < 0.2) {
        printf("W0 too small\n");
        return 0;
    }

    /***************************
     * INTERPOLATE KING VALUES *
     ***************************/

    // step size for interpolation
    h = pow(W0-rmin,2)/(1.0*M-1.0);
    // central density
    den = densty(W0);

    //x is King's W, y[0] is z**2, y[1] is 2*z*dz/dW, where z is scaled radius
    //so x(y[0]) is energy as function of radius^2 and x(y[1]) is the derivative

    xstart = 0.000001;
    ystart0 = 2.0*xstart/3.0;
    ystart1 = -2.0/3.0;

    x1 = W0 - xstart;
    x2 = 0.0;

    // integrate Poisson's eqn
    printf("Integrating Poisson's equation\n");

    odeint(ystart0, ystart1, x1, x2, den, &kount, xp, yp, M, KMAX);

    printf("No of integration steps = %i\n",kount);

    // interpolate yking
    for (k=0;k<M;k++) {
        x[k] = W0-rmin-sqrt(h*k);

        if (x[k] > xp[0]) {
            yking[k][0] = (W0-x[k])*yp[0][0]/(W0 - xp[0]);
            yking[k][1] = yp[0][1];
            printf("+");
        } else {
            i = 0;

            do {
                if ((x[k]-xp[i])*(x[k]-xp[i+1]) <= 0.0) {
                    yking[k][0] = yp[i][0] + (yp[i+1][0]-yp[i][0])*
                        (x[k]-xp[i])/(xp[i+1]-xp[i]);
                    yking[k][1] = yp[i][1] + (yp[i+1][1]-yp[i][1])*
                        (x[k]-xp[i])/(xp[i+1]-xp[i]);
                    goto jump;
                } else {
                    i++;
                }
            } while (i<kount);
        }

        jump:

        if (i >= kount) {
            yking[k][0] = yp[kount][0];
            yking[k][1] = yp[kount][1];
        }

        // get mass as a function of radius
        mass[k] = -2.0*pow(yking[k][0],1.5)/yking[k][1];


        //calculate total potential energy
        if (k == 0) {
            pot = -0.5*pow(mass[k],2)/sqrt(yking[k][0]);
        } else {
            pot -=  0.5*(pow(mass[k],2) - pow(mass[k-1],2))/(0.5
                *(sqrt(yking[k][0]) + sqrt(yking[k-1][0])));
        }
    }

    /******************************
     * DETERMINE BASIC QUANTITIES *
     ******************************/

    // Total mass
    totmas = -2.0*pow(yking[M-2][0],1.5)/yking[M-2][1];

    // Half-mass radius
    k=0;

    do {
        k++;
    } while (mass[k] < 0.5*totmas);

    zh = yking[k][0] - (yking[k][0]-yking[k-1][0])*
        (mass[k]-0.5*totmas)/(mass[k]-mass[k-1]);

    *rh = sqrt(zh);

    // Virial radius and King radius
    *rvir = -pow(totmas,2)/(2.0*pot);
    *rking = sqrt(yking[M-2][0]);

    //Central velocity dispersion**2 (3-dimensional)
    ve = sqrt(2.0*W0);
    cg1 = (-pow(ve,3)*exp(-pow(ve,2)/2.0) - 3.0*ve*exp(-pow(ve,2)/2.0) +
        3.0/2.0*sqrt(2.0*PI)*erf(ve*sqrt(2.0)/2.0) - pow(ve,5)*
        exp(-pow(ve,2)/2.0)/5.0)/(-ve*exp(-pow(ve,2)/2.0) +
        sqrt(2.0*PI)*erf(ve*sqrt(2.0)/2.0)/2.0 - pow(ve,3)*
        exp(-pow(ve,2)/2.0)/3.0);

    printf("\nTheoretical values:\n");
    printf("Total mass (King units) = %g\n", totmas);
    printf("Viriral radius (King units) = %g\n", *rvir);
    printf("Half-mass radius (King units) = %g\t(Nbody units) = %g\n",
        *rh, *rh/ *rvir);
    printf("Edge radius (King units) = %g\t(Nbody units) = %g\n",
        *rking, *rking/ *rvir);
    printf("Concentration = %g\n", log10(*rking));
    printf("Core radius (King units) = %g\t(Nbody units) = %g\n",
        1.0, 1.0/ *rvir);
    printf("3d velocity dispersion**2: %g (central)\t %g (global)\n",
        cg1, -pot/totmas);

    /***************************
     * GENERATE STAR POSITIONS *
     ***************************/

    printf("\nGenerating Stars:\n");

    if (D>=3.0) {
        for (i=0;i<N;i++) {

            if ((i/1000)*1000 == i) printf("Generating orbits #%i\n", i);

            fmass = drand48()*mass[M-1];

            if (fmass < mass[0]) {
                r2 = pow(fmass/mass[0],2.0/3.0)*yking[0][0];
            } else {
                j = 0;
                do {
                    j++;
                    if (j>M) printf("WARNING: failing iteration\n");
                } while (mass[j] <= fmass);

                r2 = yking[j-1][0] + (fmass-mass[j-1])*(yking[j][0] -
                    yking[j-1][0])/(mass[j]-mass[j-1]);
            }

            r = sqrt(r2);
            costh = 2.0*drand48()-1.0;
            sinth = sqrt(1.0-pow(costh,2));
            phi = 2.0*PI*drand48();
            xstar = r*sinth*cos(phi);
            ystar = r*sinth*sin(phi);
            zstar = r*costh;

            /****************
             * CHOOSE SPEED *
             ****************/

            if (r < *rking) {
                if (r < sqrt(yking[0][0])) {
                    w1 = x[0];
                    w = W0 - (r2/yking[0][0])*(W0 - w1);
                } else {
                    j = 0;
                    do {
                        j++;
                    } while (r > sqrt(yking[j][0]));
                    wj1 = x[j-1];
                    wj = x[j];
                    w = wj1 + (r2-yking[j-1][0])*(wj-wj1)/(yking[j][0]-yking[j-1][0]);
                }
            } else {
                printf("radius too big\n");
            }

            vmax = sqrt(2.0*w);
            do {
                speed = vmax*drand48();
                fstar = pow(speed,2)*(exp(-0.5*pow(speed,2))-exp(-1.0*w));
            } while (fstar < 2.0*drand48()/exp(1.0));

            costh = 2.0*drand48()-1.0;
            phi = 2.0*PI*drand48();
            sinth = sqrt(1.0-pow(costh,2));
            ustar = speed*sinth*cos(phi);
            vstar = speed*sinth*sin(phi);
            wstar = speed*costh;

            coord[i].rx = xstar;
            coord[i].ry = ystar;
            coord[i].rz = zstar;
            coord[i].vx = ustar;
            coord[i].vy = vstar;
            coord[i].vz = wstar;
        }

        for (i=0;i<N;i++) {
            star[i].rx = 1.0*coord[i].rx;
            star[i].ry = 1.0*coord[i].ry;
            star[i].rz = 1.0*coord[i].rz;

            star[i].vx = 1.0*coord[i].vx;
            star[i].vy = 1.0*coord[i].vy;
            star[i].vz = 1.0*coord[i].vz;
        }

    } else {

        fractalize(D, N, star, 1, symmetry);

        for (i=0;i<N;i++) {

            if ((i/1000)*1000 == i)
                printf("Generating orbits #%i\n", i);

            ri =  sqrt(pow(star[i].rx,2)+pow(star[i].ry,2)+pow(star[i].rz,2));
            fmass = ri*mass[M-1];

            if (fmass < mass[0]) {
                r2 = pow(fmass/mass[0],2.0/3.0)*yking[0][0];
            } else {
                j = 0;
                do {
                    j++;
                    if (j>M) printf("WARNING: failing iteration\n");
                } while (mass[j] <= fmass);

                r2 = yking[j-1][0] + (fmass-mass[j-1])*(yking[j][0] -
                    yking[j-1][0])/(mass[j]-mass[j-1]);
            }

            r = sqrt(r2);

            /****************
             * CHOOSE SPEED *
             ****************/

            if (r < *rking) {
                if (r < sqrt(yking[0][0])) {
                    w1 = x[0];
                    w = W0 - (r2/yking[0][0])*(W0 - w1);
                } else {
                    j = 0;
                    do {
                        j++;
                    } while (r > sqrt(yking[j][0]));
                    wj1 = x[j-1];
                    wj = x[j];
                    w = wj1 + (r2-yking[j-1][0])*(wj-wj1)/(yking[j][0]-
                        yking[j-1][0]);
                }
            } else {
                printf("radius too big\n");
            }

            vmax = sqrt(2.0*w);
            do {
                speed = vmax*drand48();
                fstar = pow(speed,2)*(exp(-0.5*pow(speed,2))-exp(-1.0*w));
            } while (fstar < 2.0*drand48()/exp(1.0));

            star[i].rx *= r/ri;
            star[i].ry *= r/ri;
            star[i].rz *= r/ri;

            star[i].vx *= speed;
            star[i].vy *= speed;
            star[i].vz *= speed;
        }

    }

    for (j=0;j<KMAX;j++)
        free (yp[j]);

    free(yp);

    return 0;
}

double densty(double z){
    double den = -sqrt(z)*(z+1.5)+0.75*sqrt(PI)*exp(z)*erf(sqrt(z));
    return den;
}

int odeint(double ystart0, double ystart1, double x1, double x2, double den,
    int *kount, double *xp, double **yp, int M, int KMAX) {

    // Minimum step size
    double HMIN = 0.0;
    // Size of first step
    double H1 = 0.0001;
    // Maximum number of steps for integration
    int MAXSTP = 100000;
    // To avoid certain numbers get zero
    double TINY = 1.0E-30;
    // Output step size for integration
    double DXSAV = 0.0001;
    // Tolerance of integration
    double TOL = 1.0E-12;

    double x;
    double h;
    int i,j;
    double y[2];
    double hdid, hnext;
    double xsav;
    double dydx[2], yscal[2];

    // King's W parameter
    x = x1;
    if (x2-x1 >= 0.0) {
        h = sqrt(pow(H1,2));
    } else {
        h = - sqrt(pow(H1,2));
    }  //step size

    // z^2
    y[0] = ystart0;
    // 2*z*dz/dW  where z is scaled radius
    y[1] = ystart1;

    xsav = x-DXSAV*2.0;

    // limit integration to MAXSTP steps
    for (i=0;i<MAXSTP;i++) {

        // find derivative
        derivs(x,y,dydx,den);

        for (j=0;j<2;j++) {
            // advance y1 and y2
            yscal[j] = sqrt(pow(y[j],2))+sqrt(h*dydx[j])+TINY;
        }

        if (sqrt(pow(x-xsav,2)) > sqrt(pow(DXSAV,2))) {

            if (*kount < KMAX-1) {
                xp[*kount] = x;

                for (j=0;j<2;j++) {
                    yp[*kount][j] = y[j];
                }

                *kount = *kount + 1;
                xsav = x;
            }
        }  //store x, y1 and y2 if the difference in x is smaller as the
           // desired output step size DXSAV

        if (((x+h-x2)*(x+h-x1)) > 0.0) h = x2-x;

        // do a Runge-Kutta step
        rkqc(y,dydx,&x,&h,den,yscal, &hdid, &hnext, TOL);

        if ((x-x2)*(x2-x1) >= 0.0) {
            ystart0 = y[0];
            ystart1 = y[1];

            xp[*kount] = x;
            for (j=0;j<2;j++) {
                yp[*kount][j] = y[j];
            }
            return 0;
            *kount = *kount +1;
        }

        if (sqrt(pow(hnext,2)) < HMIN) {
            printf("Stepsize smaller than minimum.\n");
            return 0;
        }

        h = hnext;
    }

    return 0;
}

int derivs(double x, double *y, double *dydx, double den){

    double rhox;

    if (x >= 0.0) {
        rhox =-sqrt(x)*(x+1.5)+0.75*sqrt(PI)*exp(x)*erf(sqrt(x));
    } else {
        rhox = 0.0;
    }

    dydx[0]= y[1];
    dydx[1] = 0.25*pow(y[1],2)*(6.0+9.0*y[1]*rhox/den)/y[0];

    return 0;
}

int rkqc(double *y,double *dydx, double *x, double *h, double den,
    double *yscal, double *hdid, double *hnext, double TOL) {

    double safety = 0.9;
    double fcor = 0.0666666667;
    double errcon = 6.0E-4;
    double pgrow = -0.20;
    double pshrnk = -0.25;
    double xsav;
    int i;
    double ysav[2],  dysav[2], ytemp[2];
    double errmax;
    double hh;

    xsav = *x;

    for (i=0;i<2;i++) {
        ysav[i] = y[i];
        dysav[i] = dydx[i];
    }

    do {
        hh = 0.5**h;
        rk4(xsav, ysav, dysav, hh, ytemp, den);

        *x = xsav + hh;
        // find derivative
        derivs(*x,ytemp,dydx,den);
        rk4(*x, ytemp, dydx, hh, y, den);

        *x = xsav + *h;
        if (*x  == xsav) {
            printf("ERROR: Stepsize not significant in RKQC.\n");
            return 0;
        }
        rk4(xsav, ysav, dysav, *h, ytemp, den);

        errmax = 0.0;
        for (i=0;i<2;i++) {
            ytemp[i] = y[i] - ytemp[i];
            errmax = mmax(errmax, sqrt(pow(ytemp[i]/yscal[i],2)));
        }
        errmax /= TOL;
        // if integration error is too large, decrease h
        if (errmax > 1.0)
            *h = safety**h*(pow(errmax,pshrnk));

    } while (errmax > 1.0);


    *hdid = *h;
    if (errmax > errcon) {
        // increase step size for next step
        *hnext = safety**h*(pow(errmax,pgrow));
    } else {
        // if integration error is very small increase step size significantly
        // for next step
        *hnext = 4.0**h;
    }

    for (i=0;i<2;i++) {
        y[i] += ytemp[i]*fcor;
    }

    return 0;

}

int rk4(double x, double *y, double *dydx, double h, double *yout,
    double den) {

    double hh, h6, xh;
    double yt[2], dyt[2], dym[2];
    int i;

    hh=h*0.5;
    h6=h/6.0;
    xh=x+hh;

    for (i=0;i<2;i++) {
        yt[i] = y[i] + hh*dydx[i];
    }

    // find derivative
    derivs(xh,yt,dyt,den);
    for (i=0;i<2;i++) {
        yt[i] = y[i] + hh*dyt[i];
    }

    // find derivative
    derivs(xh,yt,dym,den);
    for (i=0;i<2;i++) {
        yt[i] = y[i] + h*dym[i];
        dym[i] += dyt[i];
    }

    // find derivative
    derivs(x+h,yt,dyt,den);
    for (i=0;i<2;i++) {
        yout[i] = y[i] + h6*(dydx[i]+dyt[i]+2.0*dym[i]);
    }

    return 0;
}

int generate_subr(int N, double S, vector<struct star_data> star, double rtide, double rvir) {

    long   i, j;
    double r1, r2, tmp, rcut, rtemp;
    double U_tot, U_tot1, K_tot;
    double UT_sub, M_sub;
    double maxim, beta;
    int    N_star;
    double M_tot;
    star1 = NULL;
    star2 = NULL;

    if (S<0.5) {
        printf("Setting cut-off radius of Subr model to the approximate "
            "tidal radius\n");

        // cut-off radius for Plummer sphere = tidal radius
        rcut = rtide/rvir;
    } else {
        printf("Attention! Be aware that for about S>0.5 the Subr mass "
            "segregation may generate a virially hot system\n"
            "Setting cut-off radius of Subr model to twice the approximate "
            "tidal radius\n");
        // cut-off radius for Plummer sphere = tidal radius
        rcut = 2.0*rtide/rvir;
    }

    M_tot = 0.;
    N_star = N;
    star1 = (struct t_star1 *)realloc(star1, N_star * sizeof(struct t_star1));
    star2 = (struct t_star2 *)realloc(star2, N_star * sizeof(struct t_star2));

    for (i=0; i<N_star;i++) {
        star1[i].mass = star[i].mass_epoch;
        star1[i].m0 = star[i].mass;
        star1[i].kw = star[i].kstar;
        star1[i].epoch = star[i].epoch1;
        star1[i].spin = star[i].ospin;
        star1[i].rstar = star[i].rstar;
        star1[i].lum = star[i].lstar;
        star1[i].epochstar = star[i].epochstar;
        star1[i].zstar = star[i].zstar;
        M_tot += star1[i].mass;
    };

    quick(0, N_star-1);

    M_sub = U_tot1 = 0.;
    for (i = 0; i < N_star; i++) {

        star2[i].U = star1[i].U_tmp = star2[i].UT_add = 0.;
        star1[i].mass /= M_tot;
        M_sub += star1[i].mass;
        star2[i].M_sub = M_sub;

        if (S > 0.)
            star1[i].Ui = star1[i].mass * pow(M_sub, -S);
        else
            star1[i].Ui = star1[i].mass;

        for (j = 0; j < i; j++) {
            tmp = star1[i].Ui * star1[j].Ui;
            star2[i].UT_add += tmp;
            star1[i].U_tmp += tmp;
            star1[j].U_tmp += tmp;
        };
        U_tot1 += star2[i].UT_add;
    };

    UT_sub = U_tot = K_tot = 0.;

    for (i = 0; i < N_star; i++) {

        if ((i/1000)*1000 == i)
            printf("Generating orbit #%li\n", i);

        star2[i].UT_add *= 0.5 / U_tot1;
        star2[i].UT = star1[i].U_tmp * 0.5 / U_tot1;
        UT_sub += star2[i].UT_add;

        do {
            position(i, U_tot, UT_sub, S);
            rtemp = sqrt(pow(star1[i].r[0],2)+pow(star1[i].r[1],2)+
                pow(star1[i].r[2],2));
        } while (rtemp>rcut);
        U_tot += star1[i].mass * star2[i].U_sub;
    };

    for (i = 0; i < N_star; i++) {
        maxim = 0.25 * star1[i].mass / star2[i].UT;
        find_beta(&maxim, &beta);
        if (beta > 0.)
        {
            do
            {
                r1 = maxim * drand48();
                r2 = drand48();
                r2 *= r2;
                tmp = 1. - r2;
            }
            while (r1 > r2 * exp(beta * log(tmp)));
        }
        else
        {
            do
            {
                r1 = maxim * drand48();
                r2 = sin(M_PI_2 * drand48());
                r2 *= r2;
                tmp = sqrt(1. - r2);
            }
            while (r1 > r2 * exp((2.* beta + 1.) * log(tmp)));
        };

        star2[i].E = r2 * (star2[i].U + star2[i].U_sub);

        tmp = sqrt(2. * star2[i].E);
        isorand(star2[i].v);
        star2[i].v[0] *= tmp;
        star2[i].v[1] *= tmp;
        star2[i].v[2] *= tmp;
        K_tot += star1[i].mass * star2[i].E;
    };

    for (i = 0; i < N_star; i++){
        //printf("%.12f  %18.14f  %18.14f  %18.14f  %16.12f  %16.12f  "
        //  "%16.12f\n",star1[i].mass, star1[i].r[0], star1[i].r[1],
        //  star1[i].r[2], star2[i].v[0], star2[i].v[1], star2[i].v[2]);
        star[i].mass_epoch = star1[i].mass;
        star[i].rx = star1[i].r[0];
        star[i].ry = star1[i].r[1];
        star[i].rz = star1[i].r[2];
        star[i].vx = star2[i].v[0];
        star[i].vy = star2[i].v[1];
        star[i].vz = star2[i].v[2];
        star[i].mass = star1[i].m0;
        star[i].kstar = star1[i].kw;
        star[i].epoch1 = star1[i].epoch;
        star[i].ospin = star1[i].spin;
        star[i].rstar = star1[i].rstar;
        star[i].lstar = star1[i].lum;
        star[i].epochstar = star1[i].epochstar;
        star[i].zstar = star1[i].zstar;

    };

    return 0;
}

void quick(int start, int stop) {
    int i, j;
    double temp, median;

    median = 0.5 * (star1[start].mass + star1[stop].mass);

    i = start;
    j = stop;
    while (i <= j)
    {
        while ((i <= stop) && (star1[i].mass > median))
            i++;

        while ((j >= start) && (star1[j].mass < median))
            j--;

        if (j > i)
        {
            SWAP(i,j);
            i++;
            j--;
        }
        else if (i == j)
        {
            i++;
            j--;
        };
    };

    if (start + 1 < j)
        quick(start, j);
    else if ((start < j) && (star1[start].mass < star1[j].mass))
        SWAP(start, j);

    if (i + 1 < stop)
        quick(i, stop);
    else if ((i < stop) && (star1[i].mass < star1[stop].mass))
        SWAP(i, stop)
}

void position(int id, double U_sub, double UT_sub, double S){

    long   i;
    double rx, ry, rz, r1;
    double rfac;

    rfac = pow(star2[id].M_sub, 2. * S) / (1. - S);

    star2[id].ntry = 0;

    do
    {
        star2[id].U_sub = 0.;
        star2[id].ntry++;

        r1 = drand48();
        r1 = exp(-2. * log(r1) / 3.);
        r1 = RSCALE * sqrt(1. / (r1 - 1.));

        r1 *= rfac;

        isorand(star1[id].r);
        star1[id].r[0] *= r1;
        star1[id].r[1] *= r1;
        star1[id].r[2] *= r1;
        star2[id].rad = r1;

        for (i = 0; i < id; i++)
        {
            rx = star1[i].r[0] - star1[id].r[0];
            ry = star1[i].r[1] - star1[id].r[1];
            rz = star1[i].r[2] - star1[id].r[2];
            r1 = sqrt(rx*rx + ry*ry + rz*rz);
            star1[i].U_tmp = star1[id].mass / r1;
            star2[id].U_sub += star1[i].mass / r1;
        };
        r1 = U_sub + star1[id].mass * star2[id].U_sub;
    }
    while (fabs(r1 - UT_sub) > 0.1 * (r1 + UT_sub) / sqrt(1. + id));

    for (i = 0; i < id; i++) star2[i].U += star1[i].U_tmp;
}

void find_beta(double *avg, double *beta){

    int i;
    double a1, a2, I1, I2, J1, J2;
    double alo, ah, Ilo, Ih, Jlo, Jh;
    double tmpavg;

    if (*avg > 0.75)
    {
        *beta = -0.5;
        *avg = 1.;
        return;
    };

    Ilo = I1 = 3. * M_PI / 16.;
    Ih  = I2 = 0.2;
    Jlo = J1 = M_PI / 4;
    Jh  = J2 = 1. / 3.;
    alo = a1 = -0.5;
    ah  = a2 = 0.;

    tmpavg = I2 / J2;
    i = 0;
    while (tmpavg > *avg)
    {
        if (!(i & 1))
        {
            a1 += 1.;
            I1 /= 1. + 5. / (2. * a1);
            J1 /= 1. + 3. / (2. * a1);
            Ilo = I2; Ih = I1;
            Jlo = J2; Jh = J1;
            alo = a2; ah = a1;
        }
        else
        {
            a2 += 1.;
            I2 /= 1. + 5. / (2. * a2);
            J2 /= 1. + 3. / (2. * a2);
            Ilo = I1; Ih = I2;
            Jlo = J1; Jh = J2;
            alo = a1; ah = a2;
        };
        tmpavg = Ih / Jh;
        i++;
    };

    *beta = alo + (ah - alo) * (*avg - Ilo / Jlo) / (tmpavg - Ilo / Jlo);

    if (*beta > 0.)
        *avg = exp(*beta * log(*beta / (1. + *beta))) / (1. + *beta);
    else
        *avg = 2. * exp((2.* *beta + 1.) * log((2.* *beta + 1.) / (2.* *beta +
            3.))) / (2.* *beta + 3.);
}

void isorand(double *r) {

    double rad;

    do
    {
        r[0] = 2. * drand48() - 1.;
        r[1] = 2. * drand48() - 1.;
        r[2] = 2. * drand48() - 1.;
        rad = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    }
    while ((rad > 1.) || (rad < 0.0001));

    rad = sqrt(rad);
    r[0] /= rad;
    r[1] /= rad;
    r[2] /= rad;
}

// smallest to largest
int cmpmy(double *x1, double *x2) {
    if(*x1<*x2) return -1;
    return 1;
}

// largest to smallest
int cmpmy_reverse(double *x1, double *x2) {
    if(*x1>*x2) return -1;
    return 1;
}

double generate_profile (int N, vector<struct star_data> star, double Rmax, double Mtot,
    double *p, double *Rh, double D, int symmetry) {

    double Mnorm;
    double r, Rmin = 1E-05;
    // 51 steps is more than sufficient (in most cases)
    int steps = 51, i, h, j;
    double r_array[steps], rho_array[steps], sigma_array[steps],
        M_array[steps], X_array[steps];
    double random[7], x,y,z, ri, vx, vy, vz;
    double r_norm, v_norm;

    Mnorm = M(Rmax, p);
    p[1] = Mtot/Mnorm;

    double stepsize;
    stepsize = (log10(Rmax)-log10(Rmin))/(steps-1);

    printf("\nIntegrating profile...\n");
    for (i=0;i<steps;i++) {
        r = pow(10.0, log10(Rmin) + stepsize*i);
        //r = Rmin + pow(1.0*i/(steps-1.0),5)*(Rmax-Rmin);
        // 3D radius
        r_array[i] = r;
        // density(r)
        rho_array[i] = rho(r, p);
        // sigma3D(r)
        sigma_array[i] = sigma(r, p);
        // M(<r)
        M_array[i] = M(r, p);
        // M(<r)/Mtot
        X_array[i] = M_array[i]/Mtot;
    }

    printf("\n#r [pc]\t\trho_2D(r) [Msun/pc^2]\t\trho(r) [Msun/pc^3]\tsigma(r)"
        " [km/s]\t\tM(r) [Msun]\t\tX(r)\n");

    for (i=0;i<steps;i++) {
        printf ("%f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\n",
            r_array[i], rhoR(r_array[i],p), rho_array[i], sigma_array[i],
            M_array[i], X_array[i]);
    }

    // determine half-mass radius
    double Xi, xi;
    i = 0;
    while (X_array[i] < 0.5) {
        i++;
    };
    if (i) {
        xi = r_array[i-1];
        do {
            xi += 0.1;
            Xi = ((xi-r_array[i-1])*X_array[i]+(r_array[i]-xi)*
                X_array[i-1])/(r_array[i]-r_array[i-1]);
        } while ((Xi < 0.5) && (xi <= r_array[steps-1]));
        *Rh = 1.0*xi;
    } else {
        *Rh = r_array[0];
    }

    // draw positions & velocities
    if (D>=3.0) {

        for (i = 0; i< N; i++){

            if ((i/1000)*1000 == i)
                printf("Generating orbit #%i\n", i);

            do {
                random[0] = drand48();

                j = 0;
                while (X_array[j]<random[0]) {
                    j++;
                };
                if (j) {
                    xi = random[0];
                    Xi = ((xi-X_array[j-1])*r_array[j]+(X_array[j]-xi)*
                        r_array[j-1])/(X_array[j]-X_array[j-1]);
                    ri = 1.0*Xi;
                } else {
                    ri = r_array[0];
                }

                random[1] = drand48();
                random[2] = drand48();

                z = (1.0 - 2.0*random[1])*ri;
                x = sqrt(ri*ri-z*z)*cos(TWOPI*random[2]);
                y = sqrt(ri*ri-z*z)*sin(TWOPI*random[2]);
            //reject particles beyond Rmax
            } while (sqrt(x*x+y*y+z*z)>Rmax);

            double sigma;

            j = 0;
            while (r_array[j]<ri) {
                j++;
            };
            if (j) {
                xi = ri;
                Xi = ((xi-r_array[j-1])*sigma_array[j]+(r_array[j]-xi)*
                    sigma_array[j-1])/(r_array[j]-r_array[j-1]);
                sigma = 1.0*Xi;
            } else {
                sigma = sigma_array[0];
            }

            vx = get_gauss()*sigma;
            vy = get_gauss()*sigma;
            vz = get_gauss()*sigma;

            star[i].rx = x;
            star[i].ry = y;
            star[i].rz = z;
            star[i].vx = vx;
            star[i].vy = vy;
            star[i].vz = vz;

            //printf("%f\t%f\t%f\t%f\t%f\t%f\n", x,y,z, vx,vy,vz);
        }
    } else {
        fractalize(D, N, star, 1, symmetry);
        for (i=0;i<N;i++) {

            if ((i/1000)*1000 == i) printf("Generating orbit #%i\n", i);

            double sigma;
            ri = sqrt(pow(star[i].rx,2)+pow(star[i].ry,2)+pow(star[i].rz,2));

            j = 0;
            while (X_array[j]<ri) {
                j++;
            };
            if (j) {
                xi = ri;
                Xi = ((xi-X_array[j-1])*r_array[j]+(X_array[j]-xi)*
                    r_array[j-1])/(X_array[j]-X_array[j-1]);
                r_norm = 1.0*Xi;
            } else {
                r_norm = r_array[0];
            }

            j = 0;
            while (r_array[j]<r_norm) {
                j++;
            };
            if (j) {
                xi = r_norm;
                Xi = ((xi-r_array[j-1])*sigma_array[j]+(r_array[j]-xi)*
                    sigma_array[j-1])/(r_array[j]-r_array[j-1]);
                sigma = 1.0*Xi;
            } else {
                sigma = sigma_array[0];
            }

            v_norm = sigma;
            star[i].rx *= r_norm/ri;
            star[i].ry *= r_norm/ri;
            star[i].rz *= r_norm/ri;
            star[i].vx *= v_norm;
            star[i].vy *= v_norm;
            star[i].vz *= v_norm;
        }

    }


    return 0;

}

double dfridr(double (*func)(double, double*), double x, double h, double *err,
    double *p) {

    double CON = 1.4;
    double CON2 = (CON*CON);
    double SAFE = 2.0;

    int i,j;
    double errt,fac,hh,ans = 0;
    double a[10][10];

    double arg1, arg2;

    hh=h;
    a[1][1]=((*func)(x+hh, p)-(*func)(x-hh, p))/(2.0*hh);
    *err=BIGNUMBER;
    for (i=2;i<=10;i++) {
        hh /= CON;
        a[1][i]=((*func)(x+hh, p)-(*func)(x-hh, p))/(2.0*hh);
        fac=CON2;
        for (j=2;j<=i;j++) {
            a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
            fac=CON2*fac;
            arg1 = fabs(a[j][i]-a[j-1][i]);
            arg2 = fabs(a[j][i]-a[j-1][i-1]);
            if (arg1 > arg2) errt = arg1;
            else errt = arg2;
            if (errt <= *err) {
                *err=errt;
                ans=a[j][i];
            }
        }
        if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
    }

    return ans;
}

double midexp(double (*func)(double, double*), double aa, double bb, int n,
    double *p) {

    double x,tnm,sum,del,ddel,a,b;
    static double s;
    int it,j;

    b=exp(-aa);
    a=0.0;
    if (n == 1) {
        x = 0.5*(a+b);
        return (s=(b-a)*func(-log(x),p)/x);
    } else {
        for(it=1,j=1;j<n-1;j++)
            it *= 3;

        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;

        for (j=1;j<=it;j++) {
            sum += func(-log(x),p)/x;
            x += ddel;
            sum += func(-log(x),p)/x;
            x += del;
        }

        s=(s+(b-a)*sum/tnm)/3.0;
        return s;
    }
}

double midsql(double (*func)(double, double*), double aa, double bb, int n,
    double *p) {

    double x,tnm,sum,del,ddel,a,b;
    static double s;
    int it,j;

    b=sqrt(bb-aa);
    a=0.0;
    if (n == 1) {
        x = aa+0.5*(a+b)*0.5*(a+b);
        return (s=(b-a)*2.0*x*func(x, p));
    } else {
        for(it=1,j=1;j<n-1;j++)
            it *= 3;

        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;

        for (j=1;j<=it;j++) {
            sum += 2.0*x*func(aa+x*x, p);
            x += ddel;
            sum += 2.0*x*func(aa+x*x, p);
            x += del;
        }
        s=(s+(b-a)*sum/tnm)/3.0;
        return s;
    }
}

double midsqu(double (*func)(double, double*), double aa, double bb, int n,
    double *p) {

    double x,tnm,sum,del,ddel,a,b;
    static double s;
    int it,j;

    b=sqrt(bb-aa);
    a=0.0;
    if (n == 1) {
        x = bb-0.5*(a+b)*0.5*(a+b);
        return (s=(b-a)*2.0*x*func(x, p));
    } else {
        for(it=1,j=1;j<n-1;j++)
            it *= 3;

        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;
        for (j=1;j<=it;j++) {
            sum += 2.0*x*func(bb-x*x, p);
            x += ddel;
            sum += 2.0*x*func(bb-x*x, p);
            x += del;
        }
        s=(s+(b-a)*sum/tnm)/3.0;
        return s;
    }
}

double midinf(double (*func)(double, double*), double aa, double bb, int n,
    double *p) {

    double x,tnm,sum,del,ddel,b,a;
    static double s;
    int it,j;

    b=1.0/aa;
    a=1.0/bb;
    if (n == 1) {
        x= 0.5*(a+b);
        return (s=(b-a)*func(1.0/x,p)/x/x);
    } else {
        for(it=1,j=1;j<n-1;j++)
            it *= 3;

        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;

        for (j=1;j<=it;j++) {
            sum += func(1.0/x, p)/x/x;
            x += ddel;
            sum += func(1.0/x, p)/x/x;
            x += del;
        }

        return (s=(s+(b-a)*sum/tnm)/3.0);
    }
}

double midpnt(double (*func)(double, double*), double a, double b, int n,
    double *p) {

    double x,tnm,sum,del,ddel;
    static double s;
    int it,j;

    if (n == 1) {
        return (s=(b-a)*func(0.5*(a+b), p));
    } else {
        for(it=1,j=1;j<n-1;j++)
            it *= 3;

        tnm=it;
        del=(b-a)/(3.0*tnm);
        ddel=del+del;
        x=a+0.5*del;
        sum=0.0;

        for (j=1;j<=it;j++) {
            sum += func(x, p);
            x += ddel;
            sum += func(x, p);
            x += del;
        }

        s=(s+(b-a)*sum/tnm)/3.0;
        return s;
    }
}

double kernel (double x, double *p) {
    double h = 1.E-06;
    double result, err;
    result = dfridr(rhoR, x, h, &err, p)/sqrt(x*x-p[0]*p[0]);
    return result;
}

// user-defined 2d-density profile
double rhoR (double x, double *p) {
    // Elson, Fall & Freeman (1987)
    // return p[1]*pow(1.0+x/p[2]*x/p[2],-0.5*p[3]);

    // Nuker (Lauer et al. 1995) -> Elson, Fall & Freeman (1987)
    // for p[4] = 0 & p[5] = 2
    return p[1]*pow(2.0,0.5*(p[3]-p[4]))*pow(x/p[2],-p[4])*
        pow(1.0+pow(x/p[2],p[5]),-(p[3]-p[4])/p[5]);
}

double rho(double r, double *p) {
    int j = 0;
    double s,st,ost,os,stemp;
    if (r < p[2]) {
        ost = os = -1.0e30;
        for (j=1;j<=JMAX;j++) {
            if (p[4]) st=midinf(kernel,r,p[2],j,p);
            else st=midpnt(kernel,r,p[2],j,p);
            s=(4.0*st-ost)/3.0;
            if (fabs(s-os) < EPS*fabs(os)) break;
            if (s == 0.0 && os == 0.0 && j > 6) break;
            os=s;
            ost=st;
        }
        stemp = s;
        ost = os = -1.0e30;
        for (j=1;j<=JMAX;j++) {
            st=midinf(kernel,p[2],BIGNUMBER,j,p);
            s=(4.0*st-ost)/3.0;
            if (fabs(s-os) < EPS*fabs(os)) break;
            if (s == 0.0 && os == 0.0 && j > 6) break;
            os=s;
            ost=st;
        }
        s += stemp;
    } else {
        ost = os = -1.0e30;
        for (j=1;j<=JMAX;j++) {
            st=midinf(kernel,r,BIGNUMBER,j,p);
            s=(4.0*st-ost)/3.0;
            if (fabs(s-os) < EPS*fabs(os)) break;
            if (s == 0.0 && os == 0.0 && j > 6) break;
            os=s;
            ost=st;
        }
    }

    // no negative density!
    return mmax(s * -1.0/PI, 0.0);
}

double rho_kernel (double x, double *p) {
    double s;
    s = 4.0*PI*x*x*rho(x,p);
    return s;
}

double M(double r, double *p){
    int j = 0;
    double s,st,ost,os;

    ost = os = -1.0e30;
    for (j=1;j<=JMAX;j++) {
        st=midsql(rho_kernel,0.0,r,j,p);
        s=(4.0*st-ost)/3.0;
        if (fabs(s-os) < EPS*fabs(os)) break;
        if (s == 0.0 && os == 0.0 && j > 6) break;
        os=s;
        ost=st;
    }
    return s;

}

double sigma_kernel (double x, double *p) {
    double s;
    s = 1.0*rho(x,p)*M(x,p)/(x*x);
    return s;
}

double sigma(double r, double *p) {
    int j = 0;
    double s,st,ost,os,t;

    ost = os = -1.0e30;
    for (j=1;j<=JMAX;j++) {
        st=midexp(sigma_kernel,r,BIGNUMBER,j,p);
        s=(4.0*st-ost)/3.0;
        if (fabs(s-os) < EPS*fabs(os)) break;
        if (s == 0.0 && os == 0.0 && j > 6) break;
        os=s;
        ost=st;
    }

    if (rho(r,p)) t = sqrt(G*s/rho(r,p));
    else t = 0.0;

    return t;

}

double get_gauss(void){
    double random[2],p,q;
    do {
        random[0] = 2.0*drand48()-1.0;
        random[1] = 2.0*drand48()-1.0;
        q = random[0]*random[0]+random[1]*random[1];
    } while (q>1.0);

    p = sqrt(-2.0*log(q)/q);
    return random[0]*p;

}

double fractalize(double D, int N, vector<struct star_data> star, int radial, int symmetry) {
    int i=0, j, h, Nparent, Nparentlow;
    int Ntot = 128.0*pow(8,ceil(log(N)/log(8)));
    int Ntotorg = Ntot;
    double l = 2.0;
    double prob = pow(2.0, D-3.0);
    double scatter;
    if (radial) scatter = 0.01;
    else scatter = 0.1;
    double vx, vy, vz;
    double vscale;
    int subi;
    double morescatter = 0.1; //0.1 looks good

    printf("\nFractalizing initial conditions...\n\n");

    struct star_data star_temp[Ntot];   //temporary array for fractalized structure

    Nparent = 0; //ur-star
    Nparentlow = Nparent;
    star_temp[Nparent].rx = 0.0;//x
    star_temp[Nparent].ry = 0.0;//y
    star_temp[Nparent].rz = 0.0;//z
    star_temp[Nparent].vx = 0.0;//vx
    star_temp[Nparent].vy = 0.0;//vy
    star_temp[Nparent].vz = 0.0;//vz
    Nparent++;


    while (Nparent+i*8<Ntot) {
        l /= 2.0;
        i = 0;
        for (;Nparentlow<Nparent;Nparentlow++) {
            subi = 0;
            if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) &&
                (symmetry))) {

                star_temp[Nparent+i].mass_epoch = 1.0;
                star_temp[Nparent+i].rx = star_temp[Nparentlow].rx+l/2.0+
                    l*scatter*get_gauss();
                star_temp[Nparent+i].ry = star_temp[Nparentlow].ry+l/2.0+
                    l*scatter*get_gauss();
                star_temp[Nparent+i].rz = star_temp[Nparentlow].rz+l/2.0+
                    l*scatter*get_gauss();
                vx = get_gauss();
                vy = get_gauss();
                vz = get_gauss();
                star_temp[Nparent+i].vx = vx;
                star_temp[Nparent+i].vy = vy;
                star_temp[Nparent+i].vz = vz;
                i++;
                subi++;
            }
            if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) &&
                (symmetry))) {
                star_temp[Nparent+i].mass_epoch = 1.0;
                star_temp[Nparent+i].rx = star_temp[Nparentlow].rx+l/2.0+
                    l*scatter*get_gauss();
                star_temp[Nparent+i].ry = star_temp[Nparentlow].ry+l/2.0+
                    l*scatter*get_gauss();
                star_temp[Nparent+i].rz = star_temp[Nparentlow].rz-l/2.0+
                    l*scatter*get_gauss();
                vx = get_gauss();
                vy = get_gauss();
                vz = get_gauss();
                star_temp[Nparent+i].vx = vx;
                star_temp[Nparent+i].vy = vy;
                star_temp[Nparent+i].vz = vz;
                i++;
                subi++;
            }
            if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
                star_temp[Nparent+i].mass_epoch = 1.0;
                star_temp[Nparent+i].rx = star_temp[Nparentlow].rx+l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].ry = star_temp[Nparentlow].ry-l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].rz = star_temp[Nparentlow].rz+l/2.0+l*scatter*get_gauss();
                vx = get_gauss();
                vy = get_gauss();
                vz = get_gauss();
                star_temp[Nparent+i].vx = vx;
                star_temp[Nparent+i].vy = vy;
                star_temp[Nparent+i].vz = vz;
                i++;
                subi++;
            }
            if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
                star_temp[Nparent+i].mass_epoch = 1.0;
                star_temp[Nparent+i].rx = star_temp[Nparentlow].rx-l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].ry = star_temp[Nparentlow].ry+l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].rz = star_temp[Nparentlow].rz+l/2.0+l*scatter*get_gauss();
                vx = get_gauss();
                vy = get_gauss();
                vz = get_gauss();
                star_temp[Nparent+i].vx = vx;
                star_temp[Nparent+i].vy = vy;
                star_temp[Nparent+i].vz = vz;
                i++;
                subi++;
            }
            if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
                star_temp[Nparent+i].mass_epoch = 1.0;
                star_temp[Nparent+i].rx = star_temp[Nparentlow].rx+l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].ry = star_temp[Nparentlow].ry-l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].rz = star_temp[Nparentlow].rz-l/2.0+l*scatter*get_gauss();
                vx = get_gauss();
                vy = get_gauss();
                vz = get_gauss();
                star_temp[Nparent+i].vx = vx;
                star_temp[Nparent+i].vy = vy;
                star_temp[Nparent+i].vz = vz;
                i++;
                subi++;
            }
            if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
                star_temp[Nparent+i].mass_epoch = 1.0;
                star_temp[Nparent+i].rx = star_temp[Nparentlow].rx-l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].ry = star_temp[Nparentlow].ry+l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].rz = star_temp[Nparentlow].rz-l/2.0+l*scatter*get_gauss();
                vx = get_gauss();
                vy = get_gauss();
                vz = get_gauss();
                star_temp[Nparent+i].vx = vx;
                star_temp[Nparent+i].vy = vy;
                star_temp[Nparent+i].vz = vz;
                i++;
                subi++;
            }
            if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
                star_temp[Nparent+i].mass_epoch = 1.0;
                star_temp[Nparent+i].rx = star_temp[Nparentlow].rx-l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].ry = star_temp[Nparentlow].ry-l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].rz = star_temp[Nparentlow].rz+l/2.0+l*scatter*get_gauss();
                vx = get_gauss();
                vy = get_gauss();
                vz = get_gauss();
                star_temp[Nparent+i].vx = vx;
                star_temp[Nparent+i].vy = vy;
                star_temp[Nparent+i].vz = vz;
                i++;
                subi++;
            }
            if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
                star_temp[Nparent+i].mass_epoch = 1.0;
                star_temp[Nparent+i].rx = star_temp[Nparentlow].rx-l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].ry = star_temp[Nparentlow].ry-l/2.0+l*scatter*get_gauss();
                star_temp[Nparent+i].rz = star_temp[Nparentlow].rz-l/2.0+l*scatter*get_gauss();
                vx = get_gauss();
                vy = get_gauss();
                vz = get_gauss();
                star_temp[Nparent+i].vx = vx;
                star_temp[Nparent+i].vy = vy;
                star_temp[Nparent+i].vz = vz;
                i++;
                subi++;
            }

            //re-scaling of sub-group
            if (subi) {
                double vx, vy, vz;
                vx = 0.0; vy = 0.0; vz = 0.0;
                vscale = 0.0;
                for (j=Nparent+i-subi;j<Nparent+i;j++) {
                    vx += star_temp[j].vx;
                    vy += star_temp[j].vy;
                    vz += star_temp[j].vz;
                }
                vx /= 1.0*subi;
                vy /= 1.0*subi;
                vz /= 1.0*subi;
                for (j=Nparent+i-subi;j<Nparent+i;j++) {
                    star_temp[j].vx -= vx;
                    star_temp[j].vy -= vy;
                    star_temp[j].vz -= vz;
                }
                if (subi-1) {
                    for (j=Nparent+i-subi;j<Nparent+i;j++) {
                        vscale += star_temp[j].vx*star_temp[j].vx+
                            star_temp[j].vy*star_temp[j].vy+
                            star_temp[j].vz*star_temp[j].vz;
                    }
                    vscale = sqrt(vscale/(subi-1));
                } else {
                    vscale = 1.0;
                }
                for (j=Nparent+i-subi;j<Nparent+i;j++) {
                    //add bulk velocity of parent
                    star_temp[j].vx = star_temp[j].vx/vscale +
                        star_temp[Nparentlow].vx;
                    star_temp[j].vy = star_temp[j].vy/vscale +
                        star_temp[Nparentlow].vy;
                    star_temp[j].vz = star_temp[j].vz/vscale +
                        star_temp[Nparentlow].vz;
                }
            }

        }
        Nparent+=i;
    }
    Ntot = Nparent;

    double cmr[7];//centre-of-mass correction
    for (j=0; j<7; j++)
        cmr[j] = 0.0;

    for (j=0; j<Ntot; j++) {
        cmr[1] += star_temp[j].rx;
        cmr[2] += star_temp[j].ry;
        cmr[3] += star_temp[j].rz;
        cmr[4] += star_temp[j].vx;
        cmr[5] += star_temp[j].vy;
        cmr[6] += star_temp[j].vz;
    }

    for (j=0; j<Ntot; j++) {
        star_temp[j].rx -= cmr[1]/Ntot;
        star_temp[j].ry -= cmr[2]/Ntot;
        star_temp[j].rz -= cmr[3]/Ntot;
        star_temp[j].vx -= cmr[4]/Ntot;
        star_temp[j].vy -= cmr[5]/Ntot;
        star_temp[j].vz -= cmr[6]/Ntot;
    }

    i = 0;//randomly select stars from sample with r < 1.0
    for (i=0; i<N; i++) {
        do{
            j = drand48()*Ntot;
        } while (!(star_temp[j].mass_epoch) || (sqrt(pow(star_temp[j].rx,2)+
            pow(star_temp[j].ry,2)+pow(star_temp[j].rz,2)) > 1.0));
        star[i].rx = star_temp[j].rx;
        star[i].ry = star_temp[j].ry;
        star[i].rz = star_temp[j].rz;
        star[i].vx = star_temp[j].vx;
        star[i].vy = star_temp[j].vy;
        star[i].vz = star_temp[j].vz;
        star_temp[j].mass_epoch = 0.0;
    }

    double r, r_norm, vnorm = 0.0, start[4];
    if (radial) {
        for (i=0;i<N;i++) {
            vnorm += sqrt(pow(star[i].vx,2)+pow(star[i].vy,2)+
                pow(star[i].vz,2));
        }
        vnorm /= N;
        star[i].vx /= 0.5*vnorm;
        star[i].vy /= 0.5*vnorm;
        star[i].vz /= 0.5*vnorm;

        for (i=0;i<N;i++) {
            r = sqrt(pow(star[i].rx,2)+pow(star[i].ry,2)+pow(star[i].rz,2));
            r_norm = pow(r,3);
            do{
                start[1] = star[i].rx*r_norm/r + pow(r_norm/r,3)*morescatter*get_gauss();
                start[2] = star[i].ry*r_norm/r + pow(r_norm/r,3)*morescatter*get_gauss();
                start[3] = star[i].rz*r_norm/r + pow(r_norm/r,3)*morescatter*get_gauss();
            } while (sqrt(pow(start[1],2)+pow(start[2],2)+
                pow(start[3],2))>1.0);
            star[i].rx = start[1];
            star[i].ry = start[2];
            star[i].rz = start[3];
        }
    }

    return 0;
}

int get_binaries(struct options opt, vector<struct star_data> star, double M, double rvir,
    int *N) {

    int i, j, k;
    double m1, m2, ecc, abin;
    double eccold;
    double pmat[3][2], rop[2], vop[2], rrel[3], vrel[3];
    double ea, mm, eadot, cosi, inc, peri, node;
    double lP, P;
    double u1, u2;
    double q, p, x1;
    double lP1=0, lPmean=0, lPsigma=0;

    // metallicity parameters
    double zpars[20];
    // kick velocity for compact remnants
    double vkick[2];
    vkick[0] = 0.0;
    vkick[1] = 0.0;

    double vesc;
    if (opt.remnant) {
         vesc = sqrt(2.0*G*M/opt.Rh);
        if (opt.BSE)
            printf("Keeping only binaries with kick velocities < opt.escape "
                "velocity\n");
    } else {
        vesc = 1.0E10;
        if (opt.BSE) printf("Keeping all compact remnants\n");
    }

    for (i=0; i<20; i++)
        zpars[i] = 0;

    zcnsts_(&opt.Z,zpars);  //get metallicity parameters

    if (opt.BSE)
        printf("\nSetting up binary population with opt.Z = %.4f.\n", opt.Z);
    if (opt.epoch)
        printf("\nEvolving binary population for %.1f Myr.\n",opt.epoch);

    for (i=0; i < opt.nbin; i++) {
        do {
            //Specify component masses
            if (opt.BSE) {
                m1 = star[2*i].mass/M;
                m2 = star[2*i+1].mass/M;
            } else {
                m1 = star[2*i].mass_epoch;
                m2 = star[2*i+1].mass_epoch;
            }

            //Specify semi-major axis
            if (((m1*M>=opt.msort) || (m2*M>=opt.msort)) && (opt.OBperiods)) {
                if (opt.adis == 3) {
                    // Sana et al., (2012); Oh, S., Kroupa, P., &
                    // Pflamm-Altenburg, J. (2015)
                    double lPmin = 0.15, alpha = 0.45;
                    double ralpha = 1/alpha, eta = alpha/0.23;
                    lP = pow(pow(lPmin,alpha) + eta*drand48(),ralpha);
                } else {
                    // derive from Sana & Evans (2011) period distribution for
                    // massive binaries

                    // parameters of Sana & Evans (2011) period distribution in
                    // days (eq. 5.1)
                    double lPmin = 0.3, lPmax = 3.5, lPbreak = 1.0;
                    // fraction of binaries with periods below Pbreak
                    double Fbreak = 0.5;
                    double xperiod = drand48();
                    if (xperiod <= Fbreak)
                        lP = xperiod *(lPbreak-lPmin)/Fbreak + lPmin;
                    else
                        lP = (xperiod - Fbreak)*(lPmax-lPbreak)/(1.0-Fbreak) +
                            lPbreak;
                }

                // days
                P = pow(10.0,lP);
                // yr
                P /= 365.25;

                abin = pow((m1+m2)*M*P*P,(1.0/3.0));//AU
                abin /= 206264.806;//pc
                abin /= rvir;//Nbody units
            } else if (opt.adis == 0) {
                // flat semi-major axis distribution
                if (!i)
                    printf("\nApplying flat semi-major axis distribution with "
                        "amin = %g and amax = %g.\n", opt.amin, opt.amax);
                if (!i)
                    opt.amin /= rvir;
                if (!i)
                    opt.amax /= rvir;
                abin = opt.amin+drand48()*(opt.amax-opt.amin);
            } else if (opt.adis == 1 || opt.adis == 3) {
                // derive from Kroupa (1995) period distribution
                if (!i) {
                    printf("\nDeriving semi-major axis distribution from "
                        "Kroupa (1995) period distribution.\n");
                }
                // parameters of Kroupa (1995) period distribution
                double Pmin = 10, delta = 45, eta = 2.5;
                do {
                    lP = log10(Pmin) + sqrt(delta*(exp(2.0*drand48()/eta)-1));
                } while (lP > 8.43);

                P = pow(10.0,lP);//days
                P /= 365.25;//yr

                abin = pow((m1+m2)*M*P*P,(1.0/3.0));//AU
                abin /= 206264.806;//pc
                abin /= rvir;//Nbody units
            } else {
                // derive from Duquennoy & Mayor (1991) period distribution

                // mean of Duquennoy & Mayor (1991) period distribution
                // [log days]
                lPmean = 4.8;
                // full width half maximum of Duquennoy & Mayor (1991)
                // period distribution [log days]
                lPsigma = 2.3;
                if (!i) {
                    printf("\nDeriving semi-major axis distribution from "
                        "Duquennoy & Mayor (1991) period distribution.\n");
                }
                do {
                    //generate two random numbers in the interval [-1,1]
                    u1 = 2.0*drand48()-1.0;
                    u2 = 2.0*drand48()-1.0;

                    //combine the two random numbers
                    q = u1*u1 + u2*u2;
                } while (q > 1.0);

                p = sqrt(-2.0*log(q)/q);
                x1 = u1*p;
                lP1 = lPsigma*x1 + lPmean;

                P = pow(10,lP1);//days
                P /= 365.25;//yr

                abin = pow((m1+m2)*M*P*P,(1.0/3.0));//AU
                abin /= 206264.806;//pc
                abin /= rvir;//Nbody units
            }


            if (!i && opt.OBperiods && opt.msort) {
                if (opt.adis == 3) {
                    printf("\nDeriving semi-major axis distribution for "
                        "binaries with primary masses > %.3f Msun from Sana "
                        "et al., (2012); Oh, S., Kroupa, P., & "
                        "Pflamm-Altenburg, J. (2015) period distribution.\n",
                        opt.msort);
                }  else {
                    printf("\nDeriving semi-major axis distribution for "
                        "binaries with primary masses > %.3f Msun from Sana "
                        "& Evans (2011) period distribution.\n",
                        opt.msort);
                }
            }

            // Specify eccentricity distribution
            if (!i && opt.OBperiods && opt.msort) {
                printf("\nApplying thermal eccentricity distribution for "
                    "low-mass systems and Sana & Evans (2011) eccentricity "
                    "distribution for high-mass systems.\n");
            }
            else if (!i) {
                printf("\nApplying thermal eccentricity distribution.\n");
            }

            if (((m1*M>=opt.msort) || (m2*M>=opt.msort)) && (opt.OBperiods)) {
                double k1, k2;
                // maximum eccentricity for high-mass binaries
                double elim = 0.8;
                // fraction of circular orbits among high-mass binaries
                double fcirc = 0.3;
                k2 = (fcirc*exp(0.5*elim)-1.0)/(exp(0.5*elim)-1.0);
                k1 = fcirc-k2;
                ecc = drand48();
                if (ecc < fcirc) ecc = 0.0;
                else {
                    ecc = 2.0*log((ecc-k2)/k1);
                }
            } else {
                // Thermal distribution f(e)=2e
                ecc = sqrt(drand48());
            }
            //ecc = 0.0;   // all circular

            //Apply Kroupa (1995) eigenevolution
            eccold = ecc;

            if (opt.eigen) {
                if (!i) {
                    printf("\nApplying Kroupa (1995) eigenevolution for "
                        "short-period binaries\n");
                }
                // temporary re-scaling
                m1*=M;
                m2*=M;
                abin*=rvir;

                if (!(opt.OBperiods && ((m1>=opt.msort) || (m2>=opt.msort)))) {
                    if (m1>=m2)
                        eigenevolution(&m1, &m2, &ecc, &abin);
                    else
                        eigenevolution(&m2, &m1, &ecc, &abin);
                }

                m1/=M;
                m2/=M;
                abin/=rvir;
            }

            // Apply Binary Star Evolution (Hurley, Tout & Pols 2002)
            if (opt.BSE) {
                // stellar type
                int kw[2] = {1, 1};
                // initial mass
                double mass[2];
                // actual mass
                double mt[2];
                // radius
                double r[2] = {0.0, 0.0};
                // luminosity
                double lum[2] = {0.0, 0.0};
                // core mass
                double mc[2] = {0.0, 0.0};
                // core radius
                double rc[2] = {0.0, 0.0};
                // envelope mass
                double menv[2] = {0.0, 0.0};
                // envelope radius
                double renv[2] = {0.0, 0.0};
                // spin
                double ospin[2] = {0.0, 0.0};
                // time spent in current evolutionary state
                double epoch1[2] = {0.0, 0.0};
                // main-sequence lifetime
                double tms[2] = {0.0, 0.0};
                // initial age
                double tphys = 0.0;
                // final age
                double tphysf = opt.epoch;
                // data store value, if dtp>tphys no data will be stored
                double dtp = opt.epoch+1;
                // metallicity
                double z = opt.Z;

                // initial mass of primary
                mass[0] = m1*M;
                // initial mass of secondary
                mass[1] = m2*M;
                // actual mass
                mt[0] = mass[0];
                // actual mass
                mt[1] = mass[1];
                vkick[0] = 0.0;
                vkick[1] = 0.0;

                abin *= rvir; //pc
                abin *= 206264.806; //AU
                P = sqrt(pow(abin, 3.0)/(m1+m2));//yr
                P *= 365.25;//days

                evolv2_(kw,mass,mt,r,lum,mc,rc,menv,renv,ospin,epoch1,tms,
                    &tphys,&tphysf,&dtp,&z,zpars,&P,&eccold, vkick);

                P /= 365.25;//yr

                abin = pow((m1+m2)*P*P,(1.0/3.0));//AU
                abin /= 206264.806;//pc
                abin /= rvir;//Nbody units

                // Nbody units
                m1 = mt[0]/M;
                // Nbody units
                m2 = mt[1]/M;

                star[2*i].kstar = kw[0];
                star[2*i+1].kstar = kw[1];
                star[2*i].epoch1 = epoch1[0];
                star[2*i+1].epoch1 = epoch1[1];
                star[2*i].ospin = ospin[0];
                star[2*i+1].ospin = ospin[1];
                star[2*i].rstar = r[0];
                star[2*i+1].rstar = r[1];
                star[2*i].lstar = lum[0];
                star[2*i+1].lstar = lum[1];
            }

            //pos & vel in binary frame
            ea = rtnewt(ecc, drand48());
            rop[0] = abin*(cos(ea) - ecc);
            rop[1] = abin*sqrt(1.0-ecc*ecc)*sin(ea);

            mm = sqrt((m1+m2)/pow(abin,3));
            eadot = mm/(1.0 - ecc*cos(ea));
            vop[0] = -abin*sin(ea)*eadot;
            vop[1] = abin*sqrt(1.0-ecc*ecc)*cos(ea)*eadot;

            //Convert to cluster frame
            cosi = 2.0*drand48()-1.0;
            inc = acos(cosi);
            node = 2.0*PI*drand48();
            peri = 2.0*PI*drand48();

            pmat[0][0] = cos(peri)*cos(node) - sin(peri)*sin(node)*cosi;
            pmat[1][0] = cos(peri)*sin(node) + sin(peri)*cos(node)*cosi;
            pmat[2][0] = sin(peri)*sin(inc);
            pmat[0][1] = -sin(peri)*cos(node) - cos(peri)*sin(node)*cosi;
            pmat[1][1] = -sin(peri)*sin(node) + cos(peri)*cos(node)*cosi;
            pmat[2][1] = cos(peri)*sin(inc);

            for (j=0;j<3;j++) {
                rrel[j] = pmat[j][0]*rop[0] + pmat[j][1]*rop[1];
                vrel[j] = pmat[j][0]*vop[0] + pmat[j][1]*vop[1];
            }

            // Star2 pos
            star[2*i+1].rx = star[2*i].rx + m1/(m1+m2)*rrel[j];
            star[2*i+1].ry = star[2*i].ry + m1/(m1+m2)*rrel[j];
            star[2*i+1].rz = star[2*i].rz + m1/(m1+m2)*rrel[j];
            // Star2 vel
            star[2*i+1].vx = star[2*i].vx + m1/(m1+m2)*vrel[j];
            star[2*i+1].vy = star[2*i].vy + m1/(m1+m2)*vrel[j];
            star[2*i+1].vz = star[2*i].vz + m1/(m1+m2)*vrel[j];
            // Star1 pos
            star[2*i].rx  -= m2/(m1+m2)*rrel[0];
            star[2*i].ry  -= m2/(m1+m2)*rrel[1];
            star[2*i].rz  -= m2/(m1+m2)*rrel[2];
            // Star1 vel
            star[2*i].vx  -= m2/(m1+m2)*vrel[0];
            star[2*i].vy  -= m2/(m1+m2)*vrel[1];
            star[2*i].vz  -= m2/(m1+m2)*vrel[2];

            star[2*i].mass_epoch = m1;
            star[2*i+1].mass_epoch = m2;

        } while ((sqrt(vkick[0]) > vesc) && (sqrt(vkick[1]) > vesc));

        if (!star[2*i].mass_epoch || !star[2*i+1].mass_epoch) {
            // one binary less because one component went supernova
            opt.nbin--;

            // temporarily save surviving companion
            struct star_data startemp;
            if (!star[2*i].mass_epoch) {
                startemp = star[2*i+1];
            } else {
                startemp = star[2*i];
            }

            for (j=2*i+2;j<*N;j++) {
                //move all remaining stars two positions up in the array
                star[j-2] = star[j];
            }

            // reduce total number of stars by one, i.e. remove massless
            // supernova remnant from computations
            *N = *N-1;

            //make surviving companion the last particle in array
            star[*N-1] = startemp;

            //reduce binary counter
            i--;
        }
    }

    opt.N = *N;

    return 0;
}

// largest up
void shellsort(double **array, int N, int k) {
    int i,j,l,n;
    N = N-1;
    double swap[k];
    //guess distance n
    for (n = 1; n <= N/9; n = 3*n+1);
    for (; n > 0; n /= 3) {
        for (i = n; i <= N; i++) {
            for (l=0; l<k; l++) swap[l] = array[i][l];
            for (j = i; ((j >= n) && (array[j-n][0] < swap[0])); j -= n) {
                for (l=0; l<k; l++) array[j][l] = array[j-n][l];
            }
            for (l=0; l<k; l++) array[j][l] = swap[l];
        }
    }
}

//smallest up
void shellsort_reverse(double **array, int N, int k) {
    int i,j,l,n;
    N = N-1;
    double swap[k];
    //guess distance n
    for (n = 1; n <= N/9; n = 3*n+1);
    for (; n > 0; n /= 3) {
        for (i = n; i <= N; i++) {
            for (l=0; l<k; l++) swap[l] = array[i][l];
            for (j = i; ((j >= n) && (array[j-n][0] > swap[0])); j -= n) {
                for (l=0; l<k; l++) array[j][l] = array[j-n][l];
            }
            for (l=0; l<k; l++) array[j][l] = swap[l];
        }
    }
}

//largest up
void shellsort_1d(double *array, int N) {
    int i,j,n;
    N = N-1;
    double swap;
    //guess distance n
    for (n = 1; n <= N/9; n = 3*n+1);
    for (; n > 0; n /= 3) {
        for (i = n; i <= N; i++) {
            swap = array[i];
            for (j = i; ((j >= n) && (array[j-n] < swap)); j -= n) {
                array[j] = array[j-n];
            }
            array[j] = swap;
        }
    }
}

//smallest up
void shellsort_reverse_1d(double *array, int N) {
    int i,j,n;
    N = N-1;
    double swap;
    //guess distance n
    for (n = 1; n <= N/9; n = 3*n+1);
    for (; n > 0; n /= 3) {
        for (i = n; i <= N; i++) {
            swap = array[i];
            for (j = i; ((j >= n) && (array[j-n] > swap)); j -= n) {
                array[j] = array[j-n];
            }
            array[j] = swap;
        }
    }
}

int order(vector<struct star_data> star, int N, double M, double msort, int pairing){

    int i,j;
    int Nhighmass=0;
    int columns = 15;
    struct star_data star_temp[N];

    double **masses;
    masses = (double **)calloc(N,sizeof(double *));
    for (j=0;j<N;j++){
        masses[j] = (double *)calloc(2,sizeof(double));
        if (masses[j] == NULL) {
            printf("\nMemory allocation failed!\n");
            return 0;
        }
    }

    //temporary scaling to astrophysical units
    for (i=0;i<N;i++) {
        //star[i].mass_epoch *= M;
        //masses[i].mass_epoch = star[i].mass_epoch;
        masses[i][0] = star[i].mass;
        masses[i][1] = i;
    }

    shellsort(masses, N, 2);

    // Long Wang add pairing=3 for uniform mass ratio distribution of O type
    // binaries (Kiminki & Kobulnicky 2012; Sana et al., 2012;
    // Kobulnicky et al. 2014)
    // f(q) ~ q^n  (n = 0; 0.1<=q<=1.0)
    if (pairing == 3) {
        // Mask to show which mass is already used, 1 means used
        int *mask = (int*)calloc(N,sizeof(int));
        if (mask == NULL) {
            printf("\nMemory allocation failed!\n");
            return 0;
    }
    for (i=0;i<N;i++) {
        mask[i] = 0;
    }

    double **mmrand = (double**)calloc(N,sizeof(double *));
    for (j=0;j<N;j++) {
        mmrand[j] = (double *)calloc(2,sizeof(double));
        if (mmrand[j] == NULL) {
            printf("\nMemory allocation failed!\n");
            return 0;
        }
    }

    j = 0;
    int ileft = 0;
    for (i=0;i<N;i++) {
        if(mask[i]==0) {
            if(masses[i][0] >= msort) {
                // First star in binarie
                int m_tmp               = masses[i][1];
                star_temp[j].mass_epoch = star[m_tmp].mass_epoch;
                star_temp[j].rx         = star[m_tmp].rx;
                star_temp[j].ry         = star[m_tmp].ry;
                star_temp[j].rz         = star[m_tmp].rz;
                star_temp[j].vx         = star[m_tmp].vx;
                star_temp[j].vy         = star[m_tmp].vy;
                star_temp[j].vz         = star[m_tmp].vz;
                star_temp[j].mass       = star[m_tmp].mass;
                star_temp[j].kstar      = star[m_tmp].kstar;
                star_temp[j].epoch1     = star[m_tmp].epoch1;
                star_temp[j].ospin      = star[m_tmp].ospin;
                star_temp[j].rstar      = star[m_tmp].rstar;
                star_temp[j].lstar      = star[m_tmp].lstar;
                star_temp[j].epochstar  = star[m_tmp].epochstar;
                star_temp[j].zstar      = star[m_tmp].zstar;
                mask[i] = 1;
                j++;
                // Find the second one based on uniform mass ratio
                if (i<N-1) {
                    double mpair = (drand48()*0.9+0.1)*masses[i][0];
                    int k = i+1;
                    while (masses[k][0] > mpair) {
                        k++;
                        if (k>=N) {
                            k--;
                            break;
                        }
                    }
                    int k1 = k-1;
                    int k2 = k;
                    if (k1 == i) {
                        while(k1<N && mask[k1]==1)
                            k1++;

                        if(mask[k1]==1)
                            k = -1;
                        else
                            k = k1;
                    } else {
                        while(k1>i && mask[k1]==1)
                            k1--;

                        while(k2<N && mask[k2]==1)
                            k2++;

                        if(mask[k1]==1) {
                            if(mask[k2]==1)
                                k = -1;
                            else
                                k = k2;
                        } else {
                            if (mask[k2]==0 &&
                                mpair-masses[k2][0]<masses[k1][0]-mpair) {
                                k = k2;
                            }
                            else
                                k = k1;
                        }
                    }
                    if (k>0) {
                        if(j>=N) {
                            perror("\n Error!: Pairing binary star exceed "
                                "the maximum index!\n");
                            exit(0);
                        }
                        int k_tmp = masses[k][1];
                        // Second star in binarie
                        star_temp[j].mass_epoch = star[k_tmp].mass_epoch;
                        star_temp[j].rx         = star[k_tmp].rx;
                        star_temp[j].ry         = star[k_tmp].ry;
                        star_temp[j].rz         = star[k_tmp].rz;
                        star_temp[j].vx         = star[k_tmp].vx;
                        star_temp[j].vy         = star[k_tmp].vy;
                        star_temp[j].vz         = star[k_tmp].vz;
                        star_temp[j].mass       = star[k_tmp].mass;
                        star_temp[j].kstar      = star[k_tmp].kstar;
                        star_temp[j].epoch1     = star[k_tmp].epoch1;
                        star_temp[j].ospin      = star[k_tmp].ospin;
                        star_temp[j].rstar      = star[k_tmp].rstar;
                        star_temp[j].lstar      = star[k_tmp].lstar;
                        star_temp[j].epochstar  = star[k_tmp].epochstar;
                        star_temp[j].zstar      = star[k_tmp].zstar;
                        mask[k] = 1;
                        j++;
                    } else {
                      printf("Mass %g did not find good pair, expected pair "
                        "mass: %g\n",masses[i][0],mpair);
                    }
                }
           }
           else {
                // Store index for remaining stars
                mmrand[ileft][0] = drand48();
                mmrand[ileft][1] = masses[i][1];
                ileft++;
           }
        }
    }

    // Error check, ileft+j should be N
    if(ileft+j!=N) {
        fprintf(stderr,"\n Error! M>Msort binaries + M<Msort binaries did not "
            "match total number M>Msort: %d, M<Msort: %d, N: %d\n",j,ileft,N);
        exit(0);
    }
    if (msort) {
      Nhighmass = j;
      printf("\n Number of binaries with M>Msort: %d",j/2);
    }

    // Randomize indexs for random pairing
    shellsort(mmrand,ileft,2);

    // Random pairing remaining binaries
    for(i=0;i<ileft;i++) {
        star_temp[j+i].mass_epoch = star[(int) mmrand[i][1]].mass_epoch;
        star_temp[j+i].rx         = star[(int) mmrand[i][1]].rx;
        star_temp[j+i].ry         = star[(int) mmrand[i][1]].ry;
        star_temp[j+i].rz         = star[(int) mmrand[i][1]].rz;
        star_temp[j+i].vx         = star[(int) mmrand[i][1]].vx;
        star_temp[j+i].vy         = star[(int) mmrand[i][1]].vy;
        star_temp[j+i].vz         = star[(int) mmrand[i][1]].vz;
        star_temp[j+i].mass       = star[(int) mmrand[i][1]].mass;
        star_temp[j+i].kstar      = star[(int) mmrand[i][1]].kstar;
        star_temp[j+i].epoch1     = star[(int) mmrand[i][1]].epoch1;
        star_temp[j+i].ospin      = star[(int) mmrand[i][1]].ospin;
        star_temp[j+i].rstar      = star[(int) mmrand[i][1]].rstar;
        star_temp[j+i].lstar      = star[(int) mmrand[i][1]].lstar;
        star_temp[j+i].epochstar  = star[(int) mmrand[i][1]].epochstar;
        star_temp[j+i].zstar      = star[(int) mmrand[i][1]].zstar;
    }

    } else {

        j = 0;
        for (i=0;i<N;i++) {
            //copying to front of temporary array if massive
            if (masses[i][0] >= msort) {
                star_temp[j].mass_epoch = star[(int) masses[i][1]].mass_epoch;
                star_temp[j].rx         = star[(int) masses[i][1]].rx;
                star_temp[j].ry         = star[(int) masses[i][1]].ry;
                star_temp[j].rz         = star[(int) masses[i][1]].rz;
                star_temp[j].vx         = star[(int) masses[i][1]].vx;
                star_temp[j].vy         = star[(int) masses[i][1]].vy;
                star_temp[j].vz         = star[(int) masses[i][1]].vz;
                star_temp[j].mass       = star[(int) masses[i][1]].mass;
                star_temp[j].kstar      = star[(int) masses[i][1]].kstar;
                star_temp[j].epoch1     = star[(int) masses[i][1]].epoch1;
                star_temp[j].ospin      = star[(int) masses[i][1]].ospin;
                star_temp[j].rstar      = star[(int) masses[i][1]].rstar;
                star_temp[j].lstar      = star[(int) masses[i][1]].lstar;
                star_temp[j].epochstar  = star[(int) masses[i][1]].epochstar;
                star_temp[j].zstar      = star[(int) masses[i][1]].zstar;
                j++;
            }
        }
        if (msort)
            Nhighmass = j;

        for (i=0;i<N;i++) {
            // copying to random position in the back of temporary array
            if (masses[i][0] < msort) {
                do {
                    j = drand48()*N;
                } while (star_temp[j].mass_epoch);

                star_temp[j].mass_epoch = star[(int) masses[i][1]].mass_epoch;
                star_temp[j].rx         = star[(int) masses[i][1]].rx;
                star_temp[j].ry         = star[(int) masses[i][1]].ry;
                star_temp[j].rz         = star[(int) masses[i][1]].rz;
                star_temp[j].vx         = star[(int) masses[i][1]].vx;
                star_temp[j].vy         = star[(int) masses[i][1]].vy;
                star_temp[j].vz         = star[(int) masses[i][1]].vz;
                star_temp[j].mass       = star[(int) masses[i][1]].mass;
                star_temp[j].kstar      = star[(int) masses[i][1]].kstar;
                star_temp[j].epoch1     = star[(int) masses[i][1]].epoch1;
                star_temp[j].ospin      = star[(int) masses[i][1]].ospin;
                star_temp[j].rstar      = star[(int) masses[i][1]].rstar;
                star_temp[j].lstar      = star[(int) masses[i][1]].lstar;
                star_temp[j].epochstar  = star[(int) masses[i][1]].epochstar;
                star_temp[j].zstar      = star[(int) masses[i][1]].zstar;
            }
        }
    }

    //copying back to original array and scaling back to Nbody units
    for (i=0;i<N;i++) {
        //star[i].mass_epoch = star_temp[i].mass_epoch/M;
        star[i].mass_epoch = star_temp[i].mass_epoch;
        star[i].rx         = star_temp[i].rx;
        star[i].ry         = star_temp[i].ry;
        star[i].rz         = star_temp[i].rz;
        star[i].vx         = star_temp[i].vx;
        star[i].vy         = star_temp[i].vy;
        star[i].vz         = star_temp[i].vz;
        star[i].mass       = star_temp[i].mass;
        star[i].kstar      = star_temp[i].kstar;
        star[i].epoch1     = star_temp[i].epoch1;
        star[i].ospin      = star_temp[i].ospin;
        star[i].rstar      = star_temp[i].rstar;
        star[i].lstar      = star_temp[i].lstar;
        star[i].epochstar  = star_temp[i].epochstar;
        star[i].zstar      = star_temp[i].zstar;
    }

    for (j=0;j<N;j++)
        free (masses[j]);

    free(masses);

    if (pairing==2)
        segregate(star, Nhighmass, 0.0);

    return 0;
}

int segregate(vector<struct star_data> star, int N, double S){
    int i,j;

    int columns = 15;
    struct star_data star_temp[N];

    double **masses;
    masses = (double **)calloc(N,sizeof(double *));

    for (j=0;j<N;j++){
        masses[j] = (double *)calloc(2,sizeof(double));

        if (masses[j] == NULL) {
            printf("\nMemory allocation failed!\n");
            return 0;
        }
    }

    for (i=0;i<N;i++) {
        masses[i][0] = star[i].mass_epoch;
        masses[i][1] = i;
    }

    shellsort(masses, N, 2);

    for (i=0;i<N;i++)
        star_temp[i].mass_epoch = 0.0;

    j = 0;
    int Ntemp,l;
    Ntemp = N;

    for (i=0;i<N;i++) {
        j = 1.0*(1.0-pow(drand48(),1.0-S))*Ntemp;
        l=-1;

        do {
            l++;
            if (star_temp[l].mass_epoch) {
                j++;
            }
        } while (l<j);

        star_temp[j].mass_epoch = star[(int) masses[i][1]].mass_epoch;
        star_temp[j].rx         = star[(int) masses[i][1]].rx;
        star_temp[j].ry         = star[(int) masses[i][1]].ry;
        star_temp[j].rz         = star[(int) masses[i][1]].rz;
        star_temp[j].vx         = star[(int) masses[i][1]].vx;
        star_temp[j].vy         = star[(int) masses[i][1]].vy;
        star_temp[j].vz         = star[(int) masses[i][1]].vz;
        star_temp[j].mass       = star[(int) masses[i][1]].mass;
        star_temp[j].kstar      = star[(int) masses[i][1]].kstar;
        star_temp[j].epoch1     = star[(int) masses[i][1]].epoch1;
        star_temp[j].ospin      = star[(int) masses[i][1]].ospin;
        star_temp[j].rstar      = star[(int) masses[i][1]].rstar;
        star_temp[j].lstar      = star[(int) masses[i][1]].lstar;
        star_temp[j].epochstar  = star[(int) masses[i][1]].epochstar;
        star_temp[j].zstar      = star[(int) masses[i][1]].zstar;
        Ntemp--;
    }

    //copying back to original array
    for (i=0;i<N;i++) {
        star[i].mass_epoch = star_temp[i].mass_epoch;
        star[i].rx = star_temp[i].rx;
        star[i].ry = star_temp[i].ry;
        star[i].rz = star_temp[i].rz;
        star[i].vx = star_temp[i].vx;
        star[i].vy = star_temp[i].vy;
        star[i].vz = star_temp[i].vz;
        star[i].mass = star_temp[i].mass;
        star[i].kstar = star_temp[i].kstar;
        star[i].epoch1 = star_temp[i].epoch1;
        star[i].ospin = star_temp[i].ospin;
        star[i].rstar = star_temp[i].rstar;
        star[i].lstar = star_temp[i].lstar;
        star[i].epochstar = star_temp[i].epochstar;
        star[i].zstar = star_temp[i].zstar;
    }

    for (j=0;j<N;j++)
        free (masses[j]);

    free(masses);

    return 0;
}

int energy_order(vector<struct star_data> star, int N, int Nstars){
    int i,j;

    int columns = 15;
    struct star_data star_temp[N];

    double **energies;
    energies = (double **)calloc(N,sizeof(double *));

    for (j=0;j<N;j++){
        energies[j] = (double *)calloc(2,sizeof(double));
        energies[j][0] = 0.0;

        if (energies[j] == NULL) {
            printf("\nMemory allocation failed!\n");
            return 0;
        }
    }

    for (i=0;i<N;i++) {
        //Ekin + Epot = 0.5*1/N*v^2-1/r
        energies[i][0] = 0.5/Nstars*(star[i].vx*star[i].vx+star[i].vy*
            star[i].vy+star[i].vz*star[i].vz)-1.0/sqrt(star[i].rx*star[i].rx+
            star[i].ry*star[i].ry+star[i].rz*star[i].rz);
        energies[i][1] = i;
    }
    shellsort_reverse(energies, N, 2);

    for (i=0;i<N;i++) {//copying to random position in temporary array
        star_temp[i].mass_epoch = star[(int) energies[i][1]].mass_epoch;
        star_temp[i].rx = star[(int) energies[i][1]].rx;
        star_temp[i].ry = star[(int) energies[i][1]].ry;
        star_temp[i].rz = star[(int) energies[i][1]].rz;
        star_temp[i].vx = star[(int) energies[i][1]].vx;
        star_temp[i].vy = star[(int) energies[i][1]].vy;
        star_temp[i].vz = star[(int) energies[i][1]].vz;
        star_temp[i].mass = star[(int) energies[i][1]].mass;
        star_temp[i].kstar = star[(int) energies[i][1]].kstar;
        star_temp[i].epoch1 = star[(int) energies[i][1]].epoch1;
        star_temp[i].ospin = star[(int) energies[i][1]].ospin;
        star_temp[i].rstar = star[(int) energies[i][1]].rstar;
        star_temp[i].lstar = star[(int) energies[i][1]].lstar;
    }


    //copying back to original array
    for (i=0;i<N;i++) {
        star[i].mass_epoch = star_temp[i].mass_epoch;
        star[i].rx = star_temp[i].rx;
        star[i].ry = star_temp[i].ry;
        star[i].rz = star_temp[i].rz;
        star[i].vx = star_temp[i].vx;
        star[i].vy = star_temp[i].vy;
        star[i].vz = star_temp[i].vz;
        star[i].mass = star_temp[i].mass;
        star[i].kstar = star_temp[i].kstar;
        star[i].epoch1 = star_temp[i].epoch1;
        star[i].ospin = star_temp[i].ospin;
        star[i].rstar = star_temp[i].rstar;
        star[i].lstar = star_temp[i].lstar;
        star[i].epochstar = star_temp[i].epochstar;
        star[i].zstar = star_temp[i].zstar;
    }

    for (j=0;j<N;j++) free (energies[j]);
    free(energies);

    return 0;
}

int randomize(vector<struct star_data> star, int N){
    int i,j;

    int columns = 15;
    struct star_data star_temp[N];

    for (j=0;j<N;j++){
        star_temp[j].mass_epoch = 0.0;
    }

    for (i=0;i<N;i++) {//copying randomized to temporary array
            do {
                j = drand48()*N;
            } while (star_temp[j].mass_epoch);
            star_temp[j].mass_epoch = star[i].mass_epoch;
            star_temp[j].rx = star[i].rx;
            star_temp[j].ry = star[i].ry;
            star_temp[j].rz = star[i].rz;
            star_temp[j].vx = star[i].vx;
            star_temp[j].vy = star[i].vy;
            star_temp[j].vz = star[i].vz;
            star_temp[j].mass = star[i].mass;
            star_temp[j].kstar = star[i].kstar;
            star_temp[j].epoch1 = star[i].epoch1;
            star_temp[j].ospin = star[i].ospin;
            star_temp[j].rstar = star[i].rstar;
            star_temp[j].lstar = star[i].lstar;
            star_temp[j].epochstar = star[i].epochstar;
            star_temp[j].zstar = star[i].zstar;
    }

    //copying back to original array
    for (i=0;i<N;i++) {
        star[i].mass_epoch = star_temp[i].mass_epoch;
        star[i].rx = star_temp[i].rx;
        star[i].ry = star_temp[i].ry;
        star[i].rz = star_temp[i].rz;
        star[i].vx = star_temp[i].vx;
        star[i].vy = star_temp[i].vy;
        star[i].vz = star_temp[i].vz;
        star[i].mass = star_temp[i].mass;
        star[i].kstar = star_temp[i].kstar;
        star[i].epoch1 = star_temp[i].epoch1;
        star[i].ospin = star_temp[i].ospin;
        star[i].rstar = star_temp[i].rstar;
        star[i].lstar = star_temp[i].lstar;
        star[i].epochstar = star_temp[i].epochstar;
        star[i].zstar = star_temp[i].zstar;
    }

    return 0;
}

double rtnewt (double ecc, double ma) {

    double x1,x2,xacc,rtnewt,f,df,dx;
    int j,jmax;

    x1 = 0;
    x2 = 2*PI;
    xacc = 1E-6;
    jmax = 20;
    ma = 2*PI*ma;

    rtnewt=.5*(x1+x2);
    for (j=1;j<=jmax;j++) {
        f = ma - rtnewt + ecc*sin(rtnewt);
        df = -1 + ecc*cos(rtnewt);
        dx=f/df;
        rtnewt=rtnewt-dx;
        if ((x1-rtnewt)*(rtnewt-x2)<0)
            printf("jumped out of brackets\n");
        if(abs(dx)<xacc) return(rtnewt);
    }
    printf("RTNEWT exceeding maximum iterations\n");
    exit(-1);
}

int eigenevolution(double *m1, double *m2, double *ecc, double *abin){
    double alpha = 28.0;
    double beta = 0.75;
    double mtot,lper,lperi,ecci,mtoti,r0,rperi,qold,qnew;

    *abin *= PARSEC;
    mtot = *m1+*m2;

    lperi = sqrt(pow(*abin,3)/mtot*4.0*PI*PI/GBIN);

    ecci = *ecc;
    mtoti = mtot;

    /* Circularisation */
    r0 = alpha*RSUN;
    rperi = *abin*(1.0-*ecc);
    alpha = -pow((r0/rperi),beta);

    if (*ecc > 0) {
        *ecc = exp(alpha+log(*ecc));
    }

    /* pre-ms mass-transfer */
    qold = *m1/ *m2;
    if (qold > 1.0) qold = 1.0/qold;
    alpha = -alpha;
    if (alpha > 1.0) {
        qnew = 1.0;
    } else {
        qnew = qold + (1.0-qold)*alpha;
    }

    *m1 = mmax(*m1,*m2);
    *m2 = qnew * *m1;

    mtot = *m1 + *m2;

    lper = lperi*pow((1.0-ecci)/(1.0-*ecc),1.5);
    lper = lper*sqrt(mtoti/mtot);

    *abin = pow(mtot*lper*lper*GBIN/4.0/PI/PI,1.0/3.0);
    *abin /= PARSEC;

    return 0;
}

int radial_profile(vector<struct star_data> star, struct options opt, double rvir,
    double M, int create_radial_profile, int create_cumulative_profile,
    int *NNBMAX, double *RS0, double *Rh2D, double *Rh3D) {

    int i, j;
    *Rh2D = 0.0;
    *Rh3D = 0.0;
    double Mtemp;
    double **rarray;
    rarray = (double **)calloc(opt.N,sizeof(double *));
    for (j=0;j<opt.N;j++){
        rarray[j] = (double *)calloc(3,sizeof(double));
        if (rarray[j] == NULL) {
            printf("\nMemory allocation failed!\n");
            return 0;
        }
    }
    double **rarray2D;
    rarray2D = (double **)calloc(opt.N,sizeof(double *));
    for (j=0;j<opt.N;j++){
        rarray2D[j] = (double *)calloc(3,sizeof(double));
        if (rarray2D[j] == NULL) {
            printf("\nMemory allocation failed!\n");
            return 0;
        }
    }

    for (j=0; j<opt.N; j++) {
        rarray[j][0] = rvir*sqrt(star[j].rx*star[j].rx+star[j].ry*
            star[j].ry+ star[j].rz*star[j].rz);
        rarray[j][1] = star[j].mass_epoch*M;
        rarray[j][2] = star[j].lstar;
        rarray2D[j][0] = rvir*sqrt(star[j].rx*star[j].rx+star[j].ry*
            star[j].ry);
        rarray2D[j][1] = star[j].mass_epoch*M;
        rarray2D[j][2] = star[j].lstar;
    }
    shellsort_reverse(rarray, opt.N, 3);
    shellsort_reverse(rarray2D, opt.N, 3);

    int noofradialbins = 20;
    double rprofile[noofradialbins][9];
    double rmax = 100.0; //pc
    double rmin = 0.1;
    double stepsize;
    stepsize = (log10(rmax)-log10(rmin))/(noofradialbins-1);

    for (j=0;j<noofradialbins;j++) {
        rprofile[j][0] = pow(10.0, log10(rmin) + stepsize*j); //radius
        //rprofile[j][0] = 1.0*exp(log(rmax)*(j+1)/noofradialbins)-1.0; //radius
        rprofile[j][2] = 0; //number
        rprofile[j][3] = 0; //mass
        rprofile[j][4] = 0; //luminosity
        if (j == 0)
            rprofile[j][1] =  4.0/3.0*PI*pow(rprofile[j][0],3); //volume
        else {
            rprofile[j][1] = 4.0/3.0*PI*(pow(rprofile[j][0],3) -
                pow(rprofile[j-1][0],3)); //volume
        }
        rprofile[j][6] = 0; //number (2D)
        rprofile[j][7] = 0; //mass (2D)
        rprofile[j][8] = 0; //luminosity (2D)
        if (j == 0)
            rprofile[j][5] =  PI*pow(rprofile[j][0],2); //area
        else {
            rprofile[j][5] = PI*(pow(rprofile[j][0],2) -
                pow(rprofile[j-1][0],2)); //area
        }
    }

    j = 0; i = 0; Mtemp = 0.0;
    while ((j < noofradialbins) && (i < opt.N)) {
        if (rarray[i][0] < rprofile[j][0]) {//3D binning in astrophysical units
            rprofile[j][2] += 1.0;
            rprofile[j][3] += rarray[i][1];
            rprofile[j][4] += rarray[i][2];
            Mtemp += rarray[i][1];
            if ((Mtemp>=0.5*M) && !(*Rh3D)) *Rh3D = rarray[i][0];
            i++;
        } else {
            j++;
        }
    }

    j = 0; i = 0; Mtemp = 0.0;
    while ((j < noofradialbins) && (i < opt.N)) {
        if (rarray2D[i][0] < rprofile[j][0]) {//2D binning in astrophysical units
            rprofile[j][6] += 1.0;
            rprofile[j][7] += rarray2D[i][1];
            rprofile[j][8] += rarray2D[i][2];
            Mtemp += rarray2D[i][1];
            if ((Mtemp>=0.5*M) && !(*Rh2D)) *Rh2D = rarray2D[i][0];
            i++;
        } else {
            j++;
        }
    }

    if (create_radial_profile) {
        printf("\nRadial density profile:\n\n#  R [pc]   N [1/pc^3]   "
            "M [Msun/pc^3]   L [Lsun/pc^3]   N [1/pc^2]   M [Msun/pc^2]   "
            "L [Lsun/pc^2] \n");
        for (j=0;j<noofradialbins;j++) {
            //print radial density profile to screen
            printf("%9.4f  %11.4f  %14.4f  %14.4f  %11.4f  %14.4lf %15.4f\n",
                rprofile[j][0],rprofile[j][2]/rprofile[j][1],
                rprofile[j][3]/rprofile[j][1],rprofile[j][4]/rprofile[j][1],
                rprofile[j][6]/rprofile[j][5],rprofile[j][7]/rprofile[j][5],
                rprofile[j][8]/rprofile[j][5]);
        }
    }
    if (create_cumulative_profile) {
        double ntemp = 0.0;
        double mtemp = 0.0;
        double ltemp = 0.0;
        double ntemp2D = 0.0;
        double mtemp2D = 0.0;
        double ltemp2D = 0.0;
        double nmax, mmax, lmax, nhalf, mhalf, lhalf;
        double nhalf2D, mhalf2D, lhalf2D;
        printf("\nCumulative profile:\n\n#  R [pc]         N     M [Msun]     "
            "L [Lsun]    N (2D)    M(2D) [Msun]    L(2D) [Lsun]\n");
        for (j=0;j<noofradialbins;j++) {
            ntemp += rprofile[j][2];
            mtemp += rprofile[j][3];
            ltemp += rprofile[j][4];
            ntemp2D += rprofile[j][6];
            mtemp2D += rprofile[j][7];
            ltemp2D += rprofile[j][8];
            //print cumulative profile to screen
            printf("%9.4f  %8.0f  %11.3f  %11.3f  %8.0f  %14.4f  %14.4f\n",
                rprofile[j][0],ntemp,mtemp,ltemp,ntemp2D,mtemp2D,ltemp2D);
        }
        nmax = ntemp;
        mmax = mtemp;
        lmax = ltemp;
        nhalf = mhalf = lhalf = 0.0;
        nhalf2D = mhalf2D = lhalf2D = 0.0;
        ntemp = mtemp = ltemp = 0.0;
        ntemp2D = mtemp2D = ltemp2D = 0.0;
        for (j=0;j<opt.N;j++) {
            ntemp += 1.0;
            mtemp += rarray[j][1];
            ltemp += rarray[j][2];
            if (!(nhalf) && (ntemp>=0.5*nmax)) nhalf = rarray[j][0];
            if (!(mhalf) && (mtemp>=0.5*mmax)) mhalf = rarray[j][0];
            if (!(lhalf) && (ltemp>=0.5*lmax)) lhalf = rarray[j][0];
            ntemp2D += 1.0;
            mtemp2D += rarray2D[j][1];
            ltemp2D += rarray2D[j][2];
            if (!(nhalf2D) && (ntemp2D>=0.5*nmax)) nhalf2D = rarray2D[j][0];
            if (!(mhalf2D) && (mtemp2D>=0.5*mmax)) mhalf2D = rarray2D[j][0];
            if (!(lhalf2D) && (ltemp2D>=0.5*lmax)) lhalf2D = rarray2D[j][0];
        }
        //print radii
        printf("\nHalf-number radius = %.4f pc (%.4f pc, 2D)\nHalf-mass "
            "radius = %.4f pc (%.4f pc, 2D)\nHalf-light radius = %.4f pc "
            "(%.4f pc, 2D)\n",nhalf,nhalf2D,mhalf,mhalf2D,lhalf,lhalf2D);

    }

    //estimate NNBMAX and RS0 (Nbody6 only)
    if ((opt.code == 0) || (opt.code == 2) || opt.code == 4 || (opt.code == 5)) {
        *NNBMAX = 2.0*sqrt(opt.N);
        if (*NNBMAX < 30) *NNBMAX = 30;
        if (opt.N<=*NNBMAX) *NNBMAX = 0.5*opt.N;
        if (*NNBMAX > opt.NNBMAX_NBODY6) *NNBMAX = opt.NNBMAX_NBODY6;
        *RS0 = rarray[*NNBMAX][0];
        printf("\nEstimating appropriate NNBMAX = %i and RS0 = %f [pc]\n",
            *NNBMAX,*RS0);
    }

    for (j=0;j<opt.N;j++) free (rarray[j]);
    free(rarray);
    for (j=0;j<opt.N;j++) free (rarray2D[j]);
    free(rarray2D);

    return 0;
}

int cmd(vector<struct star_data> star, int l, double Rgal, double *abvmag, double *vmag,
    double *BV, double *Teff, double *dvmag, double *dBV) {

    double lTeff, BC, kb;
    double bvc[8], bcc[8];
    double dbmag;
    double BCsun, abvmagsun;

    // Stefan-Boltzmann constant in Lsun Rsun^-2 K^-4
    kb = 5.6704E-08*0.5*1.3914E9*0.5*1.3914E9/3.846E26;

    bvc[0] = -654597.405559323;
    bvc[1] = 1099118.61158915;
    bvc[2] = -789665.995692672;
    bvc[3] = 314714.220932623;
    bvc[4] = -75148.4728506455;
    bvc[5] = 10751.803394526;
    bvc[6] = -853.487897283685;
    bvc[7] = 28.9988730655392;

    bcc[0] = -4222907.80590972;
    bcc[1] = 7209333.13326442;
    bcc[2] = -5267167.04593882;
    bcc[3] = 2134724.55938336;
    bcc[4] = -518317.954642773;
    bcc[5] = 75392.2372207101;
    bcc[6] = -6082.7301194776;
    bcc[7] = 209.990478646363;

    BCsun = 0.11;  //sun's bolometric correction
    abvmagsun = 4.83; //sun's absolute V magnitude

    if (star[l].rstar && (star[l].kstar<14)) {
        *Teff = pow(star[l].lstar/(4.0*PI*star[l].rstar*star[l].rstar*kb),0.25);
        if ((*Teff>3000.0) && (*Teff<55000.0)) {

            lTeff = log10(*Teff);

            *BV = bvc[0] + bvc[1]*lTeff + bvc[2]*pow(lTeff,2) +
                bvc[3]*pow(lTeff,3) + bvc[4]*pow(lTeff,4) +
                bvc[5]*pow(lTeff,5) + bvc[6]*pow(lTeff,6) +
                bvc[7]*pow(lTeff,7);

            BC = bcc[0] + bcc[1]*lTeff + bcc[2]*pow(lTeff,2) +
                bcc[3]*pow(lTeff,3) + bcc[4]*pow(lTeff,4) +
                bcc[5]*pow(lTeff,5) + bcc[6]*pow(lTeff,6) +
                bcc[7]*pow(lTeff,7);

            if (star[l].lstar)
                *abvmag = -2.5*log10(star[l].lstar)-BC+BCsun+abvmagsun;

            *vmag = *abvmag + 5.0*log10(Rgal) - 5.0;


            double rand1, rand2, prand;

            do {
                rand1 = 2.0*drand48()-1.0;
                rand2 = 2.0*drand48()-1.0;
            } while (rand1*rand1+rand2*rand2 > 1.0);

            prand = sqrt(-2.0*log(rand1*rand1+rand2*rand2)/(rand1*rand1+
                rand2*rand2));
            *dvmag = rand1*prand*sqrt(pow(0.02,2) + pow(0.07*pow(10.0,
                0.4*(*vmag-25.0)),2));
            dbmag = rand2*prand*sqrt(pow(0.02,2) + pow(0.07*pow(10.0,
                0.4*(*vmag-25.0)),2));
            *dBV = *dvmag + dbmag;

            } else {
                *vmag = 9999.9;
                *abvmag = 9999.9;
                *BV = 9999.9;
                *dvmag = 0.0;
                *dBV = 0.0;
            }
        } else {
            *Teff = 0.0;
            *vmag = 9999.9;
            *abvmag = 9999.9;
            *BV = 9999.9;
            *dvmag = 0.0;
            *dBV = 0.0;
        }

    return 0;
}

int output0(struct options opt, int N, int NNBMAX, double RS0, double rvir,
    double mmean, int bin, double M, double MMAX,  double rtide,
    star_data star[], int sse) {

    //Open output files
    char PARfile[50], NBODYfile[50], SSEfile[50];
    FILE *PAR=NULL, *NBODY=NULL, *SSE12=NULL;
    sprintf(PARfile, "%s.input", opt.output);
    PAR = fopen(PARfile,"w");
    sprintf(NBODYfile, "%s.fort.10", opt.output);
    NBODY = fopen(NBODYfile,"w");

    int hrplot = 0;
    if (opt.dtplot) hrplot = 1;
    if (sse) {
        sprintf(SSEfile, "%s.fort.12", opt.output);
        SSE12 = fopen(SSEfile,"w");
        hrplot = 2;
    }

    //write to .PAR file
    fprintf(PAR,"1 5000000.0 0\n");
    fprintf(PAR,"%i 1 10 %i %i 1\n",opt.N,opt.seed,NNBMAX);
    fprintf(PAR,"0.02 0.02 %.8lf %.8lf %.8lf %.8f 1.0E-03 %.8lf %.8f\n",
        RS0,opt.dtadj,opt.dtout,opt.tcrit,rvir,mmean);
    fprintf(PAR,"2 2 1 0 1 0 2 0 0 2\n");
    fprintf(PAR,"0 %i 0 %i 2 %i %i 0 %i 3\n",
        hrplot,opt.tf,opt.regupdate,opt.etaupdate,opt.mloss);
    fprintf(PAR,"0 %i %i 0 1 2 0 1 0 1\n",bin,opt.esc);
    fprintf(PAR,"0 0 0 2 1 0 0 2 0 3\n");
    fprintf(PAR,"0 0 0 0 0 0 0 0 0 0\n");
    fprintf(PAR,"1.0E-5 1.0E-4 0.2 1.0 1.0E-06 0.001\n");
    fprintf(PAR,"2.350000 %.8lf %.8lf %i 0 %.8lf %.8lf %.8f\n",
        MMAX,mlow,opt.nbin,opt.Z,opt.epoch,opt.dtplot);
    fprintf(PAR,"%.2f 0.0 0.0 0.00000 0.125\n",opt.Q);
    if (opt.tf == 2) {
        fprintf(PAR,"%.8e %.8f\n",
            M1pointmass,sqrt(opt.RG[0]*opt.RG[0]+opt.RG[1]*opt.RG[1]+opt.RG[2]*opt.RG[2])/1000.0);
    } else if (opt.tf == 3) {
        fprintf(PAR,"%.6e %.6e %.6lf %.6lf %.6lf %.6lf %.6e %.6lf %.6f\n",
            M1allen,M2allen,a2allen,b2allen,VCIRC,RCIRC, GMB, AR, GAM);
        fprintf(PAR,"%.6lf %.6lf %.6lf %.6lf %.6lf %.6f\n",
            opt.RG[0]/1000.0,opt.RG[1]/1000.0,opt.RG[2]/1000.0,opt.VG[0],opt.VG[1],opt.VG[2]);
    }
    if (opt.tf > 2) {
        fprintf(PAR,"%.6lf %.6lf %.6lf %.6f\n", opt.extmass,opt.extrad,opt.extdecay,opt.extstart);
    }

    //write to .fort.10 file
    int j;
    for (j=0;j<opt.N;j++) {
        fprintf(NBODY,"%.16lf\t%.16lf %.16lf %.16lf\t%.16lf %.16lf %.16lf\n",
            star[j].mass_epoch,star[j].rx,star[j].ry,star[j].rz,star[j].vx,star[j].vy,
            star[j].vz);
    }

    //write to .fort.12 file
    if (sse) {
        for (j=0;j<opt.N;j++) {
            fprintf(SSE12,"%.8lf\t%.0lf %.8lf %.8lf %.8lf\n",
                star[j].mass_epoch*M,star[j].kstar,star[j].mass,star[j].epoch1,star[j].ospin);
            //,star[j].epochstar,star[j].zstar);
        }
    }

    fclose(PAR);
    fclose(NBODY);
    if (bin == 5) {
        fclose(SSE12);
        printf("\nData written to %s, %s and %s\n",
            PARfile, NBODYfile, SSEfile);
    } else {
        printf("\nData written to %s and %s\n", PARfile, NBODYfile);
    }

    return 0;

}

int output1(struct options opt, double rvir, double mmean, int bin, double M, double MMAX,
    double rtide, vector<struct star_data> star) {

    //Open output files
    char PARfile[50], NBODYfile[50];
    FILE *PAR, *NBODY;
    sprintf(PARfile, "%s.PAR",opt.output);
    PAR = fopen(PARfile,"w");
    sprintf(NBODYfile, "%s.NBODY",opt.output);
    NBODY = fopen(NBODYfile,"w");

    //write to .PAR file
    fprintf(PAR,"1 5000000.0 0\n");
    fprintf (PAR,"%i 1 10 3 8\n",opt.N);
    fprintf(PAR,"0.02 %.8lf %.8lf %.8f 100.0 1000.0 1.0E-02 %.8lf %.8f\n",
        opt.dtadj,opt.dtout,opt.tcrit,rvir,mmean);
    fprintf(PAR,"2 2 1 0 1 0 2 0 0 2\n");
    fprintf(PAR,"0 0 0 %i 2 %i %i 0 %i 2\n",opt.tf,opt.regupdate,opt.etaupdate,opt.mloss);
    fprintf(PAR,"0 %i %i 0 1 2 0 0 0 2\n",bin, opt.esc);
    fprintf(PAR,"0 0 0 0 1 0 0 2 0 3\n");
    fprintf(PAR,"0 0 0 0 0 0 0 0 0 0\n");
    fprintf(PAR,"%.8lf %.8lf %.8lf\n",M,mlow,opt.mup);
    fprintf(PAR,"1.0E-5 1.0E-4 0.01 1.0 1.0E-06 0.01\n");
    fprintf(PAR,"2.350000 %.8lf %.8lf %i %.8lf %.8f 100000.0\n",
        MMAX,mlow, opt.nbin, opt.Z, opt.epoch);
    fprintf(PAR,"%.2f 0.0 0.0 0.00000\n",opt.Q);
    if (opt.tf == 1) {
        rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*
            sqrt(opt.RG[0]*opt.RG[0]+opt.RG[1]*opt.RG[1]+opt.RG[2]*opt.RG[2]);
        fprintf(PAR,"%i %.8f\n",0,sqrt(1.0/(3.0*pow(rtide,3))));

    } else if (opt.tf == 2) {
        fprintf(PAR,"%.8e %.8f\n",M1pointmass,
            sqrt(opt.RG[0]*opt.RG[0]+opt.RG[1]*opt.RG[1]+opt.RG[2]*opt.RG[2]));
    } else if (opt.tf == 3) {
        fprintf(PAR,"%.6e %.6lf %.6e %.6lf %.6lf %.6e %.6lf %.6lf %.6lf %.6lf %.6f "
            "%.6lf %.6f\n",M1allen,b1allen,M2allen,a2allen,b2allen,M3allen,
            a3allen,opt.RG[0]/1000.0,opt.RG[1]/1000.0,opt.RG[2]/1000.0,opt.VG[0],opt.VG[1],opt.VG[2]);
    }
    if (opt.gpu) fprintf(PAR,"1.0\n");

    //write to .NBODY file
    int j;
    for (j=0;j<opt.N;j++) {
        fprintf(NBODY,"%.16lf\t%.16lf %.16lf %.16lf\t%.16lf %.16lf %.16lf\n",
            star[j].mass_epoch,star[j].rx,star[j].ry,star[j].rz,star[j].vx,star[j].vy,
            star[j].vz);
    }

    fclose(PAR);
    fclose(NBODY);
    printf("\nData written to %s and %s\n", PARfile, NBODYfile);

    return 0;

}

int output2(struct options opt, int NNBMAX, double RS0,  double rvir, double mmean,
    int bin, double M, double MMAX, double rtide, vector<struct star_data> star, int sse) {

    //Open output files
    char PARfile[50], NBODYfile[50], SSEfile[50];
    FILE *PAR=NULL, *NBODY=NULL, *SSE12=NULL;
    sprintf(PARfile, "%s.PAR",opt.output);
    PAR = fopen(PARfile,"w");
    sprintf(NBODYfile, "%s.NBODY",opt.output);
    NBODY = fopen(NBODYfile,"w");

    int hrplot = 0;
    if (opt.dtplot) hrplot = 1;
    if (sse) {
        sprintf(SSEfile, "%s.fort.12",opt.output);
        SSE12 = fopen(SSEfile,"w");
        hrplot = 2;
    }


    //write to .PAR file
    fprintf(PAR,"1 5000000.0 0\n");
    fprintf(PAR,"%i 1 10 %i %i 1\n",opt.N,opt.seed,NNBMAX);
    fprintf(PAR,"0.02 0.02 %.8lf %.8lf %.8lf %.8f 1.0E-03 %.8lf %.8f\n",
        RS0,opt.dtadj,opt.dtout,opt.tcrit,rvir,mmean);
    fprintf(PAR,"2 2 1 0 1 0 2 0 0 2\n");
    fprintf(PAR,"0 %i 0 %i 2 %i %i 0 %i 3\n",
        hrplot,opt.tf,opt.regupdate,opt.etaupdate,opt.mloss);
    fprintf(PAR,"0 %i %i 0 1 2 0 1 0 1\n",bin, opt.esc);
    fprintf(PAR,"0 0 0 2 1 0 0 2 0 3\n");
    fprintf(PAR,"0 0 0 0 0 0 0 0 0 0\n");
    fprintf(PAR,"1.0E-5 1.0E-4 0.2 1.0 1.0E-06 0.001\n");
    fprintf(PAR,"2.350000 %.8lf %.8lf %i 0 %.8lf %.8lf %.8f\n",
        MMAX,mlow,opt.nbin,opt.Z,opt.epoch,opt.dtplot);
    fprintf(PAR,"%.2f 0.0 0.0 0.00000 0.125\n",opt.Q);
    // if (opt.tf == 1) {
    //     rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*
    //          sqrt(opt.RG[0]*opt.RG[0]+opt.RG[1]*opt.RG[1]+opt.RG[2]*opt.RG[2]);
    //     fprintf(PAR,"%i %.8f\n",0,sqrt(1.0/(3.0*pow(rtide,3))));
    if (opt.tf == 2) {
        fprintf(PAR,"%.8e %.8f\n",M1pointmass,
            sqrt(opt.RG[0]*opt.RG[0]+opt.RG[1]*opt.RG[1]+opt.RG[2]*opt.RG[2])/1000.0);
    } else if (opt.tf == 3) {
        //fprintf(PAR,"%.6e %.6lf %.6e %.6lf %.6lf %.6e %.6lf %.6lf %.6lf %.6f "
        //  "%.6lf %.6lf %.6f\n",M1allen,b1allen,M2allen,a2allen,b2allen,
        //  M3allen,a3allen,opt.RG[0]/1000.0,opt.RG[1]/1000.0,opt.RG[2]/1000.0,opt.VG[0],
        //  VG[1],opt.VG[2]);
        //new version including bulge potential:
        fprintf(PAR,"%.6e %.6lf %.6e %.6lf %.6lf %.6e %.6lf %.6lf %.6lf %.6f\n",
            M1allen,b1allen,M2allen,a2allen,b2allen,M3allen,a3allen, GMB,
            AR, GAM);
        fprintf(PAR,"%.6lf %.6lf %.6lf %.6lf %.6lf %.6f\n",
            opt.RG[0]/1000.0,opt.RG[1]/1000.0,opt.RG[2]/1000.0,opt.VG[0],opt.VG[1],opt.VG[2]);
    }
    if (opt.tf > 2) {
        fprintf(PAR,"%.2lf %.2lf %.2lf %.2f\n", opt.extmass,opt.extrad,opt.extdecay,opt.extstart);
    }

    //write to .NBODY file
    int j;
    for (j=0;j<opt.N;j++) {
        fprintf(NBODY,"%.16lf\t%.16lf %.16lf %.16lf\t%.16lf %.16lf %.16lf\n",
            star[j].mass_epoch,star[j].rx,star[j].ry,star[j].rz,star[j].vx,star[j].vy,
            star[j].vz);
    }

    //write to .fort.12 file
    if (sse) {
        for (j=0;j<opt.N;j++) {
            fprintf(SSE12,"%.8lf\t%.0lf %.8lf %.8lf %.8lf\n",
                star[j].mass_epoch*M,star[j].kstar,star[j].mass,star[j].epoch1,star[j].ospin);
                    //,star[j].epochstar,star[j].zstar);
        }
    }

    fclose(PAR);
    fclose(NBODY);
    if (bin == 5) {
        fclose(SSE12);
        printf("\nData written to %s, %s and %s\n",
            PARfile, NBODYfile, SSEfile);
    } else {
        printf("\nData written to %s and %s\n", PARfile, NBODYfile);
    }
    return 0;

}

int output3(struct options opt, double rvir, double mmean, double M, double rtide,
    vector<struct star_data> star) {

    //Open output files
    char tablefile[20];
    FILE *TABLE;
    sprintf(tablefile, "%s.txt",opt.output);
    TABLE = fopen(tablefile,"w");


    //write to .txt file
    int i, j;
    //fprintf(TABLE,"# Star cluster generated by McLuster with the following "
    //  "parameters:\n");
    //fprintf(TABLE,"#\n");
    //fprintf(TABLE,"# Mass = %.3lf\t(%i stars with mean mass olf %.3f)\n",
    //  M,opt.N,mmean);
    //fprintf(TABLE,"# Half-mass radius, virial radius, tidal radius = "
    //  "%.3lf pc, %.3lf pc, %.3lf pc\n",rh, rvir, rtide);
    //fprintf(TABLE,"# Galactic orbit (opt.RGx,opt.RGy,opt.RGz),(VGx,opt.VGy,opt.VGz) = "
    //  "(%.3lf,%.3lf,%.3lf) pc, (%.3lf,%.3lf,%.3lf) km/s\n",
    //  opt.RG[0],opt.RG[1],opt.RG[2],opt.VG[0],opt.VG[1],opt.VG[2]);
    //fprintf(TABLE,"# Age = %.3lf Myr\n",opt.epoch);
    //fprintf(TABLE,"# Metallicity = %.3lf\n",opt.Z);
    //fprintf(TABLE,"#\n");


    double epot, gaspot, vscale;
    epot = 0.0;
    gaspot = 0.0;

    //rescale velocities to include effect of gas potential
    if (opt.extmass) {
    #ifndef NOOMP
    int Ntemp = opt.N;
    #pragma omp parallel shared(Ntemp, star)  private(i, j)
        {
        #pragma omp for reduction(+: epot, gaspot) schedule(dynamic)
        #endif
        for (j=0; j<opt.N; j++) {
            if (j) {
                //cluster binding energy
                for (i=0;i<j-1;i++)
                    epot -= star[i].mass_epoch*star[j].mass_epoch/sqrt((star[i].rx-
                        star[j].rx)*(star[i].rx-star[j].rx)+(star[i].ry-
                        star[j].ry)*(star[i].ry-star[j].ry)+(star[i].rz-
                        star[j].rz)*(star[i].rz-star[j].rz));
                //external gas potential binding energy
                gaspot -= star[j].mass_epoch/opt.extrad*opt.extmass/sqrt(1.0+(star[j].rx*
                    star[j].rx+star[j].ry*star[j].ry+star[j].rz*
                    star[j].rz)/(opt.extrad*opt.extrad));
            }
        }
        #ifndef NOOMP
        }
        #endif

        vscale = sqrt(1.0+0.5*gaspot/epot);

        for (j=0; j<opt.N; j++) {
            star[j].vx *= vscale;
            star[j].vy *= vscale;
            star[j].vz *= vscale;
        }

        printf("\nVelocities rescaled to virial equilibrium in Plummer "
            "background potential with mass olf %.2g Msun and a Plummer "
            "radius olf %2.f pc.\n", opt.extmass, opt.extrad);
    }

    #ifdef SSE
    double abvmag, vmag, BV, Teff, dvmag, dBV;
    fprintf(TABLE,"#Mass_[Msun]\tx_[pc]\t\t\ty_[pc]\t\t\tz_[pc]\t\t\t"
        "vx_[km/s]\t\tvy_[km/s]\t\tvz_[km/s]\t\tMass(t=0)_[Msun]\tkw\t"
        "epoch1_[Myr]\tspin\t\tR_[Rsun]\tL_[Lsun]\tepoch_[Myr]\tZ\t\t"
        "M_V_[mag]\tV_[mag]\t\tB-V_[mag]\tT_eff_[K]\tdV_[mag]\td(B-V)\n");
    for (j=0;j<opt.N;j++) {
        cmd(star, j, Rgal, &abvmag, &vmag, &BV, &Teff, &dvmag, &dBV);
        fprintf(TABLE,"%8lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t"
            "%8lf\t\t%.0lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\t"
            "%8lf\t%8lf\t%8lf\t%8lf\n",star[j].mass_epoch,star[j].rx,star[j].ry,
            star[j].rz,star[j].vx,star[j].vy,star[j].vz,star[j].mass,star[j].kstar,
            star[j].epoch1,star[j].ospin,star[j].rstar,star[j].lstar,star[j].epochstar,
            star[j].zstar, abvmag, vmag, BV, Teff, dvmag, dBV);
    }
    #else
    fprintf(TABLE,"#Mass_[Msun]\tx_[pc]\t\t\ty_[pc]\t\t\tz_[pc]\t\t\t"
        "vx_[km/s]\t\tvy_[km/s]\t\tvz_[km/s]\n");
    for (j=0;j<opt.N;j++) {
        fprintf(TABLE,"%8lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",
            star[j].mass_epoch,star[j].rx,star[j].ry,star[j].rz,star[j].vx,star[j].vy,
            star[j].vz);
    }
    #endif
    fclose(TABLE);
    printf("\nData written to %s\n", tablefile);

    return 0;

}

int output4(struct options opt, int NNBMAX, double RS0, double rvir, double mmean,
    int bin, double M, double MMAX, double rtide, vector<struct star_data> star, int sse) {

    //Open output files
    char PARfile[50], NBODYfile[50], SSEfile[50];
    FILE *PAR=NULL, *NBODY=NULL, *SSE12=NULL;
    sprintf(PARfile, "%s.input",opt.output);
    PAR = fopen(PARfile,"w");
    sprintf(NBODYfile, "%s.fort.10",opt.output);
    NBODY = fopen(NBODYfile,"w");

    int hrplot = 0;
    if (opt.dtplot) hrplot = 1;
    if (sse) {
        sprintf(SSEfile, "%s.fort.12",opt.output);
        SSE12 = fopen(SSEfile,"w");
        hrplot = 2;
    }

    //write to .PAR file
    fprintf(PAR,"1 5000000.0 0\n");
    fprintf(PAR,"%i 1 10 %i %i 1\n",opt.N,opt.seed,NNBMAX);
    fprintf(PAR,"0.02 0.02 %.8lf %.8lf %.8lf %.8f 1.0E-03 %.8lf %.8f\n",
        RS0,opt.dtadj,opt.dtout,opt.tcrit,rvir,mmean);
    fprintf(PAR,"2 2 1 0 1 0 2 0 0 2\n");
    fprintf(PAR,"-1 %i 0 %i 2 %i %i 0 %i 3\n",
        hrplot,opt.tf,opt.regupdate,opt.etaupdate,opt.mloss);
    fprintf(PAR,"0 %i %i 0 1 2 0 1 0 -1\n",bin,opt.esc);
    fprintf(PAR,"0 0 0 2 1 0 0 2 0 3\n");
    fprintf(PAR,"0 1 0 1 0 0 0 0 0 0\n");
    fprintf(PAR,"1.0E-5 1.0E-4 0.2 1.0 1.0E-06 0.001\n");
    fprintf(PAR,"2.350000 %.8lf %.8lf %i 0 %.8lf %.8lf %.8lf\n",
        MMAX, mlow, opt.nbin, opt.Z, opt.epoch, opt.dtplot);
    fprintf(PAR,"%.2f 0.0 0.0 0.00000 0.125\n",opt.Q);
    if (opt.tf == 2) {
        fprintf(PAR,"%.8e %.8f\n",M1pointmass,sqrt(opt.RG[0]*opt.RG[0]+opt.RG[1]*
            opt.RG[1]+opt.RG[2]*opt.RG[2])/1000.0);
    } else if (opt.tf == 3) {
        //old version:
        //fprintf(PAR,"%.6e %.6e %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6f "
        //  "%.6lf %.6f\n",M1allen,M2allen,a2allen,b2allen,VCIRC,RCIRC,
        //  opt.RG[0]/1000.0,opt.RG[1]/1000.0,opt.RG[2]/1000.0,opt.VG[0],opt.VG[1],opt.VG[2]);

        //new version including bulge potential:
        fprintf(PAR,"%.6e %.6e %.6lf %.6lf %.6lf %.6lf %.6e %.6lf %.6f\n",
            M1allen,M2allen,a2allen,b2allen,VCIRC,RCIRC, GMB, AR, GAM);
        fprintf(PAR,"%.6lf %.6lf %.6lf %.6lf %.6lf %.6f\n",
            opt.RG[0]/1000.0,opt.RG[1]/1000.0,opt.RG[2]/1000.0,opt.VG[0],opt.VG[1],opt.VG[2]);
    }
    if (opt.tf > 2)
        fprintf(PAR,"%.6lf %.6lf %.6lf %.6f\n",opt.extmass,opt.extrad,opt.extdecay,opt.extstart);

    fprintf(PAR,"10000.0 2 0\n");

    //write to .fort.10 file
    int j;
    for (j=0;j<opt.N;j++) {
        fprintf(NBODY,"%.16lf\t%.16lf %.16lf %.16lf\t%.16lf %.16lf %.16lf\n",
            star[j].mass_epoch,star[j].rx,star[j].ry,star[j].rz,star[j].vx,star[j].vy,
            star[j].vz);
    }

    //write to .fort.12 file
    if (sse) {
        for (j=0;j<opt.N;j++) {
            fprintf(SSE12,"%.8lf\t%.0lf %.8lf %.8lf %.8lf\n",
                star[j].mass_epoch*M,star[j].kstar,star[j].mass,star[j].epoch1,star[j].ospin);
            //,star[j].epochstar,star[j].zstar);
        }
    }

    fclose(PAR);
    fclose(NBODY);
    if (bin == 5) {
        fclose(SSE12);
        printf("\nData written to %s, %s and %s\n",
            PARfile, NBODYfile, SSEfile);
    } else {
        printf("\nData written to %s and %s\n", PARfile, NBODYfile);
    }

    return 0;

}

int output5(struct options opt, int NNBMAX, double RS0, double rvir, double mmean,
    int bin, double M, double MMAX,  double rtide, vector<struct star_data> star, int sse) {

    //Open output files
    char PARfile[50], NBODYfile[50], SSEfile[50];
    FILE *PAR=NULL, *NBODY=NULL, *SSE12=NULL;
    sprintf(PARfile, "%s.input", opt.output);
    PAR = fopen(PARfile,"w");
    sprintf(NBODYfile, "%s.dat.10", opt.output);
    NBODY = fopen(NBODYfile,"w");

    int hrplot = 0;
    if (opt.dtplot) hrplot = 1;
    if (sse) {
        sprintf(SSEfile, "%s.fort.12",opt.output);
        SSE12 = fopen(SSEfile,"w");
        hrplot = 2;
    }


    //write to .PAR file
    fprintf(PAR,"1 5000000.0 5000000.0 40 40 0\n");
    fprintf(PAR,"%i 1 10 %i %i 1\n",opt.N,opt.seed,NNBMAX);
    fprintf(PAR,"0.02 0.02 %.8lf %.8lf %.8lf %.8f 1.0E-03 %.8lf %.8f\n",
        RS0,opt.dtadj,opt.dtout,opt.tcrit,rvir,mmean);
    fprintf(PAR,"0 2 1 0 1 0 5 %i 3 2\n",(opt.nbin>0?2:0));
    fprintf(PAR,"0 %i 0 %i 2 %i %i 0 %i 6\n",
        hrplot,opt.tf,opt.regupdate,opt.etaupdate,opt.mloss);
    fprintf(PAR,"0 6 %i 0 1 2 2 0 0 1\n", opt.esc);
    fprintf(PAR,"1 0 3 2 1 0 0 2 0 0\n");
    fprintf(PAR,"0 0 0 0 0 1 2 0 0 0\n");
    fprintf(PAR,"1.0E-5 1.0E-4 0.2 1.0 1.0E-06 0.001 0.125\n");
    fprintf(PAR,"2.350000 %.8lf %.8lf %i 0 %.8lf %.8lf %.8f\n",
        MMAX,opt.mlow,opt.nbin,opt.Z,opt.epoch,opt.dtplot);
    fprintf(PAR,"%.2f 0.0 0.0 0.00000\n",opt.Q);
//    if (opt.tf == 1) {
//        rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*
//          sqrt(opt.RG[0]*opt.RG[0]+opt.RG[1]*opt.RG[1]+opt.RG[2]*opt.RG[2]);
//        fprintf(PAR,"%i %.8f\n",0,sqrt(1.0/(3.0*pow(rtide,3))));
    if (opt.tf == 2) {
        fprintf(PAR,"%.8e %.8f\n",M1pointmass,
            sqrt(opt.RG[0]*opt.RG[0]+opt.RG[1]*opt.RG[1]+opt.RG[2]*opt.RG[2])/1000.0);
    } else if (opt.tf == 3) {
        //fprintf(PAR,"%.6e %.6lf %.6e %.6lf %.6lf %.6e %.6lf %.6lf %.6lf %.6f "
        //  "%.6lf %.6lf %.6f\n",M1allen,b1allen,M2allen,a2allen,b2allen,M3allen,
        //  a3allen,opt.RG[0]/1000.0,opt.RG[1]/1000.0,opt.RG[2]/1000.0,opt.VG[0],opt.VG[1],opt.VG[2]);
        //new version including bulge potential:
        fprintf(PAR,"%.6e %.6lf %.6e %.6lf %.6lf %.6e %.6lf %.6lf %.6lf %.6f\n",
            M1allen,b1allen,M2allen,a2allen,b2allen,M3allen,a3allen,
            GMB, AR, GAM);
        fprintf(PAR,"%.6lf %.6lf %.6lf %.6lf %.6lf %.6f\n",
            opt.RG[0]/1000.0,opt.RG[1]/1000.0,opt.RG[2]/1000.0,opt.VG[0],opt.VG[1],opt.VG[2]);
    }
    if (opt.tf > 2)
        fprintf(PAR,"%.2lf %.2lf %.2lf %.2f\n",opt.extmass,opt.extrad,opt.extdecay,opt.extstart);

    //write to .NBODY file
    int j;
    for (j=0;j<opt.N;j++) {
        fprintf(NBODY,"%.16lf\t%.16lf %.16lf %.16lf\t%.16lf %.16lf %.16lf\n",
            star[j].mass_epoch,star[j].rx,star[j].ry,star[j].rz,star[j].vx,star[j].vy,
            star[j].vz);
    }

    //write to .fort.12 file
    if (sse) {
        for (j=0;j<opt.N;j++) {
            fprintf(SSE12,"%.8lf\t%.0lf %.8lf %.8lf %.8lf\n",
                star[j].mass_epoch*M,star[j].kstar,star[j].mass,star[j].epoch1,star[j].ospin);
                //,star[j].epochstar,star[j].zstar);
        }
    }

    fclose(PAR);
    fclose(NBODY);
    if (bin == 5) {
        fclose(SSE12);
        printf("\nData written to %s, %s and %s\n", PARfile, NBODYfile, SSEfile);
    } else {
        printf("\nData written to %s and %s\n", PARfile, NBODYfile);
    }
    return 0;

}

void info(options opt) {

    int i;
    char INFOfile[50];
    FILE *INFO;
    sprintf(INFOfile, "%s.info",opt.output);
    INFO = fopen(INFOfile,"w");

    time_t timep;
    time(&timep);
    fprintf(INFO, "\nMcLuster model generated: %s\n\n",ctime(&timep));

    fprintf(INFO, "N = %i\n",opt.N);
    fprintf(INFO, "Mcl = %g\n",opt.Mcl);
    fprintf(INFO, "profile = %i\n", opt.profile);
    fprintf(INFO, "W0 = %g\n",opt.W0);
    fprintf(INFO, "S = %g\n",opt.S);
    fprintf(INFO, "D = %g\n",opt.D);
    fprintf(INFO, "Q = %g\n",opt.Q);
    fprintf(INFO, "Rh = %g\n",opt.Rh);
    fprintf(INFO, "gammas = %g\t%g\t%g\n",opt.gamma[0],opt.gamma[1],opt.gamma[2]);
    fprintf(INFO, "a = %g\n",opt.a);
    fprintf(INFO, "Rmax = %g\n",opt.Rmax);
    fprintf(INFO, "tcrit = %g\n",opt.tcrit);
    fprintf(INFO, "tf = %i\n",opt.tf);
    fprintf(INFO, "RG = %g\t%g\t%g\n",opt.RG[0],opt.RG[1],opt.RG[2]);
    fprintf(INFO, "VG = %g\t%g\t%g\n",opt.VG[0],opt.VG[1],opt.VG[2]);
    fprintf(INFO, "mfunc = %i\n",opt.mfunc);
    fprintf(INFO, "single_mass = %g\n",opt.single_mass);
    fprintf(INFO, "mlow = %g\n",opt.mlow);
    fprintf(INFO, "mup = %g\n",opt.mup);
    fprintf(INFO, "alpha = ");
    for(i=0;i<MAX_AN;i++) fprintf(INFO, "%g\t",opt.alpha[i]);
    fprintf(INFO, "\n");
    fprintf(INFO, "mlim = ");
    for(i=0;i<MAX_MN;i++) fprintf(INFO, "%g\t",opt.mlim[i]);
    fprintf(INFO, "\n");
    fprintf(INFO, "alpha_L3 = %g\n",opt.alpha_L3);
    fprintf(INFO, "beta_L3 = %g\n",opt.beta_L3);
    fprintf(INFO, "mu_L3 = %g\n",opt.mu_L3);
    fprintf(INFO, "weidner = %i\n",opt.weidner);
    fprintf(INFO, "mloss = %i\n",opt.mloss);
    fprintf(INFO, "remnant = %i\n",opt.remnant);
    fprintf(INFO, "epoch = %g\n",opt.epoch);
    fprintf(INFO, "Z = %g\n",opt.Z);
    fprintf(INFO, "prantzos = %i\n",opt.prantzos);
    fprintf(INFO, "nbin = %i\n",opt.nbin);
    fprintf(INFO, "fbin = %g\n",opt.fbin);
    fprintf(INFO, "pairing = %i\n",opt.pairing);
    fprintf(INFO, "msort = %g\n",opt.msort);
    fprintf(INFO, "adis = %i\n",opt.adis);
    fprintf(INFO, "OBperiods = %i\n",opt.OBperiods);
    fprintf(INFO, "amin = %g\n",opt.amin);
    fprintf(INFO, "amax = %g\n",opt.amax);
    fprintf(INFO, "eigen = %i\n",opt.eigen);
    fprintf(INFO, "BSE = %i\n",opt.BSE);
    fprintf(INFO, "extmass = %g\n",opt.extmass);
    fprintf(INFO, "extrad = %g\n",opt.extrad);
    fprintf(INFO, "extdecay = %g\n",opt.extdecay);
    fprintf(INFO, "extstart = %g\n",opt.extstart);
    fprintf(INFO, "code = %i\n",opt.code);
    fprintf(INFO, "seed = %i\n",opt.seed);
    fprintf(INFO, "dtadj = %g\n",opt.dtadj);
    fprintf(INFO, "dtout = %g\n",opt.dtout);
    fprintf(INFO, "dtplot = %g\n",opt.dtplot);
    fprintf(INFO, "gpu = %i\n",opt.gpu);
    fprintf(INFO, "regupdate = %i\n",opt.regupdate);
    fprintf(INFO, "etaupdate = %i\n",opt.etaupdate);
    fprintf(INFO, "esc = %i\n",opt.esc);
    fprintf(INFO, "units = %i\n",opt.units);
    fprintf(INFO, "match = %i\n",opt.match);
    fprintf(INFO, "symmetry = %i\n",opt.symmetry);

    fclose(INFO);

}

void help(double msort) {
    printf("\n Usage: mcluster -[N|M|P|W|R|r|c|g|S|D|T|Q|C|A|O|G|o|f|a|m|B|b|p|s|t|e|Z|X|V|x|u|h|?]\n");
    printf("                                                                     \n");
    printf("       -N <number> (number of stars)                                 \n");
    printf("       -M <value> (mass of cluster; specify either N or M)           \n");
    printf("       -P <0|1|2|3|-1> (density profile; 0= Plummer, 1= King (1966), \n");
    printf("                   2= Subr et al. (2007) mass-segregated,            \n");
    printf("                   3= 2-dimensional EFF/Nuker template,              \n");
    printf("                   -1= no density gradient)                          \n");
    printf("       -W <1-12> (W0 parameter for King model)                       \n");
    printf("       -R <value> (half-mass radius [pc], ignored for P = 3)         \n");
    printf("       -r <value> (scale radius of EFF/Nuker template [pc])          \n");
    printf("       -c <value> (cut-off radius of EFF/Nuker template [pc])        \n");
    printf("       -g <value> (power-law slope(s) of EFF/Nuker template; use     \n");
    printf("                   once for EFF template; use three times for Nuker  \n");
    printf("                   template (outer slope, inner slope, transition)   \n");
    printf("       -S <0.0-1.0> (degree of mass segregation; 0.0= no segregation)\n");
    printf("       -D <1.6-3.0> (fractal dimension; 3.0= no fractality)          \n");
    printf("       -T <value> (tcrit in N-body units)                            \n");
    printf("       -Q <value> (virial ratio)                                     \n");
    printf("       -C <0|1|3|5> (code; 0= Nbody6, 1= Nbody4, 3= table of stars, 5= Nbody6++)    \n");
    printf("       -A <value> (dtadj in N-body units)                            \n");
    printf("       -O <value> (deltat in N-body units)                           \n");
    printf("       -G <0|1> (gpu usage; 0= no GPU, 1= use GPU)                   \n");
    printf("       -o <name> (output name of cluster model)                      \n");
    printf("       -f <0|1|2|3|4> (IMF; 0= no IMF, 1= Kroupa (2001),             \n");
    printf("             2= user defined, 3= Kroupa (2001) with optimal sampling,\n");
    printf("             4= L3 IMF (Maschberger 2012))                           \n");
    printf("       -a <value> (IMF slope; for user defined IMF, may be used      \n");
    printf("                   multiple times, from low mass to high mass;       \n");
    printf("                   for L3 IMF use three times for alpha, beta and mu)\n");
    printf("       -m <value> (IMF mass limits, has to be used multiple times    \n");
    printf("                 (at least twice), from low mass to high mass [Msun])\n");
    printf("       -B <number> (number of binary systems)                        \n");
    printf("       -b <value> (binary fraction, specify either B or b)           \n");
    printf("       -p <0|1|2|3> (binary pairing, 0= random, 1= ordered for M>%.1f Msun,\n",msort);
    printf("                   2= random but separate pairing for M>%.1f Msun)\n",msort);
    printf("                   3= random but uniform distribution of mass ratio (0.1<q<1.0) for M>%.1f Msun)\n",msort);
    printf("       -s <number> (seed for randomization; 0= randomize by timer)   \n");
    printf("       -t <0|1|2|3> (tidal field; 0= no tidal field, 1= near-field,  \n");
    printf("                    2= point-mass, 3= Milky-Way potential)           \n");
    printf("       -e <value> (epoch for stellar evolution [Myr])                \n");
    printf("       -Z <value> (metallicity [0.0001-0.03, 0.02 = solar])          \n");
    printf("       -X <value> (galactocentric radius vector, use 3x, [pc])       \n");
    printf("       -V <value> (cluster velocity vector, use 3x, [km/s])          \n");
    printf("       -x <value> (specify external (gas) Plummer potential, use 4x, \n");
    printf("                  1. gas mass [Msun], 2. Plummer radius [pc]         \n");
    printf("                  3. decay time for gas expulsion [Myr], 4. delay    \n");
    printf("                  time for start of gas expulsion [Myr])             \n");
    printf("       -u <0|1> (output units; 0= Nbody, 1= astrophysical)           \n");
    printf("       -h (display this help)                                        \n");
    printf("       -? (display this help)                                        \n");
    printf("                                                                     \n");
    printf(" Examples: mcluster -N 1000 -R 0.8 -P 1 -W 3.0 -f 1 -B 100 -o test1  \n");
    printf("           mcluster -f 2 -m 0.08 -a -1.35 -m 0.5 -a -2.7 -m 100.0    \n");
    printf("           mcluster -t 3 -X 8500 -X 0 -X 0 -V 0 -V 220 -V 0          \n");
    printf("           mcluster -D 1.6 -Q 0.4 -P -1                              \n");
    printf("                                                                     \n");

}
