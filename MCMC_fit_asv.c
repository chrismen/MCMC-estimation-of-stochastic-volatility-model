

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <time.h>
#include <math_recipes.h>
#include "normsinv.h"
#include "normsinv.cpp"

double Z(const double x);
double CDF(const double x);
double min(double x, double y);
double max(double x, double y);
double normpdf(double x, double mu, double stdv);

//----------------------------------------------
int main ()
{
    FILE * outputfile;
    FILE * input;
    FILE *para_output;
    const int N =5000;
    const int T =1500;
    const int n=N/5;
// the burn in number of the measurement dynamic to be stable
    const int k1=50000;
// simulated values of the measurement dynamics
    float y [T] = {};
    float phi  ;
// store the variances for the latent disturbances
    float sigma;
    float tau;
    float psi;
    float rho;
    float mu;
// some temp variables
    float a ,b ,a1, b1,a2,b2,a3,b3,a4,b4,a11,a1_right,b11,b1_left,c,c1,c1_left,c1_right,rnd1, latent, bounder,left,right;
    float u1, u2,u3,d;
    float z[T]={};

    float alpha,alpha_t, phi_t, alpha_star, u, rnd , phi_temp, chi_temp,temp1,temp;
    int j,k,h,i,g,burnin,m,count_phi,sign;

// the followings are used for the final calculation
    float phi_mean=0.0;
    float tau_mean=0.0;
    float psi_mean=0.0;
    float mu_mean=0;
    float rho_mean=0.0;
    float sigma_mean=0.0;
    float norm_temp, norm1_temp;
    float unif_temp;
    float phi_para, sigma_para,rho_para,mu_para;
    float mean;
    float mu1, sig,p;
    float alpha_mu, sig_mu, alpha_phi, sig_phi, alpha_psi, sig_psi, sig_tau;
    int alpha_tau;

    alpha_mu=0.0 ;
    sig_mu=10.0;
    alpha_phi=0.0;
    sig_phi=10.0;
    alpha_psi=0.0;
    sig_psi=10.0;
    alpha_tau=2  ;
    sig_tau=0.0025;


// initialized the random generation seed
    long a_time, b_time;
    time_t t1;

    (void) time(&t1);
    a_time=(long) t1;
    b_time=-(a_time+1);
// the initial values for the start of the Gibbs
    phi=0.7;
    sigma=1.5; //rooted
    rho=-0.2;
    mu=0.8;
    phi_para=0.5;
    rho_para=-0.5;
    sigma_para=1.0;//rooted
    mu_para=0.5;

    para_output = fopen ( "MCMC-time-series.txt", "w");

//-------------------------------------------------------------
    input = fopen ("ibm.txt", "r");
    h=0;// the time of the repeating sampling//
    while (h< T)
    {
        fscanf (input, "%f", y+h);
        h ++;
    }

    fclose (input);

    norm_temp=gasdev(&b_time);
    z[0]=mu + sigma/sqrt(1.0-phi*phi)*norm_temp;

    for ( j=1; j< T;  j++)
    {
        norm_temp=gasdev(&b_time);
        z[j]=mu + phi * (z[j-1]-mu)+ sigma*norm_temp;
    }

    tau=sigma*sigma*(1.0-pow(rho,2.0));
    psi=rho*sigma;
// sample for all random variable in the model
    for ( k=1;k<=N; k++)
    {
        cout<< k<<"\n";
        //generating the first value of x(1) for the Gibbs
        u1=ran1(&b_time);
        a1=u1*exp(-z[0]/2.0);
        a1_right=-2.0*log(a1);
        u2=ran1(&b_time);
        b1=u2*exp(-y[0]*y[0]/2.0*exp(-z[0]));

        if (y[0]!=0.0)
        {
            b1_left=-log(-2.0/y[0]/y[0]*log(b1));
        }
        u3=ran1(&b_time);
        d=mu ;
        c1=u3*exp(-pow(z[0]-d,2.0)/2.0/tau/(1.0-pow(phi,2.0)));
        c1_left=d-sqrt(-2.0*tau*log(c1)/(1.0-phi*phi));
        c1_right=d+sqrt(-2.0*tau*log(c1)/(1.0-phi*phi));
        if (y[0]!=0.0)
        {
            left=max(b1_left,c1_left);
            right=min(c1_right,a1_right);
        }
        else

        {
            left=c1_left;
            right=min(c1_right,a1_right);
        }
        u=ran1(&b_time);
        rnd=left + (right-left)*u;
        a1=exp(-pow((z[1] -mu -phi*(rnd-mu)-psi*exp(-rnd/2.0)*y[1]),2.0)/2.0/tau);
        u=ran1(&b_time);
        if (u<=a1)
        {
            z[0]=rnd;
        }
// generate the random number for the random variabes x(2:T-1)
        for (j=1; j<T; j++)
        {
            u1=ran1(&b_time);
            a1=u1*exp(-z[j]/2.0);
            a1_right=-2.0*log(a1);
            if (y[j]!=0.0)
            {
                u2=ran1(&b_time);
                b1=u2*exp(-y[j]*y[j]/2.0*exp(-z[j]));
                b1_left=-log(-2.0/y[j]/y[j]*log(b1));
            }
            u3=ran1(&b_time);
            d=mu + phi*(z[j-1]-mu) + psi*exp(-z[j-1]/2.0)*y[j-1];
            c1=u3*exp(-pow(z[j]-d,2.0)/2.0/tau);
            c1_left=d-sqrt(-2.0*tau*log(c1));
            c1_right=d+sqrt(-2.0*tau*log(c1));
            if (y[j]!=0.0)
            {
                left=max(b1_left,c1_left);
                right=min(c1_right,a1_right);
            }
            else
            {
                left=c1_left;
                right=min(c1_right,a1_right);
            }
            u=ran1(&b_time);
            rnd=left + (right-left)*u;
            a1=exp(-pow((z[j+1] -mu -phi*(rnd-mu)-psi*exp(-rnd/2.0)*y[j]),2.0)/2.0/tau);
            u=ran1(&b_time);
            if (u<=a1)
            {
                z[j]=rnd;
            }
        }  // of j

// we generate the random number for x[T]
        u1=ran1(&b_time);
        a1=u1*exp(-z[T-1]/2.0);
        a1_right=-2.0*log(a1);
        if (y[T-1]!=0.0)
        {
            u2=ran1(&b_time);
            b1=u2*exp(-y[T-1]*y[T-1]/2.0*exp(-z[T-1]));
            b1_left=-log(-2.0/y[T-1]/y[T-1]*log(b1));
        }
        u3=ran1(&b_time);
        d=mu + phi*(z[T-2]-mu) + psi*exp(-z[T-2]/2.0)*y[T-2];
        c1=u3*exp(-pow(z[T-1]-d,2.0)/2.0/tau);
        c1_left=d-sqrt(-2.0*tau*log(c1));
        c1_right=d+sqrt(-2.0*tau*log(c1));
        if (y[T-1]!=0.0)
        {
            left=max(b1_left,c1_left);
            right=min(c1_right,a1_right);
        }
        else
        {
            left=c1_left;
            right=min(c1_right,a1_right);
        }
        u=ran1(&b_time);
        z[T-1]=left + (right-left)*u;

/// generate the random number for the variable mu
        a1=float(T-1)*pow(1.0-phi, 2.0) + (1.0-phi*phi);
        a1=a1/tau;
        a1=a1+1.0/pow(sig_mu, 2.0);
        b1=0.0;
        for (h=0; h<T-1; h++)
        {
            b1=b1+(1.0-phi)*(z[h+1]-phi*z[h]-psi*exp(-z[h]/2.0)*y[h])     ;
        }
        b1=b1+(1.0-phi*phi)*z[0];
        b1=b1/tau;
        b1=b1+ alpha_mu/sig_mu/sig_mu ;
        temp1=gasdev(&b_time);
        mu=b1/a1+1.0/sqrt(a1)*temp1;
/// generate the random number for the variable phi
        u=ran1(&b_time);
        latent=sqrt(1.0-phi*phi)*u;
        bounder=min(1.0, sqrt(1.0-latent*latent));
        a1=0.0;
        for ( h=0;h<=T-2; h++)
            a1=a1+(z[h]-mu)*(z[h]-mu);
        a1=a1-(z[0]-mu)*(z[0]-mu);
        a1=a1/tau;
        a1=a1+1.0/sig_phi/sig_phi;
        b1=0.0;
        for ( h=0; h<T-1; h++)
            b1=b1+(z[h]-mu)*(z[h+1]- mu -psi*exp( -z[h]/2.0 )*y[h]   );
        b1=b1/tau;
        b1=b1+alpha_phi/sig_phi/sig_phi;
        mu1=b1/a1;
        sig=1.0/sqrt(a1);
        u=ran1(&b_time);
        p=CDF((-bounder-mu1)/sig)+(CDF((bounder-mu1)/sig)-CDF((-bounder-mu1)/sig))*u;
        phi=mu1+sig* normsinv(p);
        count_phi=count_phi+1;
/// we generate the sample for disturbance variance of psi
        a1=0.0;
        for ( h=0;h<=T-2; h++)
            a1=a1+pow(exp(-z[h]/2.0)*y[h],2.0);
        a1=a1/tau;
        a1=a1+  1.0/sig_psi/sig_psi;
        b1=0.0;
        for ( h=0; h<T-1; h++)
            b1=b1+exp(-z[h]/2.0)*y[h]*(z[h+1]-mu-phi*(z[h]-mu));
        b1=b1/tau;
        b1=b1+alpha_psi/sig_psi/sig_psi;
        mu1=b1/a1;
        sig=1.0/sqrt(a1);
        temp1=gasdev(&b_time);
        psi=mu1+sig* temp1;
// generate the random numbers for tao

        a1=0.0;
        for ( h=0;h<T-1; h++)
            a1=a1+pow((z[h+1]-mu - phi*(z[h]-mu)-psi*exp(-z[h]/2.0)*y[h]),2.0);

        a1=a1+pow(z[0]-mu, 2.0) * (1.0-phi*phi );
        a1=a1/2.0;
        a1=a1+sig_tau;
        tau=a1/gamdev(T/2+alpha_tau,&b_time);
        sigma=sqrt(psi*psi+ tau);
        rho=psi/sigma;
        fprintf (para_output, " %f %f %f %f \n",mu,phi, rho,sigma);
    } /// the end of the loop k with the length N for gererating samples

    fclose (para_output);
    return 0;
}


///-----------------------------------------------------------------------------------

double f_star(double y, double x)
{
    double fvalue;
    double pi=3.14159265358979323846;
    fvalue=1.0/sqrt(2.0*pi*exp(x))* exp(-y*y /2.0/exp(x));
    return fvalue;
}
///------------------------------------------------------------------------------
double min(double x, double y)
{
    if (x<y) return x;
    else
        return y;
}
///----------------------------------------
double max(double x, double y)
{
    if (x<y) return y;
    else
        return x;
}
//-----------------------------------------------------
double normpdf(double x, double mu, double stdv)
{
    double temp2;
    double pi=3.14159265358979323846;
    temp2=1.0/sqrt(2.0*pi)/stdv*exp(-1.0/2.0/stdv/stdv*(x-mu)*(x-mu));
    return temp2;
}
// Start of the  Abromowitz and Stegun approximation function
double CDF(const double x)
{
    const double b1 =  0.319381530;
    const double b2 = -0.356563782;
    const double b3 =  1.781477937;
    const double b4 = -1.821255978;
    const double b5 =  1.330274429;
    const double p  =  0.2316419;
    if (x >= 0.0)
    {
        double t = 1.0 / (1.0 + p*x);
        return (1.0 - Z(x)*t*
                (t*(t*(t*(t*b5 + b4) + b3) + b2) + b1));
    }
    else
    {
        double t = 1.0 / ( 1.0 - p * x );
        return ( Z(x)*t*
                 (t*(t*(t*(t*b5 + b4) + b3) + b2) + b1));
    }
}
// Functions
// The Gaussian p.d.f with mean = 0 and stddev = 1
double Z(const double x)
{
    double PI = 3.14159265358979323846;
    return (1.0/sqrt(2.0*PI))*exp(-x*x/2.0 );
}
