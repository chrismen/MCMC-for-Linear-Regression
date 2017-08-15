/* -----------------------------------------------------------------------

        This c file do slice samplig for ASV mdoels

     the slice sampler within Gibbs and acceptance-rejection method
       with normal nosies

                            January, 2009
Aug 5, 2009
--------------------------------------------------------------------------*/


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
double f_star(double y,  double x);
double g_star1(double y,  double x);
double g_star2(double y, double x, double z);
double min(double x, double y);
double max(double x, double y);
double normpdf(double x, double mu, double stdv);


double gamma(double x);
//----------------------------------------------
int main ()
{
    FILE * outputfile;
    FILE * input;
    //char     outputfile [] = "sv lev.txt";
    FILE *para_output; // for each g_th sample   "sv lev.txt"
    FILE *para_result; //for g=1 to G         "sv lv estimates.txt"
    FILE *y_output; // the y dynamics
    //  FILE *x_output; // the x dynamics
    // FILE *z_output; // the last record of the Gth  sample
    // FILE *z1_output;



    const int G=1;

// the length of the sampling for each of the individula variable
    const int N =20000;

    const int rows=10000;
    const int cols=3 ;


// the burn in number of the sample
    const int n=N/5;

// some temp variables
    float a1 ,b1 ,c ,d, psi,temp1, temp_current1,temp_new1, temp_current2,temp_new2, u,mu_r,sigma_r;

    float   mu1, sigma1, mu2, sigma2, rho1, r, r_star, low, high, likelihood1, likelihood2,  likelihood_curre1, likelihood_new1,temp_curr1,temp_curr2;

    float x[rows+1][cols+1]={};  //

    // int   y1[rows+1]={};

// remember the states for each loop of Giblls
    float y[rows+1]={};  // the response vector

    float beta1[cols+1]={};
    //  float beta_new1[cols+1]={};

    float beta2[cols+1]={};
    //  float beta_new2[cols+1]={};



    //float beta_star;


    float yy[cols+1]={};

    int j,k,h,i,g,   rows1, rows2,;

// the followings are used for the final calculation


    double pi=3.14159265358979323846;
// Constants and global variables
    double PI = 3.14159265358979323846;

    /*-----------------------------------------------------------------------------

    // initialized the random generation seed
    int seed = (int)time(0);            /// random seed
    StochasticLib1 sto(seed);           /// make instance of random library

    // choose one of the random number generators:
    CRandomMersenne RanGen(seed);       /// make instance of random number generator
    ---------------------------------------------------------------------------------*/

    //  float dic_temp, predict_dic;



// initialized the random generation seed
    long a_time, b_time;
    time_t t1;

    (void) time(&t1);
    a_time=(long) t1;

    b_time=-(a_time+1);



    temp1=gasdev(&b_time);


    /*--------------the code of the body----------------*/

// the initial values for the start of the Gibbs

    //  sigma1=0.15; //rooted
    //  mu1=-0.75;



    mu1=0.0;
    sigma1=0.05;


    mu2=0.75;
    sigma2=0.6;


    mu_r=0.0;
    sigma_r=1.005;




    r=0.0;
    r_star=0.5;

    low=-0.67;
    high=0.67;

    para_output = fopen ( "C-MCMC threshold.txt", "w");


//-------------------------------------------------------------

    input = fopen ("C-simulated-tlinear-data-x.txt", "r");

    h=1;


    for ( j=1; j<=3;  j++)
    {
        beta1[j]=2.0;


        beta2[j]=2.0;


        for ( h=1; h<=rows;  h++)
        {
            fscanf (input, "%f", &rho1);
            x[h][j]=rho1;

            //  cout<< x[h][j]<<" ";
            //   cout<< "\n";
        }

        //


    }

    fclose (input);
//
//    for ( h=1; h<=rows;  h++)
//    {
//
//        for ( j=1; j<=3;  j++)
//        {
//
//
//            cout<< x[h][j]<<" ";
//        }
//
//        cout<< "\n";
//    }
//


// the initial value for the state x_0




    input = fopen ("C-simulated-tlinear-data-y.txt", "r");
    h=1;// the time of the repeating sampling//

    for ( h=1; h<=rows;  h++)
    {


        // fscanf (para_output, " %f %f %f %f \n",&aa,&bb, &cc,&dd);
        fscanf (input, "%f \n", &rho1);
        y[h]=rho1;
//         if ( y[h] ==1.0 )
//            {
//
//
        //   cout<< y[h]<<" ";

        // cout<< "\n";
//            }
    }

    fclose (input);



//
//    for ( h=1; h<=rows;  h++)
//    {
//        //y[h] = y1[h] * 1.0f;
//
//     //  y[h] =y[h];
//
//          cout<< y[h]<<" ";
//
//    }



// sample for all random variable in the model
    for ( g=1;g<=N; g++)
    {

        cout<< g<<"\n ";


        for ( k=1;k<=cols; k++)

        {
            c = 0.0;
            d = 0.0;


            for (j=2; j<=rows; j++)

            {
                if ( x[j-1][2]>r )

                {

                    c = c + pow(x[j][k],2.0) / sigma1/sigma1;
                    psi=0.0;

                    for ( i=1;i<=cols; i++)

                    {
                        if ( i!=k)
                        {
                            psi=psi+beta1[i]*x[j][i];
                        }
                    }


                    d = d +  (y[j] -psi - mu1) * x[j][k] / sigma1/sigma1;
                }
            }


            temp1=gasdev(&b_time);

            beta1[k]=d/c+ 1.0/sqrt(c)*temp1;

        }
//        beta[1]=2.0;
//        beta[2]=3.0;
//        beta[3]=4.0;

/// generate the random number for the variable mu1


        b1=0.0;
        rows1=0;


        for (j=2; j<=rows; j++)

        {
            if ( x[j-1][2]>r )

            {
                rows1=rows1+1;
                psi=0.0;

                for ( k=1;k<=cols; k++)

                {
                    psi=psi+beta1[k]*x[j][k];
                }

                b1 = b1 +  (y[j] - psi) / sigma1/sigma1;


            }
        }
        a1= (float)rows1/sigma1/sigma1;

        temp1=gasdev(&b_time);

        mu1=b1/a1+ 1.0/sqrt(a1)*temp1;


        //   mu1=-0.75;

/// generate the random numbers for sigma1^2


        a1= rows;

        b1=0.0;

        for (j=2; j<=rows; j++)

        {


            if ( x[j-1][2]>r )

            {

                psi=0.0;

                for ( k=1;k<=cols; k++)

                {
                    psi=psi+beta1[k]*x[j][k];
                }

                b1 = b1 +  pow( y[j] - psi - mu1, 2.0) / 2.0;


            }
        }


        sigma1=pow(b1/gamdev(rows1/2,&b_time),0.5);

        //    sigma1=0.2;


///   -----------------------------------------------------------------------------------

        for ( k=1;k<=cols; k++)

        {
            c = 0;
            d = 0;


            for (j=2; j<=rows; j++)

            {
                if ( x[j-1][2] <= r )

                {

                    c = c + pow(x[j][k],2.0) / sigma2/sigma2;
                    psi=0.0;

                    for ( i=1;i<=cols; i++)

                    {
                        if ( i!=k)
                        {
                            psi=psi+beta2[i]*x[j][i];
                        }
                    }


                    d = d +  (y[j] -psi - mu2) * x[j][k] / sigma2/sigma2;
                }
            }
            temp1=gasdev(&b_time);

            beta2[k]=d/c+ 1.0/sqrt(c)*temp1;

        }
////        beta[1]=2.0;
////        beta[2]=3.0;
////        beta[3]=4.0;
//
///// generate the random number for the variable mu1


        b1=0.0;
        rows2=0;


        for (j=2; j<=rows; j++)

        {
            if ( x[j-1][2]  <= r )

            {
                rows2=rows2+1;
                psi=0.0;

                for ( k=1;k<=cols; k++)

                {
                    psi=psi+beta2[k]*x[j][k];
                }

                b1 = b1 +  (y[j] - psi) / sigma2/sigma2;


            }
        }
        a1= (float)rows2/sigma2/sigma2;

        temp1=gasdev(&b_time);

        mu2=b1/a1+ 1.0/sqrt(a1)*temp1;



//
///// generate the random numbers for sigma1^2
//

        //  a1= rows;

        b1=0.0;

        for (j=2; j<=rows; j++)

        {


            if ( x[j-1][2] <=r)

            {

                psi=0.0;

                for ( k=1;k<=cols; k++)

                {
                    psi=psi+beta2[k]*x[j][k];
                }

                b1 = b1 +  pow( y[j] - psi - mu2, 2.0) / 2.0;


            }
        }


        sigma2=pow(b1/gamdev(rows2/2,&b_time),0.5);

//        //  sigma1=0.2,
//






        // ##  simulate r

        // r=0.3;
      //  r_star=low +(high-low)*ran1(&b_time);

        r_star=r +sigma_r*gasdev(&b_time)/8.0;
        while (r_star<low || r_star>high)

        {
            r_star=r +sigma_r*gasdev(&b_time)/8.0;

        }

        likelihood_new1=0.0;
        likelihood_curre1=0.0;

        for (h=2; h<=rows; h++)

        {
            //# for current
//            if (x[h-1][2]>r)
//            {
//                temp_curr1=mu1 +beta1[1]*x[h][1] +beta1[2]*x[h][2] + beta1[3]*x[h][3];
//                likelihood_curre1 = likelihood_curre1 - temp_curr1*temp_curr1/2.0/sigma1/sigma1;
//            }
//
//            if (x[h-1][2]<=r)
//            {
//                temp_curr2=mu2 +beta2[1]*x[h][1] +beta2[2]*x[h][2] + beta2[3]*x[h][3];
//                likelihood_curre1 = likelihood_curre1 - temp_curr2*temp_curr2/2.0/sigma2/sigma2;
//
//
//            }
//            // # for r_star
//            if (x[h-1][2]>r_star)
//            {
//                temp_curr1=mu1 +beta1[1]*x[h][1] +beta1[2]*x[h][2] + beta1[3]*x[h][3];
//                likelihood_new1 = likelihood_new1 - temp_curr1*temp_curr1/2.0/sigma1/sigma1;
//
//
//            }
//
//            if (x[h-1][2]<=r_star)
//            {
//
//                temp_curr2=mu2 +beta2[1]*x[h][1] +beta2[2]*x[h][2] + beta2[3]*x[h][3];
//                likelihood_new1 = likelihood_new1 - temp_curr2*temp_curr2/2.0/sigma2/sigma2;
//
//
//            }




            if (x[h-1][2]>r)
            {
                temp_curr1=mu1 +beta1[1]*x[h][1] +beta1[2]*x[h][2] + beta1[3]*x[h][3];
                likelihood_curre1 = likelihood_curre1 +log(normpdf(y[h], temp_curr1, sigma1));
            }

            if (x[h-1][2]<=r)
            {
                temp_curr2=mu2 +beta2[1]*x[h][1] +beta2[2]*x[h][2] + beta2[3]*x[h][3];
                likelihood_curre1 = likelihood_curre1  +log(normpdf(y[h], temp_curr2, sigma2));


            }
            // # for r_star
            if (x[h-1][2]>r_star)
            {
                temp_curr1=mu1 +beta1[1]*x[h][1] +beta1[2]*x[h][2] + beta1[3]*x[h][3];
                likelihood_new1 = likelihood_new1 +log(normpdf(y[h], temp_curr1, sigma1));

            }

            if (x[h-1][2]<=r_star)
            {

                temp_curr2=mu2 +beta2[1]*x[h][1] +beta2[2]*x[h][2] + beta2[3]*x[h][3];
                likelihood_new1 = likelihood_new1 +log(normpdf(y[h], temp_curr2, sigma2));


            }
        }
        u=ran1(&b_time);

        if ( u<min(1.0, exp(likelihood_new1-likelihood_curre1)))
        {
            r=r_star;
        }




        if ( g>100)
        {
            // temp_new=exp(likelihood);
            fprintf (para_output, "%d  %f %f  %f %f  %f %f %f  %f %f %f %f \n",g,  r, mu1, beta1[1],  beta1[2], beta1[3], sigma1, mu2,beta2[1],  beta2[2], beta2[3], sigma2);
        }


    } /// the end of the loop k with the length N for gererating samples

    fclose (para_output);

//
//
//    para_output = fopen ( "test y.txt", "w");
//
//    h=1;// the time of the repeating sampling//
//
//    for ( h=1; h<=rows;  h++)
//    {
//
//        fprintf (para_output, "%f \n",y[h] );
//    }
//
//    fclose (input);
//



//    cout <<"mu="<< mu_mean/float(G)<< "    phi="<<phi_mean/float(G) <<"   rho="<<rho_mean/float(G)<<"   sigma="<<sigma_mean/float(G)<< "\n";
// cout << endl;
// system("press any key to terminate the program");
//----------------------------------------
//
//       The end of the code
//
//----------------------------------------

    return 0;


}
///-----------------------------------------------------------------------------------
///  the following are some functions that are used in the main boday

double f_star(double y, double x)
{
    double fvalue;
    double pi=3.14159265358979323846;

    fvalue=1.0/sqrt(2.0*pi*exp(x))* exp(-y*y /2.0/exp(x));

    return fvalue;
}
///-----------------------------------------------------------------------------------

double g_star1(double y, double x)
{

/// this fuction is by me
    double fvalue;
    double pi=3.14159265358979323846;

    fvalue=exp(-1.0/2.0*log(2.0*pi)-1.0/2.0*x-y*y /2.0*(1-x));

    return fvalue;
}

///----------------------------------------------------------------------------
double g_star2 (double y, double x, double z)
{
/// this fuction is from the papaer
    double fvalue;
    double pi=3.14159265358979323846;
    fvalue=exp(-1.0/2.0*log(2.0*pi)-1.0/2.0*x-y*y /2.0*exp(-z)*(1+z-x));

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




double gamma(double x)
{
    int i,k,m;
    double ga,gr,r,z;

    static double g[] =
    {
        1.0,
        0.5772156649015329,
        -0.6558780715202538,
        -0.420026350340952e-1,
        0.1665386113822915,
        -0.421977345555443e-1,
        -0.9621971527877e-2,
        0.7218943246663e-2,
        -0.11651675918591e-2,
        -0.2152416741149e-3,
        0.1280502823882e-3,
        -0.201348547807e-4,
        -0.12504934821e-5,
        0.1133027232e-5,
        -0.2056338417e-6,
        0.6116095e-8,
        0.50020075e-8,
        -0.11812746e-8,
        0.1043427e-9,
        0.77823e-11,
        -0.36968e-11,
        0.51e-12,
        -0.206e-13,
        -0.54e-14,
        0.14e-14
    };

    if (x > 171.0) return 1e308;    // This value is an overflow flag.
    if (x == (int)x)
    {
        if (x > 0.0)
        {
            ga = 1.0;               // use factorial
            for (i=2;i<x;i++)
            {
                ga *= i;
            }
        }
        else
            ga = 1e308;
    }
    else
    {
        if (fabs(x) > 1.0)
        {
            z = fabs(x);
            m = (int)z;
            r = 1.0;
            for (k=1;k<=m;k++)
            {
                r *= (z-k);
            }
            z -= m;
        }
        else
            z = x;
        gr = g[24];
        for (k=23;k>=0;k--)
        {
            gr = gr*z+g[k];
        }
        ga = 1.0/(gr*z);
        if (fabs(x) > 1.0)
        {
            ga *= r;
            if (x < 0.0)
            {
                ga = -M_PI/(x*ga*sin(M_PI*x));
            }
        }
    }
    return ga;
}
