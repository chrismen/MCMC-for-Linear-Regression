/// use Code::Blocks 8.02
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <time.h>
#include <math_recipes.h>

double min(double x, double y);
double max(double x, double y);

//----------------------------------------------
int main ()
{
    FILE * outputfile;
    FILE * input;

    FILE *para_output;
    FILE *para_result;

    const int G=1;
    const int N =5000;
    const int rows=10000;
    const int cols=3 ;
    const int n=N/5;
    float a ,b ,c ,d, psi,temp1,a1,b1;
    float mu1, sigma1, mu2, sigma2, rho1;
    float x[rows][cols]={};
    float y[rows]={};
    float beta[cols+1]={};
    float yy[cols+1]={};
    int j,k,h,i,g;
    long a_time, b_time;
    time_t t1;

    (void) time(&t1);
    a_time=(long) t1;
    b_time=-(a_time+1);
    sigma1=0.15; //rooted
    mu1=-0.75;
    beta[1]=2.0;
    beta[2]=3.0;
    beta[3]=4.0;
    para_output = fopen ( "MCMC-time-seris.txt", "w");
//-------------------------------------------------------------
    input = fopen ("C-simulated-logist-data-x.txt", "r");

    for ( h=0; h<rows;  h++)
    {
        for ( j=0; j<(cols-1);  j++)
        {
            fscanf (input, "%f ", &temp1);
            x[h][j]=temp1;
        }
        fscanf (input, "%f \n", &temp1);
        x[h][(cols-1)]=temp1;
    }
    fclose (input);

    input = fopen ("C-simulated-logist-data-y.txt", "r");
    for ( h=0; h<rows;  h++)
    {
        fscanf (input, "%f", y+h);
    }
    fclose (input);

    for ( g=1;g<=N; g++)
    {
        cout<< g<<"\n ";
/// sample the coefficients beta
        for ( k=0;k<cols; k++)
        {
            c = 0;
            d = 0;
            for (j=0; j<rows; j++)
            {
                c = c + pow(x[j][k],2.0) / sigma1/sigma1;
                psi=0.0;
                for ( i=0;i<cols; i++)
                {
                    if ( i!=k)
                    {
                        psi=psi+beta[i]*x[j][i];
                    }
                }
                d = d +  (y[j] -psi - mu1) * x[j][k] / sigma1/sigma1;
            }
            temp1=gasdev(&b_time);
            beta[k]=d/c+ 1.0/sqrt(c)*temp1;
        }
/// generate the random number for the variable mu1
        a1= rows/sigma1/sigma1;
        b1=0.0;
        for (j=0; j<rows; j++)
        {
            psi=0.0;
            for ( k=0;k<cols; k++)
            {
                psi=psi+beta[k]*x[j][k];
            }
            b1 = b1 +  (y[j] - psi) / sigma1/sigma1;
        }
        temp1=gasdev(&b_time);
        mu1=b1/a1+ 1.0/sqrt(a1)*temp1;
/// generate the random numbers for sigma1^2
        a1= rows;
        b1=0.0;
        for (j=0; j<rows; j++)
        {
            psi=0.0;
            for ( k=0;k<cols; k++)
            {
                psi=psi+beta[k]*x[j][k];
            }
            b1 = b1 +  pow( y[j] - psi - mu1, 2.0) / 2.0;
        }
        sigma1=pow(b1/gamdev(rows/2,&b_time),0.5);

        fprintf (para_output, "%d  %f ",g, mu1);
        for ( j=0; j<cols;  j++)
        {
            fprintf (para_output, " %f", beta[j]);
        }
        fprintf (para_output, " %f", sigma1);
        fprintf (para_output, " \n");

    } /// the end of the loop k with the length N for gererating samples
    fclose (para_output);
    cout << endl;
    system("press any key to terminate the program");
    return 0;
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
