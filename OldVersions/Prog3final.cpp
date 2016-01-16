#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <list>
#include <set>
#include <vector>

/******************************************************************************
* Class: CSC 372 - Analysis of Algorithms
* Professor: Dr. Antonette Logar
* Teacher's Assistant: Dr. Edward Corwin
* Author: Benjamin Kaiser
* Due: December 1, 2015
* Dates Modified:  
* 	11 - 18 - 15: Started/Copied Teacher's code
* 	11 - 25 - 15: Spun my wheels
* 	11 - 27 - 15: Called FFT
* 	11 - 30 - 15: Things made sense. Mostly written ;) (Thank you for
* 		the notes on it!!!!! :D
* 	12 - 1 - 15: Finished and commented
* Description: This program takes points which simulate a digital signal
* 	and then it calcluates the  
******************************************************************************/
using namespace std;

const int MAX = 15000;
const double EPSILON = 1.0e-12;

const complex<double> I (0, 1);

////////////////////////////////////////////////////////////////////
// complex_round - just set to zero values that are very small    //
//             makes output easier to read                        //
////////////////////////////////////////////////////////////////////
void complex_round(complex <double> a[], int n)
  {
  int    i;
  float  x;
  float  y;

  for (i = 0; i < n; i++)
    {
    x = a[i].real();
    y = a[i].imag();
    if (fabs(x) < EPSILON)
      x = 0;
    if (fabs(y) < EPSILON)
      y = 0;
    a[i] = complex<double>(x, y);
    }
  }

////////////////////////////////////////////////////////////////////
//                     print a complex polynomial                 //
////////////////////////////////////////////////////////////////////
void print_polynomial(complex <double> a[], int n)
  {
  int i;
    for (i = 0; i < n; i++)
    cout << a[i] << "  " << endl;
  cout << endl;
  }

////////////////////////////////////////////////////////////////////
//                        FFT based on CLRS page 911              //
////////////////////////////////////////////////////////////////////
void fft(complex <double> a[], int n, complex <double> y[],
  complex<double> omegas[])
  {
  complex  <double> even[n];
  complex  <double> even_fft[n];
  int      i;
  int      j;
  int      n2;
  complex  <double> odd[n];
  complex  <double> odd_fft[n];
  complex  <double> omega;
  complex  <double> omega_power;

  if (n == 1)
    {
    y[0] = a[0];
    return;
    }
  
  n2 = n / 2;
  j = 0;
  for (i = 0; i < n; i += 2)
    {
    even[j] = a[i];
    odd[j] = a[i + 1];
    j ++;
    }
    
  fft(even, n2, even_fft, omegas);
  fft(odd, n2, odd_fft, omegas);

    omega = cos(-2.0 * M_PI / n) + I * sin(-2.0 * M_PI / n);
 omega_power = 1;
  for (i = 0; i < n2; i ++)
    {
	omegas[i] = omega_power;
    y[i] = even_fft[i] + omega_power * odd_fft[i];
    y[i + n2] = even_fft[i] - omega_power * odd_fft[i];
    omega_power = omega * omega_power;
    }
  }

/******************************************************************************
* Author: Benjamin Kaiser
* 
******************************************************************************/
void calcQ(complex<double> y[], complex<double> omegas[], complex<double> a[], int n, int it, double ttl[])
  {
  for (int i = 0; i < n/2; i++)
    {
    y[i] = ((y[i] - a[it]) + a[it + n]) * (1.0/omegas[i]);
    ttl[i] += abs(y[i]);
    }
  }

/******************************************************************************
* Author: Benjamin Kaiser
*
******************************************************************************/
void findMaxes(double a[], int n, double maxes[], int indices[])
{
	double* copy1  = new double[n];
	double* copy2 = new double[n];
	for (int i = 0; i < n/2; i++)	
	{
		copy1[i] = a[i];
		copy2[i] = a[i];
	}
	/*for (int i = 2; i < n/2 - 3; i++)
	{
		copy1[i] = copy2[i] + copy2[i - 2] + copy2[i - 1] + copy2[i + 1] + copy2[i + 2];
	}*/
	//double* uniq = unique(copy1, copy1 + n);
	//sort(uniq, uniq + n);
	sort(copy1, copy1 + (n));
	int index = 0;
		for (int j = n  - 1; j > n - 6; j--)
		{
			maxes[index] = copy1[j];
			for (int l = 0; l < n/2; l++)
			{
				if (maxes[index] == a[l])
				{
					maxes[index] = a[l] + a[l-2] + a[l-1] + a[l+1] + a[l+2];
					indices[index] = l;
					cout << l << endl;
				}
			}
			index++;
		}
	delete [] copy1;
	delete [] copy2;
		
}

/******************************************************************************
* Author: Benjamin Kaiser
*
******************************************************************************/
int main(int argc, char** argv)
  {
  ifstream fin;
  ofstream fout;

  cout << setprecision(2) << fixed;
  fout << setprecision(2) << fixed;
	
	
  string temps = argv[1];
  int dotind = temps.rfind('.');
  temps = temps.replace(dotind, 5, ".out");


  fin.open(argv[1]);
  fout.open(temps);

  int numinfft;
  fin >> numinfft;
  int numoffft;
  fin >> numoffft;	

  complex<double>* aval = new complex<double>[numinfft+numoffft];
  complex<double>* yval = new complex<double>[numinfft+numoffft];
  complex<double>* omeg = new complex<double>[numinfft+numoffft];
  complex<double>* total = new complex<double>[numinfft+numoffft];
  double* avg = new double[numinfft+numoffft];
  double maxes[5] = {0};
  int indices[5] = {0};

  for (int i = 0; i < numinfft + numoffft; i++)
    {
    avg[i] = 0;
    }

  for (int i = 0; i < numinfft + numoffft - 1; i++)
    {
    fin >> aval[i];
    }

  //calculate the first fft
  fft(aval, numinfft, yval, omeg);
	
  for (int i = 0; i < numinfft/2; i++)
    {
    avg[i] += abs(yval[i]);
    }

  for (int j = 0; j < numoffft - 1; j++)
    {
    calcQ(yval, omeg, aval, numinfft, j, avg);
    }

  for (int i = 0; i < numinfft/2; i++)
    {
    avg[i] = avg[i]/numoffft;
    }

  findMaxes(avg, numinfft, maxes, indices); 

  for (int i = 0; i < 5; i++)
    {
    fout << indices[i] << " " << maxes[i] << endl;
    }

  delete [] aval;
  delete [] yval;
  delete [] omeg;
  delete [] total;
  delete [] avg;
  fin.close();
  fout.close();
	

  return 0;
  }
