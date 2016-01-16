#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <list> //not sure if necessary
#include <set> //not sure if necessary
#include <vector> //not sure if necessary

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
* Name: CalcQ
* Author: Benjamin Kaiser
* Description: This function calculates the fft's in linear time.  This is the
* 	"sliding" window function.  This was computed using Dr. Logar's notes.  
* Parameters:
*	y: the y values from the previous function
*	omegas: the omega values calculated from the fft
*	a: the a values from the file
*	n: the number of values in the fft
*	it:	this is k.  Stands for iteration
*	ttl: the total array that needs to be added to.
******************************************************************************/
void calcQ(complex<double> y[], complex<double> omegas[], complex<double> a[], int n, int it, double ttl[])
  {
  for (int i = 0; i < n/2; i++)
    {
    y[i] = ((y[i] - a[it]) + a[it + n]) * (1.0/omegas[i]); //computer new y's
    ttl[i] += abs(y[i]); //add the running total
    }
  }

/******************************************************************************
* Name: findMaxes
* Author: Benjamin Kaiser
* Description: This function finds the top 5 peaks.  This was modified from
*	Dr. Logar's notes because the algorithm in her notes did not quite get
*	the values from the output file because it added before finding peaks.
*	I add after finding max values.  
* Parameters:
*	a: is an array in which we are trying to find the max values
*	n: this is the number of points in the fft
*	maxes: these are the max values which were found
*	indices: this is the array containing the indices at which the maxes
*		were found
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
	// this commented code is leftovers from using Dr. Logar's notes
	/*for (int i = 2; i < n/2 - 3; i++)
	{
		copy1[i] = copy2[i] + copy2[i - 2] + copy2[i - 1] + copy2[i + 1] + copy2[i + 2];
	}*/
	//double* uniq = unique(copy1, copy1 + n);
	//sort(uniq, uniq + n);
	//sort these so we can find the top values
	sort(copy1, copy1 + (n));
	int index = 0;
		for (int j = n  - 1; j > n - 6; j--)
		{
			maxes[index] = copy1[j];
			//check for this in the list to find index and do the add stuff
			for (int l = 0; l < n/2; l++)
			{
				if (maxes[index] == a[l])
				{
					maxes[index] = a[l] + a[l-2] + a[l-1] + a[l+1] + a[l+2];
					indices[index] = l;
					//cout << l << endl;
				}
			}
			index++;
		}
	delete [] copy1;
	delete [] copy2;
		
}

/******************************************************************************
* Name: Main
* Author: Benjamin Kaiser
* Description: This function does the main processing of the function.  It
*	calls other functions.  It also computes the averages of the totals.  
*	it also outputs to the file.  
* Parameters:
*	argc: this is the number of command line arguments
*	argv: this is the command line arguments
******************************************************************************/
int main(int argc, char** argv)
  {
	//output files and set precision
  ifstream fin;
  ofstream fout;
  cout << setprecision(2) << fixed;
  fout << setprecision(2) << fixed;
	
	//declaring variables
  string temps = argv[1];
  int dotind = temps.rfind('.');
  temps = temps.replace(dotind, 5, ".out");


  fin.open(argv[1]);
  fout.open(temps);

  int numinfft;
  fin >> numinfft;
  int numoffft;
  fin >> numoffft;	


	//declare all the arrays that we need
  complex<double>* aval = new complex<double>[numinfft+numoffft];
  complex<double>* yval = new complex<double>[numinfft+numoffft];
  complex<double>* omeg = new complex<double>[numinfft+numoffft];
  complex<double>* total = new complex<double>[numinfft+numoffft];
  double* avg = new double[numinfft+numoffft];
  double maxes[5] = {0};
  int indices[5] = {0};
//initialize
  for (int i = 0; i < numinfft + numoffft; i++)
    {
    avg[i] = 0;
    }
//read in values
  for (int i = 0; i < numinfft + numoffft - 1; i++)
    {
    fin >> aval[i];
    }

  //calculate the first fft
  fft(aval, numinfft, yval, omeg);

//running total from fft function call	
  for (int i = 0; i < numinfft/2; i++)
    {
    avg[i] += abs(yval[i]);
    }

//calculate the next k - 1 ffts
  for (int j = 0; j < numoffft - 1; j++)
    {
    calcQ(yval, omeg, aval, numinfft, j, avg);
    }

//get the average
  for (int i = 0; i < numinfft/2; i++)
    {
    avg[i] = avg[i]/numoffft;
    }

//get the five max points
  findMaxes(avg, numinfft, maxes, indices); 

//print to file
  for (int i = 0; i < 5; i++)
    {
    fout << indices[i] << " " << maxes[i] << endl;
    }

//clean up our mess
  delete [] aval;
  delete [] yval;
  delete [] omeg;
  delete [] total;
  delete [] avg;
  fin.close();
  fout.close();
	

  return 0;
  }
