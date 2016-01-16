/***************************************************************************************************************************************
Program Name- Prog 3

Author - Akshay Singh

Class - CSC 372

Course Instructor - Dr. A. Logar

Due Date - 1st December , 2015

Written on - 1st December , 2015

Program Description and usage - The program reads in the real coefficients of a polynomial from an input file and calculates the fft.At the
                                end of the program ,the program produces an output file with 5 values.
                                
Compilation Instructions-       g++ Prog3.C -o [executable name]

Issues and Bugs:                None that I have come across yet
                                
Files - Input file - argv[1]-command line argument        
****************************************************************************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h> 
#include <vector>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <queue>
#include <iomanip>

using namespace std;

void fft(complex <double> a[], int n, complex <double> y[]);

const int MAX = 4096;
const double EPSILON = 1.0e-12;
const complex<double> I (0, 1);

/***************************************************************************************************************************************
Function Name- main()

Author - Akshay Singh

Usage - The function takes an input file and reads the values and stores it in an array.After that the function calls the file and fft functions
        and at the end prints out the top 5 values with their corresponding indexes.
         
Parameters- ifstream fin             -  file to read in values from pairs.in
            ofstream fout            -  file to write all the results to pairs.out
****************************************************************************************************************************************/

int main(int argc, char **argv)
{
	
	
	double read;
	int n, k, x;
	
	//Dynamic Array initialisation
	double *total =  NULL;  
	complex<double> *real = NULL;
	complex<double> *y = NULL;
	priority_queue < pair<double, int> > q;
	complex<double> omega, omega_power;
	

	

        ifstream fin;	
	ofstream fout;
	string fname=argv[ 1 ];
	int index=fname.find_last_of( "." );
	fname=fname.substr( 0 , index);
	string outfile=fname+".out";
	
	fin.open(argv[1]);
	fout.open(outfile);

	//read in n and k
	fin >> n >> k;

	//Allocating storage to the dynamic array
	real = new complex<double>[n];
	y = new complex<double>[n];
	total =  new double [n];

	//read in real number coefficients
	for( int i = 0; i < n+k; i++)
	{
		fin >> read;
		real[i] = read;
	}

	fft(real, n, y);

	for(int i = 0; i < n; i++)
	{
		total[i] = abs(y[i]);
	}

	for(int i = 0; i < k-1; i++)
	{
		omega = cos(-2.0 * M_PI / n) + I * sin(-2.0 * M_PI / n);
 		omega_power = 1;
		for (int j = 0; j < n/2; j++)
		{
			y[j] = (y[j] - real[i] + real[n+i])/omega_power;
			omega_power *= omega;
			total[j]+= abs(y[j]);
		}
	}

	//find top 5 val in the array
	for( int i = 0; i < n/2; i++)
	{
		q.push(std::pair<double, int>(total[i], i));
	}
	for(int i = 0; i < 5; i++)
	{
		x = q.top().second;
		read = (total[x-2] + total[x-1] + total[x] + total[x+1] + total[x+2]) / k;
		fout << x << " " << fixed << setprecision(2) << read << endl;
		q.pop();
	}

	//free the array 
	delete [] y;
	delete [] real;
	delete [] total;
	fin.close();
	fout.close();
	return 0;
}



/***************************************************************************************************************************************
Function Name- fft()

Author - Professor Logar

Usage - The function is given two complex arrays and n number of points.the function would compute fft and output it to the y array.The function also computes
        omega.
         
Parameters- a             - complex array with the real number coefficients
            n             -  number of pouints in the array a
            y             -complex array which stores all the computed y
****************************************************************************************************************************************/
void fft(complex <double> a[], int n, complex <double> y[])
{
  	complex  <double> even[MAX];
  	complex  <double> even_fft[MAX];
  	int      i;
  	int      j;
  	int      n2;
  	complex  <double> odd[MAX];
  	complex  <double> odd_fft[MAX];
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
    
  	fft(even, n2, even_fft);
  	fft(odd, n2, odd_fft);

  	
    omega = cos(-2.0 * M_PI / n) + I * sin(-2.0 * M_PI / n);
 	omega_power = 1;
  	for (i = 0; i < n2; i ++)
    {
    	y[i] = even_fft[i] + omega_power * odd_fft[i];
    	y[i + n2] = even_fft[i] - omega_power * odd_fft[i];
    	omega_power = omega * omega_power;
    }
    return;
}
