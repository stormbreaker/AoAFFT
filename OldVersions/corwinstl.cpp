////////////////////////////////////////////////////////////////////
//           FFT based on code from CLRS Algorithms Text          //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <complex>
using namespace std;

const int MAX = 8096;
const double EPSILON = 1.0e-12;
const complex<double> I (0, 1);
enum  direction {FORWARD, INVERSE};

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
  direction dir, complex<double> omegas[])
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
    
  fft(even, n2, even_fft, dir, omegas);
  fft(odd, n2, odd_fft, dir, omegas);

  if (dir == FORWARD)
    omega = cos(-2.0 * M_PI / n) + I * sin(-2.0 * M_PI / n);
  else
    omega = cos(2.0 * M_PI / n) + I * sin(2.0 * M_PI / n);
  omega_power = 1;
  for (i = 0; i < n2; i ++)
    {
	  omegas[i] = omega_power;
    y[i] = even_fft[i] + omega_power * odd_fft[i];
    y[i + n2] = even_fft[i] - omega_power * odd_fft[i];
	omega_power = omega * omega_power;
    }
  }

////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
  {
  double   a_real[MAX];
  complex  <double> a[MAX];
  complex  <double> b[MAX] = {0};
  int      i;
  ifstream inf;
  int      n;
  complex  <double> y[MAX] = {0};
	complex<double> o[MAX] = {0};

  inf.open(argv[1]);
  if (inf.fail())
    {
    cout << "Unable to open fft_2011.in" << endl;
    exit(1);
    }

  inf >> n;
  cout << "n = " << n << endl;
  for (i = 0; i < n; i++)
    {
    inf >> a_real[i];
    a[i] = a_real[i];
    }

  cout << "Original:       ";
  print_polynomial(a, n);

  // Do forward FFT

  fft(a, n, y, FORWARD, o);
	/*for (int i = n/2; i < n; i++)
	{
		o[i] = conj(o[i - n/2]);
	}*/

  // Clean up result by setting small values to zero

  complex_round(y, n);

  cout << "Forward FFT:    ";
  print_polynomial(y, n);
	complex_round(o,n);
	print_polynomial(o, n);

  return 0;   // don't do inverse for comparisons
  // Do inverse FFT

  fft(y, n, b, INVERSE, o);

  // Clean up result by setting small values to zero

  complex_round(b, n);
  cout << "Inverse FFT:    ";
  print_polynomial(b, n);

  // Reconstruct original polynomial by dividing by n

  for (i = 0; i < n; i++)
    {
    b[i] /= n;
    }
  cout << "Reconstructed:  ";
  print_polynomial(b, n);

  }

////////////////////////////////////////////////////////////////////

