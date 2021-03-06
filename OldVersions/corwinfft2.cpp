#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <complex>


using namespace std;

const int MAX = 1024;
const double EPSILON = 1.0e-12;

void round(complex<double> a[], int n)
{
	int i;
	float x;
	float y;

	for (i = 0; i < n; i++)
	{
		x = a[i].real();
		y = a[i]. imag();
		if (fabs(x) < EPSILON)
			x = 0;
		if (fabs(y) < EPSILON)
			y = 0;
		a[i] = complex<double>(x,y);
	}
}

void swap(complex<double> &x, complex<double> &y)
{
	complex<double> old_x;
	complex<double> old_y;

	old_x = x;
	old_y = y;
	x = old_y;
	y = old_x;
}

void reverse(complex<double> a[], int n)
{
	int i;
	for (i = 0; i < n/2; i++)
		swap(a[i], a[n - i - 1]);
}

void print_polynomial(complex<double> a[], int n)
{
	int i;
	for (i = 0; i < n; i++)
		cout << a[i] << " ";
	cout << endl;
}

void forward_fft(complex<double> a[], int n, complex<double> y[])
{
	complex<double> even[MAX];
	complex<double> even_fft[MAX];
	int i;
	int j;
	int n2;
	complex<double> odd[MAX];
	complex<double> odd_fft[MAX];
	complex<double> omega;
	complex<double> omega_power;

	if (n == 1)
	{
		y[0] = a[0];
		return;
	}

	n2 = n/2;
	j = 0;
	for (i = 0; i < n; i +=2)
	{
		even[j] = a[i];
		odd[j] = a[i+1];
		j++;
	}
	forward_fft(even, n2, even_fft);
	forward_fft(odd, n2, odd_fft);

	omega = complex<double>(cos(2.0 * M_PI /n), sin(2.0 * M_PI /n));
	omega_power = 1;
	for (i = 0; i < n2; i++)
	{
		omega_power = omega * omega_power;
		y[i] = even_fft[i] + omega_power * odd_fft[i];
		y[i + n2] = even_fft[i] - omega_power * odd_fft[i];
	}
}

void inverse_fft(complex<double> y[], int n, complex<double> a[])
{
	complex<double> even[MAX];
	complex<double> even_fft[MAX];
	int i;
	int j;
	int n2;
	complex<double> odd[MAX];
	complex<double> odd_fft[MAX];
	complex<double> omega;
	complex<double> omega_power;

	if(n == 1)
	{
		a[0] = y[0];
		return;
	}

	n2 = n/2;
	j = 0;
	for (i = 0; i < n; i+=2)
	{
		even[j] = y[i];
		odd[j] = y[i+1];
		j++;
	}

	inverse_fft(even, n2, even_fft);
	inverse_fft(odd, n2, odd_fft);

	omega = complex<double>(cos(-2.0 * M_PI /n), sin(-2.0 * M_PI /n));
	omega_power = 1;
	for (i = 0; i < n2; i++)
	{
		omega_power = omega * omega_power;
		a[i] = even_fft[i] + omega_power * odd_fft[i];
		a[i+n2] = even_fft[i] - omega_power * odd_fft[i];
	}
}

int main(int argc, char** argv)
{
	complex<double> a[MAX];
	complex<double> b[MAX] = {0};
	ifstream fin;
	int i;
	int n;
	complex<double> y[MAX] = {0};
	
	fin.open(argv[1]);
	if (fin.fail())
	{
		cout << "Unable to open " << argv[1] << endl;
	}

	fin >> n;
	cout << "n = " << n << endl;
	for (i = 0; i < n; i++)
	{
		fin >> a[i];
	}

	cout << "Original: ";
	print_polynomial(a, n);

	//forward fft

	forward_fft(a, n, y);

	//clean up result by setting small values to zero
	
	round(y, n);

	cout << "Forward FFT:  ";
	print_polynomial(y, n);

	//do inverse fft, first reverse the y vector
	
	reverse(y, n);

	inverse_fft(y, n, b);

	//clean up
	
	round(b, n);
	cout << "Inverse FFT: ";
	print_polynomial(b, n);

	// reconstruct original polynomial by reversing and dividing by n
	

	reverse(b, n);
	for (i = 0; i < n; i++)
	{
		b[i] = b[i]/n;
	}
	cout << "Reconstructed: ";
	print_polynomial(b, n);
}
