#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char** argv)
{
	ifstream fin;
	ofstream fout;
	fin.open("kaiser.in");
	fout.open("kaiser.out");

	int numsinwaves;
	double amp[50];
	double freq[50];
	double shift[50];
	
	int evaltimes;
	int pointsfft;

	double evaled;

	double ang;

	//cout << "Yo, bro.  How many sin functions be ya wantin'? (We can serve ya up to 50) " << endl;
	//cin >> numsinwaves;
	fin >> numsinwaves;

	fin >> pointsfft;
	fout << pointsfft << endl;

	int numfft;
	fin >> numfft;
	fout << numfft << endl;
	

	//cout << "Dude, ya gotta tell me how many points I should crunch at! " << endl;

	//cin >> evaltimes;
	evaltimes = pointsfft + numfft - 1;

	for (int i = 0; i < numsinwaves; i++)
	{	
		//cout << "That's rad man! Whatcha wantin' that ampleetude to be? " << endl;
		//cin >> amp[i];
		fin >> amp[i];
		//cout << "Coolio! Now I kinder need ta know dat frequency " << endl;
		//cin >> freq[i];
		fin >> freq[i];
		//cout << "Brah!  Looking sharp man!  You phase shiftin'?" << endl;
		fin >> shift[i];
	}


	cout << "Right.  Let's get dis bad boy rollin'" << endl;

	//cout << "Oh yeah, man, forgot.  I also need the starting point " << endl;
	//cin >> ang;
	ang = .1;

	double stepval = .06221995;

	for (double i = 0; i < evaltimes; i++)
	{
		evaled = 0;
		for (int j = 0; j < numsinwaves; j++)
		{
			evaled +=  amp[j] * sin(freq[j] * ang + shift[j]);
		}
		ang += stepval;
		fout << evaled << endl;
	}
	fin.close();
	fout.close();
	
	return 0;
}
