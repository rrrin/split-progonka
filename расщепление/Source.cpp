#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include<string.h>
#include <fstream>
#include <conio.h>
#include <vector>
#include <math.h>
#include <iomanip>
#include <algorithm>

using namespace std;
typedef vector<vector<double>> mat_t;
mat_t create_mat(size_t nX, size_t nY) //создаю матрицу 
{
	ifstream file_in;
	file_in.open("text.txt");
	mat_t ans;
	ans.resize(nY);
	double x;
	for (size_t i = 0; i < nY; i++)
	{
		ans[i].resize(nX);
	}

	return ans;
}

double func_u(double x)
{
	return pow(x, 4)/8. - 7./24.*x*x;
}

vector<double> progonka(const mat_t &coeff, int n, double a0, double b0, double an, double bn, double hx)
{
	vector<double> y(n+1), p(n+1), x(n+1);

	y[1] = a0;
	p[1] = b0;
	for (int i = 1; i <=n-1; i++)
	{
		y[i+1]=-coeff[2][i] /(coeff[0][i] *y[i]+coeff[1][i]);
		p[i+1]=(coeff[3][i] - coeff[0][i] * p[i]) / (coeff[0][i] * y[i] + coeff[1][i]);
	}
	
	x[n] = (bn + p[n] * an) / (1. - y[n] * an);
	for (int i = n-1; i >=0; i--)
	{
		x[i] = y[i + 1] * x[i + 1] + p[i + 1];
	}
	return x;
}

void decision(double a, double b, int n)
{
	double hx = (b - a) / (n);
	vector<double> xi(n+1), f(n+1);
	for (int i = 0; i <= n; i++)
	{
		xi[i] = a + i * hx;
		f[i] = func_u(xi[i]);
	}
	mat_t coeff = create_mat(n+1, 4);
	for (int i = 0; i < n; i++)
	{
		coeff[0][i] = 1./(hx*hx);
		coeff[1][i] = -2./(hx*hx);
		coeff[2][i] = 1. / (hx*hx);
		coeff[3][i] = 12./8.*xi[i] * xi[i]-14./24.;
	}
	ofstream file;
	file.open("split.txt");
	file.setf(ios::fixed);
	vector<double> y= progonka(coeff, n, 0., -1./6., 0., 5./6., hx);
	for (int i = 0; i <= n; i++)
	{
		file << setprecision(15) << y[i] << " " << f[i] << " " << abs(f[i] - y[i]) << endl;
	}

}

int main()
{
	double a = 1.;
	double b = 2.;
	int n = 50;
	decision(a, b, n);
	return 0;
}