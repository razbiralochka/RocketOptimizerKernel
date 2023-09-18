#include <iostream>
#include <math.h>
#include <Eigen/Dense>


using namespace std;

double const Vx = 9529;
double const w1 = 3300, w2 = 3700;
double const s1 = 11, s2 = 9;

double p1 = (s1 / (s1 - 1));
double p2 = (s2 / (s2 - 1));

double border(double x1, double x2)
{
	
	double ogr;

	ogr = Vx - w1 * log((1 + p1 * x1 + p2 * x2) / (1 + p1 * x1 + p2 * x2 - x1)) - w2 * log((1 + p2 * x2) / (1 + p2 * x2 - x2));

	return ogr;

}

double func(double x1, double x2, double lamda)
{

	double f = 1 + p1 * x1 + p2 * x2 + lamda * border(x1,x2);
	return f;
}

Eigen::Vector3d equals(Eigen::Vector3d arg)
{
	Eigen::Vector3d out;
	
	double x1 = arg(0);
	double x2 = arg(1);
	double lamda = arg(2);
	
	double h = 1e-6;
	
	out(0) = (func(x1 + h, x2, lamda) - func(x1 - h, x2, lamda)) / (2 * h);
	out(1) = (func(x1, x2 + h, lamda) - func(x1, x2 - h, lamda)) / (2 * h);
	out(2) = border(x1, x2);


	return out;
}

int main()
{
	
	Eigen::Vector3d X(20, 10, 0.1), Y, dX;
	
	Eigen::Vector3d d1(0.001, 0, 0);
	Eigen::Vector3d d2(0, 0.001, 0);
	Eigen::Vector3d d3(0, 0, 0.001);

	Eigen::Matrix3d Jac;

	for (int i = 0; i <= 10; i++) 
	{

		Jac.col(0) = (equals(X + d1) - equals(X - d1)) / 0.002;
		Jac.col(1) = (equals(X + d2) - equals(X - d2)) / 0.002;
		Jac.col(2) = (equals(X + d3) - equals(X - d3)) / 0.002;

		Y = -equals(X);

		dX = Jac.lu().solve(Y);

		X += dX;

	}
	cout << X << endl;
	return 0;
}