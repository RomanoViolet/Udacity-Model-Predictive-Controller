/*
 * utils.cpp
 *
 *  Created on: Dec 5, 2017
 *      Author: dumbledore
 */
#include <iostream>
#include "utils.hpp"
#include <Eigen/QR>
#include <assert.h>
#include <limits>
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
Utils::Utils() {
  // nothing for now
}

Utils::~Utils() {

  this->previousPsi = 0.0;
  this->previousCTE = 0.0;
}

bool Utils::Compare(double a, double b)
{
	//https://stackoverflow.com/a/17341
	std::cout << __FILE__ << ": " << __LINE__ << "\t Comparing: " << a << " vs " << b << std::endl;
    return fabs(a - b) < std::numeric_limits<double>::epsilon();
}

Eigen::VectorXd Utils::polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order)
{
	  assert(xvals.size() == yvals.size());
	  assert(order >= 1 && order <= xvals.size() - 1);
	  Eigen::MatrixXd A(xvals.size(), order + 1);

	  for (int i = 0; i < xvals.size(); i++) {
	    A(i, 0) = 1.0;
	  }

	  for (int j = 0; j < xvals.size(); j++) {
	    for (int i = 0; i < order; i++) {
	      A(j, i + 1) = A(j, i) * xvals(j);
	    }
	  }

	  auto Q = A.householderQr();
	  auto result = Q.solve(yvals);
	  return result;
}


// Evaluate a polynomial.
double Utils::polyeval(Eigen::VectorXd coeffs, double x)
{
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}


double Utils::velocityInMetersPerSecondFromMilesPerHour(const double v)
{
	return(v*0.44704);
}


void Utils::coordinatesInVehicleReference(std::vector<double>& wayPoints_ptsx, std::vector<double>& wayPoints_ptsy, double& location_px, double& location_py, double& psi)
{
	/*
	 * psi is 0 degrees in the direction of the vehicle, and increases counter-clockwise
	 */

	assert (wayPoints_ptsx.size() == wayPoints_ptsy.size());

	/*double check_x;
	double check_y;*/
	double newX;
	double newY;

	for (size_t i = 0; i< wayPoints_ptsx.size(); i++)
	{
		newX = ((wayPoints_ptsx[i] - location_px) * std::cos(psi)) + ((wayPoints_ptsy[i] - location_py) * std::sin(psi));
		newY = ((location_px - wayPoints_ptsx[i]) * std::sin(psi)) - ((location_py - wayPoints_ptsy[i]) * std::cos(psi));

		/*newX = ((wayPoints_ptsx[i] - location_px) * std::cos(-psi)) - ((wayPoints_ptsy[i] - location_py) * std::sin(-psi));
		newY = ((wayPoints_ptsx[i] - location_px) * std::sin(-psi)) + ((wayPoints_ptsy[i] - location_py) * std::cos(-psi));*/

/*		double x = wayPoints_ptsx[i] - location_px;
		double y = wayPoints_ptsy[i] - location_py;*/

		wayPoints_ptsx[i] = newX;
		wayPoints_ptsy[i] = newY;

/*

		check_x = x * cos(-psi) - y * sin(-psi);
		check_y = x * sin(-psi) + y * cos(-psi);

		assert(this->Compare(check_x, wayPoints_ptsx[i]));
		assert(this->Compare(check_y, wayPoints_ptsy[i]));
*/

		/*wayPoints_ptsx[i] = (wayPoints_ptsx[i] - location_px)*cos(-psi) - (wayPoints_ptsy[i] - location_py)*sin(-psi);
		wayPoints_ptsy[i] = (wayPoints_ptsx[i] - location_px)*sin(-psi) + (wayPoints_ptsy[i] - location_py)*cos(-psi);*/

	}

	/* in vehicle's reference, the vehicle is always at (0,0), and heading 0 degrees*/
	location_px = 0.0;
	location_py = 0.0;
	psi = 0.0;
}
