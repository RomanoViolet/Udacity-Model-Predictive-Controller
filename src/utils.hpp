/*
 * utils.hpp
 *
 *  Created on: December 5, 2017
 *      Author: dumbledore
 */

#ifndef SRC_UTILS_HPP_
#define SRC_UTILS_HPP_

#include <array>
#include <Eigen/Core>

class Utils {
 public:
  // default constructor
  Utils();

  // default destructor
  ~Utils();

  double previousPsi;
  double previousCTE;

  Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals, int order);
  double polyeval(Eigen::VectorXd coeffs, double x);
  void coordinatesInVehicleReference(std::vector<double>& wayPoints_ptsx, std::vector<double>& wayPoints_ptsy, double& location_px, double& location_py, double& psi);
  double velocityInMetersPerSecondFromMilesPerHour(const double v);
  bool Compare(double a, double b);

 private:

};


/*#include "templatedUtils.cpp"*/
#endif /* SRC_UTILS_HPP_ */
