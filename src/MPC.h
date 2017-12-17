#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "utils.hpp"

using namespace std;

class MPC {
public:

	// default constructor
	MPC();

	// default destructor
	~MPC();

	// the utilities object
	Utils utils;

	// Number of steps each of duration dt
	size_t N;

	/* the duration of each step: Each tuple of (steering angle, throttle) is computed for time dt into the future.
	 * The value chosen for dt is 20% more than the MPC-to-actuator latency as this was empiritcally seen to reduce the
	 * strong corrections at the beginning of each new iteration thereby causing significant oscillations.
	 * The longer time-window allows for more gentle maneuvers as the correction is now spread over a larger time window.
	 * However, a very long time-window will result in vehicle coming off the track since rapid changes in directions are not accounted for.
	 */
	double dt;

	// Chosen to predict 100ms into the future -- the latency between the time instant a control vector is computed, and the time when the actuator executes it.
	double latency;


	// This is the length from front to CoG that has a similar radius.
	double Lf;

	// setup target velocity
	double ref_v;

	// Various helper indices
	size_t x_start;
	size_t y_start;
	size_t psi_start;
	size_t v_start;
	size_t cte_start;
	size_t epsi_start;
	size_t delta_start;
	size_t a_start;

	// private member variables
	// TODO: When stuff works, cleanup here.

	std::vector<double> setupSolver(const std::vector<double>& ptsx,
			const std::vector<double>& ptsy, const double px, const double py,
			const double psi, const double v, const double currentSteeringAngle,
			const double currentAcceleration, std::vector<double>& planned_x,
			std::vector<double>& planned_y);

	void moveStateSingleLatencyForward(const double px, const double py,
			const double psi, const double v, const double cte,
			const double epsi, const double currentSteeringAngle,
			const double currentAcceleration, const double dt,
			Eigen::VectorXd& newState);

	// Solve the model given an initial state and polynomial coefficients.
	// Return the first actuations.
	vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);

	unsigned getIterationLimit();

private:
	const unsigned iters = 10;
};

#endif /* MPC_H */
