#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include <assert.h>



using CppAD::AD;


class FG_eval {
 public:

  /*Fitted polynomial coefficients*/
  Eigen::VectorXd coeffs;
  MPC* mpc;
	FG_eval(Eigen::VectorXd coeffs, MPC* mpc) {
		this->coeffs = coeffs;
		this->mpc = mpc;

	}

	typedef CPPAD_TESTVECTOR(AD<double>)ADvector;
	void operator()(ADvector& fg, const ADvector& vars) {

		/* `fg` a vector of the cost constraints,
		 * `vars` is a vector of variable values (state & actuators)
		 */


		/*https://www.coin-or.org/CppAD/Doc/ipopt_solve_get_started.cpp.htm*/
		/*
		 * cost function always goes into fg[0], and is expressed in terms of tunables
		 * We will try to minimize:
		 * 1. Error between predicted and actual trajectory (i.e., the cross-track error), yaw-error, and velocity error
		 * 2. The jerk (i.e., rapid change of steering angle, and also rapid changes of acceleration) [Not Done Yet]
		 */
		fg[0] = 0; /* initialize*/
		for (size_t t = 0; t < this->mpc->N; t++)
		{
			// high multiplicative factor for CTE, Yaw and velocity errors forces the solver to reduce errors in these components as a first class optimization objective.
			fg[0] += 2000.0 * CppAD::pow(vars[this->mpc->cte_start + t], 2) /*Minimize CTE*/;
			fg[0] += 10000.0 * CppAD::pow(vars[this->mpc->epsi_start + t], 2) /*Minimize Yaw Error*/;
			fg[0] += 10000.0 * CppAD::pow(vars[this->mpc->v_start + t] - this->mpc->ref_v, 2) /*Minimize Velocity Error*/;
		}

		/* 2. The jerk (i.e., rapid change of steering angle, and also rapid changes of acceleration)*/
		for (size_t t = 0; t < this->mpc->N - 2/* changes to actuations are 1 less than number of actuations*/; t++)
		{
			// A high multiplicative factor for steering jerk prevents rapid changes in the steering angle.
			fg[0] += 10000.0 * CppAD::pow(vars[this->mpc->delta_start + t + 1] - vars[this->mpc->delta_start + t], 2) /*Minimize steering jerk*/;

			// The optimizer can work with throttle as the allowed degree of freedom to acheive the objectives related to CTE, Yaw Angle, and Velocity.
			fg[0] += 100.0 * CppAD::pow(vars[this->mpc->a_start + t + 1] - vars[this->mpc->a_start + t], 2) /*Minimize acceleration jerk*/;
		}

		/*Minimize the use of actuators.*/
		for (size_t t = 0; t < this->mpc->N - 1; t++) {
			/* Reduce the number of times a steering angle command is changed.
			 * A very high number may lead to an unresponsive steering.
			*/
			fg[0] += 10000.0 * CppAD::pow(vars[this->mpc->delta_start + t], 2);

			/*
			 * The optimizer can use the throttle/brake much more freely.
			 */
			fg[0] += 10.0 * CppAD::pow(vars[this->mpc->a_start + t], 2);
		}

		// adjust the indices since the solver will always put the overall cost function at index 0.
		fg[1 + this->mpc->x_start] = vars[this->mpc->x_start];
		fg[1 + this->mpc->y_start] = vars[this->mpc->y_start];
		fg[1 + this->mpc->psi_start] = vars[this->mpc->psi_start];
		fg[1 + this->mpc->v_start] = vars[this->mpc->v_start];
		fg[1 + this->mpc->cte_start] = vars[this->mpc->cte_start];
		fg[1 + this->mpc->epsi_start] = vars[this->mpc->epsi_start];


		// The rest of the constraints
		for (size_t t = 1; t < this->mpc->N; t++) {
			AD<double> x1 = vars[this->mpc->x_start + t];
			AD<double> y1 = vars[this->mpc->y_start + t];
			AD<double> psi1 = vars[this->mpc->psi_start + t];
			AD<double> v1 = vars[this->mpc->v_start + t];
			AD<double> cte1 = vars[this->mpc->cte_start + t]/* cross track error */;
			AD<double> epsi1 = vars[this->mpc->epsi_start + t ]/* acceleration */;

			AD<double> x0 = vars[this->mpc->x_start + t - 1];
			AD<double> y0 = vars[this->mpc->y_start + t - 1];
			AD<double> psi0 = vars[this->mpc->psi_start + t - 1];
			AD<double> v0 = vars[this->mpc->v_start + t - 1];
			AD<double> cte0 = vars[this->mpc->cte_start + t - 1];
			AD<double> epsi0 = vars[this->mpc->epsi_start + t - 1]/* acceleration */;

			AD<double> delta0 = vars[this->mpc->delta_start + t - 1];
			AD<double> a0 = vars[this->mpc->a_start + t - 1]/* acceleration */;


			AD<double> refY0 =  coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2)/* reference y coordinate, assuming that x coordinate is always correct*/;
			AD<double> psiDesired0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 )/* d(refY0)/dx, evaluated at x=x0, and then converted to radians*/;
			// Here's `x` to get you started.
			// The idea here is to constraint this value to be 0.
			//
			// NOTE: The use of `AD<double>` and use of `CppAD`!
			// This is also CppAD can compute derivatives and pass
			// these to the solver.
			fg[1 + this->mpc->x_start + t] 	= x1 - (x0 + v0 * CppAD::cos(psi0) * this->mpc->dt)/* how to go from x[t] to x[t+1]*/;
			fg[1 + this->mpc->y_start + t] 	= y1 - (y0 + v0 * CppAD::sin(psi0) * this->mpc->dt)/* how to go from y[t] to y[t+1]*/;
			fg[1 + this->mpc->psi_start + t] 	= psi1 - (psi0 + ((v0/this->mpc->Lf) * -1.0/* delta>0 turns the steering right instead of assumed left*/ * delta0 * this->mpc->dt))/* how to go from psi[t] to psi[t+1]*/;
			fg[1 + this->mpc->v_start + t] 	= v1 - (v0 + a0 * this->mpc->dt)/* how to go from v[t] to v[t+1]*/;
			fg[1 + this->mpc->cte_start + t] 	= cte1 - ((refY0 - y0)/* CTE at time t*/ + /*Change in CTE*/(v0 * CppAD::sin(epsi0) * this->mpc->dt)) /* how to go from cte[t] to cte[t+1]*/;
			fg[1 + this->mpc->epsi_start + t] 	= epsi1 - ((psi0 - psiDesired0)/* epsi at time t*/ + /*Change in epsi*/(((v0/this->mpc->Lf) * -1.0/* delta>0 turns the steering right instead of assumed left*/ * delta0 * this->mpc->dt))) /* how to go from epsi[t] to epsi[t+1]*/;
		}

	}
};

//
// MPC class definition implementation.
//
MPC::MPC() {

	// Number of steps each of duration dt
	N = 5;

	/* the duration of each step: Each tuple of (steering angle, throttle) is computed for time dt into the future.
	 * The value chosen for dt is 20% more than the MPC-to-actuator latency as this was empirically seen to reduce the
	 * strong corrections at the beginning of each new iteration thereby causing significant oscillations.
	 * The longer time-window allows for more gentle maneuvers as the correction is now spread over a larger time window.
	 * However, a very long time-window will result in vehicle coming off the track since rapid changes in directions are not accounted for.
	 *
	 * 0.15 also works, and 0.1 (equal to latency) also works.
	 */
	dt = 0.1;

	// Chosen to predict 100ms into the future -- the latency between the time instant a control vector is computed, and the time when the actuator executes it.
	latency = 0.1;

	// This value assumes the model presented in the classroom is used.
	//
	// It was obtained by measuring the radius formed by running the vehicle in the
	// simulator around in a circle with a constant steering angle and velocity on a
	// flat terrain.
	//
	// Lf was tuned until the the radius formed by the simulating the model
	// presented in the classroom matched the previous radius.
	//
	// This is the length from front to CoG that has a similar radius.
	Lf = 2.67;

	// setup target velocity
	ref_v = 1.2*17.8816 /*in m/s, equal to 48 miles per hour*/;

	// The solver takes all the state variables and actuator
	// variables in a singular vector. Thus, we should to establish
	// when one variable starts and another ends to make our lifes easier.
	x_start = 0;
	y_start = x_start + N;
	psi_start = y_start + N;
	v_start = psi_start + N;
	cte_start = v_start + N;
	epsi_start = cte_start + N;
	delta_start = epsi_start + N;
	a_start = delta_start + N - 1;
}
MPC::~MPC() {}

unsigned MPC::getIterationLimit()
{
	return (this->iters);
}

void MPC::moveStateSingleLatencyForward(const double px, const double py,
		const double psi, const double v, const double cte, const double epsi,
		const double currentSteeringAngle, const double currentAcceleration,
		const double dt, Eigen::VectorXd& newState) {
	// new x coordinate
	newState[0] = (px + v * std::cos(psi) * dt);

	// new y coordinate
	newState[1] = (py + v * std::sin(psi) * dt);

	// new yaw angle
	newState[2] = (psi
			+ ((v / Lf)
					* (-1.0/* since delta>0 turns the car to the right instead of normally assumed left*/
							* currentSteeringAngle) * dt));

	// new velocity
	newState[3] = (v + currentAcceleration * dt);

	// new cte
	newState[4] = (cte + (v * std::sin(epsi) * dt));

	// new yaw error
	newState[5] = (epsi
			+ ((v / Lf)
					* (-1.0/* since delta>0 turns the car to the right instead of normally assumed left*/
							* currentSteeringAngle) * dt));

}

std::vector<double> MPC::setupSolver(const std::vector<double>& ptsx, const std::vector<double>& ptsy, const double px, const double py, const double psi, const double v, const double currentSteeringAngle, const double currentAcceleration, std::vector<double>& planned_x, std::vector<double>& planned_y)
{
	Eigen::VectorXd state = Eigen::VectorXd::Zero(6);
	Eigen::VectorXd ptsx_;
	Eigen::VectorXd ptsy_;

	/* fit a polynomial to the above x and y coordinates
	 * We throw away the last three waypoints since it fitting the polynomial to all give waypoints
	 * resulted in a large CTE and yaw-angle errors (epsi)
	 */
	ptsx_ = Eigen::VectorXd::Map(ptsx.data(), ptsx.size());
	ptsy_ = Eigen::VectorXd::Map(ptsy.data(), ptsy.size());
	auto coeffs = utils.polyfit(ptsx_.head(ptsx_.size() - 3), ptsy_.head(ptsy_.size() - 3), 2)/* Let's start with a quadratic */;

	// calculate the cross track error
	double cte = utils.polyeval(coeffs, px)/* actual x coordinate of the vehicle*/ - /*actual y coordinate of the vehicle*/py;

	// calculate the orientation error
	double epsi = /* required yaw: tangent to the trajectory*/-1* std::atan(coeffs[1] + (2 * coeffs[2] * px))  - /* actual yaw of the vehicle */ psi /* since all waypoints have been rearranged to have 0 degrees pointing in the direction of the car*/;


	/*
	 * The sensor inputs (i.e., current steering angle, throttle) are already committed to for
	 * actuation at the next time-step (i.e., current time + 100ms (latency).
	 * Therefore, the optimization step starts with predicting the state of the vehicle at the
	 * next time step when currently committed actuation commands have been executed, and from that point
	 * (i.e., current time + 100ms (latency) as a reference, start the optimization that will compute the
	 * actuation commands to be executed at (i.e., current time + 100ms (latency) + 100ms (horizon)).
	 *
	 * Therefore, we first move the state of the vehicle forward by 1 latency.
	 * Use state transformation equations
	 *
	 * https://discussions.udacity.com/t/how-to-incorporate-latency-into-the-model/257391/3
	 */
	moveStateSingleLatencyForward(px, py, psi, v, cte, epsi, currentSteeringAngle, currentAcceleration, latency, state);

	//state << px, py, psi, v, cte, epsi;
	std::vector<double> x_vals = { state[0] };
	std::vector<double> y_vals = { state[1] };
	std::vector<double> psi_vals = { state[2] };
	std::vector<double> v_vals = { state[3] };
	std::vector<double> cte_vals = { state[4] };
	std::vector<double> epsi_vals = { state[5] };
	std::vector<double> delta_vals = { };
	std::vector<double> a_vals = { };
	std::vector<double> vars;
	vars.reserve(8);

	/*
	 * Compute the actuation commands, and well as the vehicle state
	 * at current time + 100ms (latency) + 100ms (horizon).
	 */
	vars = Solve(state, coeffs);

	/*
	 * Clean the vectors
	 */
	planned_x.clear();
	planned_y.clear();
	for (size_t j = 0; j < iters; ++j)
	{
		planned_x.push_back(vars[2 + 2*j]);
		planned_y.push_back(vars[2 + 2*j + 1]);
	}

	return(vars);

}


vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
	typedef CPPAD_TESTVECTOR(double)Dvector;

	double x = state[0];
	double y = state[1];
	double psi = state[2];
	double v = state[3];
	double cte = state[4];
	double epsi = state[5];

	// number of independent variables
	// N timesteps == N - 1 actuations
	size_t n_vars = (N * 6)/* each of 6 state variables [x,y,ψ,v,cte,eψ]for each time step*/ + ((N - 1) * 2) /* actuator controls [δ,a] between successive time steps*/;

	// Number of constraints
	size_t n_constraints = (N * 6) /* each of [x,y,ψ,v,cte,eψ] has a constraint associated to it, indicating how to go from step [t] to [t+1]*/;

	// Initial value of the independent variables.
	// Should be 0 except for the initial values.
	Dvector vars(n_vars);
	for (size_t i = 0; i < n_vars; i++) {
		vars[i] = 0.0;
	}

	// Set the initial variable values
	vars[x_start] = x;
	vars[y_start] = y;
	vars[psi_start] = psi;
	vars[v_start] = v;
	vars[cte_start] = cte;
	vars[epsi_start] = epsi;

	// Lower and upper limits for x
	/*
	 * An upper and a lower bound exists for each state variable and actuator variables for each of the predicted timesteps.
	 * Notice that non-actuator variables exist for all N timesteps, whereas actuator control variables exist for N-1 timesteps.
	 */
	Dvector vars_lowerbound(n_vars);
	Dvector vars_upperbound(n_vars);

	// Set all non-actuators upper and lowerlimits
	// to the max negative and positive values.
	/*
	 * All non-actuator variables are unconstrained
	 */
	for (size_t i = 0; i < delta_start; i++) {
		vars_lowerbound[i] = -1.0e19;
		vars_upperbound[i] = 1.0e19;
	}

	// The upper and lower limits of delta are set to -25 and 25
	// degrees (values in radians).
	// NOTE: Feel free to change this to something else.
	/*
	 * Constrain the value of the steering angle for all time steps to be between -25 and +25 degrees (converted to radians)
	 */
	for (size_t i = delta_start; i < a_start; i++) {
		vars_lowerbound[i] = -0.436332;
		vars_upperbound[i] = 0.436332;
	}

	// Acceleration/deceleration upper and lower limits.
	// NOTE: Feel free to change this to something else.
	/*
	 * Constrain the value of the braking/acceleration for all time steps to be between -1 and +1.
	 */
	for (size_t i = a_start; i < n_vars; i++) {
		vars_lowerbound[i] = -1.0;
		vars_upperbound[i] = 1.0;
	}



	/*
	 * Each of the six constraints for all N steps will have upper and lower bounds associated.
	 * Since we model the constraints as errors (e.g., v[t+1]-(v[t]+ increment) = 0), all constraints
	 * are set to zero, except for initial states.
	 */

	// Lower and upper limits for constraints
	// All of these should be 0 except the initial
	// state indices.
	Dvector constraints_lowerbound(n_constraints);
	Dvector constraints_upperbound(n_constraints);
	for (size_t i = 0; i < n_constraints; i++) {
		constraints_lowerbound[i] = 0;
		constraints_upperbound[i] = 0;
	}
	constraints_lowerbound[x_start] = x;
	constraints_lowerbound[y_start] = y;
	constraints_lowerbound[psi_start] = psi;
	constraints_lowerbound[v_start] = v;
	constraints_lowerbound[cte_start] = cte;
	constraints_lowerbound[epsi_start] = epsi;

	constraints_upperbound[x_start] = x;
	constraints_upperbound[y_start] = y;
	constraints_upperbound[psi_start] = psi;
	constraints_upperbound[v_start] = v;
	constraints_upperbound[cte_start] = cte;
	constraints_upperbound[epsi_start] = epsi;


  // object that computes objective and constraints
  FG_eval fg_eval(coeffs, this);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;

  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";

  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;
  std::vector<double> result;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
  assert(ok==true);

  // Cost: commented out to reduce unnecessary screen prints
  //auto cost = solution.obj_value;


  result.push_back(solution.x[delta_start]);
  result.push_back(solution.x[a_start]);

  /*
   * Also return the array of predicted x-values and predicted y-values
   * for all time steps until the chosen horizon (t, t+N*dt] so that the
   * path can be drawn in the simulator
   */

  for (size_t i = 0; i < this->getIterationLimit(); ++i)
  {
	  result.push_back(solution.x[x_start + i]);
	  result.push_back(solution.x[y_start + i]);
  }

  return(result);

}
