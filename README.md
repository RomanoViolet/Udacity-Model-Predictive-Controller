# Model Predictive Controller
An implementation of a Model Predictive Controller based on lectures from the Udacity's Self Driving Care Nano Degree (SDCND) Program 

## About
This project implements a Model Predictive Controller (MPC) as required by the Udacity's [Self Driving Car Nano-Degree program.](https://www.udacity.com/course/self-driving-car-engineer-nanodegree--nd013)

The input data to the Model Predictive Controller is the cross-track error (CTE), current vehicle velocity, and current steering angle, current location of the vehicle, and a sequence of waypoints.
The MPC aims to maneuver the car so that it passes as close to the waypoints as possible, and as close to the set speed as possible (the optimization objectives), while avoiding steering or throttle jerks (cnstraints).
The waypoints are assumed to be generated by the trajectory planning algorithms implemented separately.

In the provided implementation, the speed set for the MPC to aim for is ~21 m/s (75 km/h, or 47 miles/h).
It is possible to tune the MPC controller to aim for higher vehicle speed.


### How it Looks
This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases). Select "MPC Controller" option. 

A sample frame from one of the videos as the vehicle negotiates the circuit inside the simulator is ![shown](https://github.com/RomanoViolet/Udacity-Model-Predictive-Controller/blob/master/Results/screenshot.png)
The yellow line represents the second order polynomial fitted to the provided waypoints, whereas the green line represents the path computed by the optimizer, and is to be followed by the vehicle.

### Prerequisites
_Ad-verbatim from Udacity's Instructions:_

uWebSocketIO Starter Guide

All of the projects in Term 2 and some in Term 3 involve using an open source package called [uWebSocketIO](https://github.com/uNetworking/uWebSockets). This package facilitates the same connection between the simulator and code that was used in the Term 1 Behavioral Cloning Project, but now with C++. The package does this by setting up a web socket server connection from the C++ program to the simulator, which acts as the host. In the project repository there are two scripts for installing uWebSocketIO - one for Linux and the other for macOS.

Note: Only uWebSocketIO branch e94b6e1, which the scripts reference, is compatible with the package installation.
Linux Installation:

From the project repository directory run the script: install-ubuntu.sh

### Additonal Dependencies
The MPC project requires installation of following addition dependencies:
- C++ Algorithmic Differentiation ([CPPAD]) Library and Development Headers (https://www.coin-or.org/CppAD/Doc/install.htm)
- Interior Point Optimization Library, [IPOpt](https://www.coin-or.org/download.html). Before building the library, please also build and install additional libraries contained in the folder "Third Party", especially, LAPACK, BLAS, and HSL libraries. [PDF Instructions](https://projects.coin-or.org/Ipopt/browser/stable/3.10/Ipopt/doc/documentation.pdf?format=raw)
- [MA27 Solver](http://www.hsl.rl.ac.uk/ipopt/)

Once all libraries have been built and installed, the `CMakeLists.txt` needs to be suitably adapted. A reference `CMakeLists.txt.reference` has been provided here which contains the working configuration used for developing the MPC project. Pay special attention to line `include_directories` and `link_directories` parameters. 


### Structure of the Project
The project is structured as follows:
- Folder "src": Contains the core logic required to build the MPC.
- Folder "Results": Contains a sample video and a screenshot when testing the MPC using the simulator supplied by Udacity. The yellow line represents the polynomial fitted to the provided waypoints, whereas the green line representst the trajectory computed by the optimizer, and is the path to be followed by the vehicle.
- `build.sh`: A shell script to build the MPC project.
- `install-ubuntu.sh`: A script required for installing dependencies for running the simulator. See the section above on "Prerequisites"
- `run.sh`: A shell script to run the binary built using `build.sh`, and connects to the Udacity's Term 2 simulator. See the section above "How it Looks".


### Running the Model Predictive Controller
- `build.sh`, followed by
- `run.sh`.
- Start the Udacity Term 2 Simulator, and choose the option "MPC Contoller"

### Pending Improvements
- The speed of the vehicle can be increased even further. Suggestions welcome.

### Credits
- Udacity: Lecturers, and mentors;
- Internet: for examples and samples.

### Disclaimer
Some of the ideas are borrowed and adapted from other people's work.

