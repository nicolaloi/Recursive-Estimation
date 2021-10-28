# Recursive Estimation
Implementation of an Hybrid Extended Kalman Filter to estimate the states of a boat and implementation of a Particle Filter to localize a robot in a partially known environment.

##

Project rated with full marks.

## Tasks

### Hybrid Extended Kalman Filter
Design a Hybrid EKF to estimate the full state of a boat during its trip. The position, orientation and linear velocity of the boat, the wind direction, and the gyroscope sensor drift are the estimator states. 
  
Since the sensors do not provide direct measurements for all states and each sensor has a different variance, the estimates of some states may be worse than the estimates of other states.  
  
Examples of some of the estimated states during a random generated boat trajectory:

![EKF_test](https://user-images.githubusercontent.com/79461707/139303647-5f1d54b5-f994-4f92-a5cc-a3290ff10eed.png)

### Particle Filter
Design a Particle Filter that tracks the position and orientation of a mobile robot, which is moving in a closed room with a partially known contour; the *x* position of one wall is uncertain. A distance sensor is mounted on the robot, pointing in the same direction as the robot is heading. This sensor measures the distance between the robot and the first opposing wall straight in front of the robot. The robot is controlled with a known input, which prevents the robot from driving into any wall.

The objective is to design a PF that estimates the location and heading of the robot, as well as the offset of the left wall from its nominal position. Video example: 


https://user-images.githubusercontent.com/79461707/139303731-3c80d083-d6ab-4e71-9203-876dad169f7f.mp4

## Setup

Further info on the programming exercises can be found in [Assignment_Hybrid_Extended_Kalman_Filter.pdf](Assignment_Hybrid_Extended_Kalman_Filter.pdf) and [Assignment_Particle_Filter.pdf](Assignment_Particle_Filter.pdf).  
  
Run the Hybrid Extended Kalman Filter from [Hybrid_Extended_Kalman_Filter/run_HEKF.m](Hybrid_Extended_Kalman_Filter/run_HEKF.m), and the Particle Filter from [Particle_Filter/run_PF.m](Particle_Filter/run_PF.m)
