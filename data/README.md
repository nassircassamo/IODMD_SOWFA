## Data repository

Two different folders exist in the link, corresponding to two different datasets where two distinct wind farm control strategies are tested:

### Pitch Control 

1.  ***pitch_control***: axial induction control is used as a wind farm control strategy. Wind turbine blades are collectively pitched. In this directory, two additional folders and two additional files are present:
* 1.1. ***steps_theta_col_new***: folder containing all turbine relevant information saved during simulation to generate idenfitication data set.
* 1.2. ***steps_theta_col_new_val***: folder containing all turbine relevant information saved during simulation to generate validation data set.
* 1.3. ***U_data_complete_vec_pulse.mat***: MATLAB structure containing flow field information of identification data set.
* 1.4. ***U_data_complete_vec_pulse_val.mat***: MATLAB structure containing flow field information of validation data set.

### Yaw Control 

2. ***yaw_control***: a yaw control strategy is used: wake redirection control is used as a wind farm control strategy. The wind turbine tower nacelle is yawed. Similarly to pitch control, in this directory two additional folders and two additional files are present:
* 1.1. ***steps_yaw_20deg_10offset***: folder containing all turbine relevant information saved during simulation to generate identification data set.
* 1.2. ***steps_yaw_20deg_10offset_val***: folder containing all turbine relevant information saved during simulation to generate validation data set.
* 1.3. ***U_data_complete_vec_yaw_off.mat***: MATLAB structure containing flow field information of identification data set.
* 1.4. ***U_data_complete_vec_yaw_off_val.mat***: MATLAB structure containing flow field information of validation data set.

The MATLAB structures containing flow field information have the following variables:
* ***Decimate***: Number of spacing used to resample original grid information
* ***Nbegin***: Time sample from which to start saving simulation information
* ***Nend***: Time sample to end simulation information
* ***QQ_u***: Reasampled (in space) and resized streamwise velocity, for all time instants
* ***QQ_v***: Reasampled (in space) and resized spanwise velocity, for all time instants
* ***QQ_w***: Reasampled (in space) and resized vertical velocity, for all time instants
* ***dN***: flow field sampling period 
* ***x***: streamwise grid points
* ***y***: spanwise grid points
* ***z***: vertical grid points

The data were collected in a high fidelity using SOWFA (Simulator for Offshore Wind Farm Applications), a high fidelity simulator developed by the National Renewable Energy Laboratory (NREL).