# Koopman Dynamic Mode Decomposition for Wind Farm Control

## Thesis to obtain the Master of Science Degree in Mechanical Engineering

## Repository
This repository contains all the work developed in the context of the above mentioned Master Thesis dissertation. The main outcomes of this work have been segmented as follows:
* **Thesis.pdf**: corresponds to the final dissertation, with an introduction to wind turbine control, wind farm control, data driven modelling within fluid dynamics, dynamic mode decomposition and variant algorithms suited for control. All results are also contained in this document.
* **Thesis_presentation.pdf**: a presentation of the *Thesis.pdf* contents.
* **ExtendedAbstract.pdf**: a 10 page summary, in two column format, of the *Thesis.pdf*.
* **poster_thesis.pdf**: a poster format summary of the objectives and results in the *Thesis.pdf*.
* **KOPMAN_IODMD_1.0**: the source code, developed in Matlab and making use of existing functions developed by the National Renewable Energy laboratory (NREL), used to obtain the results in the *Thesis.pdf*.
* **animations**: animations where the datasets can be visualised.
* **articles**: articles published based on the work developed within the Thesis.
* **data**: contains the datasets used to test the proposed algorithms in *Thesis.pdf*.

## Thesis Abstract
Sitting wind turbines together in wind farms is economically advantageous. However, as the first turbine extracts energy from the wind, less power is available for downstream turbines. Current industry practices neglect the aerodynamic interaction, optimizing only at the individual turbine level, which leads to suboptimal behaviour of the total wind farm. Controlling wind farms as a whole is becoming increasingly important. Nevertheless, due to the fact that wind farms are high order systems whose dynamics are governed by nonlinear partial differential equations with no known analytical solution, the design and implementation of numerical optimal controllers in high fidelity simulators becomes computationally expensive and unsuitable for real time usage. Reduce order state models provide a possible route to the design and implementation of practical cooperative wind farm controllers. This thesis makes use of an innovative algorithm in the context of wind farm modelling - Input Output Dynamic Mode Decomposition - and the ideas of non linear dynamical system theory - the Koopman Operator - to find suitable reduced order models to be used for model predictive control. The wind farm control strategy of wake redirection control is studied. It is shown that a reduced state space model with 37 states can accurately reproduce the downstream turbine generator power dynamics with a variance accounted for of 88%, rebuild the upstream turbine wake with an average normalized root mean squared error of 4% and that controllers can be designed and implemented in a high fidelity simulator for a collective power reference tracking problem.

## Source code
The source code is contained in *KOPMAN_IODMD_1.0*. For simplicity and organisational purposes, the code has been segmented in different folders, according to the final use of a certain function.

Nevertheless, all files should be added to the same directory (either Matlab's or a user created directory) to ensure the program runs smoothly.

## Data sets
Two different data sets are contained in the drive mentioned in *data*, concerning simulations where axial induction control by collectively pitching the blades is used and another where wake redirection control by misaligning the rotor is used. Each folder contains flow field and turbine related information. A more detailed description is present in the folder *data*.

## Application 
The Matlab script *Main.m* contained in *KOPMAN_IODMD_1.0* sets the relevant parameters and calls the appropriate functions. 

The script functionalities are further described in *KOPMAN_IODMD_1.0* folder.

#### Author: Nassir Rodrigues Cassamo: nassir.cassamo@tecnico.ulisboa.pt
#### October 2020




