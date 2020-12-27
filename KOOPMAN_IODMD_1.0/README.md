# Source Code

The *Main.m* script calls the relevant functions. 
The main functionalities of the program include:

## 0/1. Initialise/Assess data
The data sets can be easily visualised, *i.e*., snapshots and animations showing the upstream turbine wake behaviour as wind turbine control parameters change and plot showing the wind turbines variables changing during the simulations (control variables, forces, generator power, rotor speed). 

It goes through the following steps:
* **1. Create directory**
* **2. Load identification (turbine related) data set**
* **3. Load validation (turbine related) data set**
* **4. Define program parameters**
* **5. Define simulation parameters**
* **6. Visualise turbine related data sets**
* **7. Visualise turbine and wake related data sets**

## 2. Dynamic Mode Decomposition
The relevant information is preprocessed and used to identify/compute the Reduce Order Models by using a suitable Dynamic Mode Decomposition algorithm, such as Input Output Dynamic Mode Decomposition.

It goes through the following steps:
* **1. Define time interval**
* **2. Read and pre process input-output identification data**
* **3. Read and pre process input-output validation data**
* **4. Define high order states for identification and validation**
* **5. Use Dynamic Mode Decomposition to compute models**

## 3. Validation

## 4. Dynamical Analysis

## 5. Rebuild

## 6. Model Predictive Control

## Functions