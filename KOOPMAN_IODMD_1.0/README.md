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
All computed Reduced Order Models are validated. Each one is simulated and using the Variance Accounted For criteria their fitness is assessed. A figure displaying the fitness of each model, with increasing number of modes used, given the identification and validation data set is then presented. The relevant variables are all saved to a specified directory, so that the models can be loaded afterwards.

* **1. Validate Reduced Order Models**
* **2. Present results**
* **3. Save results**

## 4. Dynamical Analysis
The dynamical characteristics of each model can then be analysed by evaluating the natural frequencies of each mode (given a certain Reduced Order Model), the location of the system's poles within the complex plane, the relative importante of each mode, etc. Other functions which attempt to compute a state space representation of each mode have been developed (based on the modal state space representation approach), although their applications and results have not been fully exploited and are not present in the existing documentation. The functions developed also allow for the 2D and 3D representation of the modes.

* **1. Select best performing model to analyse**
* **2. Compute dynamical properties**
* **3. Visualise dynamical modes**

## 5. Rebuild
The capabilities of the Reduced Order Models in terms of reconstruction the full wake are then assessed. The state space trajectory can be represented in the high dimensional space  and the functions here present allow to compare the modes reconstruction to the identification and validation data retrieved in SOWFA, as well as to quantify deviations between the two.

* **1. Compute high order representation of state space trajectory**
* **2. Compare linear simulation results with validation data**
* **3. Evaluate deviations**

## 6. Model Predictive Control
The actuator dynamics are modelled as a low pass filter and joined in series to the chosen  Reduced Order Model, reducing the feedthrough term of the state space model to zero. Several types of Model Predictive Controllers can be designed, such as Input Output (IO) - where the optimal control action is computed at each time - or Incremental Input Output (IIO) - where each optimal incremental control action is computed at each time. The problem is formulated as a Quadratic Programming problem and solved with MATLAB inbuilt tools.

* **1. Model actuator dynamics**
* **2. Build Model Predictive Control framework**

## Functions
These include several functions needed by the program, such as:
* **1. Export Figure**: figures are continuously exported automatically to pre specified directories
* **2. cbrewer**: specified colours for wake snapshots
* **3. LTIToolbox**: functions from the Delft Center for Systems and Control (DCSC) regarding Linear Times Invariant systems
* **4. postProcessingTools**: function developed by the National Renewable Energy Laboratory (NREL) to read turbine information computed when using SOWFA.
