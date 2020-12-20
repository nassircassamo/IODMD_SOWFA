
%% Functions -- By Category
%
%% Optimization interface routines
%   doptlti    - Time-domain optimization of discrete-time LTI state
%                space models.
%   foptlti    - Frequency-domain optimization of discrete-time and
%                continuous-time LTI state space models.
%
%% Optimization functions
%   lmmore     - More-Hebden trust-region based implementation of the
%                Levenberg-Marquardt gradient-search optimization 
%                algorithm.
%
%% Cost-functions
%   dfunlti    - Cost-function for doptlti.
%   ffunlti    - Cost-function for foptlti
%
%% (De-)parametrization functions
%   dss2th     - Parametrization of discrete-time models
%   css2th     - Parametrization of continuous-time models
%   dth2ss     - De-parametrization of discrete-time models
%   cth2ss     - De-parametrization of continuous-time models
%
%% Miscellaneous support-functions
%   destmar    - Estimation of a multivariable AR process.
%   cholicm    - Cholesky-factor of the Inverse Covariance Matrix
%                of a multivariable AR process.
%   mkoptstruc - An empty MATLAB 6 optimization options-structure.
%   optim5to6  - MATLAB 5 to 6 optimization options conversion.
%   vaf        - Variance Accounted For between two signals.
%   prbn       - Pseudo random binary noise.
%   gbn        - Generalized binary noise test-signal.
%   shave      - Remove spikes from measured signals.
%
%% Low-level calculation routines
%   ltiitr     - State-trajectory calculation for LTI systems.
%   ltifrf     - Frequency Response Function calculation for LTI systems.
%   simlns     - Similarity-map left null-space for the full parametrization.
%   dltisim    - Simulates a discrete-time LTI state-space system. 
%   ss2frf     - This function calculates the frequency-response function
%
%% Time-domain subspace identification functions
%   dordpo     - Compress data and estimate order for PO-MOESP.
%   dordpi     - Compress data and estimate order for PI-MOESP.
%   dordrs     - Compress data and estimate order for RS-MOESP.
%   dmodpo     - Calculate A and C using PO-MOESP.
%   dmodpi     - Calculate A and C using PI-MOESP.
%   dmodrs     - Calculate A and C using RS-MOESP.
%   dac2b      - Calculate B matrix (assume D is zero).
%   dac2bd     - Calculate B and D matrices.
%   dinit      - Calculate initial condition.
%
%% Frequency-domain subspace identification functions
%   fdordom    - Compress data and estimate order for discrete-time models
%   fcordom    - Compress data and estimate order for continuous-time models
%   fdmodom    - Calculate A and C for discrete-time models
%   fcmodom    - Calculate A and C for continuous-time models
%   fac2b      - Calculate B matrix (assume D is zero)
%   fac2bd     - Calculate B and D matrices.
