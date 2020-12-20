%% Functions -- By Category

%% Data Processing
%   fit             - Model fit [%].
%   idmultisine     - Generation of a Multi Sine Signal.
%   idprbs          - Generation of a Pseudo-Random Binary Signal (PRBS).
%   pec             - Prediction error cost.
%   sigscale        - Scaling of identifications signals.
%   snr             - Signal to Noise Ratio (SNR).
%   vaf             - Variance Accounted For [%].

%% LTI Model Identification
%   abc2idss        - Create IDSS model structure.
%   abcd2idss       - Create IDSS model structure.
%   abcdk2idss      - Create IDSS model structure.
%   abck2idss       - Create IDSS model structure.
%   dmodx           - Closed-loop LTI system identification using the PBSIDopt method.
%   dordfir         - Open-loop LTI system identification using the PBSIDopt method.
%   dordvarmax      - Closed-loop LTI system identification using the PBSIDopt method.
%   dordvarx        - Closed-loop LTI system identification using the PBSIDopt method.
%   dvar2eig        - Eigenvalues and its covariance estimation.
%   dvar2frd        - System and its covariance estimation in frequency domain.
%   dvar4abcdk      - Asymptotic variance of the PBSIDopt (VARX only) estimation.
%   dvar4abck       - Asymptotic variance of the PBSIDopt (VARX only) estimation.
%   dvar4varx       - Asymptotic variance of the VARX estimation.
%   dx2abc          - Estimates the matrices A, B, and C of the state space model.
%   dx2abcd         - Estimates the matrices A, B, C, and D of the state space model.
%   dx2abcdk        - Estimates the matrices A, B, C, D, and K of the state space model.
%   dx2abck         - Estimates the matrices A, B, C, and K of the state space model.
%   spaavf          - Spectral analysis with frequency averaging.

%% Hammerstein/Wiener Model Identification
%   hmodx           - Closed-loop Hammerstein system identification using the PBSIDopt method.
%   hordvarx        - Closed-loop Hammerstein system identification using the PBSIDopt method.
%   hwmodx          - Closed-loop Hammerstein-Wiener system identification using the PBSIDopt method.
%   hwordvarx       - Closed-loop Hammerstein-Wiener system identification using the PBSIDopt method.
%   hwx2abcdk       - Estimates the matrices A, B, C, D, and K of the state space model.
%   hwx2abck        - Estimates the matrices A, B, C, and K of the state space model.
%   hx2abcdk        - Estimates the matrices A, B, C, D, and K of the state space model.
%   hx2abck         - Estimates the matrices A, B, C, and K of the state space model.
%   wmodx           - Closed-loop Wiener system identification using the PBSIDopt method.
%   wordvarx        - Closed-loop Wiener system identification using the PBSIDopt method.
%   wx2abcdk        - Estimates the matrices A, B, C, D, and K of the state space model.
%   wx2abck         - Estimates the matrices A, B, C, and K of the state space model.

%% LPV Model Identification
%   lmodx           - Closed-loop LPV system identification using the PBSIDopt method.
%   lordvarx        - Closed-loop LPV system identification using the PBSIDopt method.
%   lpv2ltv         - Transforms a LPV system matrices to a LTV system matrices.
%   lx2abcdk        - Estimates the matrices A, B, C, D, and K of the state space model.
%   lx2abck         - Estimates the matrices A, B, C, and K of the state space model.

%% Periodic-LPV Model Identification
%   plpv2ltv        - Transforms a PLPV system matrices to a LTV system matrices.
%   pmodx           - Closed-loop PLPV system identification using the PBSIDopt method.
%   pordvarx        - Closed-loop PLPV system identification using the PBSIDopt method.
%   pschedclust     - Periodic clustering of scheduling sequence
%   px2abcdk        - Estimates the matrices A, B, C, D, and K of the state space model.
%   px2abck         - Estimates the matrices A, B, C, and K of the state space model.

%% Recursive-LTI Model Identification
%   rpbsid          - Recursive Predictor-Based Subspace Identification.

%% Plot Functions
%   dbode           - Bode diagram.
%   dbodemag        - Bode diagram (magnitude only).
%   dbodemagpatch   - Bode diagram with given error bounds (magnitude only).
%   dbodemagsd      - Bode diagram with probalistic error bounds (magnitude only).
%   dbodepatch      - Bode diagram with given error bounds.
%   dbodesd         - Bode diagram with probalistic error bounds.
%   deigen          - Eigenvalues plot.
%   deigensd        - Eigenvalues plot with probalistic error bounds.
%   dnyquist        - Nyquist diagram.
%   dnyquistdet     - Nyquist diagram (det(G)).
%   dnyquistdetsd   - Nyquist diagram with probalistic error bounds (det(G)).
%   dnyquistsd      - Nyquist diagram with probalistic error bounds.

%% IDAFFLPV Model Functions
%   bodemag         - Mae a Bode magnitude plot for the IDAFFLPV model
%   c2d             - Continuous-to-discrete conversion of state-space models.
%   d2d             - Discrete-to-continuous conversion of state-space models.
%   findstates      - Estimate initial states of the model for a given data set.
%   idafflpv2idss   - Convert IDAFFLPV model to multiple IDSS models.
%   idafflpv2ss     - Convert IDAFFLPV model to multiple SS models.
%   idafflpv        - Construct an IDAFFLPV model structure.
%   idss2idafflpvA  - Convert single IDSS model to IDAFFLPV (constant A) model.
%   idafflpvA2idss  - Convert IDAFFLPV (constant A) model to single IDSS model.
%   iosize          - Returns I/O size of dynamical systems.
%   isct            - True for continuous-time IDAFFLPV models.
%   isdt            - True for discrete-time IDAFFLPV models.
%   isempty         - True for empty IDAFFLPV models.
%   isstatic        - Returns TRUE if model is a pure gain.
%   issiso          - True for SISO IDAFFLPV models.
%   order           - Computes order of IDAFFLPV models
%   pe              - Computes the prediction errors of an IDAFFLPV.
%   pem             - Computes the prediction error estimate of an IDAFFLPV.
%   predict         - Linear response simulation of affine LPV state-space predictor.
%   predstab        - Assess stability of the predictor form of the IDAFFLPV model.
%   resid           - Compute the residuals associated with an IDAFFLPV.
%   sim             - Linear response simulation of affine LPV state-space model.
%   size            - Size of affine LPV state-space system.
%   ssdata          - Returns state-space matrices for IDAFFLPV models.