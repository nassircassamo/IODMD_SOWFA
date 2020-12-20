echo off
clc
close all;
format compact;
echo on
% This a demo of the LTI Identification Systems Toolbox.
% It illustrates the use of iterative LTI model refinement, both in a
% non-linear least squares sense and in a maximum-likelihood sense.

% A critical step in iterative refinement is finding an accurate
% initial estimate. Subspace identification will be used for this
% estimate.
%
% The first step is loading the measured input/output data.
%
load examples/ltidemo.mat u y;
%
% The data is detrended prior to starting the identification.
% The resulting sequences are plotted in Figures 1 and 2.
%
u=detrend(u);y=detrend(y);
%
echo off
% Plot the data
figure(1); clf;
for i=1:size(u,2),
    subplot(size(u,2),1,i);
    plot(1:size(u,1),u(:,i));
end;
subplot(size(u,2),1,1); title('Input signals');
figure(2); clf;
for i=1:size(y,2),
    subplot(size(y,2),1,i);
    plot(1:size(y,1),y(:,i));
end;
subplot(size(y,2),1,1); title('Output signals');

echo on
pause; % Press any key to start the identification
echo off
% Clear screen any remove figures
clc
close all;

echo on
%
% PO-MOESP subspace identification will be used to generate an
% initial estimate. A block-size of 12 will be used to estimate
% a forth-order model.
%
[S,R]=dordpo(u,y,12);
[Ai,Ci]=dmodpo(R,4);
[Bi,Di]=dac2bd(Ai,Ci,u,y);
%
pause; % Press any key to continue.
echo off
clc;
echo on
%
% In order to determine to quality of the model estimated by PO-MOESP,
% the Variance Accounted For is calculated. This figure of merit indicates,
% for each output, the "match" between the measured output and the output
% simulated using the estimated model. A VAF of 100% indicates a perfect
% match.
%
% Simulate the output
xi=ltiitr(Ai,Bi,u,[],[]);
yi=xi*Ci'+u*Di';
%
% Calculate the VAF
vaf_po=vaf(y,yi)
%
% It is clear that the second output is predicted poorly by the estimated
% model. Figures 1 shows the actual outputs and the outputs predicted by
% PO-MOESP.
%
% Press any key to continue ...
echo off
figure(1); clf;
for i=1:size(y,2),
    subplot(size(y,2),1,i);
    plot(1:size(y,1),y(:,i),1:size(y,1),yi(:,i));
end;
subplot(size(y,2),1,1); title('Measured output (blue) and output predicted by PO-MOESP (green)');

pause;

clc
close all;
echo on
%
% The PO-MOESP model will now be refined using a nonlinear least squares
% algorithm. This may take some time....
%
[Ao,Bo,Co,Do]=doptlti(u,y,Ai,Bi,Ci,Di);
%
% An optimized estimate (Ao,Bo,Co,Do) for the system model has now been
% obtained. Its VAF is calculated and compared to the PO-MOESP VAF.
%
xo=ltiitr(Ao,Bo,u,[],[]);
yo=xo*Co'+u*Do';
vaf_opt=vaf(y,yo);
%
echo off;
fprintf('\nOutput | PO-MOESP VAF | Refined VAF\n');
fprintf('-----------------------------------\n');
fprintf('   %2d  |  %5f   |  %5f\n',1,vaf_po(1),vaf_opt(1));
fprintf('   %2d  |  %5f   |  %5f\n',2,vaf_po(2),vaf_opt(2));
fprintf('-----------------------------------\n');
% Plot signals
figure(1); clf;
for i=1:size(y,2),
    subplot(size(y,2),1,i);
    plot(1:size(y,1),y(:,i),1:size(y,1),yo(:,i));
end;
subplot(size(y,2),1,1); title('Measured output (blue) and output predicted by refined model (green)');

echo on;
%
% It is clear that the second output is predicted much better now. The
% signals are plotted in Figure 1.
%
pause; % Press any key to continue
echo off;

clc;
close all;

echo on;
%
% The model can be refined even further using a maximum-likelihood
% estimation. Such an estimation requires a noise-model, which is
% obtained from the residual of the model obtained by nonlinear
% least squares. A third-order noise-model is estimated:
%
w=y-yo;
[Af,Ab,Sf,Sb]=destmar(w,3);
%
% This model is put into the nonlinear optimization function to
% obtain the maximum likelihood estimate:
%
[Am,Bm,Cm,Dm]=doptlti(u,y,Ao,Bo,Co,Do,[],[],[],[],[Sf,Sb],[Af;Ab]);
%
% VAF results are given in the following table:
%
echo off;
xm=ltiitr(Am,Bm,u,[],[]);
ym=xm*Cm'+u*Dm';
vaf_ml=vaf(y,ym);
fprintf('\nOutput | PO-MOESP VAF | Refined VAF |  Maximum Likelihood VAF\n');
fprintf('---------------------------------------------------------------\n');
fprintf('   %2d  |  %5f   |  %5f   |  %5f\n',1,vaf_po(1),vaf_opt(1),vaf_ml(1));
fprintf('   %2d  |  %5f   |  %5f   |  %5f\n',2,vaf_po(2),vaf_opt(2),vaf_ml(2));
fprintf('---------------------------------------------------------------\n');
echo on
%
% Surprisingly, the maximum likelihood VAF is not as high as the nonlinear
% least squares VAF. However...
%
pause; % Press any key to continue

echo off;
clc;
close all;
echo on;
%
% Although the VAF of the model obtained by maximum likelihood is lower than
% that of the model obtained using nonlinear least squares, the maximum
% likelihood model is not necessarily worse.
%
% Indeed, there are other ways to assess a model. An important aspect of LTI
% models is the location of their poles. These poles can be compared to those
% of the actual model (which will now be loaded from the data file).
%
% The table below shows how accurately the estimated poles match the actual
% system poles. The number is calculated by obtaining the eigenvalues of
% the A-matrices of the estimated models, subtracting the actual system poles,
% and summing the absolute values of these differences.
%
echo off;
load examples/ltidemo.mat A B C D;
pd_po=sum(abs(eig(Ai)-eig(A)));
pd_opt=sum(abs(eig(Ao)-eig(A)));
pd_ml=sum(abs(eig(Am)-eig(A)));
fprintf('\n PO-MOESP   | Nonlinear LS |  Maximum Likelihood\n');
fprintf('-----------------------------------------------\n');
fprintf('  %5f   |  %5f   |  %5f\n',pd_po,pd_opt,pd_ml);
fprintf('-----------------------------------------------\n');
echo on;
%
% It is clear that the model obtained using maximum likelihood yields
% the best estimates for the system's poles.
%
pause; % Press any key to quit
