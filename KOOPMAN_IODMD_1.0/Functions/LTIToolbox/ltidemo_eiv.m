clc
echo on 
% This is errors-in-variables of the LTItoolbox
% This file shows how you may use the errors-in-variables (EIV) functions 
% of the SMI toolbox to perform closed-loop identification 
% 
% The model structure for the EIV problem is an innovations model:
% 
%     x(k+1)  = Ax(k) + Bu~(k) + f(k)
%     y~(k)   = Cx(k) + Du~(k) 
%     
% with noisy input u(k) and noisy output y(k)
% 
%     u(k) = u~(k) + w(k)
%     y(k) = y~(k) + v(k)
%     
% available for identification. 
% Here f(k), w(k) and v(k) are zero-mean white noise sequences
% independent of the noise-free input u~(k).
%
% Press any key to continue
pause 

clc 
% The plant is a mechanical system consisting of 2 rotating discs
% and an integrator. Its transfer function (order = 5) is 
Pnum = [0 0.98 12.99 18.59 3.30 -0.02]*1e-3; 
Pden = [1 -4.4 8.09 -7.83 4 -0.86]; 
P = tf(Pnum,Pden,1);

% A controller has been designed based on H_infinity method
Cnum = [0.61 -2.03 2.76 -1.83 0.49];
Cden = [1 -2.65 3.11 -1.75 0.39];
C = tf(Cnum,Cden,1);

% The process noise is shaped by the filter 
Fnum = [0 2.89 11.13 2.74]*0.01;
Fden = [1 -2.7 2.61 -0.9];
F = tf(Fnum,Fden,1);

% An external input, a white noise of standard deviation (std) 1, 
% is used to excite the system. The white noise driving the process
% noise filter $F(z)$ has a std 1/9. The input and output noise are
% both of std 0.01. This results in SNRs of approximately 20dB at 
% the plant output and 5dB at the plant input. 

% Press any key to continue
pause 

clc 
% Form closed-loop system 

% Build the combined deterministic and stochastic system 
% with 4 inputs and 2 outputs where 
% Input 1: Plant input
% Input 2: Process noise
% Input 3: Input measurement noise 
% Input 4: Ouput measurement noise 
% Output 1: Noisy plant output
% Output 2: Noisy plant input 
PF = parallel(P,F,[],[],1,1);
[apf,bpf,cpf,dpf] = ssdata(ss(PF));
ac = apf;
bc = [bpf zeros(8,2)];
cc = [cpf ; zeros(1,8)];
dc = [dpf 0 1; 1 0 1 0];
PFE = ss(ac,bc,cc,dc,1);

% Build the feedback system 
CL = feedback(PFE,C,1,1); 

% Press any key to continue
pause 

clc 
% Perform simulation
N = 1200; 
r = randn(N,1);               % external reference input 
pn = (1/3)*randn(N,1);        % process noise 
in = 0.01*randn(N,1);         % input measurement noise
on = 0.01*randn(N,1);         % output measurement noise 
ycl = lsim(CL,[r pn in on]);
u = ycl(:,2);
y = ycl(:,1); 

% Plots of simulated input and output
figure 
subplot(211)
plot(u)
xlabel('sample index')
title('plant input u')
subplot(212)
plot(y)
xlabel('sample index')
title('plant output y')

% Now we have the data needed for identification

% Press any key to continue
pause 

clc 
% The data available for identification is plant input (u), 
% plant output (y) and the external excitation (r). 
% 
% The first step in identification is to find out about the order 
% of the system 
% 

[S,R] = dordeiv(u,y,r,20);

% An indicator of the order is given by the singular values in S
clf 
semilogy(S,'x')
title('singular values') 

% The singular value plot indicates that the order of the system
% is 7. The expected order is 8 which is the sum of the orders
% of the plant and the process noise shaping filter. It means that 
% a pole cannot be recovered from identification and this pole 
% belongs to the shaping filter. The information contents on the 
% process noise shaping filter is usually very little in 
% closed-loop as the controller is designed to reject this signal.

% Press any key to continue
pause 

clc 
% We proceed to estimate the system matrices 
% First the A and C matrices 

[A,C] = dmodeiv(R,7);

% and then B and D matrices

% Note: SMI toolbox has 2 functions to estimate (B,D) from closed-loop
%       data, eiv_bd works best for reference signal whose first half
%       is exactly the same as its second half. Otherwise, it may be 
%       best to use the routine cl_bd to estimate (B,D). This is more
%       elaborate and is explained in the example in the manual. 

[B,D] = dac2bd_eiv(A,C,u,y,r);

% Press any key to continue
pause 

clc 
% Finally, we compare the frequency responses of the estimated 
% and true systems 
[mt,pt,w] = bode(P);
[me,pe] = bode(ss(A,B,C,D,1),w);
[xt,yt] = pol2cart(pt(:)/180*pi,mt(:));
[xe,ye] = pol2cart(pe(:)/180*pi,me(:));
[ep,em] = cart2pol(xe-xt,ye-yt);
clf 
loglog(w,mt(:),w,me(:),w,em)
xlabel('frequency')
legend('True FRF','Estimated FRF','Error')

% End of LTIdemo_eiv
echo off 


