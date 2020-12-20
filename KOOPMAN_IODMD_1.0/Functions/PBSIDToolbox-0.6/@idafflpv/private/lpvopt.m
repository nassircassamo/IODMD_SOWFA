function [tho,OPTIONS,CostFunction,JACOB] = lpvopt(th,thc,d,u,y,p,parain)
%LPVOPT     Optimize LPV system represented by a parameter vector.
%           Optimize the following LPV system represented by a 
%           parameter vector created by lpv2par with respect to an 
%           output error criterium, using an iterative Gauss-Newton
%           or Levenberg-Marquardt algorithm.
%
%           x(k+1) = A(:,1:n)*x(k) + A(:,n+1:n*(q+1))*kron(p(k),x(k))
%                  + B(:,1:m)*u(k) + B(:,m+1:m*(q+1))*kron(p(k),u(k))
%           y(k)   = C(:,1:n)*x(k) + C(:,n+1:n*(q+1))*kron(p(k),x(k))
%                  + D(:,1:m)*u(k) + D(:,m+1:m*(q+1))*kron(p(k),u(k))
%
% Syntax: 
%           tho = lpvopt(th,thc,dim,u,y,p)
%           [tho,options,CostFun,Jacob] = lpvopt(th,thc,dim,u,y,p,options)
% 
% Input:
%           th       Parameter vector created by lpv2par containing
%                    the initial starting point.
%           thc      Parameter vector created by lpv2par containing
%                    zeros for the parameters that are optimized 
%                    and ones for the parameters that are kept constant.
%           dim      Dimension vector created by lpv2par. 
%           y        N x l matrix containing N data points of the l
%                    outputs.
%           u        N x m matrix containing N data points of the m
%                    inputs.
%           p        N x q matrix containing N data points of the q
%                    time varying parameters.
%           options  Matlab 5 compatible foptions vector.
%
% Output:
%           tho      The optimal parameter vector.
%           options  Matlab 5 compatible foptions vector.
%           CostFun  The value of the output error criterium at 
%                    the solution th. 
%           Jacob    The Jacobian of the output error criterium at 
%                    the solution TH.
% 
% Remarks:
%           This a modified version of the NLSQ function of the
%           Matlab 5 optimization toolbox to solve nonlinear
%           least squares problems. It uses a projected gradient to
%           deal with the nonuniqueness created by similarity
%           transformations in the identification of LPV sytems.
%           The default algorithm is the Levenberg-Marquardt method
%           with a  mixed quadratic and cubic line search procedure.
%
% See also: lpv2par, par2lpv.
%
% Written by Vincent Verdult, February 2001.

% ------------Initialization----------------

if nargin<7
  parain=[];
end
if size(th,2)>size(th,1)
  th=th';
end
if size(thc,2)>size(thc,1)
  thc=thc';
end
if size(d,2)>size(d,1)
  d=d';
end
if size(th,2)~=1
  error('th must be a vector.')
end
if size(thc,2)~=1
  error('thc must be a vector.')
end
if isempty(thc);
  thc=zeros(size(th,1),1);
elseif size(th,1)~=size(thc,1)
  error('th and thc must have the same number of rows.')
end
if size(d,2)~=1
  error('dim must be a vector.')
end
n=d(1);
m=d(2);
l=d(3);
s=d(4);
if size(th,1)~=(n+l)*(n+m)*(s+1)
  error('tn and dim are not compatible.')
end

N=size(u,1);
if size(u,2)~=m
  error('Number of colums in u conflicts with dim.')
end
if size(y,1)~=N
  error('u and y must have the same number of rows.')
end
if size(y,2)~=l
  error('Number of colums in y conflicts with dim.')
end
if size(p,1)~=N
  error('u and p must have the same number of rows.')
end
if size(p,2)~=s
  error('Number of colums in p conflicts with dim.')
end

% thc contains ones for the constant parameters
% i.e. the parameters that are not optimized

% convert original parameter vector th
% ind contains the index of the non-constant parameters
% only the non-constant parameters are copied into x.
ind=find(thc<1);
x=th(ind);

% optimization routine
XOUT = x(:);
[nvars] = length(XOUT);
how = [];
nfun = N*d(3);  % N*l

% Set default options
sizep=length(parain);
OPTIONS=zeros(1,18);
OPTIONS(1:sizep)=parain(1:sizep);
default_options=[0,1e-4,1e-4,1e-6,0,0,0,0,0,0,0,0,0,200,0,1e-8,0.1,0];
for i=1:18
  if OPTIONS(i)==0
    OPTIONS(i)=OPTIONS(i)+default_options(i);
  end
end
if OPTIONS(14)==0 
   OPTIONS(14)=length(XOUT)*100;
end 

% Determine display style
if OPTIONS(1)>0
  disp('')
  if isinf(OPTIONS(1))
     disp('f-COUNT      RESID    STEP-SIZE      GRAD/SD        LAMBDA LINE-SEARCH')
  else
     disp('f-COUNT      RESID    STEP-SIZE      GRAD/SD        LAMBDA')
  end
end
PCNT = 0;
EstSum=0.5;
GradFactor=1; 
CHG = 1e-7*abs(XOUT)+1e-7*ones(nvars,1);

OPTIONS(10)=1;
status=-1;

% --- Main loop ---   
while (status~=1)
   OPTIONS(11)=OPTIONS(11)+1;
   if OPTIONS(9)
      % numerical approximation of gradient
      GRAD=zeros(nvars,nfun);
      [CostFunction,Vperp]=lpvfun(XOUT,th,d,ind,u,y,p,1);
      OLDF=CostFunction;
      CHG = sign(CHG+eps).*min(max(abs(CHG),OPTIONS(16)),OPTIONS(17));
      for gcnt=1:nvars
         temp = XOUT(gcnt);
         XOUT(gcnt) = temp +CHG(gcnt);
         x(:) = XOUT;
         CostFunction = lpvfun(XOUT,th,d,ind,u,y,p,2);
         CostFunction = CostFunction(:);
         GRAD(gcnt,:)=(CostFunction-OLDF)'/(CHG(gcnt));
         XOUT(gcnt) = temp;
      end
      CostFunction = OLDF;
      if nfun==1
         if GRAD==0
            CHG = inf;
         else
            CHG = nfun*1e-8./GRAD;
         end
      else
         sumabsGRAD = sum(abs(GRAD),2);
         ii = (sumabsGRAD == 0);
         CHG(ii) = inf;
         CHG(~ii) = nfun*1e-8./sumabsGRAD(~ii);
      end
   else
      % use analytical gradient
      [CostFunction,Vperp,GRAD]=lpvfun(XOUT,th,d,ind,u,y,p,0);
   end
   Vnvars=size(Vperp,2);
   VGRAD=Vperp'*GRAD;
   GradF=2*VGRAD*CostFunction;
   HessF=VGRAD*VGRAD';
   NewF = CostFunction'*CostFunction;
   % --- Initial step and search direction ---
   if status==-1
      OLDX=XOUT;
      MATX=zeros(3,1);
      MATL=[NewF;0;0];
      if cond(VGRAD)>1e8
         GradFactor=norm(VGRAD)+1; 
         how='COND';
      end
      SD=-pinv(HessF+GradFactor*(eye(Vnvars,Vnvars)))*GradF;
      FIRSTF=NewF;
      GDOLD=GradF'*SD;
      % OPTIONS(18) controls the initial starting step-size.
      if OPTIONS(18) == 0, 
         OPTIONS(18)=1; 
      end
      if OPTIONS(1)>0
         disp([sprintf('%5.0f %12.6g %12.3g ',OPTIONS(10),NewF,OPTIONS(18)),sprintf('%12.3g  ',GDOLD)]);
      end
      XOUT=XOUT+OPTIONS(18)*Vperp*SD;   % first step
      estf=VGRAD'*SD+CostFunction;
      GradFactor=estf'*estf;
      SD=-pinv(HessF+GradFactor*(eye(Vnvars,Vnvars)))*GradF;
      XOUT=XOUT+OPTIONS(18)*Vperp*SD;   % second step
      estf=VGRAD'*SD+CostFunction;
      EstSum=estf'*estf;
      status=0;
      if OPTIONS(7)==0; 
         PCNT=1; 
      end
      % --- End of initial step and search direction ---
   else
      % --- Step and search direction update ---
      gdnew=GradF'*SD;
      if OPTIONS(1)>0, 
         num=[sprintf('%5.0f %12.6g %12.3g ',OPTIONS(10),NewF,OPTIONS(18)),sprintf('%12.3g  ',gdnew)];
      end
      % --- Case 1: Interpolation ---
      if gdnew>0 && NewF>FIRSTF
         how='inter';  % interpolation
         [stepsize]=cubici1(NewF,FIRSTF,gdnew,GDOLD,OPTIONS(18));
         OPTIONS(18)=0.9*stepsize;
      % --- Case 2: Update stepsize and search direction ---
      elseif NewF<FIRSTF
         [newstep,fbest] =cubici3(NewF,FIRSTF,gdnew,GDOLD,OPTIONS(18));
         if fbest>NewF,
            fbest=0.9*NewF; 
         end 
         if gdnew<0
            how='incstep';
            if newstep<OPTIONS(18),  
               newstep=(2*OPTIONS(18)+1e-4); 
               how=horzcat(how,'IF'); 
            end
            OPTIONS(18)=abs(newstep);
         else 
            if OPTIONS(18)>0.9
               how='int_step';
               OPTIONS(18)=min([1,abs(newstep)]);
            end
         end
         % --- Levenberg-Marquardt method --- 
         % EstSum is the estimated sum of squares.
         % GradFactor is the value of lambda.
         % Estimated Residual:
         if EstSum>fbest
            GradFactor=GradFactor/(1+OPTIONS(18));
         else
            GradFactor=GradFactor+(fbest-EstSum)/(OPTIONS(18)+eps);
         end
         SD=-(HessF+GradFactor*(eye(Vnvars,Vnvars)))\GradF; 
         OPTIONS(18)=1; 
         estf=VGRAD'*SD+CostFunction;
         EstSum=estf'*estf;
         if OPTIONS(1)>0, 
            num= horzcat(num,sprintf('%12.6g ',GradFactor)); 
         end
         % --- End of Levenberg-Marquardt method ---
         gdnew=GradF'*SD;
         
         OLDX=XOUT;
         % Save Variables
         FIRSTF=NewF;
         GDOLD=gdnew;    
         
         % If quadratic interpolation set PCNT
         if OPTIONS(7)==0, 
            PCNT=1; MATX=zeros(3,1); MATL(1)=NewF; 
         end
      % --- Case 3: Reduce stepsize ---
      else 
         how='Red_Step';
         if NewF==FIRSTF,
            if OPTIONS(1)>0
               disp('No improvement in search direction: Terminating')
            end
            status=1;
         else
            OPTIONS(18)=OPTIONS(18)/8;
            if OPTIONS(18)<1e-8
               OPTIONS(18)=-OPTIONS(18);
            end
         end
      end
      % --- End of cases ---
      XOUT=OLDX+OPTIONS(18)*Vperp*SD;
      if isinf(OPTIONS(1))
         disp([num,how])
      elseif OPTIONS(1)>0
         disp(num)
      end
   end 
   % --- End of step and search direction update ---
   if OPTIONS(7)==0, 
      PCNT=1; MATX=zeros(3,1);  MATL(1)=NewF; 
   end
   % --- Check termination --- 
   if max(abs(Vperp*SD))< OPTIONS(2) && (GradF'*SD) < OPTIONS(3)
      if OPTIONS(1) > 0
         disp('Optimization Terminated Successfully')  
		 if max(abs(GradF)) > 10*(OPTIONS(3)+OPTIONS(2))
			disp('WARNING: Gradient > 10*(TolX+TolFun)')
		 end
      end
      status=1; 
   elseif OPTIONS(10)>OPTIONS(14)
      disp('maximum number of iterations has been exceeded');
      if OPTIONS(1)>0
         disp('Increase OPTIONS(14)')
      end
      status=1;
   else
      % --- Line search --- 
      % mixed polynomial interpolation and extrapolation.
      if PCNT~=0
         % initialize OX and OLDF 
         OX = XOUT; OLDF = CostFunction;
         while PCNT > 0 
            CostFunction=lpvfun(XOUT,th,d,ind,u,y,p,2);
            OPTIONS(10)=OPTIONS(10)+1;
            NewF = CostFunction'*CostFunction;
            % <= used in case when no improvement found.
            if NewF <= OLDF'*OLDF, 
               OX = XOUT; 
               OLDF=CostFunction; 
            end
            [PCNT,MATL,MATX,steplen,NewF,how]=searchq(PCNT,NewF,OLDX,MATL,MATX,Vperp*SD,GDOLD,OPTIONS(18),how);
            OPTIONS(18)=steplen;
            XOUT=OLDX+steplen*Vperp*SD;
            if NewF==FIRSTF,  
               PCNT=0; 
            end
         end
         XOUT = OX;
         CostFunction=OLDF;
      else
         CostFunction=lpvfun(XOUT,th,d,ind,u,y,p,2);
         OPTIONS(10)=OPTIONS(10)+1;
      end
      % --- End of line search ---
   end
   % --- End of check termination ---
end
% --- End of main loop ---
OPTIONS(8) = NewF;
XOUT=OLDX;
JACOB = GRAD.';

% convert optimized parameters XOUT to full parameter vector
% including the constant parameters

tho=th;
tho(ind)=XOUT;

% --- End ---



