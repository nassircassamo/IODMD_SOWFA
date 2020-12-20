function inov = abcd2idss(data,A,B,C,D)
%ABCD2IDSS Create IDSS model structure
%  M = ABCD2IDSS(IDDATA,A,B,C,D) returnes as a model structure object
%  describing the discrete-time model:
%  
%     x[k+1] = A x[k] + B u[k];             x[0] = X0
%       y[k] = C x[k] + D u[k] + e[k]
%
%  See also: IDSS.


% get sizes
n = size(A,1);
r = size(B,2);
l = size(C,1);
if isempty(D)
    D = zeros(l,r);
end

% store state-space matrices in inovation system object
if iscell(data.Ts)
    Ts = data.Ts{1};
else
    Ts = data.Ts;
end
inov = idss(A,B,C,D,zeros(l,r),zeros(n,1),Ts,'InputName',data.Inputname,...
    'OutputName',data.Outputname,'DisturbanceModel','Zero',...
    'InitialState','Auto','SSParameterization','Free');

% estimate the prediction error and initial state
[e,x0] = pe(inov,data,'e');
inov.x0 = x0(:,1);
lambda = reshape(covf(e,1),l,l);
inov = pvset(inov,'NoiseVariance',lambda);

% store information about estimation 
est = pvget(inov,'EstimationInfo');
if iscell(e.OutputData)
    N = 0;
    for i = 1:length(e.OutputData)
        N = N + length(data.InputData{i});
    end
else
    N = length(data.InputData);
end
est.DataLength = N;
est.DataTs = data.Ts;
est.LossFcn = det(lambda);
npar = n^2 + n*r + n*l + l*n;
est.FPE = det(lambda)*(1+2*npar/N);
dn = data.Name;
if isempty(dn)
    dn = inputname(1);
end
est.DataName = dn;
est.DataDomain = 'Time';
est.Status = 'Estimated model (PBSIDopt)';
est.Method = 'PBSIDopt';
if ~isempty(wtxt)
    est.Warning = wtxt;
end
est.InitialState = 'Estimate';
inov = pvset(inov,'EstimationInfo',est);



