function sys = idafflpv(varargin)
%IDAFFLPV Construct an IDAFFLPV model structure
%
%   M = IDAFFLPV(A,B,C,D)
%   M = IDAFFLPV(A,B,C,D,K,X0,Ts)
%   M = IDAFFLPV(A,B,C,D,K,X0,Ts,'Property',Value,..)
%
%   M: returned as a model structure object describing the linear parameter
%   varying discrete-time model
%
%     x(k+1) = A kron(mu(k),x(k)) + B kron(mu(k),u(k)) + K kron(mu(k),e(k))
%     y(k)   = C kron(mu(k),x(k)) + D kron(mu(k),u(k)) + e(k); x(0) = X0.
%
%   A,B,C,D and K are the LPV state-space matrices. X0 is the initial
%   condition, if any, and Ts is the sampling time. For Ts == 0, a
%   continuous-time model is constructed. The outputs are the linear
%   parameter-varying matrices A, B, C, D,  K, where matrices have the form
%   of A=[A(1) A(2) ... A(m)].

ni = nargin;
if ni > 0 && isa(varargin{1},'idafflpv'),
    % Quick exit for IDAFFLPV(SYS) with SYS of class IDAFFLPV
    if ni == 1
        sys = varargin{1};
    else
        error('Use SET to modify the properties of IDAFFLPV objects.');
    end
    return
end

% Dissect input list
DataInputs = 0;
for i=1:ni
    nextarg = varargin{i};
    if ~(ischar(nextarg) || iscell(nextarg))
        DataInputs = DataInputs+1;
    end
end

% Process numerical data
switch DataInputs,
    case 0
        if ni,
            error('Too many LPV arguments or missing numerical data.');
        else
            % Empty model
            a = [];  b = [];  c = [];  d = []; k = [];
        end

    case 1
        % Gain matrix
        a = [];  b = [];  c = []; k = [];
        d = checkMatrixData(varargin{1},'D');

    case 2
        error('Undefined C and D matrices.');

    case 3
        error('Undefined D matrix.');
        
    case 4
        % A,B,C,D specified: validate data
        a = checkMatrixData(varargin{1},'A');
        b = checkMatrixData(varargin{2},'B');
        c = checkMatrixData(varargin{3},'C');
        d = checkMatrixData(varargin{4},'D');
        k = zeros(size(a,1),size(a,2)/size(a,1)*size(c,1));
        x0 = zeros(size(a,1),1);
    case 5
        % A,B,C,D,K specified: validate data
        a = checkMatrixData(varargin{1},'A');
        b = checkMatrixData(varargin{2},'B');
        c = checkMatrixData(varargin{3},'C');
        d = checkMatrixData(varargin{4},'D');
        k = checkMatrixData(varargin{5},'K');
        x0 = zeros(size(a,1),1);
        if isempty(k)
            k = zeros(size(a,1),size(a,2)/size(a,1)*size(c,1));
        end
    otherwise
        % A,B,C,D,K specified: validate data
        a = checkMatrixData(varargin{1},'A');
        b = checkMatrixData(varargin{2},'B');
        c = checkMatrixData(varargin{3},'C');
        d = checkMatrixData(varargin{4},'D');
        k = checkMatrixData(varargin{5},'K');
        x0 = varargin{6};
        if isempty(k)
            k = zeros(size(a,1),size(a,2)/size(a,1)*size(c,1));
        end
        if isempty(x0)
            x0 = zeros(size(a,1),1);
        end
end

% Sample time
if DataInputs==7
    % Discrete SS
    try
        Ts = ltipack.utValidateTs(varargin{7});
    catch
        rethrow(lasterror);
    end
else
    Ts = 0;
end

% Array size
if ni > 0
    [Ns,Ny,Nu,Np] = checkMatrixDims(a,b,c,d,k,x0);
else
    Ny = 0;
    Nu = 0;
    Np = 0;
end

% Process char data
sa = '';
sn(1:Ns,1) = {''};
sp(1:Np,1) = {''};
su(1:Nu,1) = {''};
sy(1:Ny,1) = {''};
for i = DataInputs+1:2:ni
    if strcmpi(varargin{i},'StateName')
        nextarg = varargin{i+1};
        if iscell(nextarg)
            if length(sn) ~= Ns
                error('Wrong number of state names!')
            end
            if size(nextarg,1) < size(nextarg,2)
                sn = nextarg';
            else
                sn = nextarg;
            end
        elseif ischar(nextarg)
           sn{1} = nextarg; 
        else
            error('Type not recognized!')
        end
    elseif strcmpi(varargin{i},'SchedulingName')
        nextarg = varargin{i+1};
        if iscell(nextarg)
            if length(sp) ~= Np
                error('Wrong number of state names!')
            end
            if size(nextarg,1) < size(nextarg,2)
                sp = nextarg';
            else
                sp = nextarg;
            end
        elseif ischar(nextarg)
           sp{1} = nextarg; 
        else
            error('Type not recognized!')
        end
    elseif strcmpi(varargin{i},'InputName')
        nextarg = varargin{i+1};
        if iscell(nextarg)
            if length(su) ~= Nu
                error('Wrong number of input names!')
            end
            if size(nextarg,1) < size(nextarg,2)
                su = nextarg';
            else
                su = nextarg;
            end
        elseif ischar(nextarg)
           su{1} = nextarg; 
        else
            error('Type not recognized!')
        end
    elseif strcmpi(varargin{i},'OutputName')
        nextarg = varargin{i+1};
        if iscell(nextarg)
            if length(sy) ~= Ny
                error('Wrong number of output names!')
            end
            if size(nextarg,1) < size(nextarg,2)
                sy = nextarg';
            else
                sy = nextarg;
            end
        elseif ischar(nextarg)
           sy{1} = nextarg; 
        else
            error('Type not recognized!')
        end
    elseif strcmpi(varargin{i},'Name')
        nextarg = varargin{i+1};
        if iscell(nextarg)
            sa = nextarg{1};
        elseif ischar(nextarg)
            sa = nextarg;
        else
            error('Type not recognized!')
        end
    end
end

% Define class-specific properties
% RE: a,b,c,d,StateName,SchedulingName all virtual
superiorto('double')
sys = struct;

% Add afflpv system data to sys structure
sys.a = a;
sys.b = b;
sys.c = c;
sys.d = d;
sys.k = k;
sys.Ts = Ts;
sys.StateName = sn;
sys.SchedulingName = sp;
sys.InputName = su;
sys.OutputName = sy;
sys.Name = sa;
sys.NoiseVariance = eye(Ny,Ny);
sys.x0 = x0;

% Create AFFLPV system
sys = class(sys,'idafflpv');

