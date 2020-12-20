function pind = pschedclust(mu,f,p,tol1,tol2)
%PSCHEDCLUST Periodic clustering of scheduling sequence
%  pind=pschedclust(mu,f,p,tol1,tol2) clusters the (quasi) scheduling sequence
%  mu for periodic repetitions of p sized windows. It is intended as a
%  preprocessor for pordvarx. The algortihm uses a divide and conqeur
%  approach for sorting the p windowed vectors. After sorting the vectors
%  are bined in group up to a tolerance tol.
%
%  See also: pordvarx, pmodx, px2abck, px2abcdk.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% assign default values to unspecified parameters
if (nargin < 5) || isempty(tol2)
    tol1 = 1e-3;
end
if (nargin < 5) || isempty(tol1)
    tol2 = 1e-3;
end

% check dimensions of inputs
if size(mu,2) < size(mu,1)
    mu = mu';
end
N = size(mu,2);
s = size(mu,1);

% check the size of the windows
if f > p
    error('Future window size f must equal or smaller then past window p. (f <= p)')
end

% clustering of scheduling periods
k = 0;
pvec = [];
MU = zeros((p+f)*s,N-p-f+1);
for i = 1:N-p-f+1
    d = mu(:,i:i+p+f-1);
    MU(:,i) = d(:);
end

% cluster by divide and conquer
msum = sum(MU,1)./(p+f);
for i = 1:N-p-f+1
    if isempty(pvec) || ~any(pvec == i)
        k = k + 1;
        pvec = [pvec i];
        pind{k,1} = i;
        mind = find((msum >= (msum(i) - tol1)) & (msum <= (msum(i) + tol1)));
        for j = 1:length(mind)
            if i ~= mind(j) && ~any(pvec == mind(j))
                if norm(MU(:,i) - MU(:,mind(j)))/sqrt(p+f) <= tol2
                    pvec = [pvec mind(j)];
                    pind{k,1} = [pind{k,1} mind(j)];
                end
            end
        end
    end
end
MU = [];
end

