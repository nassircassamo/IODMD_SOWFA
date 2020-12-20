function display(sys)
%IDAFFLPV/DISPLAY   Pretty-print for IDAFFLPV models.

% Get size
[Ny,Nu] = size(sys);

% Use ISSTATIC to account for delays
StaticFlag = isstatic(sys);

% Handle various types
if ((Ny==0 || Nu==0) && StaticFlag)
   disp(xlate('Empty affine LPV state-space model.'))
else
   % Display name
   disp(xlate(sys.Name))

   % Single IDAFFLPV model
   dispsys(sys,'')
   
   % Display IDAFFLPV properties (sample times)
   if ~StaticFlag,
       if sys.Ts<0,
           disp(xlate('Sampling time: unspecified'))
       elseif sys.Ts>0,
           disp(sprintf('Sampling time: %0.5g',sys.Ts))
       end
   end
   
   % Last line
   if StaticFlag,
      disp(xlate('Static gain.'))
   elseif sys.Ts==0,
      disp(xlate('Continuous-time affine LPV state-space model.'))
   else
      disp(xlate('Discrete-time affine LPV state-space model.'));
   end
end
end

function dispsys(sys,LeftMargin)
%DISPLAY  Pretty-print for affine LPV state-space models.

% Print matrices
printsys(sys.a,sys.b,sys.c,sys.d,sys.k,sys.InputName,sys.OutputName,sys.StateName,sys.SchedulingName,LeftMargin);
end

function printsys(a,b,c,d,k,ulabels,ylabels,xlabels,plabels,LeftMargin)
%PRINTSYS  Print system in pretty format.
%
%   PRINTSYS is used to print state space systems with labels to the
%   right and above the system matrices.
%
%   PRINTSYS(A,B,C,D,K,E,ULABELS,YLABELS,XLABELS) prints the state-space
%   system with the input, output and state labels contained in the
%   cellarrays ULABELS, YLABELS, and XLABELS, respectively.
%
%   PRINTSYS(A,B,C,D,K) prints the system with numerical labels.
%
%   See also: PRINTMAT

nx = size(a,1);
np = size(a,2)/nx;
[ny,nu] = size(d);
nu = nu/np;

if ((isempty(ulabels) || isequal('',ulabels{:})) && (isempty(plabels) || isequal('',plabels{:})))
    for j=1:np
        for i=1:nu,
            if j == 1
                ulabels{i} = sprintf('u%d',i);
            else
                ulabels{(j-1)*nu+i} = sprintf('u%d*p%d',i,j-1);
            end
        end
    end
elseif (isempty(ulabels) || isequal('',ulabels{:}))
    for j=1:np
        if j ~= 1
            if isempty(plabels{j-1})
                plabels{j-1} = '?';
            end
        end
        for i=1:nu,
            if j == 1
                ulabels{i} = sprintf('u%d',i);
            else
                ulabels{(j-1)*nu+i} = strcat(sprintf('u%d',i),'*',plabels{j-1});
            end
        end
    end
else
    for j=1:np
        if j ~= 1
            if isempty(plabels{j-1})
                plabels{j-1} = '?';
            end
        end
        for i=1:nu,
            if j == 1
                if isempty(ulabels{i})
                    ulabels{i} = '?';
                end
            else
                ulabels{(j-1)*nu+i} = strcat(ulabels{i},'*',plabels{j-1});
            end
        end
    end
end

if ((isempty(ylabels) || isequal('',ylabels{:})) && (isempty(plabels) || isequal('',plabels{:})))
    for j=1:np
        for i=1:ny,
            if j == 1
                ylabels{i} = sprintf('y%d',i);
            else
                ylabels{(j-1)*ny+i} = sprintf('y%d*p%d',i,j-1);
            end
        end
    end
elseif (isempty(ylabels) || isequal('',ylabels{:}))
    for j=1:np
        if j ~= 1
            if isempty(plabels{j-1})
                plabels{j-1} = '?';
            end
        end
        for i=1:ny,
            if j == 1
                ylabels{i} = sprintf('y%d',i);
            else
                ylabels{(j-1)*ny+i} = strcat(sprintf('y%d',i),'*',plabels{j-1});
            end
        end
    end
else
    for j=1:np
        if j ~= 1
            if isempty(plabels{j-1})
                plabels{j-1} = '?';
            end
        end
        for i=1:ny,
            if j == 1
                if isempty(ylabels{i})
                    ylabels{i} = '?';
                end
            else
                ylabels{(j-1)*ny+i} = strcat(ylabels{i},'*',plabels{j-1});
            end
        end
    end
end

if ((isempty(xlabels) || isequal('',xlabels{:})) && (isempty(plabels) || isequal('',plabels{:})))
    for j=1:np
        for i=1:nx,
            if j == 1
                xlabels{i} = sprintf('x%d',i);
            else
                xlabels{(j-1)*nx+i} = sprintf('x%d*p%d',i,j-1);
            end
        end
    end
elseif (isempty(xlabels) || isequal('',xlabels{:}))
    for j=1:np
        if j ~= 1
            if isempty(plabels{j-1})
                plabels{j-1} = '?';
            end
        end
        for i=1:nx,
            if j == 1
                xlabels{i} = sprintf('x%d',i);
            else
                xlabels{(j-1)*nx+i} = strcat(sprintf('x%d',i),'*',plabels{j-1});
            end
        end
    end
else
    for j=1:np
        if j ~= 1
            if isempty(plabels{j-1})
                plabels{j-1} = '?';
            end
        end
        for i=1:nx,
            if j == 1
                if isempty(xlabels{i})
                    xlabels{i} = '?';
                end
            else
                xlabels{(j-1)*nx+i} = strcat(xlabels{i},'*',plabels{j-1});
            end
        end
    end
end

disp(' ')
if isempty(a),
    % Gain matrix
    printmat(d,[LeftMargin 'd'],ylabels,ulabels);
else
    printmat(a,[LeftMargin 'a'],xlabels,xlabels);
    printmat(b,[LeftMargin 'b'],xlabels,ulabels);
    printmat(c,[LeftMargin 'c'],ylabels,xlabels);
    printmat(d,[LeftMargin 'd'],ylabels,ulabels);
    printmat(k,[LeftMargin 'k'],xlabels,ylabels);
end
end

function printmat(a,name,rlab,clab)
%PRINTMAT Print matrix with labels.
%   PRINTMAT(A,NAME,RLAB,CLAB) prints the matrix A with the row labels
%   RLAB and column labels CLAB.  NAME is a string used to name the
%   matrix.  RLAB and CLAB are cell vectors of strings.
%
%   See also  PRINTSYS.

CWS = get(0,'CommandWindowSize');      % max number of char. per line
MaxLength = round(.9*CWS(1));

[nrows,ncols] = size(a);
len = 12;    % Max length of labels
Space = ' ';

% Print name
%disp(' ')
if ~isempty(name),
    disp([name,' = ']),
end

% Empty case
if (nrows==0) || (ncols==0),
    if (nrows==0) && (ncols==0),
        disp('     []')
    else
        disp(sprintf('     Empty matrix: %d-by-%d',nrows,ncols));
    end
    disp(' ')
    return
end

% Row labels
RowLabels = strjust(strvcat(' ',rlab{1:nrows}),'left');
RowLabels = RowLabels(:,1:min(len,end));
RowLabels = [Space(ones(nrows+1,1),ones(3,1)),RowLabels];

% Construct matrix display
Columns = cell(1,ncols);
prec = 3 + isreal(a);
for ct=1:ncols,
    clab{ct} = clab{ct}(:,1:min(end,len));
    col = [clab(ct); cellstr(deblank(num2str(a(:,ct),prec)))];
    col = strrep(col,'+0i','');  % xx+0i->xx
    Columns{ct} = strjust(strvcat(col{:}),'right');
end

% Equalize column width
lc = cellfun('size',Columns,2);
lcmax = max(lc)+2;
for ct=1:ncols,
    Columns{ct} = [Space(ones(nrows+1,1),ones(lcmax-lc(ct),1)) , Columns{ct}];
end

% Display MAXCOL columns at a time
maxcol = max(1,round((MaxLength-size(RowLabels,2))/lcmax));
for ct=1:ceil(ncols/maxcol)
    disp([RowLabels Columns{(ct-1)*maxcol+1:min(ct*maxcol,ncols)}]);
    disp(' ');
end
end