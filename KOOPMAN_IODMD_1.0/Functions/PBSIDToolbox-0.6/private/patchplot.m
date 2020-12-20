function [Pa,Li,t] = patchplot(varargin)
%PATCHPLOT plots errorbars around a given curve
%
%     [Pa,Li,t] = PATCHPLOT(x,y,L,U,'r','g')
%     Plots a gray Jackknife around the line displayed in black very useful for
%     funky nature style error bars which are shaded.
%
%         Pa is a patch object for more help on patch objects see below
%         Li is a line object, more help on line object is available in MATLAB
%
%     USAGE :
%              1)   [Pa,Li] = patchplot(x,y,E)
%                    Calculates the Lower and upper errorbars as
%                    L = Y-E and U = Y+E. It then takes a default gray color
%                    as patch color, and the line color as black and plots
%                    it around the line using a patch object.
%
%              2)   [Pa,Li] = patchplot(x,y,E,LineColor,PatchColor)
%                    Calculates the Lower and upper errorbars as
%                    L = Y-E and U = Y+E. It then takes PatchColor
%                    as patch color, and the Line Color from the LineColor
%                    variable. It then plots it around the line
%                    using a patch object.
%
%              3)   [Pa,Li] = patchplot(x,y,L,U)
%                    User Supplied bounds are taken as L and U, It then takes
%                    a default gray color as patch color, and the line color
%                    as black and plots it around the line using a patch object.
%
%              4)   [Pa,Li] = patchplot(x,y,L,U,LineColor,PatchColor)
%                    User Supplied bounds are taken as L and U, It then takes
%                    PatchColor as patch color, and the Line Color from the LineColor
%                    variable. It then plots it around the line using a
%                    patch object.
%      CAVEATS
%                 1) Can be Slow sometimes for length(Array) > 10000,
%                 2) Needs better vectorization
%      EXAMPLE
%                         t = [-5:0.05:5];
%                         Y = sin(t);
%                         E = 0.4*rand(1,length(t));
%                         [Pa,Li] = patchplot(t,Y,E);
%                         xlabel('time');
%                         ylabel('Amplitude');
%                         title('Using Errors alone');
% 
%                         figure;
%                         L = Y - E;
%                         U = Y + E;
%                         [Pa,Li] = patchplot(t,Y,L,U);
%                         xlabel('time');
%                         ylabel('Amplitude');
%                         title('Using Lower and Upper Confidence Intervals');
%                         hold on;
% 
%                         Y1 = 2*Y;
%                         L = Y1 - 0.2;
%                         U = Y1 + 0.2;
%                         [Pa,Li] = patchplot(t,Y1,L,U,[255 51 51]./255,[255 153 102]./255);
%                         hold on;
%                         [Pa,Li] = patchplot(t,Y1*2,E,[51 51 153]./255,[102 153 204]./255);
%                         [Pa,Li] = patchplot(t,Y1*2,E,'r','g');
% 
% See also ERRORBAR, PATCH, LINE

% Version 0.001 Chandramouli Chandrasekaran (Chandt) - 13 April 2006.

switch(nargin)
    case 3,
        % If there are 3 inputs it means its just the errors
        x = varargin{1};
        y = varargin{2};
        E = varargin{3};
        L = y - E;
        U = y + E;
        LineColor = 'k';
        PatchColor = [0.85 0.85 0.85];
    case 4,
        % If there are 4 inputs it means they entered the Lower and upper
        % bounds
        x = varargin{1};
        y = varargin{2};
        L = varargin{3};
        U = varargin{4};
        LineColor = 'k';
        PatchColor = [0.85 0.85 0.85];
    case 5,
        x = varargin{1};
        y = varargin{2};
        E = varargin{3};
        L = y - E;
        U = y + E;
        LineColor =  varargin{4};
        PatchColor = varargin{5};
    case 6,
        % If there are 6 inputs then
        x = varargin{1};
        y = varargin{2};
        L = varargin{3};
        U = varargin{4};
        LineColor =  varargin{5};
        PatchColor = varargin{6};
end
tic
Xcoords = [x x(end:-1:1)];
Ycoords = [U L(end:-1:1)];

Pa = patch(Xcoords,Ycoords,PatchColor);
set(Pa,'linestyle','none','linewidth',2);
hold on;
Li = plot(x,y,'color',LineColor,'linewidth',2);
hold on;
t = toc;
end




