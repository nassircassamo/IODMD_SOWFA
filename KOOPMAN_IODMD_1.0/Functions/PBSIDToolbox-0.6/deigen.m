function deigen(E,Er)
%DEIGEN Eigenvalues plot
%  deigen(E) plot the system eigenvalues in E.
%
%  deigen(E,Er) also plot the system eigenvalues in Er.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% Plot eigenvalues
n = size(E,1);
N = size(E,2);
[cx,cy] = pol2cart(linspace(0,2*pi,1000),ones(1,1000));
hold on
plot(cx,cy,'Color',[0.4 0.4 0.4]);
if nargin == 4
plot(real(Er),imag(Er),'k+','LineWidth',0.1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',30);
end
for i = 1:N
    plot(real(E(:,i)),imag(E(:,i)),'x','MarkerEdgeColor',[0.7 0.7 0.7],'MarkerFaceColor',[0.7 0.7 0.7],'MarkerSize',10);
end
axis equal
axis([-1 1 -1 1]);
xlabel('Real axis')
ylabel('Imaginary axis')
box on;
hold off
