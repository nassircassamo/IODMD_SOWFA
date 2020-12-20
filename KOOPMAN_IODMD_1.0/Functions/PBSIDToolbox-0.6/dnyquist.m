function dnyquist(G)
%DNYQUIST Nyquist diagram
%  dnyquist(G) plots the nyquist diagram using the estimated
%  frequency response given in G.


%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

l = size(G,1);
r = size(G,2);
for i = 1:l
    for j = 1:r
        subplot(l,r,j+(i-1)*r);
        hold on;
        [sdreal,sdimag] = pol2cart(unwrap(squeeze(angle(G(i,j,:)))),max(squeeze(abs(G(i,j,:))),1e-5));
        plot(sdreal',sdimag','k','Linewidth',2);
        %axis equal
        vline(0,'k:');
        hline(0,'k:');
        xlabel('Real axis');
        ylabel('Imaginary axis');
        box on;
        hold off;
    end
end
end
