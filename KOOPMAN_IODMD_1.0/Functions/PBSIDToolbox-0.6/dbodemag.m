function dbodemag(frd,w,h)
%DBODEMAG Bode diagram (magnitude only)
%  dbodemag(G,w,h) plots the bode diagram (magnitude only) using the
%  estimated frequency response given in G, the frequencies w, and sample
%  time h.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

l = size(frd,1);
r = size(frd,2);
for i = 1:l
    for j = 1:r
        subplot(l,r,j+(i-1)*r)
        hold on;
        plot(w./(2*pi),mag2db(squeeze(abs(frd(i,j,:))))','k','linewidth',2);
        vline((pi/h)/(2*pi),'k');
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        logx;
        box on;
        hold off;
    end
end
end
