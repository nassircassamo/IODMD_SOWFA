function dbode(frd,w,h)
%DBODE Bode diagram
%  dbode(G,w,h) plots the bode diagram (magnitude and phase) using the
%  estimated frequency resonse given in G, the frequencies w, and sample
%  time h.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

% frd = unwrap(frd);
l = size(frd,1);
r = size(frd,2);
for i = 1:l
    for j = 1:r
        subplot(l*2,r,j+(i-1)*2*r);
        hold on;
        plot(w./(2*pi),mag2db(squeeze(abs(frd(i,j,:)))),'k','linewidth',2);
        vline((pi/h)/(2*pi),'k');
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        logx;
        box on;
        hold off;
        
        subplot(l*2,r,j+r+(i-1)*2*r);
        hold on;
        plot(w./(2*pi),unwrap(rad2deg(squeeze(angle(frd(i,j,:))))),'k','linewidth',2);
        vline((pi/h)/(2*pi),'k');
        xlabel('Frequency (Hz)');
        ylabel('Phase (degrees)');
        logx;
        box on;
        hold off;
    end
end
end
