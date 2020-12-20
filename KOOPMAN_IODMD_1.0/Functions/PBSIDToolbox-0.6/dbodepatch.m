function dbodepatch(G,G_min,G_max,w,h,G_real)
%DBODEPATCH Bode diagram with given error bounds
%  dbodepatch(G,Gmin,Gmax,w,h) plots the bode diagram (magnitude and phase)
%  using the estimated frequency response given in G and the frequency
%  bounds given in Gmin and Gmax, the frequencies w, and sample time h.
%
%  dbodepatch(G,Gmin,Gmax,w,h,Greal) includes the plot of the frequency
%  response given in Greal.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

if nargin > 5
   if ~isequal(size(G,3),size(G_min,3),size(G_max,3),size(G_real,3),length(w))
       error('Number of frequencies points should be the same!')
   end
else
   if ~isequal(size(G,3),size(G_min,3),size(G_max,3),length(w))
       error('Number of frequencies points should be the same!')
   end    
end

l = size(G,1);
r = size(G,2);
for i = 1:l
    for j = 1:r
        subplot(l*2,r,j+(i-1)*2*r);
        hold on;
        patchplot(w./(2*pi),mag2db(squeeze(abs(G(i,j,:))))',mag2db(squeeze(abs(G_min(i,j,:))))',mag2db(squeeze(abs(G_max(i,j,:))))');
        if nargin > 5
            plot(w./(2*pi),mag2db(squeeze(abs(G_real(i,j,:))))','k--');
        end
        vline((pi/h)/(2*pi),'k');
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        logx;
        box on;
        hold off;

        subplot(l*2,r,j+r+(i-1)*2*r);
        hold on;
        patchplot(w./(2*pi),rad2deg(unwrap(squeeze(angle(G(i,j,:)))))',rad2deg(unwrap(squeeze(angle(G_min(i,j,:)))))',rad2deg(unwrap(squeeze(angle(G_max(i,j,:)))))');
        if nargin > 5
            plot(w./(2*pi),rad2deg(unwrap(squeeze(angle(G_real(i,j,:)))))','k--');
        end
        vline((pi/h)/(2*pi),'k');
        xlabel('Frequency (Hz)');
        ylabel('Phase (degrees)');
        logx;
        box on;
        hold off;
    end
end
end





