function dbodemagsd(G,covG,sd,w,h,G_real)
%DBODEMAGSD Bode diagram with probalistic error bounds (magnitude only)
%  dbodemagsd(G,cogG,sd,w,h) plots the bode diagram (magnitude only)
%  using the estimated frequency response given in G and the frequency
%  bounds given in covG and sd, the frequencies w, and sample time h. The
%  value sd is the standard deviation and is larger than zero.
%
%  dbodemagsd(G,cogG,sd,w,h,Greal) includes the plot of the frequency
%  response given in Greal.

if nargin > 5
   if ~isequal(size(G,3),size(covG,3),size(G_real,3),length(w))
       error('Number of frequencies points should be the same!')
   end
else
   if ~isequal(size(G,3),size(covG,3),length(w))
       error('Number of frequencies points should be the same!')
   end    
end

% from Lennart
sda = real(sqrt((real(G).^2).*covG(:,:,:,1,1) + 2*(real(G).*imag(G)).*covG(:,:,:,2,1) + (imag(G).^2).*covG(:,:,:,2,2)))./abs(G);

l = size(G,1);
r = size(G,2);
for i = 1:l
    for j = 1:r
        subplot(l,r,j+(i-1)*r);
        hold on;
        patchplot(w./(2*pi),mag2db(squeeze(abs(G(i,j,:))))',mag2db(max(squeeze(abs(G(i,j,:)))-sd.*squeeze(sda(i,j,:)),1e-5))',mag2db(max(squeeze(abs(G(i,j,:)))+sd.*squeeze(sda(i,j,:)),1e-5))');
        if nargin > 5
            plot(w./(2*pi),mag2db(squeeze(abs(G_real(i,j,:))))','k--');
        end
        vline((pi/h)/(2*pi),'k');
        xlabel('Frequency (Hz)');
        ylabel('Magnitude (dB)');
        logx;
        box on;
        hold off;
    end
end
end
