function dnyquistsd(G,covG,sd,G_real)
%DNYQUISTSD Nyquist diagram with probalistic error bounds
%  dnyquistsd(G,covG,sd) plots the nyquist diagram using the estimated
%  frequency response given in G and the frequency bounds given in covG and
%  sd. The value sd is the standard deviation and is larger than zero.
%
%  dnyquistsd(G,covG,sd,Greal) includes the plot of the frequency response
%  given in Greal.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

if nargin > 3
   if ~isequal(size(G,3),size(covG,3),size(G_real,3))
       error('Number of frequencies points should be the same!')
   end
else
   if ~isequal(size(G,3),size(covG,3))
       error('Number of frequencies points should be the same!')
   end    
end

l = size(G,1);
r = size(G,2);
for i = 1:l
    for j = 1:r
        subplot(l,r,j+(i-1)*r);
        hold on;
        [sdreal,sdimag] = pol2cart(unwrap(squeeze(angle(G(i,j,:)))),max(squeeze(abs(G(i,j,:))),1e-5));
        for k = 1:size(G,3)
            ellipsebnd(squeeze(covG(i,j,k,:,:)),[sdreal(k); sdimag(k)],'conf',erf(sd/sqrt(2)),'style','k')
        end
        plot(sdreal',sdimag','k','Linewidth',2);
        if nargin > 3
            [sdreal,sdimag] = pol2cart(unwrap(squeeze(angle(G_real(i,j,:)))),max(squeeze(abs(G_real(i,j,:))),1e-5));
            plot(sdreal',sdimag','k--','Linewidth',1);
        end
        %axis equal;
        vline(0,'k:');
        hline(0,'k:');
        xlabel('Real axis');
        ylabel('Imaginary axis');
        box on;
        hold off;
    end
end
end
