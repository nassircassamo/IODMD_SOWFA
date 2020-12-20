function dnyquistdet(G)
%DNYQUISTDET Nyquist diagram (det(G))
%  dnyquistdet(G) plots the nyquist diagram using the determinant of the
%  estimated frequency response given in G.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

l = size(G,1);
r = size(G,2);
Gdet = zeros(1,size(G,3));
for i = 1:size(G,3)
   Gdet(i) = det(G(:,:,i)); 
end

if l == r
    figure
    hold on;
    [sdreal,sdimag] = pol2cart(unwrap(squeeze(angle(Gdet))),max(squeeze(abs(Gdet)),1e-5));
    plot(sdreal',sdimag','k','Linewidth',2);
    %axis equal;
    vline(0,'k:');
    hline(0,'k:');
    xlabel('Real axis');
    ylabel('Imaginary axis');
    box on;
    
    hold off;
else
    error('Only for square MIMO systems!')
end
end


