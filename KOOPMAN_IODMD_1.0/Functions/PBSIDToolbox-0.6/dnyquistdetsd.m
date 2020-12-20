function dnyquistdetsd(G,covG,sd,Greal)
%DNYQUISTDETSD Nyquist diagram with probalistic error bounds (det(G))
%  dnyquistdetsd(G,covG,sd) plots the nyquist diagram using the determinant
%  of the estimated frequency response given in G and the frequency bounds
%  given in covG and sd. The value sd is the standard deviation and is
%  larger than zero.
%
%  dnyquistdetsd(G,covG,sd,Greal) includes the plot of the determinant of
%  the frequency response given in Greal

if nargin > 3
   if ~isequal(size(G,3),size(covG,3),size(Greal,3))
       error('Number of frequencies points should be the same!')
   end
else
   if ~isequal(size(G,3),size(covG,3))
       error('Number of frequencies points should be the same!')
   end    
end

l = size(G,1);
r = size(G,2);
Gdet = zeros(1,size(G,3));
covGdet = zeros(2,2,size(G,3));
if nargin > 3
    Gdetreal = zeros(1,size(G,3));
end
for k = 1:size(G,3)
   Gdet(1,k) = det(G(:,:,k)); 
   if nargin > 3
       Gdetreal(1,k) = det(Greal(:,:,k));
   end
       
   X = zeros(2*l,r);
   X(1:2:end,:) = real(G(:,:,k));
   X(2:2:end,:) = imag(G(:,:,k));
   J = jacobianest(@(x) deter(x,l,r),X(:));
   P = zeros(l*r*2);
   for i = 1:l
       for j = 1:r
           P((i-1)*r*2+(j-1)*2+(1:2),(i-1)*r*2+(j-1)*2+(1:2)) = squeeze(covG(i,j,k,:,:));
       end
   end
   covGdet(:,:,k) = J*P*J';
end

if l == r
    figure
    hold on;   
    [sdreal,sdimag] = pol2cart(unwrap(angle(Gdet)),max(squeeze(abs(Gdet)),1e-5));
    for k = 1:size(G,3)
        ellipsebnd(covGdet(:,:,k),[sdreal(k); sdimag(k)],'conf',erf(sd/sqrt(2)),'style','k')
    end
    plot(sdreal',sdimag','k','Linewidth',2);
    if nargin > 3
        [sdreal,sdimag] = pol2cart(unwrap(angle(Gdetreal)),max(squeeze(abs(Gdetreal)),1e-5));
        plot(sdreal',sdimag','k--','Linewidth',1);
    end
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

function dD = deter(G,l,r)
Gr = G(1:2:end,:);
Gi = G(2:2:end,:);
G = reshape(Gr,l,r) + 1i.*reshape(Gi,l,r);
D = det(G);
dD = zeros(2,1);
dD(1,1) = real(D);
dD(2,1) = imag(D);
end    
