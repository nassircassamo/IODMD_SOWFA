function P = predictor(M,type)
%PREDICTOR returns the predictor form of the LPV model M

%  [E,X0] = PE(M,U,Y,T,MU,'Type') specifies the type of LPV predictor. 'K'
%  specifies that K is not dependent on scheduling, or 'CD' specifies
%  that C and D is not dependent on scheduling. Default is Type = 'CD'.

[A,B,C,D,K] = getABCDK(M);
[Ny,Nu,Nx,Np] = size(M);

if strcmpi(type,'CD')
    for i = 1:Np+1
        A(:,(i-1)*Nx+1:i*Nx) = A(:,(i-1)*Nx+1:i*Nx) - K(:,(i-1)*Ny+1:i*Ny)*C(:,1:Nx);
        B(:,(i-1)*Nu+1:i*Nu) = B(:,(i-1)*Nu+1:i*Nu) - K(:,(i-1)*Ny+1:i*Ny)*D(:,1:Nu);
    end
elseif strcmpi(type,'K')
    for i = 1:Np+1
        A(:,(i-1)*Nx+1:i*Nx) = A(:,(i-1)*Nx+1:i*Nx) - K(:,1:Ny)*C(:,(i-1)*Nx+1:i*Nx);
        B(:,(i-1)*Nu+1:i*Nu) = B(:,(i-1)*Nu+1:i*Nu) - K(:,1:Ny)*D(:,(i-1)*Nu+1:i*Nu);
    end
else
    error('Type not recognized!')
end

P = idafflpv(A,B,C,D,K,M.x0,M.Ts);

end

