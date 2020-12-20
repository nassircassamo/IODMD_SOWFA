function [G,covG] = dvar2frd(P,w,h,varargin)
%DVAR2FRD System and its covariance estimation in frequency domain
%  [G,covG]=dvar2frd(P,w,h,p,VARX) returns the system and covariance in the
%  frequency domain for the frequencies w and sample time h.
%
%  [G,covG]=dvar2frd(P,w,h,p,A,B,C,K) returns the system and covariance in
%  the frequency domain for the frequencies w and sample time h.
%
%  [G,covG]=dvar2frd(P,w,h,p,A,B,C,D,K) returns the system and covariance
%  in the frequency domain for the frequencies w and sample time h.

%  Ivo Houtzager
%  Delft Center of Systems and Control
%  Delft University of Technology 
%  The Netherlands, 2010

if nargin == 5
    p = varargin{1};
    VARX = varargin{2};
    m = floor(size(VARX,2)/p);
    l = size(VARX,1);
    r = m-l;
    noD = ~(size(VARX,2)/m > p);
    N = length(w);
    G = zeros(l,r+l,N);
    covG = zeros(l,r+l,N,2,2);

    for q = 1:N
        for j = 0:p
            if j == 0
                if noD ~= 1
                    G(:,:,q) = [VARX(:,p*m+1:p*m+r) eye(l)];
                else
                    G(:,:,q) = [zeros(l,r) eye(l)];
                end
             else
                G(:,:,q) = G(:,:,q) + [exp(-j*w(q)*1i*h).*VARX(:,(p-j)*m+1:(p-j)*m+r) -exp(-j*w(q)*1i*h).*VARX(:,(p-j)*m+r+1:(p-j)*m+r+l)];
            end
        end
        G(:,:,q) = [pinv(G(:,r+1:r+l,q))*G(:,1:r,q) pinv(G(:,r+1:r+l,q))];
        
        
        %T = kron([eye(r) zeros(r,l); -G(:,r+1:r+l,q)*G(:,r+1:r+l,q)*G(:,1:r,q) -G(:,r+1:r+l,q)].',G(:,r+1:r+l,q));
        T = kron([eye(r) zeros(r,l); zeros(l,r) -G(:,r+1:r+l,q)].',G(:,r+1:r+l,q));
        gradVARX = zeros(l*(l+r),p*l*(l+r)+~noD*r*l);
        for c = 1:l
            for v = 1:(r+l)
                gradc = zeros(l,1);
                gradc(c) = 1;
                gradv = zeros(l*(l+r),1); 
                gradv((v-1)*l+1:v*l,1) = gradc;
                gradVARX((v-1)*l+c,1:p*l*(l+r)) = kron(exp(-(p:-1:1).*w(q)*1i*h)',gradv)';
                if ~noD
                    if v <= r
                        gradv = zeros(l*r,1);
                        gradv((v-1)*l+1:v*l,1) = gradc;
                        gradVARX((v-1)*l+c,p*l*(l+r)+(1:l*r)) = gradv(1:l*r);
                    else
                        gradVARX((v-1)*l+c,p*l*(l+r)+(1:l*r)) = zeros(l*r,1);
                    end
                end

            end
        end
        gradVARX = T*gradVARX;
        for c = 1:l
            for v = 1:(r+l)
                covG(c,v,q,1,1) = real(gradVARX((v-1)*l+c,:))*P*real(gradVARX((v-1)*l+c,:))';
                covG(c,v,q,1,2) = real(gradVARX((v-1)*l+c,:))*P*imag(gradVARX((v-1)*l+c,:))';
                covG(c,v,q,2,1) = imag(gradVARX((v-1)*l+c,:))*P*real(gradVARX((v-1)*l+c,:))';
                covG(c,v,q,2,2) = imag(gradVARX((v-1)*l+c,:))*P*imag(gradVARX((v-1)*l+c,:))';
            end
        end
    end
        
elseif nargin == 7
    A = varargin{1};
    B = varargin{2};
    C = varargin{3};
    K = varargin{4};
    l = size(C,1);
    n = size(C,2);
    r = size(B,2);
    B = [B K];
    D = [zeros(l,r) eye(l)];
    s = exp(w*1i*h);
    N = length(s);
    gradABKC = zeros(n^2+n+n,N,l,r);
    G = zeros(l,r+l,N);
    covG = zeros(l,r+l,N,2,2);
    
    for q = 1:N
        Ainv = (s(q).*eye(n) - A)^(-1);
        dCr = Ainv*B;
        dBr = (C*Ainv)';
        dAr = zeros(n^2,l,l+r);
        for j = 1:n
            for b = 1:n
                sigma = zeros(n,n);
                sigma(b,j) = 1;
                for c = 1:l
                    for v = 1:(r+l)
                        dAr(b+(j-1)*n,c,v) = C(c,:)*Ainv*sigma*Ainv*B(:,v);
                    end
                end
            end
        end
        
        for c = 1:l
            for v = 1:(r+l)
                gradABKC(:,q,c,v) = conj([dAr(:,c,v); dBr(:,c); dCr(:,v)]);
                ind = [1:n^2 n^2+(v-1)*n+1:n^2+n*v n^2+(r+l)*n+c:l:n^2+(r+l)*n+n*l];
                covG(c,v,q,1,1) = real(gradABKC(:,q,c,v))'*P(ind,ind)*real(gradABKC(:,q,c,v));
                covG(c,v,q,1,2) = real(gradABKC(:,q,c,v))'*P(ind,ind)*imag(gradABKC(:,q,c,v));
                covG(c,v,q,2,1) = imag(gradABKC(:,q,c,v))'*P(ind,ind)*real(gradABKC(:,q,c,v));
                covG(c,v,q,2,2) = imag(gradABKC(:,q,c,v))'*P(ind,ind)*imag(gradABKC(:,q,c,v));
                G(c,v,q) = D(c,v) + C(c,:)*Ainv*B(:,v);
            end
        end
    end
    
    
elseif nargin == 8
    A = varargin{1};
    B = varargin{2};
    C = varargin{3};
    D = varargin{4};
    K = varargin{5};
    l = size(C,1);
    n = size(C,2);
    r = size(B,2);
    B = [B K];
    D = [D eye(l)];
    s = exp(w*1i*h);
    N = length(s);
    gradABKCD = zeros(n^2+n+n+1,N,l,r);
    G = zeros(l,r+l,N);
    covG = zeros(l,r+l,N,2,2);
    
    for q = 1:N
        Ainv = (s(q).*eye(n) - A)^(-1);
        dCr = Ainv*B;
        dBr = (C*Ainv)';
        dAr = zeros(n^2,l,l+r);
        for j = 1:n
            for b = 1:n
                sigma = zeros(n,n);
                sigma(b,j) = 1;
                for c = 1:l
                    for v = 1:(r+l)
                        dAr(b+(j-1)*n,c,v) = C(c,:)*Ainv*sigma*Ainv*B(:,v);
                    end
                end
            end
        end
        
        for c = 1:l
            for v = 1:(r+l)
                if v > r
                    gradABKCD(:,q,c,v) = conj([dAr(:,c,v); dBr(:,c); dCr(:,v); 0]);
                    ind = [1:n^2 n^2+(v-1)*n+1:n^2+n*v n^2+(r+l)*n+c:l:n^2+(r+l)*n+n*l];
                    covG(c,v,q,1,1) = real(gradABKCD(:,q,c,v))'*blkdiag(P(ind,ind),0)*real(gradABKCD(:,q,c,v));
                    covG(c,v,q,1,2) = real(gradABKCD(:,q,c,v))'*blkdiag(P(ind,ind),0)*imag(gradABKCD(:,q,c,v));
                    covG(c,v,q,2,1) = imag(gradABKCD(:,q,c,v))'*blkdiag(P(ind,ind),0)*real(gradABKCD(:,q,c,v));
                    covG(c,v,q,2,2) = imag(gradABKCD(:,q,c,v))'*blkdiag(P(ind,ind),0)*imag(gradABKCD(:,q,c,v));
                else
                    gradABKCD(:,q,c,v) = conj([dAr(:,c,v); dBr(:,c); dCr(:,v); 1]);
                    ind = [1:n^2 n^2+(v-1)*n+1:n^2+n*v n^2+(r+l)*n+c:l:n^2+(r+l)*n+n*l n^2+(r+l)*n+n*l+(v-1)*l+c];
                    covG(c,v,q,1,1) = real(gradABKCD(:,q,c,v))'*P(ind,ind)*real(gradABKCD(:,q,c,v));
                    covG(c,v,q,1,2) = real(gradABKCD(:,q,c,v))'*P(ind,ind)*imag(gradABKCD(:,q,c,v));
                    covG(c,v,q,2,1) = imag(gradABKCD(:,q,c,v))'*P(ind,ind)*real(gradABKCD(:,q,c,v));
                    covG(c,v,q,2,2) = imag(gradABKCD(:,q,c,v))'*P(ind,ind)*imag(gradABKCD(:,q,c,v));
                end

                G(c,v,q) = D(c,v) + C(c,:)*Ainv*B(:,v);
            end
        end
    end
else
    error('DVAR2RFD requires five, seven or eight input arguments.')
end

end


















