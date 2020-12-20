function H = bhankel(R,N)
% hankel matrix construction
% R is a block-row vector and 
% N the number of block rows

% Rufus Fraanje, April 2005, visit ISVR


[m,NSMP]=size(R);

L = NSMP - N + 1;

H = zeros(m*N,L);

for i=1:N,
  H((i-1)*m+1:i*m,:) = R(:,i:i+L-1);
end;

