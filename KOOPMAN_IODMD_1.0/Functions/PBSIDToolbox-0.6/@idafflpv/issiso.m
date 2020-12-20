function boo = issiso(sys)
%IDAFFLPV/ISSISO  True for SISO IDAFFLPV models.

sizes = size(sys);
boo = all(sizes(1:2) == 1);

