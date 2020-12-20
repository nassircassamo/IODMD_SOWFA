function boo = isempty(sys)
%IDAFFLPV/ISEMPTY  True for empty IDAFFLPV models.

boo = any(size(sys)==0);

