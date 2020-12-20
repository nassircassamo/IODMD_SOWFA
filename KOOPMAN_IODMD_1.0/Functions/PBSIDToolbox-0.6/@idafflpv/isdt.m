function boo = isdt(sys)
%IDAFFLPV/ISDT  True for discrete-time IDAFFLPV models.

% SYS is discrete if Ts ~= 0
boo = (isempty(sys.Ts) || sys.Ts~=0 || isstatic(sys));
