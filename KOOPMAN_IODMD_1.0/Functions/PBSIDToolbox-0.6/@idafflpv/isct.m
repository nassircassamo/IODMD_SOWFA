function boo = isct(sys)
%IDAFFLPV/ISCT  True for continuous-time IDAFFLPV models.

% SYS is continuous if Ts = 0
boo = isempty(sys.Ts) || sys.Ts==0;
