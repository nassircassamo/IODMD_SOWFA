function boo = isstatic(sys)
%IDAFFLPV/ISSTATIC Returns TRUE if model is a pure gain.

boo = isempty(sys.a);
