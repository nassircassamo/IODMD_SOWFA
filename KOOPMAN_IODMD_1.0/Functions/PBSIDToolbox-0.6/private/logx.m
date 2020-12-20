function logx
%LOGX  Turn the X axis to LOG
%   LOGX turns the X axis of the current plot to log coordinates.
%
%   See also LINX, LINY, LOGY.

set(gca,'XScale','log');
end