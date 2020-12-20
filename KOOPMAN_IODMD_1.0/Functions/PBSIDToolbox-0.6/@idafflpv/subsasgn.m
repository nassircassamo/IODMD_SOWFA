function S = subsasgn(S,subs,A)
%IDAFFLPV/SUBSASGN  Define index assignment for AFFLPV system

sub = subs(1);
type = sub.type;
index = sub.subs;

switch type
    case '.'
        switch lower(index)
            case 'a'
                S.a = A;
            case 'b'
                S.b = A;
            case 'c'
                S.c = A;
            case 'd'
                S.d = A;
            case 'k'
                S.k = A;
            case 'x0'
                S.x0 = A;
            case 'ts'
                S.Ts = A;
            case 'statename'
                S.StateName = A;
            case 'schedulingname'
                S.SchedulingName = A;
            case 'noisevariance'
                S.NoiseVariance = A;
            case 'inputname'
                S.InputName = A;
            case 'outputname'
                S.Outputname = A;
            case 'name'
                S.Name = A;
            otherwise
                error('Variable not known')
        end
        
    otherwise
        error('Invalid field name')
end

end
