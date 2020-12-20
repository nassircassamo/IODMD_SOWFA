function A = subsref(S,subs)
%IDAFFLPV/SUBSREF  Define subscripted reference for IDAFFLPV models.

sub = subs(1);
type = sub.type;
index = sub.subs;

switch type
    case '()'
        switch length(index)
            case 2
                [Ny,Nu,Nx,Np] = size(S);
                if strcmp(index{1},':')
                    i = 1:Ny;
                elseif strcmp(index{1},'end')
                    i = Ny;
                else
                    if index{1} > Ny
                        error('Index out of bounds')
                    end
                    i = index{1};
                end
                
                if strcmp(index{2},':')
                    j = 1:Nu;
                elseif strcmp(index{2},'end')
                    j = Nu;
                else
                    if index{2} > Nu
                        error('Index out of bounds')
                    end
                    j = index{2};
                end
                A = S;
                T = zeros(Nx,length(j)*(Np+1));
                for h = 1:(Np+1)
                    T(:,(h-1)*Nu+1:(h-1)*Nu+length(j)) = S.b(:,(h-1)*Nu+j);
                end
                A.b = T;
                A.c = S.c(i,:);
                T = zeros(length(i),length(j)*(Np+1));
                for h = 1:Np
                    T(:,(h-1)*Nu+1:(h-1)*Nu+length(j)) = S.d(i,(h-1)*Nu+j);
                end
                A.d = T;
                T = zeros(Nx,length(i)*(Np+1));
                for h = 1:Np
                    T(:,(h-1)*Ny+1:(h-1)*Ny+length(i)) = S.k(:,(h-1)*Ny+i);
                end
                A.k = T;
                A.NoiseVariance = A.NoiseVariance(i,i);
                A.InputName = S.InputName(j);
                A.OutputName = S.OutputName(i);
                
            otherwise
                error('Index out of bounds')
        end
        
    case '.'
        switch lower(index)
            case 'a'
                A = S.a;
            case 'b'
                A = S.b;
            case 'c'
                A = S.c;
            case 'd'
                A = S.d;
            case 'k'
                A = S.k;
            case 'x0'
                A = S.x0;
            case 'ts'
                A = S.Ts;
            case 'statename'
                A = S.StateName;
            case 'noisevariance'
                A = S.NoiseVariance;
            case 'schedulingname'
                A = S.SchedulingName;
            case 'inputname'
                A = S.InputName;
            case 'outputname'
                A = S.Outputname;
            case 'name'
                A = S.Name;
            otherwise
                error('Variable not known')
        end
        
end
end


