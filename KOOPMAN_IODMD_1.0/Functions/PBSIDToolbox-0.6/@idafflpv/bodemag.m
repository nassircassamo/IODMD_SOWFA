function [h] = bodemag(sys,p,w,displaytype,hr,Y,U)
%  BODEMAG(SYS,P,W,DISPLAYTYPE,HR,U,Y) can be used to make a Bode magnitude 
%  plot of the IDAFFLPV model sys for a frequency range W and a scheduling 
%  parameter range P.
%  P is an N-by-m data matrix where N is the number of samples, m
%  is the number of scheduling parameters. W is the frequency range.
%  displaytype is 'colors', 'greyscale' or 'lines'
%  HR are the units for the frequency: 'hz' for Hertz, default is rad/s
%  In U,Y we can select the inputs and output channels, by entering a
%  vector of input and output channel indices. default:  Y = 1:Ny; U =
%  1:Nu+Ny.
if nargin<5 || isempty(hr)
    hr = 'rads';
end

[Ny,Nu,Ns,Np] = size(sys);
if nargin<6
    Y = 1:Ny;
    U = 1:Nu+Ny;
end  

if Np~=size(p,2)
    error('P should be an N-by-m data matrix where N is the number of samples, m is the number of scheduling parameters.')
end

sys = idafflpv2ss(sys,p);
N = size(p,1);

if strcmpi(displaytype,'lines')
    leg = cell(N,1);
    count = 0;
    for u = U
        for y = Y
            count = count + 1;
            subplot(length(Y),length(U),count)
            hold on
            opts = bodeoptions;
            opts.PhaseVisible = 'off';
            if strcmpi(hr,'hz') 
                opts.FreqUnits = 'Hz';
            end
            for i = 1:N        
                if strcmpi(hr,'hz')
                    bodeplot(sys(y,u,i,1),w*2*pi,opts);
                else
                    bodeplot(sys(y,u,i,1),w,opts);
                end
                leg{i} = ['\mu = ',num2str(p(i,:))];
            end
            legend(leg,'Location','Best')
        end
    end
    h = gcf;
elseif strcmpi(displaytype,'colors') || strcmpi(displaytype,'greyscale')
    
    if size(p,2)>1
        warning('bodemag:multiplescheduling','You entered scheduling vectors in BODEMAG with displaytype=colors. The magnitude is plotted against frequency and the FIRST scheduling parameter. This makes sense if other schedulding parameters are functions of the first scheduling parameter.');
    end
    
    if strcmpi(hr,'hz')
        H = freqresp(sys,w*2*pi);
    else
        H = freqresp(sys,w);
    end
    count = 0;
    for u = U
        for y = Y
            count = count + 1;
            subplot(length(Y),length(U),count)
            if strcmpi(displaytype,'colors')
                colormap jet
            else strcmpi(displaytype,'greyscale')
                colormap gray
            end
            handle = surf(w,p(:,1),db(abs(squeeze(H(y,u,:,:))))');
            shading interp
            view(90,-90)
            axis tight
            set(handle,'edgecolor','none');
            ylabel('Scheduling \mu', 'Interpreter', 'Tex')
            if strcmpi(hr,'hz')
                xlabel('Frequency (Hz)')
            else
                xlabel('Frequency (rad/s)')
            end
            if isempty(sys.Outputname{y})
                outputname = ['Out(',num2str(y),')'];
            else
                outputname = sys.Outputname{y};
            end
            if isempty(sys.Inputname{u})
                inputname = ['In(',num2str(u),')'];
            else
                inputname = sys.Inputname{u};
            end
            title([inputname,' to ', outputname])
            t = colorbar('peer',gca);
            set(get(t,'ylabel'),'String', 'Magnitude [dB]');
            h = gcf;
       end
    end
else
    error 'type not recognized'
end

end
