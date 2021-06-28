function [f,LambdaDiag, P, phi,damping,b,fnatural]=dynamicalanalysis(sys_red, U, S,V, dt,X_p,X, method,mn,scaling,D,Uups,Deterministic,r,dirdmd,n,Xd,best)

%% DMD STANDARD MODAL ANALYSIS

if method==2 %ioDMD
    
    [W,Lambda]=eig(sys_red{mn}.A);
    LambdaDiag=diag(Lambda);
    omega=log(LambdaDiag)/dt/2/pi;
    f=abs(imag(omega))*D/Uups; %convertion from Hz to Strouhal number
    fnatural=sqrt(abs(imag(omega)).^2+abs(real(omega)).^2);
    damping=-cos(atan2(imag(log(LambdaDiag)),real(log(LambdaDiag))));
    
    % do the same for best validation model 
    [WBest,LambdaBest]=eig(sys_red{best}.A);
    LambdaDiagBest=diag(LambdaBest);
    omegaBest=log(LambdaDiagBest)/dt/2/pi;
    fBest=abs(imag(omegaBest))*D/Uups; %convertion from Hz to Strouhal number
    fnaturalBest=sqrt(abs(imag(omegaBest)).^2+abs(real(omegaBest)).^2);
    dampingBest=-cos(atan2(imag(log(LambdaDiagBest)),real(log(LambdaDiagBest))));
    
    if scaling==0
        %phi=X_p*V*inv(S)*W;
        phi=U(:,1:mn)*W;
        b=phi(1:n,:)\X(1:n,1);
        P=abs(b);
        
        phiBest=U(:,1:best)*WBest;
        bBest=phiBest(1:n,:)\X(1:n,1);
        PBest=abs(bBest);
        
    elseif scaling==1 %common scaling by eigenvalues
        Ahat=(S(1:mn,1:mn)^(-1/2)) * sys_red{mn}.A * (S(1:mn,1:mn)^(1/2));
        [What,Dhat]=eig(Ahat);
        W_r=(S(1:mn,1:mn)^(1/2))*What;
        %Phi=U*W_r;
        phi=X_p*V(:,1:mn)/S(1:mn,1:mn)*W_r;
        P=(diag(phi'*phi));
    end

%% STATES HAVE BEEN EXTENDED TO ACCOUNT FOR DETERMINISTIC VARIABLES
elseif method==3 %extioDMD
    fromhere=size(Xd,1)+1;
    [W,Lambda]=eig(sys_red{mn}.A);
    LambdaDiag=diag(Lambda);
    omega=log(LambdaDiag)/dt/2/pi;
    ft=abs(imag(omega))*D/Uups; %convertion from Hz to Strouhal number
    dampingt=-cos(atan2(imag(log(LambdaDiag)),real(log(LambdaDiag))));
    f=ft(fromhere:end);
    damping=dampingt(fromhere:end);
    
    %do the same for best validation model
    [WBest,LambdaBest]=eig(sys_red{best}.A);
    LambdaDiagBest=diag(LambdaBest);
    omegaBest=log(LambdaDiagBest)/dt/2/pi;
    ftBest=abs(imag(omegaBest))*D/Uups; %convertion from Hz to Strouhal number
    dampingtBest=-cos(atan2(imag(log(LambdaDiagBest)),real(log(LambdaDiagBest))));
    fBest=ftBest(fromhere:end);
    dampingBest=dampingtBest(fromhere:end);


    if scaling==0
        %phi=X_p*V*inv(S)*W;
        phii=blkdiag(eye(size(Deterministic,1),size(Deterministic,1)), U(1:n,1:mn))*W;
        bt=phii\[Deterministic(:,1);X(1:n,1)];
        phi=phii(fromhere:end,fromhere:end);
        b=bt(fromhere:end);
        P=abs(b);
        
        %do the same for best validation model 
        phiiBest=blkdiag(eye(size(Deterministic,1),size(Deterministic,1)), U(1:n,1:best))*WBest;
        btBest=phiiBest\[Deterministic(:,1);X(1:n,1)];
        phiBest=phiiBest(fromhere:end,fromhere:end);
        bBest=btBest(fromhere:end);
        PBest=abs(bBest);
        
    elseif scaling==1
        Ahat=(S(1:mn,1:mn)^(-1/2)) * sys_red{mn}.A(fromhere:end,fromhere:end) * (S(1:mn,1:mn)^(1/2));
        [What,Dhat]=eig(Ahat);
        W_r=(S(1:mn,1:mn)^(1/2))*What;
        %Phi=U*W_r;
        phi=X_p(1:n,:)*V(:,1:mn)/S(1:mn,1:mn)*W_r;
        P=(diag(phi'*phi));  
        
        %do same for best model 
        AhatBest=(S(1:best,1:best)^(-1/2)) * sys_red{best}.A(fromhere:end,fromhere:end) * (S(1:best,1:best)^(1/2));
        [WhatBest,Dhat]=eig(AhatBest);
        W_rBest=(S(1:best,1:best)^(1/2))*WhatBest;
        %Phi=U*W_r;
        phiBest=X_p(1:n,:)*V(:,1:best)/S(1:best,1:best)*W_rBest;
        PBest=(diag(phiBest'*phiBest));  
        
    end

end

%% EIGENVALUE VISUALISATION 

% IN COMPLEX PLANE
figure460=figure('Position', [100 100 600 300],'Visible','off');
set(gcf,'color','w','Position', get(0, 'Screensize'));
subplot(1,2,1)
p=plot(LambdaDiag, 'o');
p.Color=[0.2 0.2 1];
p.MarkerSize=10;
p.MarkerEdgeColor=[0 0 0];
p.MarkerFaceColor=[.2 .2 1];
rectangle('Position', [-1 -1 2 2], 'Curvature', 1,'EdgeColor', 'k', 'LineStyle', '--');
axis(1.2*[-1 1 -1 1])
axis square
xlabel('Real axis \Re')
ylabel('Imaginary axis \Im')
title('Eigenvalue \lambda visualisation on the complex plane')
set(gca, 'FontSize', 14)
grid on
grid minor
hold on
%plot the best model eigenvalues ontop
p2=plot(LambdaDiagBest, 'o');
p2.Color=[1 0.2 .2];
p2.MarkerSize=10;
p2.MarkerEdgeColor=[0 0 0];
p2.MarkerFaceColor=[1 .2 .2];
legend('Highest Order Model','Highest FIT Model','Location','bestoutside','Orientation','horizontal')
legend('boxoff')

% PER FREQUENCY OF OSCILLATION
subplot(1,2,2)
p3=plot(omega, 'o');
p3.Color=[0.2 0.2 1];
p3.MarkerSize=10;
p3.MarkerEdgeColor=[0 0 0];
p3.MarkerFaceColor=[.2 .2 1];
line([0 0], 0.3*[-1 1], 'Color', 'k','LineStyle', '--'); 
xlabel('Frequency \Omega [Hz]')
ylabel('Imaginary axis \Im')
title('Eigenvalue visualisation in frequencies per oscillation')
axis square
set(gca, 'FontSize', 14)
grid on
grid minor
hold on
p4=plot(omegaBest, 'o');
p4.Color=[1 0.2 .2];
p4.MarkerSize=10;
p4.MarkerEdgeColor=[0 0 0];
p4.MarkerFaceColor=[1 .2 .2];
legend([p3 p4],'Highest Order Model','Highest FIT Model','Location','bestoutside','Orientation','horizontal')
legend('boxoff')
export_fig(figure460,strcat(dirdmd,'/image','polescomplexplane'),'-nocrop','-m2');



%% POWER SPECTRUM 
 figure470=figure('Position', [100 100 600 300],'Visible','off');
 set(gcf,'color','w','Position', get(0, 'Screensize'));
 s=stem(f,P,'k','filled');
 grid on
 xlabel('Adimensionalised frequency St=f*D/U')
 ylabel('Power of mode || \phi ||')
 title('Dynamical Mode Decomposition Power Spectrum');
 set(gca, 'FontSize', 14)
 hold on
 s2=stem(fBest,PBest,'r','filled');
 legend('Highest Order Model','Highest FIT Model','Location','bestoutside','Orientation','horizontal')
 legend('boxoff')
 export_fig(figure470,strcat(dirdmd,'/image','dmdpowerspectrum'),'-nocrop','-m2');
 
 %% POWER ACCOUNTED FOR BY EACH MODE
 totalpower=sum(P);
 pfraction=P/totalpower;
 
 incpower=0;
 
 for i=1:length(P)
     incpower=incpower+pfraction(i);
     pp(i)=incpower;
 end
 
 fig403=figure(403);
 set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');
 a=1:1:length(f);
 sid=scatter(a,pp,'o');
 hold on
 sid.MarkerFaceColor = [0 0.8 0.2];
 sid.MarkerEdgeColor = [0 0.8 0.2];
 powerplot=plot(a,pp,'LineWidth',1.5','color','red');
 powerplot.LineStyle='- -';
 powerplot.Color=[0 0.8 0.2];
 xlabel(' Number of singular values');
 ylabel('Total Energy maesured as || \phi || (%)')
 title('Energy accounted for by number of singular values used');
 set(gca, 'FontSize', 14)
 grid on
 grid minor
 hold off
 export_fig(fig403,strcat(dirdmd,'/image','poweraccountedforwithsv'),'-nocrop','-m2');
 
printresults(imag(omega),P,LambdaDiag,damping,method,dirdmd)
