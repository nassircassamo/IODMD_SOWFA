function [FITje,OMEGA,DAMPING,fig1,x]=evaluatemodel(sys_red,si,Inputs, Outputs,FITje,OMEGA,DAMPING,purpose,x,states,U,Deterministic,method)
dN=2;
%Estimate the initial state, given the estimated system matrices, and a  set of input/output data.

a=strcmp(purpose, 'identification');

if a
    if method ==2
       % [xo]=dinit(sys_red{si}.A,sys_red{si}.B,sys_red{si}.C,sys_red{si}.D,[Inputs]',[Outputs]');
         xo=U(1:size(states,1),1:si)'*states(:,1);
    elseif method==3
        xo=U(1:size(states,1),1:si)'*states(:,1);
        xo=[Deterministic(:,1);xo];
    end
    
    
    %% VARIANCE ACCOUNTED FOR IN THE MODEL FOR IDENTIFICATION DATA
    %time response of dynamic system based on (1) derived state space
    %model (2) system inputs (3) computed initial conditions 
    [ysim, t, xout]=lsim(sys_red{si}, [Inputs]',[],xo);  
    x{si}=xout;

    %compute the variance accounted for in the model based on simualton
    %of model and real output (for each time instant
    %FITje(:,si)=vaf(ysim,[Outputs]');
    FITje(:,si)=vaf([Outputs]',ysim);
  

    %% GRAPHICAL VISUALISAITON OF MODEL PREDICTION AND TRUE SIMULATION (IDENTIFICATION) RESULTS
    fig1=figure(1000+si);
    set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');
    subplot(2,1,1)
    plot([Outputs(1,1:end-1);]','b','LineWidth',1.6'); %1
    hold on; 
    plot(ysim(:,1),'Color',[0, 0.5, 0.0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    xlabel('Time instant [ ]')
    ylabel(' Power [MW]')
    title(['Model fitness: generator power for first turbine. VAF of ',num2str(round(FITje(1,si),2)),' %. '])
    legend({'Real Generator Power','Model Output'},'Location','bestoutside','Orientation','horizontal')
    legend('boxoff')
    set(gca,'fontsize', 14)

    subplot(2,1,2)
    plot([Outputs(2,1:end-1)]','b','LineWidth',1.6'); %1
    hold on; 
    plot(ysim(:,2),'Color',[0, 0.5, 0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    xlabel('Time instant [ ]')
    ylabel(' Power [MW]')
    title(['Model fitness: generator power for second turbine. VAF of ',num2str(round(FITje(2,si),2)),' %. '])
    legend({'Real Generator Power','Model Output'},'Location','bestoutside','Orientation','horizontal')
    legend('boxoff')
    set(gca,'fontsize', 14)

    %% GRAPHICAL VISUALISAITON OF MODEL PREDICTION AND TRUE SIMULATION (VALIDATION) RESULTS
else
    
     %[xo_val]=dinit(sys_red{si}.A,sys_red{si}.B,sys_red{si}.C,sys_red{si}.D,[Inputs]',[Outputs]');
     xo_val=U(1:size(states,1),1:si)'*states(:,1);
     
     if method==3
         xo_val=[Deterministic(:,1);xo_val];
     end
     
     [ysim_val, t, xout]=lsim(sys_red{si}, [Inputs]',[],xo_val);  

     %FITje(:,si)=vaf(ysim_val,[Outputs]');  
     FITje(:,si)=vaf([Outputs]',ysim_val);

     x{si}=xout;
    
    fig1=figure(2000+si);
    fig1.Visible='off';
    set(gcf,'color','w','Position', get(0, 'Screensize'));
    subplot(2,1,1)
    plot([Outputs(1,1:end-1);]','b','LineWidth',1.6'); %1
    hold on; 
    plot(ysim_val(:,1),'Color',[0, 0.5, 0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    xlabel('Time sample k [-]')
    ylabel(' \delta Power [MW]')
    title(['Model fitness: generator power for first turbine. VAF of ',num2str(round(FITje(1,si),2)),' %.'])
    legend({'Real Generator Power','Model Output'},'Location','bestoutside','Orientation','horizontal')
    legend('boxoff')
    set(gca,'fontsize', 14)
             
    subplot(2,1,2,'Visible','off')
    plot([Outputs(2,1:end-1)]','LineWidth',1.6'); %1
    hold on; 
    plot(ysim_val(:,2),'Color',[0, 0.5, 0],'LineWidth',1.6,'LineStyle','--') %3
    grid on
    xlabel('Time instant [ ]')
    ylabel(' Power [MW]')
    title(['Model fitness: generator power for second turbine. VAF of ',num2str(round(FITje(2,si),2)),' %. '])
    legend({'Real Generator Power','Model Output'},'Location','bestoutside','Orientation','horizontal')
    legend('boxoff')
    set(gca,'fontsize', 14)

end

