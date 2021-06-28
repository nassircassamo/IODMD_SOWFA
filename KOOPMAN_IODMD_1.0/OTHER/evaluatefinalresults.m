function evaluatefinalresults(Preffinal,dirName)

for j=1:1:length(dirName)
        cd(dirName{j})  
        cd ..

        [nTurbine,time1,dt,nVal,powerRotor] = readTurbineOutputGlobal(dirName{j},'rotorPower');
        [time2,pitch] = readPitchData(strcat(dirName{j},'/1000/bladePitch')); %pitch=pitch{:};
        [nTurbine,time3,dt,nVal,powerGenerator] = readTurbineOutputGlobal(dirName{j},'generatorPower');
        [nTurbine,time4,dt,nVal,thrust] = readTurbineOutputGlobal(dirName{j},'rotorAxialForce');
        [nTurbine,time4,dt,nVal,thrustv] = readTurbineOutputGlobal(dirName{j},'rotorVerticalForce');
        [nTurbine,time4,dt,nVal,thrusth] = readTurbineOutputGlobal(dirName{j},'rotorHorizontalForce');
        [nTurbine,time5,dt,nVal,torqueGen] = readTurbineOutputGlobal(dirName{j},'generatorTorque');
        [nTurbine,time6,dt,nVal,torqueRotor] = readTurbineOutputGlobal(dirName{j},'generatorTorque');
        [nTurbine,time7,dt,nVal,rotSpeed] = readTurbineOutputGlobal(dirName{j},'rotorSpeed');

        %rotor speed F ?
        [nTurbine,time7,dt,nVal,rotSpeedF] = readTurbineOutputGlobal(dirName{j},'rotorSpeedFiltered');
        [nTurbine,time8,dt,nVal,nacYaw] = readTurbineOutputGlobal(dirName{j},'nacelleYaw');
        [nTurbine,time8,dt,nVal,rotorAzimuth] = readTurbineOutputGlobal(dirName{j},'rotorAzimuth');
        %[nTurbine2,bla,time9,dt9,nVal9,bladeV] = readTurbineOutputSectional(dirName{j},'bladePointVmag');

        %yaw data
        [nTurbine,timeyaw,dt,nVal,nacelleYaw] = readTurbineOutputGlobal(dirName{j}, 'nacelleYaw');
end
    
        %% YAW AND DELTA YAW
        yaw=nacelleYaw(:,1);
        [yaw] = resampleedgeeffect(yaw,10);
        yaw=yaw';
        
        fig1=figure;
        fig1.Visible='off';
        set(gcf,'color','w','Position', get(0, 'Screensize'));
        subplot(2,1,1);
        p=plot((timeyaw(1:10:end)-1000)./60,yaw(:,1),'LineWidth',1.5);
        p.Color=[0.8 0.1 0.1];
        ax=gca;
        set(gca,'XTick',0:1:34)
        grid on
        hold on
        p=plot([50/60 50/60],[min(yaw(:,1))-5 max(yaw(:,1)) ]);
        p.Color=[0 0 0];
        p.LineWidth=1.5;
        xaxis([0 max((timeyaw(1:10:end)-1000)./60)])
        yaxis([min(yaw(:,1))-5 max(yaw(:,1)) ])
        title('Open loop optimal contral action','FontName','Helvetica')
        xlabel('Time [ min ]','FontName','Helvetica')
        ylabel('Yaw angle [deg]','FontName','Helvetica')
        set(gca,'fontsize', 14)
        text(50/60,max(yaw)-2.00,'\leftarrow Controller implementation starts','FontSize',14)
        
        subplot(2,1,2);
        deltau=diff(yaw(:,1));
        p=plot((timeyaw(2:10:end)-1000)./60,deltau,'LineWidth',1.5);
        p.Color=[0.8 0.1 0.1];
        ax=gca;
        set(gca,'XTick',0:1:34)
        grid on
        hold on
        p=plot([50/60 50/60],[-0.3 0.3 ]);
        p.Color=[0 0 0];
        p.LineWidth=1.5;
        xaxis([0 max((timeyaw(1:10:end)-1000)./60)])
        yaxis([-0.3 0.3 ])
        ylabel('\Delta \gamma [º s^{-1}]');
        xlabel('Time [ min ]','FontName','Helvetica')
        title('Optimal control action  rate for generator power tracking');
        grid on
        set(gca,'fontsize', 14)
        text(50/60,max(deltau)-0.05,'\leftarrow Controller implementation starts','FontSize',14)
        shg
        %% POWER
        fig2=figure;
        fig2.Visible='off';
        set(gcf,'color','w','Position', get(0, 'Screensize'));
        subplot(2,1,1);
        p=plot((timeyaw(1:10:end)-1000)./60,yaw(:,1),'LineWidth',1.5);
        p.Color=[0.8 0.1 0.1];
        ax=gca;
        set(gca,'XTick',0:1:34)
        grid on
        hold on
        p=plot([50/60 50/60],[min(yaw(:,1))-5 max(yaw(:,1)) ]);
        p.Color=[0 0 0];
        p.LineWidth=1.5;
        xaxis([0 max((timeyaw(1:10:end)-1000)./60)])
        yaxis([min(yaw(:,1))-5 max(yaw(:,1)) ])
        title('Open loop optimal contral action','FontName','Helvetica')
        xlabel('Time [ min ]','FontName','Helvetica')
        ylabel('Yaw angle [deg]','FontName','Helvetica')
        set(gca,'fontsize', 14)
        text(50/60,max(yaw)-2.00,'\leftarrow Controller implementation starts','FontSize',14)
        
        subplot(2,1,2);
        set(gcf,'color','w','Position', get(0, 'Screensize'));
        a=powerGenerator(1:end,1)./1e6/rho;
        b=powerGenerator(1:end,2)./1e6/rho;
        c=a+b; 
        p=plot((timeyaw-1000)./60, c(1:end),'LineWidth',1.5);
        p.Color=[0.1 0.1 0.9];
        hold on
        p=plot([50/60 50/60],[6 8]);
        p.Color=[0 0 0];
        p.LineWidth=1.5;
        grid on
        ax=gca;
        set(gca,'XTick',0:1:34)
        xaxis([0 max((timeyaw(1:10:end)-1000)./60)])
        title('Wind farm total power output reference tracking','FontName', 'Helvetica')
        xlabel('Time [ min ]','FontName', 'Helvetica')
        ylabel('Power [MW]','FontName', 'Helvetica')
        set(gca,'fontsize', 14)
        yaxis([6 8])
        text(50/60,max(c)-0.10,'\leftarrow Controller implementation starts','FontSize',14)
        p=plot((timeyaw(1:10:end)-1000)./60,Pref(1:1001)+meanY1*10^-6+meanY2*10^-6);
        shg
        
        %% PREDICTIONS AND REALITY TURBINE PWOER
        %simulations with yaw angle
        xo=dinit(model.A,model.B,model.C,model.D,Inputs_test, Outputs_test);
        [ysim, t, xout]=lsim(model, (yaw-260),[],xo);  ysim=ysim';
           
        fig3=figure;
        fig3.Visible='off';
        set(gcf,'color','w','Position', get(0, 'Screensize'));
        
        subplot(3,1,1);
        p=plot((timeyaw(1:10:end)-1000)./60,ysim(1,:)+5.352);
        p.Color=[0.1 0.1 0.9];
        p.LineWidth=1.5;
        p.LineStyle='--';
        a=powerGenerator(1:end,1)./1e6/rho;
        [a] = resampleedgeeffect(a,10);
        hold on
        p=plot((timeyaw(1:10:end)-1000)./60,a);
        p.Color=[0.1 0.1 0.9];
        p.LineWidth=1.5;
        p=plot([50/60 50/60],[4.4 5.5 ]);
        p.Color=[0 0 0];
        p.LineWidth=1.5;
        grid on
        title('Upstream turbine generator power','FontName', 'Helvetica')
        xlabel('Time [ min ]','FontName', 'Helvetica')
        ylabel('Power [MW]','FontName', 'Helvetica')
        set(gca,'fontsize', 14)
        text(50/60,4.6,'\leftarrow Controller implementation starts','FontSize',14)
        xaxis([0 max((timeyaw(1:10:end)-1000)./60)])
        ax=gca;
        set(gca,'XTick',0:1:34)
        yaxis([4.4 5.5])
        legend('Model prediciton','SOWFA','Location','best','Orientation','vertical')
        legend('boxoff')   
        
        subplot(3,1,2);
        p=plot((timeyaw(1:10:end)-1000)./60,ysim(2,:)+0.9512);
        p.Color=[0.1 0.1 0.9];
        p.LineWidth=1.5;
        p.LineStyle='--';
        b=powerGenerator(1:end,2)./1e6/rho;
        [b] = resampleedgeeffect(b,10);
        hold on
        p=plot((timeyaw(1:10:end)-1000)./60,b);
        p.Color=[0.1 0.1 0.9];
        p.LineWidth=1.5;
        p=plot([50/60 50/60],[0 10 ]);
        p.Color=[0 0 0];
        p.LineWidth=1.5;
        grid on
        title('Downstream turbine generator power','FontName', 'Helvetica')
        xlabel('Time [ min ]','FontName', 'Helvetica')
        ylabel('Power [MW]','FontName', 'Helvetica')
        set(gca,'fontsize', 14)
        text(50/60,2,'\leftarrow Controller implementation starts','FontSize',14)
        xaxis([0 max((timeyaw(1:10:end)-1000)./60)])
        ax=gca;
        set(gca,'XTick',0:1:34)
        yaxis([0 3])
        legend('Model prediciton','SOWFA','Location','best','Orientation','vertical')
        legend('boxoff')   
        
        subplot(3,1,3);
        p=plot((timeyaw(1:10:end)-1000)./60,ysim(1,:)+ysim(2,:)+0.9512+5.352);
        p.Color=[0.1 0.8 0.1];
        p.LineStyle='--';
        p.LineWidth=1.5;
        c=powerGenerator(1:end,1)./1e6/rho+powerGenerator(1:end,2)./1e6/rho;
        [c] = resampleedgeeffect(c,10);
        hold on
        p=plot((timeyaw(1:10:end)-1000)./60,c);
        p.Color=[0.1 0.8 0.1];
        p.LineWidth=1.5;
        p=plot([50/60 50/60],[0 10]);
        p.Color=[0 0 0];
        p.LineWidth=1.5;
        grid on
        title('Wind farm total generator power','FontName', 'Helvetica')
        xlabel('Time [ min ]','FontName', 'Helvetica')
        ylabel('Power [MW]','FontName', 'Helvetica')
        set(gca,'fontsize', 14)
        text(50/60,7,'\leftarrow Controller implementation starts','FontSize',14)
        xaxis([0 max((timeyaw(1:10:end)-1000)./60)])
        ax=gca;
        set(gca,'XTick',0:1:34)
        yaxis([6 7.5])
        legend('Model prediciton','SOWFA','Location','best','Orientation','vertical')
        legend('boxoff')   
        shg
