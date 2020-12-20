function [modelVAF_val]=idvaloverview(FITje,FITje_val,dirdmd,name)

figure350= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
%% TURBINE 1
subplot(2,1,1)
turbine1=[FITje(1,:);FITje_val(1,:)]';
x=turbine1(:,1);
y=turbine1(:,2);
a=size(turbine1,1);
modes=1:1:a;
b1 = x\y;
yexp = b1*x;
sid=scatter(modes,x,'o');
hold on
sid.MarkerFaceColor = [1 0 0];
sid.MarkerEdgeColor = [1 0 0];
pid=plot(x,'LineWidth',1.1','color','red');
ylim([min(x)-20 100])
yticks=0:10:100;
pid.LineStyle='- -';
pid.Color=[1 0.7 0.7];
hold on
sval=scatter(modes,y,'o');
sval.MarkerFaceColor = [0 0 1];
sval.MarkerEdgeColor = [0 0 1];
hold on
pval=plot(y,'LineWidth',1.1','color','blue');
pval.LineStyle='- -';
pval.Color=[0.7 0.7 1];
%sval=plot(x,yexp);
warning off
lgd=legend([sid sval],'VAF: Identification Models','VAF: Validation Models',...
    'VAF: Validation Models','VAF: Validation Models','Location',...
    'best','Orientation','horizontal');
warning on
lgd.FontSize=14;
legend('boxoff')
xlabel('Number of modes in model');
ylb=ylabel('Variance Accounted For (VAF) [%]');
title('VAF for tubine 1 by different models for idenitfication and validation');
set(gca, 'FontSize', 14);
ax=gca;
ax.YAxis.TickValues=[0 10 20 30 40 50 60 70 80 90 100];
ax.XAxis.TickValues=0:10:200;
grid on
hold off

%% TURBINE 2
subplot(2,1,2)
turbine2=[FITje(2,:);FITje_val(2,:)]';
x=turbine2(:,1);
y=turbine2(:,2);
a=size(turbine2,1);
modes=1:1:a;
b1 = x\y;
yexp = b1*x;
sid=scatter(modes,x,'o');
hold on
sid.MarkerFaceColor = [1 0 0];
sid.MarkerEdgeColor = [1 0 0];
pid=plot(x,'LineWidth',1.1','color','red');
yticks=0:10:100;
pid.LineStyle='- -';
pid.Color=[1 0.7 0.7];
hold on
sval=scatter(modes,y,'o');
sval.MarkerFaceColor = [0 0 1];
sval.MarkerEdgeColor = [0 0 1];
hold on
pval=plot(y,'LineWidth',1.1','color','blue');
yticks=0:10:100;
pval.LineStyle='- -';
pval.Color=[0.7 0.7 1];
%sval=plot(x,yexp);
warning off
lgd2=legend([sid sval],'VAF: Identification Models','VAF: Validation Models',...
    'VAF: Validation Models','VAF: Validation Models','Location',...
    'bestoutside','Orientation','horizontal');
warning on
lgd2.FontSize=14;
legend('boxoff')
xlabel('Number of modes in model');
ylb=ylabel('Variance Accounted For (VAF) [%]');
title('VAF for tubine 2 by different models for idenitfication and validation');
set(gca, 'FontSize', 14);
ax=gca;
ax.YAxis.TickValues=[0 10 20 30 40 50 60 70 80 90 100];
ax.XAxis.TickValues=0:10:200;
grid on
hold off

%% overall results, giving a weighted average VAF for a model 
% subplot(3,1,3)
 turbine1=[FITje(1,:);FITje_val(1,:)]';
 turbine2=[FITje(2,:);FITje_val(2,:)]';
% 
 alpha=0.5; %weight for VAF of turbine 1
 beta=0.5; %wight of VAF for turbine 2
% 
 modelVAF_id=turbine1(:,1)*alpha+turbine2(:,1)*beta;
 modelVAF_val=turbine1(:,2)*alpha+turbine2(:,2)*beta;
% 
% x=modelVAF_id;
% y=modelVAF_val;
% a=size(turbine2,1);
% modes=1:1:a;
% b1 = x\y;
% yexp = b1*x;
% sid=scatter(modes,x,'o');
% hold on
% sid.MarkerFaceColor = [1 0 0];
% sid.MarkerEdgeColor = [1 0 0];
% pid=plot(x,'LineWidth',1.1','color','red');
% yticks=0:10:100;
% pid.LineStyle='- -';
% pid.Color=[1 0.7 0.7];
% hold on
% sval=scatter(modes,y,'o');
% sval.MarkerFaceColor = [0 0 1];
% sval.MarkerEdgeColor = [0 0 1];
% hold on
% pval=plot(y,'LineWidth',1.1','color','blue');
% yticks=0:10:100;
% pval.LineStyle='- -';
% pval.Color=[0.7 0.7 1];
% %sval=plot(x,yexp);
% legend('VAF: Identification Models','VAF: Identification Models',...
%     'VAF: Validation Models','VAF: Validation Models','Location',...
%      'bestoutside','Orientation','horizontal');
% legend('boxoff')
% xlabel('Number of modes in model');
% ylb=ylabel('Variance Accounted For (VAF)');
% title('Average VAF for both turbines by different models for identification and validation');
% set(gca, 'FontSize', 12);
% ax=gca;
% ax.YAxis.TickValues=[0 10 20 30 40 50 60 70 80 90 100];
% ax.XAxis.TickValues=0:10:200;
% grid on
hold off
set(gcf,'color','w','Position', get(0, 'Screensize'));
export_fig(figure350,strcat(dirdmd,'/image',name),'-nocrop','-m2');
