function [fig200]=VAFpermodes(FITje,r,Xd)

a=size(Xd,1);
lastmode=r+a;

vectormodes=1:1:lastmode;
h = zeros(1,2);

fig200=figure(200);
fig200.Visible='off';
set(gcf,'color','w','Position', get(0, 'Screensize'));
h(1)=plot(FITje(1,:),'LineWidth',1.6,'color','blue'); %1
hold on; 
s=scatter(vectormodes, FITje(1,:),'o');
h(2)=plot(FITje(2,:),'LineWidth',1.6','color','red'); %3
s.MarkerFaceColor = [0 0 1];
s.MarkerEdgeColor = [0 0 1];
s2=scatter(vectormodes, FITje(2,:),'o');
s2.MarkerFaceColor = [1 0 0];
s2.MarkerEdgeColor = [1 0 0];
grid on
xlabel('Number of modes used','FontSize',18,'FontName','Times')
ylabel(' VAF [%]','FontSize',18,'FontName','Times')
title('Variance Accounted For (VAF) with number of modes for both turbines','FontSize',22,'FontName','Times','Fontweight','Normal')
legend(h,{'Turbine 1','Turbine 2'},'Location','southeast','FontSize',18,'FontName','Times')
get(gca,'fontname')  % shows you what you are using.
set(gca,'fontname','times')  % Set it to times