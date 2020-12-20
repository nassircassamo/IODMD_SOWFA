function [nrmse, nrmsevalid]=evaluatetimevaryingerror(X,statesrebuild,Xvalid, statesrebuildvalid, dirdmd,filename)

timevecc=size(X,2);
timevec=1:1:timevecc;

%conversion of time instant to continuous time in simulation
for l=1:length(timevec)
    k=timevec(l)+250;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    timeveccont(l)=integ;
end

% %calculate percentual deviations
% for t=1:length(timevec)
%     deltastate(:,t)=abs(X(:,t)-real(statesrebuild(:,t)));
%     deltastatep(:,t)=abs(X(:,t)-real(statesrebuild(:,t)))./abs(X(:,t));
%     meandeltainst(:,t)=mean(deltastate(:,t));
%     meandeltainstp(:,t)=mean(deltastatep(:,t))*100;
%     %deltainst(:,t)=norm(deltastate(:,t));
% end

%% Calculate normalized root mean squared error
for t=1:length(timevec)
    diff=X(:,t)-statesrebuild(:,t);
    diffsquared=diff.^2;
    rmse=sqrt(sum(diffsquared)/length(diff));
    scale=max(X(:,t))-min(X(:,t));
    nrmse(t)=rmse/scale*100;
end

%% NRMSE for validaiton
%% Calculate normalized root mean squared error
for t=1:length(timevec)
    diff=Xvalid(:,t)-statesrebuildvalid(:,t);
    diffsquared=diff.^2;
    rmse=sqrt(sum(diffsquared)/length(diff));
    scale=max(Xvalid(:,t))-min(Xvalid(:,t));
    nrmsevalid(t)=rmse/scale*100;
end


figure550=figure;
figure550.Visible='off';
%82 is index right before frst yaw on identification data
set(gcf,'color','w','Position', get(0, 'Screensize')); 
%sid=scatter(timevec(82:end),meandeltainstp(82:end),'o');
sid=scatter(timevec(82:end),nrmse(82:end),'o');
hold on
sid.MarkerFaceColor = [0.8 0.1 0.1];
sid.MarkerEdgeColor = [0.8 0.1 0.1];
pid=plot(nrmse,'LineWidth',1.1','color','red');
pid.LineStyle='- -';
pid.Color=[1 0.8 0.8];

hold on
sidvalid=scatter(timevec(82:end),nrmsevalid(82:end),'o');
sidvalid.MarkerFaceColor = [0.1 0.1 0.9];
sidvalid.MarkerEdgeColor = [0.1 0.1 0.9];
pidvalid=plot(nrmsevalid,'LineWidth',1.1','color','green');
pidvalid.LineStyle='- -';
pidvalid.Color=[0.8 0.8 1];

xlabel('Time [minutes]');
axis([-0.3667 750 0 100])
ax=gca;
set(gca,'XTick',(0-0.3667:60:max(timevec))+1)
ax.XTickLabel = {'8','10','12','14','16','18','20','22','24','26','28','30','32','34'};
set(gca,'YTick',0:5:100)
ylabel('Normalized Root Mean Squared Error (%) ');
title('Dynamic Mode Decomposition State Reconstruction: Normalized Root Mean Squared Error through time');
set(gca, 'FontSize', 14);
grid on
hold off
ytickformat('percentage')

legend([sid sidvalid],'NRMSE for identification',...
            'NRMSE for validation ',...
            'Location','bestoutside','Orientation','horizontal')
        legend('boxoff')
   
export_fig(figure550,strcat(dirdmd,'/image',filename),'-nocrop','-m2'); 
