clear all;
close all;
clc;

Voltage=7;          %CHANGE
NumberOfFiles=10;   %Number of Trajectories (CHANGE)
cal=0.95;           %pixels to micron
fps=2;              %Frames per second
Dv=4;               %Number of frames used to calculate the magnitude of the instantaneous velocity (V=dr/Dv) 
D=2;                %Number of frames used to calculate the orientation of the instantaneous velocity
LimCorrFit1=fps*5;
LimCorrFit2=fps*25;
Tv=0.5;               %Velocity threshold (if wanted)

for i=1:1:NumberOfFiles
   A = importdata([num2str(Voltage) 'v\' num2str(Voltage) 'v' num2str(i) '.txt']);   %CHANGE!
   P(i).ID=(A.data);
   P(i).FileName=[num2str(Voltage) 'v\' num2str(Voltage) 'v' num2str(i) '.txt'];     %CHANGE!
end

% 
% % Plot trajectory
% 
% for i=1:1:length(P)
% F=figure(i);
% plot(P(i).ID(:,2),P(i).ID(:,3),'Color','b','LineWidth',2.5);
% xlim([0 700])
% ylim([0 500])
% title(P(i).FileName)
% saveas(F,[P(i).FileName(1:end-4) '_Traj.png'])
% close(F)
% end

%% Instant velocity as distance/time, where time=Dv/2;

V=zeros(1,length(P));

for i=1:1:length(P)

    TB = P(i).ID(:,1);
    XB = P(i).ID(:,2);
    YB = P(i).ID(:,3);        
    tic
    diff2 = NaN(max(TB)-1,max(TB)-1);                          %define a "place-holder"
    for k = 2:length(TB)                                       %"from"
        for l = 1:k-1                                          %"to"
            T1 = TB(k)-TB(1);                                  %"from"
            T2 = TB(l)-TB(1);                                  %"to"
            disp1 = sqrt((XB(k)-XB(l))^2+(YB(k)-YB(l))^2);
            diff2(T1,T1-T2) = disp1;
        end
    end
    diff2 = diff2*cal;
    deltaT(1:length(TB)-1)=[1:1:length(TB)-1]./fps;
    InstantVelocity=diff2(:,D)./deltaT(D);

val=InstantVelocity>Tv;
v_mean = nanmean(InstantVelocity(val)); %you can select a velocity threshold
v_std = nanstd(InstantVelocity);
  
G=figure(i);
plot([0:1:length(InstantVelocity(2:end))-1]./fps,InstantVelocity(2:end),'Color','b','LineWidth',2.5);
hold on
plot([0 (length(InstantVelocity(2:end))-1)/fps],[v_mean v_mean],'--','Color','r','LineWidth',2.5)
xlabel('t [s]','FontName','Times','FontSize',50)
ylabel('V [\mum/s]','FontName','Times','FontSize',50)
set(gca,'linewidth',3,'FontName','Times','FontSize',40,'box','on','TickLength',[0.02,0.02])
set(G,'Position',[10, 100, 1000, 800]);
title(['dt =' num2str(D/fps) 's, <V> =' num2str(v_mean)])
%saveas(G,[P(i).FileName(1:end-4) '_Velocity.png'])
close(G)

V(i)=v_mean;

end



%% Velocity correlation

diff=NaN(1000,20);

for i=1:1:length(P)
TB = P(i).ID(:,1);
XB = P(i).ID(:,2);
YB = P(i).ID(:,3);        
tic
diff2 = NaN(max(TB)-1,max(TB)-1);                        
  for k = D+1:length(TB)                                     
        for l = D+1:k                                        
            T1 = TB(k)-TB(1);                                 
            T2 = TB(l)-TB(1);                                 
            disp1 = dot([XB(k)-XB(k-D) YB(k)-YB(k-D)],[XB(l)-XB(l-D) YB(l)-YB(l-D)]);
            diff2(T1,T1-T2+1) = disp1;
        end
 end
diff(1:length(nanmean(diff2)),i)=nanmean(diff2);
Deltat=[1:1:length(diff(:,i))]./fps;

H=figure(i)
plot(Deltat,diff(:,i),'Color','b','LineWidth',2.5);
xlabel('\Deltat [s]','FontName','Times','FontSize',50)
ylabel('<V(t+\Deltat)\cdotV(t)>','FontName','Times','FontSize',50)
xlim([0 50])
ylim([0 max(diff(:,i))])
set(gca,'linewidth',3,'FontName','Times','FontSize',40,'box','on','TickLength',[0.02,0.02])
set(H,'Position',[10, 100, 1000, 800]);
title(['dt =' num2str(D/fps) 's'])
%saveas(H,[P(i).FileName(1:end-4) '_Vcorr.png'])
close(H)

end

%% Mean Velocity Correlation

diff3=diff';
corr=nanmean(diff3);

I=figure;
plot(Deltat,corr,'Color','b','LineWidth',2.5);
xlabel('\Deltat [s]','FontName','Times','FontSize',50)
ylabel('<V(t+\Deltat)\cdotV(t)>','FontName','Times','FontSize',50)
xlim([0 50])
ylim([0 max(corr)])
set(gca,'linewidth',3,'FontName','Times','FontSize',40,'box','on','TickLength',[0.02,0.02])
set(I,'Position',[10, 100, 1000, 800]);
hold on
%Fitting
x=Deltat(LimCorrFit1:LimCorrFit2);
y=corr(LimCorrFit1:LimCorrFit2);
fitfun = fittype( @(a,b,x) a-b*x );
x0 = [1 1];
[fitted_curve,gof] = fit(x',y',fitfun,'StartPoint',x0);
coeffvals = coeffvalues(fitted_curve);
tau=mean(V)^2/coeffvals(2);
plot(x,fitted_curve(x),'--','Color','r','LineWidth',2.5);
title(['\tau=' num2str(tau) 's'])

% saveas(I,['Vcorr_ALL_' num2str(Voltage) 'v.png'])
% save(['Velocities_' num2str(Voltage) 'v.mat'],'V');
