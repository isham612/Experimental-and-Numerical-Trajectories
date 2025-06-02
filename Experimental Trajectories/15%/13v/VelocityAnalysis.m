clear all;
close all;
clc;

[DataT.File,DataT.Path] = uigetfile('*.txt');
A=readtable([DataT.Path,DataT.File]);
ID=table2array(A);


cal=0.95;                                                     %pixels to micron
fps=2;                                                      %frames per second

%% Plot trajectory

figure(1)
plot(ID(:,2),ID(:,3),'Color','b','LineWidth',2.5);
xlim([0 700])
ylim([0 500])
title([DataT.File])


% We calculate the instant velocity as distance/time, where time=1 second (2 frames)


    TB = ID(:,1);
    XB = ID(:,2);
    YB = ID(:,3);        
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
    InstantVelocity=diff2(:,2)./deltaT(2);
    
    F=figure(2)
    plot([0:1:length(InstantVelocity(2:end))-1]./fps,InstantVelocity(2:end),'Color','b','LineWidth',2.5);
    xlabel('t [s]','FontName','Times','FontSize',50)
    ylabel('V [\mum/s]','FontName','Times','FontSize',50)
    set(gca,'linewidth',3,'FontName','Times','FontSize',40,'box','on','TickLength',[0.02,0.02])
    set(F,'Position',[10, 100, 1000, 800]);
    
    v_mean = nanmean(InstantVelocity); %if parts are wrong.. nanmean(InstantVelocity(A*2:B*2) (300:400)  (200*2:400*2)
    v_std = nanstd(InstantVelocity);



