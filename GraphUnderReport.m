%----------------------------Introduction------------------------------%
% Visuliaztion of All 54 African Countries Fitted Parameter Values
% Input: 'meanpara' from DataFit54Country.m
% Output: Grahphs + Matrices meanpara & meanstd into Latex matrix form
%----------------------------------------------------------------------%

global color meanpara meanstd;

%--------------Input----------------%
parameters=meanpara;
%-----------------------------------%



%-----------------Regional Meta Mean and STD ---------------------%
npara=length(parameters(1,:));
mean1=parameters([6 7 9 10 14 17 21 12 41],:);% central africa
mean2=parameters([11 15 18 20 26 30 34 40 43 45 47 48 49 52],:);% eastern africa
mean3=parameters([1 16 29 33 35 51],:);% northern africa
mean4=parameters([2 4 19 27 31 36 37 46 53 54],:);% southern africa
mean5=parameters([3 5 8 13 22 23 24 25 28 32 38 39 42 44 50],:);% western africa
table1=zeros(2,npara);% central africa
table2=zeros(2,npara);% eastern africa
table3=zeros(2,npara);% northern africa
table4=zeros(2,npara);% southern africa
talbe5=zeros(2,npara);% western africa
for i=1:npara
    table1(1,i)=mean(mean1(:,i));
    table1(2,i)=std(mean1(:,i));
    table2(1,i)=mean(mean2(:,i));
    table2(2,i)=std(mean2(:,i));
    table3(1,i)=mean(mean3(:,i));
    table3(2,i)=std(mean3(:,i));
    table4(1,i)=mean(mean4(:,i));
    table4(2,i)=std(mean4(:,i));
    table5(1,i)=mean(mean5(:,i));
    table5(2,i)=std(mean5(:,i));
end

regionaltable=[table1(1,[1 2 8 9]);table2(1,[1 2 8 9]);table3(1,[1 2 8 9]);table4(1,[1 2 8 9]);table5(1,[1 2 8 9])];
regionaltable(:,[3 4])=regionaltable(:,[3 4])*4;
regionalstring={'Central','Eastern','Northern','Southern','Western'};


%colors
color=[252 41 30;250 200 205;219 249 244;54 195 201;0 70 222;255 255 199;255 170 50]./255;


%-------------------------------Figures------------------------------------%
%-----------R0-----------%
clf
close all
figure(1)
[a,b]=sort(parameters(:,1),'ascend');
bb=barh(parameters(b,[2 1]),1.7);
set(bb(1),'FaceColor',color(1,:));
set(bb(2),'FaceColor',color(4,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
grid on;
legendstr1=['R_0 for reported infections (Mean ',num2str(roundn(meanstd(1,2),-2)),', SD ',num2str(roundn(meanstd(2,2),-2)),')'];
legendstr2=['R_0 for all infections (Mean ',num2str(roundn(meanstd(1,1),-2)),', SD ',num2str(roundn(meanstd(2,1),-2)),')'];
legend({legendstr1,legendstr2},'FontSize',14,'Location','SouthOutside');
title('Basic reproduction numbers','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40]);
set(gca,'XAxisLocation','top');
saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','R0.png']);


%-----------Report Rate-----------%
figure(2)
[a,b]=sort(parameters(:,9),'ascend');
bb=barh(parameters(b,[8 9]),1.7);%r_m and overall r
set(bb(1),'FaceColor',color(2,:));
set(bb(2),'FaceColor',color(1,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
grid on;
legendstr1=['Report rate of mild infections (Mean ',num2str(roundn(meanstd(1,8),-2)),'%, SD ',num2str(roundn(meanstd(2,8),-2)),'%)'];
legendstr2=['Report rate of all infections (Mean ',num2str(roundn(meanstd(1,9),-2)),'%, SD ',num2str(roundn(meanstd(2,9),-2)),'%)'];
legend({legendstr1,legendstr2},'FontSize',14,'Location','SouthOutside');
title('Report rates (%)','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','r.png']);


%---------------Proportion of Severe, Mild and Asympt. Cases----------------------------------%
figure(3)
[a,b]=sort(parameters(:,5),'ascend');
proportion=[parameters(b,5) parameters(b,6) parameters(b,7)];% p_s*(1-p_a),(1-p_a)*(1-p_s),p_a
bb=barh(proportion,'stacked');
set(bb(1),'FaceColor',color(1,:));
set(bb(2),'FaceColor',color(2,:));
set(bb(3),'FaceColor',color(3,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
axis([0 100 0.4 54.6]);
legendstr1=['Severe infections (Mean ',num2str(roundn(meanstd(1,5),-2)),'%, SD ',num2str(roundn(meanstd(2,5),-2)),'%)'];
legendstr2=['Mild infections (Mean ',num2str(roundn(meanstd(1,6),-2)),'%, SD ',num2str(roundn(meanstd(2,6),-2)),'%)'];
legendstr3=['Asymptomatic infections (Mean ',num2str(roundn(meanstd(1,7),-2)),'%, SD ',num2str(roundn(meanstd(2,7),-2)),'%)'];
legend({legendstr1,legendstr2,legendstr3},'FontSize',14,'Location','SouthOutside');
title('Porportion of severe, mild and asymptomatic infections (%)','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','p.png']);


%---------------Relative infectiousness of severe and asymptomatic cases----------------------%
figure(4)
[a,b]=sort(parameters(:,4),'ascend');
bb=barh(parameters(b,[4 3]),1.7);% rho_a,rho_s
set(bb(1),'FaceColor',color(6,:));
set(bb(2),'FaceColor',color(7,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
axis([0 Inf 0.4 54.6]);
grid on;
legendstr1=['Relative infectiousness of asymptomtaic infection (Mean ',num2str(roundn(meanstd(1,4),-2)),', SD ',num2str(roundn(meanstd(2,4),-2)),')'];
legendstr2=['Relative infectiousness of severe infection (Mean ',num2str(roundn(meanstd(1,3),-2)),', SD ',num2str(roundn(meanstd(2,3),-2)),')'];
legend({legendstr1,legendstr2},'FontSize',14,'Location','SouthOutside');
title('Relative infectiousness of severe and asymptomatic infections','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','rho.png']);


%-----------Number of Severe, Mild, Asympt., Recoverd and Exposed at Day 0------------------%
figure(5)
[a,b]=sort(parameters(:,14),'ascend');
initial=[parameters(b,[14 13 11 12 10])];% R,A_0,I_m0,I_s0 and E
bb=barh(initial(1:end,:),'stacked');
set(bb(1),'FaceColor',color(4,:));
set(bb(2),'FaceColor',color(3,:));
set(bb(3),'FaceColor',color(2,:));
set(bb(4),'FaceColor',color(1,:));
set(bb(5),'FaceColor',color(7,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
grid on;
legendstr1=['Recovered (Mean ',num2str(roundn(meanstd(1,14),-2)),', SD ',num2str(roundn(meanstd(2,14),-2)),')'];
legendstr2=['Asymptomatic infections (Mean ',num2str(roundn(meanstd(1,13),-2)),', SD ',num2str(roundn(meanstd(2,13),-2)),')'];
legendstr3=['Mild infections (Mean ',num2str(roundn(meanstd(1,11),-2)),', SD ',num2str(roundn(meanstd(2,11),-2)),')'];
legendstr4=['Severe infections (Mean ',num2str(roundn(meanstd(1,12),-2)),', SD ',num2str(roundn(meanstd(2,12),-2)),')'];
legendstr5=['Exposed (Mean ',num2str(roundn(meanstd(1,10),-2)),', SD ',num2str(roundn(meanstd(2,10),-2)),')'];
legend({legendstr1,legendstr2,legendstr3,legendstr4,legendstr5},'FontSize',14,'Location','SouthOutside');
title('Number of initial recovered, asymptomatic, mild, severe and exposed infections','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','0.png']);


%-----------------Regional Bar Plot------------------------------%
%figure(6)
%bb=bar(regionaltable);
%set(bb(1),'FaceColor',color(4,:));
%set(bb(2),'FaceColor',color(3,:));
%set(bb(3),'FaceColor',color(2,:));
%set(bb(4),'FaceColor',color(1,:));
%ax=gca;
%ax.XTick=1:5; 
%ax.XTickLabels=regionalstring;
%yyaxis left
%ylabel('Basic reproduction number','FontSize',22)
%yyaxis right
%ylable('Report rate (%)','FontSize',22)
%grid on;
%legend({'R_0','R_K0','r_m','r_overall'},'FontSize',14,'Location','northeast')
%title('Regional Mean Estimates','FontSize',22);
%set(gcf,'unit','centimeters','position',[3 5 30 15]);

   
%----------Output all parameters mean values & meta mean std----------------%
latex2(parameters,2,AfricaCountry)
latex2(meanstd,2,{'Mean' 'SD'})


