%----------------------------Introduction------------------------------%
% Visuliaztion of All 54 African Countries Fitted Parameter Values
% Input: 'meanpara' from DataFitUnderReport54Country.m or 
%        'meanparaE' from DataFitUnderReport54CountryWithE.m
% Output matrix meanpara into Latex matrix form
%----------------------------------------------------------------------%

global color;

%--------------Input----------------%
input='meanparaE';
%-----------------------------------%

parameters=eval(input);

%-----------------Mean and SD of all Parameters Mean Values---------------------%
npara=length(parameters(1,:));
mean1=parameters([6 7 9 10 14 17 21 12 41],:);
mean2=parameters([11 15 18 20 26 30 34 40 43 45 47 48 49 52],:);
mean3=parameters([1 16 29 33 35 51],:);
mean4=parameters([2 4 19 27 31 36 37 46 53 54],:);
mean5=parameters([3 5 8 13 22 23 24 25 28 32 38 39 42 44 50],:);
table=zeros(2,npara);
table1=zeros(2,npara);
table2=zeros(2,npara);
table3=zeros(2,npara);
table4=zeros(2,npara);
talbe5=zeros(2,npara);
for i=1:npara
    table(1,i)=mean(parameters(:,i));
    table(2,i)=std(parameters(:,i));
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

regionaltable=[table1(1,[1 2 7 8]);table2(1,[1 2 7 8]);table3(1,[1 2 7 8]);table4(1,[1 2 7 8]);table5(1,[1 2 7 8])];
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
legend({'R_0 for reported cases (Mean 0.92 SD 0.65)','R_0 (Mean 2.20 SD 0.72)'},'FontSize',14,'Location','SouthOutside')
title('Basic reproduction number','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40]);
set(gca,'XAxisLocation','top');
if strcmp(input,'meanpara')
   saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Latex/','R0.png']);
end
if strcmp(input,'meanparaE')
   saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','R0WithE.png']);
end


%-----------Report Rate-----------%
figure(2)
[a,b]=sort(parameters(:,8),'ascend');
bb=barh(parameters(b,[7 8])*100,1.7);%r_m and overall r
set(bb(1),'FaceColor',color(2,:));
set(bb(2),'FaceColor',color(1,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
grid on;
legend({'Mild case report rate (Mean 49% SD 23%)','Overall report rate (Mean 35% SD 22%)'},'FontSize',14,'Location','SouthOutside')
title('Report rate (%)','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
if strcmp(input,'meanpara')
   saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Latex/','r.png']);
end
if strcmp(input,'meanparaE')
   saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','rWithE.png']);
end


%---------------Proportion of Severe, Mild and Asympt. Cases----------------------------------%
figure(3)
[a,b]=sort(parameters(:,5),'ascend');
proportion=[parameters(b,5) ones(54,1)-parameters(b,5)-parameters(b,6) parameters(b,6)];% p_s,p_m,p_a
bb=barh(proportion*100,'stacked');
set(bb(1),'FaceColor',color(1,:));
set(bb(2),'FaceColor',color(2,:));
set(bb(3),'FaceColor',color(3,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
axis([0 100 0.4 54.6]);
legend({'Severe cases (Mean 19% SD 14%)','Mild cases','Asympt. cases (Mean 50% SD 23%)'},'FontSize',14,'Location','SouthOutside')
title('Porportion of severe, mild and asymptomatic cases (%)','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
if strcmp(input,'meanpara')
   saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Latex/','p.png']);
end
if strcmp(input,'meanparaE')
   saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','pWithE.png']);
end

%---------------Relative infectiousness of severe and asymptomatic cases----------------------%
figure(4)
[a,b]=sort(parameters(:,4),'ascend');
bb=barh(parameters(b,[4 3]),1.7);% rho_a,rho_s
set(bb(1),'FaceColor',color(6,:));
set(bb(2),'FaceColor',color(7,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
axis([0 1.5 0.4 54.6]);
grid on;
legend({'Relative infectiousness of asympt. cases (Mean 0.56 SD 0.11)','Relative infectiousness of severe cases (Mean 1.25 SD 0.05)'},'FontSize',14,'Location','SouthOutside')
title('Relative infectiousness of severe and asymptomatic cases','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
if strcmp(input,'meanpara')
   saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Latex/','p.png']);
end
if strcmp(input,'meanparaE')
   saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','rhoWithE.png']);
end

%-----------Number of Severe, Mild, Asympt. and Recoverd Cases (Exposed) at Day 0------------------%
figure(5)

if strcmp(input,'meanpara')
[a,b]=sort(parameters(:,11),'ascend');
proportion=[parameters(b,[11 10 8 9])];% R,A_0,I_m0,I_s0
bb=barh(proportion(1:end,:),'stacked');
set(bb(1),'FaceColor',color(4,:));
set(bb(2),'FaceColor',color(3,:));
set(bb(3),'FaceColor',color(2,:));
set(bb(4),'FaceColor',color(1,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
grid on;
legend({'Recovered','Asympt. cases','Mild cases','Severe cases'},'FontSize',14,'Location','southeast')
title('Number of initial recovered, asymptomatic, mild and severe cases','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top','bottom');
%saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Latex/','0.png']);
end

if strcmp(input,'meanparaE')
[a,b]=sort(parameters(:,13),'ascend');
initial=[parameters(b,[13 12 10 11 9])];% R,A_0,I_m0,I_s0 and E
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
legend({'Recovered','Asympt. cases','Mild cases','Severe cases','Exposed'},'FontSize',14,'Location','southeast')
title('Number of initial recovered, asymptomatic, mild, severe cases and exposed','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','0WithE.png']);
end

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

   
%-----------------Output all parameters mean values------------------------------%
latex2(parameters,2,AfricaCountry)

latex2noname(table,2)
latex2noname(table1,2) 
latex2noname(table2,2)
latex2noname(table3,2)
latex2noname(table4,2)
latex2noname(table5,2)



%------------------------Matlab Matrix to Latex-----------------------------------%
function f = latex2(A, precision, name)
% if no precision is given, the default precision is 4
if nargin == 1
    precision = '4';
else
    precision = int2str(precision);
end
% 定义单一元素输出格式
out_num = [' %0.' precision 'f &'];
% 用作整数输出判断
z = zeros(1, str2num(precision) + 1);
z(1) = '.';
z(2 : end) = '0';
z = char(z);
% size of matrix A
[r c] = size(A);
nc = zeros(1, c);
nc(:)= 108;  % character l
% first Latex sentence
out = sprintf('\\left(\n\t\\begin{array}{%s}', char(nc));
% 二重循环，用作生成整个矩阵的Latex语句
for i = 1 : r
    out = [out sprintf('\n\t') name{i} '&']; % 换行
    for j = 1 : c
        temp = sprintf(out_num, A(i, j));
        % 小数位皆为零时，把数取整。如1.0001取为1
        dot_position = find(temp == '.');
        if temp(dot_position : end - 2) == z
            temp = temp(1 : dot_position - 1);
            temp = [temp ' &'];
            % 要取整时，如有负号，则必须丢掉
            if temp(2) == '-'
               temp = [temp(1) temp(3 : end)]; 
            end
        end
        out = [out temp];
    end
    % discard '&'
    out = out(1 : end - 1);
    % '\\' in the end of each row
    out = [out '\\'];
end
% last Latex sentence
out = [out sprintf('\n\t\\end{array}\n\\right)')];
f = out;
end

%------------------------Matlab Matrix to Latex-----------------------------------%
function f = latex2noname(A, precision)
% if no precision is given, the default precision is 4
if nargin == 1
    precision = '4';
else
    precision = int2str(precision);
end
% 定义单一元素输出格式
out_num = [' %0.' precision 'f &'];
% 用作整数输出判断
z = zeros(1, str2num(precision) + 1);
z(1) = '.';
z(2 : end) = '0';
z = char(z);
% size of matrix A
[r c] = size(A);
nc = zeros(1, c);
nc(:)= 108;  % character l
% first Latex sentence
out = sprintf('\\left(\n\t\\begin{array}{%s}', char(nc));
% 二重循环，用作生成整个矩阵的Latex语句
for i = 1 : r
    out = [out sprintf('\n\t')  '&']; % 换行
    for j = 1 : c
        temp = sprintf(out_num, A(i, j));
        % 小数位皆为零时，把数取整。如1.0001取为1
        dot_position = find(temp == '.');
        if temp(dot_position : end - 2) == z
            temp = temp(1 : dot_position - 1);
            temp = [temp ' &'];
            % 要取整时，如有负号，则必须丢掉
            if temp(2) == '-'
               temp = [temp(1) temp(3 : end)]; 
            end
        end
        out = [out temp];
    end
    % discard '&'
    out = out(1 : end - 1);
    % '\\' in the end of each row
    out = [out '\\'];
end
% last Latex sentence
out = [out sprintf('\n\t\\end{array}\n\\right)')];
f = out;
end