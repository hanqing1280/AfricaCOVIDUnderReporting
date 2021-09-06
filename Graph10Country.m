%----------------------------Introduction------------------------------%
% Visuliaztion of Selected African Countries Fitted Parameter Values
% Input: 'meanpara' from DataFitUnderReport54Country.m or 
%        'meanparaE' from DataFitUnderReport54CountryWithE.m
% Output matrix meanpara into Latex matrix form
%----------------------------------------------------------------------%

global color;

parameters=meanparaE([4 7 19 36 37 39 40 46 53 54],:);
SelectedCountries={'Botswana','Cameroon','Eswatini','Mozambique','Namibia','Nigeria','Rwanda','South Africa','Zambia','Zimbabwe'};

%-----------------Mean and SD of all Parameters Mean Values---------------------%
npara=length(parameters(1,:));
table=zeros(2,npara);
for i=1:npara
    table(1,i)=mean(parameters(:,i));
    table(2,i)=std(parameters(:,i));
end


%colors
color=[252 41 30;250 200 205;219 249 244;54 195 201;0 70 222;255 255 199;255 170 50]./255;

%-------------------------------Figures------------------------------------%
%-----------R0-----------%
clf
close all
figure(1)
[a,b]=sort(parameters(:,1),'ascend');
bb=barh(parameters(b,[2 1]),1.5);
set(bb(1),'FaceColor',color(1,:));
set(bb(2),'FaceColor',color(4,:));
ax=gca;
ax.YTick=1:10; 
ax.YTickLabels=SelectedCountries(b);
grid on;
legend({'R_0 for reported cases (Mean 0.89 SD 0.73)','R_0 (Mean 2.24 SD 0.76)'},'FontSize',18,'Location','SouthOutside')
set(gca,'FontSize',20);
title('Basic reproduction number','FontSize',24);
set(gcf,'unit','centimeters','position',[3 5 30 40]);
set(gca,'XAxisLocation','top'); 
saveas(gcf,['/Users/Qing/Desktop/','R0.png']);


%-----------Report Rate-----------%
figure(2)
[a,b]=sort(parameters(:,8),'ascend');
bb=barh(parameters(b,[7 8])*100,1.5);%r_m and overall r
set(bb(1),'FaceColor',color(2,:));
set(bb(2),'FaceColor',color(1,:));
ax=gca;
ax.YTick=1:10; 
ax.YTickLabels=SelectedCountries(b);
grid on;
legend({'Mild case report rate (Mean 48% SD 25%)','Overall report rate (Mean 33% SD 23%)'},'FontSize',18,'Location','SouthOutside')
set(gca,'FontSize',20);
title('Report rate (%)','FontSize',24);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
saveas(gcf,['/Users/Qing/Desktop/','r.png']);



%---------------Proportion of Severe, Mild and Asympt. Cases----------------------------------%
figure(3)
[a,b]=sort(parameters(:,5),'ascend');
proportion=[parameters(b,5) ones(10,1)-parameters(b,5)-parameters(b,6) parameters(b,6)];% p_s,p_m,p_a
bb=barh(proportion*100,'stacked');
set(bb(1),'FaceColor',color(1,:));
set(bb(2),'FaceColor',color(2,:));
set(bb(3),'FaceColor',color(3,:));
ax=gca;
ax.YTick=1:10; 
ax.YTickLabels=SelectedCountries(b);
axis([0 100 0.4 10.6]);
legend({'Severe cases (Mean 19% SD 18%)','Mild cases','Asympt. cases (Mean 48% SD 24%)'},'FontSize',18,'Location','SouthOutside')
set(gca,'FontSize',20);
title('Porportion of severe, mild and asymptomatic cases (%)','FontSize',24);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top'); 
saveas(gcf,['/Users/Qing/Desktop/','p.png']);

%---------------Relative infectiousness of severe and asymptomatic cases----------------------%
figure(4)
[a,b]=sort(parameters(:,4),'ascend');
bb=barh(parameters(b,[4 3]),1.5);% rho_a,rho_s
set(bb(1),'FaceColor',color(6,:));
set(bb(2),'FaceColor',color(7,:));
ax=gca;
ax.YTick=1:10; 
ax.YTickLabels=SelectedCountries(b);
axis([0 1.4 0.4 10.6]);
grid on;
legend({'Asympt. cases (Mean 0.53 SD 0.13)','Severe cases (Mean 1.26 SD 0.03)'},'FontSize',18,'Location','SouthOutside')
set(gca,'FontSize',20);
title('Relative infectiousness of severe and asymptomatic cases','FontSize',24);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top'); 
saveas(gcf,['/Users/Qing/Desktop/','rho.png']);

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
ax.YTick=1:10; 
ax.YTickLabels=SelectedCountries(b);
grid on;
legend({'Recovered','Asympt. cases','Mild cases','Severe cases'},'FontSize',14,'Location','SouthOutside')
title('Number of initial recovered, asymptomatic, mild and severe cases','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top','bottom');
%saveas(gcf,['/Users/Qing/Desktop/','0.png']);
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
ax.YTick=1:10; 
ax.YTickLabels=SelectedCountries(b);
grid on;
legend({'Recovered','Asympt. cases','Mild cases','Severe cases','Exposed'},'FontSize',14,'Location','SouthOutside')
title('Number of initial recovered, asymptomatic, mild, severe cases and exposed','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
%saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','0WithE.png']);
end



   
%-----------------Output all parameters mean values------------------------------%
latex2(parameters,2,SelectedCountries)

latex2noname(table,2)





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