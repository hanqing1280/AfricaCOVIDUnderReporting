%----------------------------Introduction---------------------------------%
% All 54 African Countries Data Fit with EXPOSED STAGE
% Method: MCMC
% Data Set: Reported Cumulative Cases (C_K) + Reported New Cases (I_K)
% Input: mindays (minimum days of data fitting; default=10)
% Countries: Country Name given in AfricaCountry in AfricaCountryData.mat
%-------------------------------------------------------------------------%

clear model data params options
load('AfricaCovidData');
load('AfricaCountryData');

global gammas gammam gammaa sigma nsimu meanparaE np ck0 paramstr ncountry nparam color mindays;

%--------------Input----------------%
mindays=10;% minmun days of data fitting
%-----------------------------------%

ncountry=length(AfricaCountry);
sigma=1./5.1;
gammaa=1./5.88;
gammas=1./4.2;
gammam=1./4.2;
nsimu=5000;% number of MCMC simulations
paramstr={'R_0';'\rho_s';'\rho_a';'p_s';'p_a';'r_m';'E_0';'I_{m0}';'I_{s0}';'A_0';'R'};
nparam=length(paramstr);

%color
color=[252 41 30;250 200 205;219 249 244;54 195 201;0 70 222;255 187 0;255 170 50]./255;


meanparaE=zeros(ncountry,nparam+2);% stores mean of all parameter values & r_overall & R_0K


for i=1:ncountry
    global np ck0;
    %countrywise parameter values
    country=AfricaCountry{i};
    name=AfricaCountryName{i};
    np=AfricaPopulationSize(i);
    ck0=CK0(i);
    xtick=XTick(i,1:(length(XTick(XTick(i,:)>0))+1));
    xticklabel=XTickLabel{i};
    
    legendon=1; %legend is on
    if (i>1)
        legendon=0; %legend is off
    end
    
    if (Days(i)>mindays)
        days=Days(i);
    else
        days=mindays;
    end
    
    %data of cumulative cases and new cases
    countrydata=[eval([country,'CumCase']) eval([country,'NewCase'])];
    data.ydata=countrydata((1:days),:);
    data.time=(0:length(data.ydata(:,1))-1)';
    data.ydata=[data.time data.ydata];
    
    theta00=[2.5,1,0.55,0.16,0.46,0.5,ck0,ck0,ck0,ck0,0]';
    [theta0,ss0]=fminsearch(@urmodelss,theta00,[],data);
    mse=ss0/(length(data.ydata)-nparam);

    params = {
   %  name,                 init,       min,    max
    {'R_0',                theta00(1),   1,     3.9}
    {'\rho_s',             theta00(2),   1,     1.5}
    {'\rho_a',             theta00(3),   0,     1}
    {'p_s',                theta00(4),   0,     0.9}
    {'p_a',                theta00(5),   0,     0.9}
    {'r_m' ,               theta00(6),   0,     1}
    {'E_0' ,               theta00(7),  ck0,    np}
    {'I_{m0}' ,            theta00(7),  ck0,    np}
    {'I_{s0}' ,            theta00(8),  ck0,    np}
    {'A_0' ,               theta00(9),  ck0,    np}
    {'R' ,                 theta00(10),  0,     np}
         };
         
    model.ssfun = @urmodelss;
    model.sigma2 = mse;

    options.nsimu = nsimu;
    options.updatesigma = 1; 
    options.method= 'dram'; % adaptation method, 'mh', 'dr', 'am', or 'dram'
    options.adaptint= 500;    % how often to adapt the proposal

    [results,chain,s2chain] = mcmcrun(model,data,params,options);
    mR0=mean(chain(:,1));
    mrhos=mean(chain(:,2));
    mrhoa=mean(chain(:,3));
    mps=mean(chain(:,4));
    mpa=mean(chain(:,5));
    mrm=mean(chain(:,6));
    r=mrm*(1-mps-mpa)+mps;% overall report rate
    mbeta=mR0./(mps*mrhos./gammas+(1-mps-mpa)./gammam+mpa*mrhoa./gammaa);% beta
    R0k=mps*mrhos*mbeta./gammas+(1-mps-mpa)*mrm*mbeta./gammam;% R0 for reported cases
    meanparaE(i,[1 2 8])=[mR0 R0k r];
    for j=3:7
        meanparaE(i,j)=mean(chain(:,j-1));
    end
    for j=9:nparam+2
        meanparaE(i,j)=mean(chain(:,j-2));
    end
    
    %-----------------Mean and SD of all Parameters Mean Values---------------------%
    npara=length(meanparaE(1,:));
    table=zeros(2,npara);
    for i=1:npara
        table(1,i)=mean(meanparaE(:,i));
        table(2,i)=std(meanparaE(:,i));
    end

    
    %---------------------------Countrywise Figures------------------------------------%
    %-----------Parameter Chain Plot-----------%
    clf
    close all
    figure(1);
    mcmcplot(chain,[],results,'chainpanel')
    %mcmcplot(sqrt(s2chain),[],[],'dens',2)
    %title('Error Std')
    sgtitle(country,'FontSize',30);
    %saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Latex/',[country,'ParaPanelWithE.png']]);

    %-----------Parameter Density Plot-----------%
    figure(2);
    for j=1:nparam
        subplot(1,nparam,j)
        histogram(chain(:,j),nsimu);
        title(paramstr{j},'FontSize',20);
    end
    sgtitle(name,'FontSize',30);
    set(gcf,'unit','centimeters','position',[3 5 55 10])
    saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/',[country,'ParaDensityWithE.png']]);
    
    %------------------Fit Plot (NO CI)-------------------%
    figure(3);
    [t,y]=urmodel(data.time,mean(chain));
    plot(data.time,data.ydata(:,2),'o','Color',color(1,:),'LineWidth',2);
    hold on
    plot(t,y(:,2),'-','Color',color(1,:),'LineWidth',2);
    hold on
    plot(data.time,data.ydata(:,3),'o','Color',color(7,:),'LineWidth',2);
    hold on
    plot(t,y(:,7),'-','Color',color(7,:),'LineWidth',2);
    hold on
    plot(t,y(:,1),'-','Color',color(4,:),'LineWidth',2);
    if (legendon==1)
       legend('Data C_K','Modeled C_K','Data I_K','Modeled I_K','Modeled C_I','FontSize',24,'Location','northwest');
    end  
    set(gca,'XTick',xtick,'XTickLabel',xticklabel,'XTickLabelRotation',45,'FontSize',24);
    xlabel('Date','FontSize',32);ylabel('Cases','FontSize',32);
    title(name,'FontSize',46);
    set(gcf,'unit','centimeters','position',[10 15 20 15]);
    %saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/',[country,'FitWithE.png']]);
    
    %------------Fit Plot (WITH CI)--------------%
    figure(4);
    [t,y]=urmodel(data.time,mean(chain));
    out = mcmcpred(results,chain,s2chain,data,@urmodelik,500);
    mcmcpredplot(out,6,7);alpha(0.5);
    hold on
    out = mcmcpred(results,chain,s2chain,data,@urmodelck,500);
    mcmcpredplot(out,2,1);alpha(0.5);
    hold on
    plot(data.time,data.ydata(:,2),'o','Color',color(1,:),'LineWidth',2);
    %hold on
    %plot(t,y(:,2),'-','Color',color(1,:),'LineWidth',2);
    hold on
    plot(data.time,data.ydata(:,3),'o','Color',color(7,:),'LineWidth',2);
    %hold on
    %plot(t,y(:,7),'-','Color',color(4,:),'LineWidth',2);
    axis([0 inf 0 inf]);
    if (legendon==1)
       legend('95% HDI','Modeled I_K','95% HDI','Modeled C_K','Data C_K','Data I_K','FontSize',24,'Location','northwest');
    end
    set(gca,'XTick',xtick,'XTickLabel',xticklabel,'XTickLabelRotation',45,'FontSize',24);
    xlabel('Date','FontSize',32);ylabel('Cases','FontSize',32);
    title(name,'FontSize',46);
    set(gcf,'unit','centimeters','position',[10 15 20 15]);
    saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/',[country,'FitWithE.png']]);

    %---------Predictive Envelop Plot for Cumulative Cases----------%
    figure(5);
    out = mcmcpred(results,chain,s2chain,data,@urmodelci,500);
    mcmcpredplot(out,3,4);
    hold on
    %plot(t,y(:,1),'-','Color',color(1,:),'LineWidth',2);
    %hold on
    plot(data.time,data.ydata(:,2),'o','Color',color(1,:),'LineWidth',2);
    hold on
    plot(t,y(:,2),'-','Color',color(1,:),'LineWidth',2);
    %hold on
    %plot(data.time,data.ydata(:,3),'o','Color',color(7,:),'LineWidth',2);
    %hold on
    %plot(t,y(:,7),'-','Color',color(7,:),'LineWidth',2);
    axis([0 inf 0 inf]);
    if (legendon==1)
       legend('95% HDI','Modeled C_I','Data C_K','Modeled C_K','FontSize',24,'Location','north west');
    end
    set(gca,'XTick',xtick,'XTickLabel',xticklabel,'XTickLabelRotation',45,'FontSize',24);
    xlabel('Date','FontSize',32);ylabel('Cumulative cases','FontSize',32);
    title(name,'FontSize',46);
    set(gcf,'unit','centimeters','position',[10 15 20 15]);
    saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/',[country,'C_IWithE.png']]);
    
end
close all

%---------------------------Summarization Figures------------------------------------%
%-----------R0-----------%
figure(1)
[a,b]=sort(meanparaE(:,1),'ascend');
bb=barh(meanparaE(b,[2 1]),1.7);
set(bb(1),'FaceColor',color(1,:));
set(bb(2),'FaceColor',color(4,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
grid on;
legend({'R_0 for reported cases','R_0'},'FontSize',14,'Location','SouthOutside')
title('Basic reproduction number','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40]);
set(gca,'XAxisLocation','top');
%saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','R0WithE.png']);


%-----------Report Rate-----------%
figure(2)
[a,b]=sort(meanparaE(:,8),'ascend');
bb=barh(meanparaE(b,[7 8])*100,1.7);%r_m and overall r
set(bb(1),'FaceColor',color(2,:));
set(bb(2),'FaceColor',color(1,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
grid on;
legend({'Mild case report rate','Overall report rate'},'FontSize',14,'Location','SouthOutside')
title('Report rate (%)','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
%saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','rWithE.png']);

%---------------Proportion of Severe, Mild and Asympt. Cases----------------------------------%
figure(3)
[a,b]=sort(meanparaE(:,5),'ascend');
initial=[meanparaE(b,5) ones(54,1)-meanparaE(b,5)-meanparaE(b,6) meanparaE(b,6)];% p_s,p_m,p_a
bb=barh(initial,'stacked');
set(bb(1),'FaceColor',color(1,:));
set(bb(2),'FaceColor',color(2,:));
set(bb(3),'FaceColor',color(3,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
axis([0 1 0.4 54.6]);
legend({'Severe cases','Mild cases','Asympt. cases'},'FontSize',14,'Location','SouthOutside')
title('Porportion of severe, mild and asymptomatic cases','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
%saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','pWithE.png']);

%---------------Relative infectiousness of severe and asymptomatic cases----------------------%
figure(4)
[a,b]=sort(meanparaE(:,4),'ascend');
bb=barh(meanparaE(b,[4 3]),1.7);% rho_a,rho_s
set(bb(1),'FaceColor',color(3,:));
set(bb(2),'FaceColor',color(4,:));
ax=gca;
ax.YTick=1:54; 
ax.YTickLabels=AfricaCountryName(b);
axis([0 inf 0.4 54.6]);
grid on;
legend({'Relative infectiousness of asympt. cases','Relative infectiousness of severe cases'},'FontSize',14,'Location','SouthOutside')
title('Relative infectiousness of severe and asymptomatic cases','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
%saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','rhoWithE.png']);


%-----------Number of Severe, Mild, Asympt. and Recoverd Cases (Exposed) at Day 0------------------%
figure(5)
[a,b]=sort(meanparaE(:,13),'ascend');
initial=[meanparaE(b,[13 12 10 11 9])];% R,A_0,I_m0,I_s0 and E
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
legend({'Recovered','Asympt. cases','Mild cases','Severe cases','Exposed'},'FontSize',14,'Location','SouthOutside')
title('Number of initial recovered, asymptomatic, mild, severe cases and exposed','FontSize',22);
set(gcf,'unit','centimeters','position',[3 5 30 40])
set(gca,'XAxisLocation','top');
%saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/','0WithE.png']);



latex2(meanparaE,2,AfricaCountry)



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

%-------------------------------Functions------------------------------------%
function ss=urmodelss(theta,data)
time=data.time; 
[t,y]=urmodel(time,theta); 
ss=sum((data.ydata(:,2)-y(:,2)).^2)+sum((data.ydata(:,3)-y(:,7)).^2); % the total SS
end

function ckmod=urmodelck(data,theta)
time=data.time;
[t,y]=urmodel(time,theta);
ckmod=y(:,2);
end

function cimod=urmodelci(data,theta)
time=data.time;
[t,y]=urmodel(time,theta);
cimod=y(:,1);
end

function ikmod=urmodelik(data,theta)
time=data.time;
[t,y]=urmodel(time,theta);
ikmod=y(:,7);
end

function [t,y]=urmodel(time,theta)
global ck0 np gammas gammam gammaa sigma;
y0(3)=theta(7);
y0(4)=theta(8);
y0(5)=theta(9);
y0(6)=theta(10);
y0(2)=ck0;
y0(1)=sum(theta(8:11));
[t,yy]=ode45(@urode,time,y0,[],theta);
ik=(theta(6)*(1-theta(4)-theta(5))+theta(4))*sigma*yy(:,3);
y=[yy ik];
end


function dy=urode(t,y,theta,np)
% known parameter values
global np gammas gammam gammaa sigma;

% take parameters and components out from y and theta
R0=theta(1);
rhos=theta(2);
rhoa=theta(3);
ps=theta(4);
pa=theta(5);
rm=theta(6);
beta=R0./((ps*rhos)./gammas+(1-ps-pa)./gammam+(pa*rhoa)./gammaa);


% variables
ci=y(1);
ck=y(2);
e=y(3);
im=y(4);
is=y(5);
a=y(6);

% define the ODE
dy(1)=sigma*e;
dy(2)=rm*(1-ps-pa)*sigma*e+ps*sigma*e;
dy(3)=beta*(1-(e+ci)./np)*(rhos*is+im+rhoa*a)-sigma*e;
dy(4)=(1-ps-pa)*sigma*e-gammam*im;
dy(5)=ps*sigma*e-gammas*is;
dy(6)=pa*sigma*e-gammaa*a;

dy=dy(:);
end% make sure that we return a column vector