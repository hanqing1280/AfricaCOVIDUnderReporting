%----------------------------Introduction---------------------------------%
% All 54 African Countries Data Fit
% Method: MCMC
% Data Set: Reported Cumulative Cases (C_K) + Reported New Cases (I_K)
% Input: mindays (minimum days of data fitting; default=9)
% Output: Country-wise Grahphs + Matrices meanpara & meanstd into Latex matrix form
% Countries: Country Name given in AfricaCountry in AfricaCountryData.mat
%-------------------------------------------------------------------------%
addpath('/Users/Qing/Desktop/IDRC/Under-report/Under-reporting MATLAB Codes/mcmcstat-master')

clear model data params options
load('AfricaCovidData');
load('AfricaCountryData');

global gammas gammam gammaa sigma nsimu meanpara np ck0 paramstr ncountry nparam color mindays;

%--------------Input----------------%
mindays=9; % minmun days of data fitting
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


meanpara=zeros(ncountry,nparam+3);% stores mean of all parameter values & r_overall & R_0K & (1-p_a)*p_s (instead of p_s) & (1-p_a)*(1-p_s)


for i=1:ncountry
    global np ck0;
    %countrywise parameter values
    country=AfricaCountry{i};
    name=AfricaCountryName{i};
    np=AfricaPopulationSize(i);
    ck0=CK0(i);
    xtick=XTick(i,1:(length(XTick(XTick(i,:)>0))+1));
    xticklabel=XTickLabel{i};% dates considered for each country
    
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
    
    theta00=[2.5,1,0.55,0.202,0.46,0.5,ck0,ck0,ck0,ck0,0]';
    [theta0,ss0]=fminsearch(@urmodelss,theta00,[],data);
    mse=ss0/(length(data.ydata)-nparam);

    params = {
   %  name,                 init,       min,    max
    {'R_0',                theta00(1),   1,     3.9}
    {'\rho_s',             theta00(2),   1,     1.5}
    {'\rho_a',             theta00(3),   0,     1}
    {'p_s',                theta00(4),   0,     1}
    {'p_a',                theta00(5),   0,     1}
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
    options.adaptint= 100;    % how often to adapt the proposal

    %resulting chain of parameters
    [results,chain,s2chain] = mcmcrun(model,data,params,options);
    stat=chainstats(chain,results)
    str_var=[country,'ChainStat'];
    eval([str_var,'=stat']);
    if (i==1)
       save('ChainStatistics.mat',[country,'ChainStat']);% save chain statistics for each country
    end
    save('ChainStatistics.mat',[country,'ChainStat'],'-append');% save chain statistics for each country
    
    mR0=mean(chain(:,1));
    mrhos=mean(chain(:,2));
    mrhoa=mean(chain(:,3));
    mps=mean(chain(:,4));
    mpa=mean(chain(:,5));
    mrm=mean(chain(:,6));
    r=mrm*(1-mpa)*(1-mps)+(1-mpa)*mps;% overall report rate
    mbeta=mR0./((1-mpa)*mps*mrhos./gammas+(1-mpa)*(1-mps)./gammam+mpa*mrhoa./gammaa);% beta
    R0k=(1-mpa)*mps*mrhos*mbeta./gammas+(1-mpa)*(1-mps)*mrm*mbeta./gammam;% R0 for reported cases
    meanpara(i,[1 2 3 4 5 6 7 8 9])=[mR0 R0k mrhos mrhoa (1-mpa)*mps*100 (1-mpa)*(1-mps)*100 mpa*100 mrm*100 r*100];
    for j=10:nparam+3
        meanpara(i,j)=mean(chain(:,j-3));
    end
    
    %-----------------Mean and SD of all Parameters Mean Values---------------------%
    npara=length(meanpara(1,:));
    meanstd=zeros(2,npara);
    for i=1:npara
        meanstd(1,i)=mean(meanpara(:,i));
        meanstd(2,i)=std(meanpara(:,i));
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
    saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/',[country,'ParaPanel.png']]);

    %-----------Parameter Density Plot-----------%
    figure(2);
    for j=1:nparam
        subplot(1,nparam,j)
        [counts,centers] = hist(chain(:,j),50);
        b=bar(centers,counts/sum(counts),'k');
        threshold=mean(chain(:,j));
        yrange=get(gca,'ylim');
        hold on;
        plot([threshold threshold],yrange,'r-','LineWidth',1);
        title(paramstr{j},'FontSize',20);
    end
    sgtitle(name,'FontSize',30);
    set(gcf,'unit','centimeters','position',[3 5 55 10])
    saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/',[country,'ParaDensity.png']]);
    
    %------------------Fit Plot (NO HDI)-------------------%
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
    
    
    %------------Fit Plot (WITH HDI)--------------%
    figure(4);
    [t,y]=urmodel(data.time,mean(chain));
    out = mcmcpred(results,chain,s2chain,data,@urmodelik,500);
    mcmcpredplotwithcolor(out,6,7);alpha(0.5);
    hold on
    out = mcmcpred(results,chain,s2chain,data,@urmodelck,500);
    mcmcpredplotwithcolor(out,2,1);alpha(0.5);
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
    saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/',[country,'Fit.png']]);

    %---------Predictive Envelop Plot for Cumulative Cases----------%
    figure(5);
    out = mcmcpred(results,chain,s2chain,data,@urmodelci,500);
    mcmcpredplotwithcolor(out,3,4);
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
    saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Manuscript/',[country,'C_I.png']]);
    
end
close all

%----------Output all parameters mean values & meta mean std----------------%
latex2(meanpara,2,AfricaCountry)
latex2noname(meanstd,2,{'Mean' 'SD'})
save('MeanPara.mat','meanpara','meanstd');% save meanpara meanstd


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
global ck0 np sigma;
y0(3)=theta(7)/np;
y0(4)=theta(8)/np;
y0(5)=theta(9)/np;
y0(6)=theta(10)/np;
y0(2)=ck0/np;
y0(1)=sum(theta(8:11))/np;
[t,yy]=ode45(@urode,time,y0,[],theta);
ik=(theta(6)*(1-theta(4))*(1-theta(5))+(1-theta(5))*theta(4))*sigma*yy(:,3);
y=[yy ik]*np;
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
beta=R0./(((1-pa)*ps*rhos)./gammas+(1-ps)*(1-pa)./gammam+(pa*rhoa)./gammaa);


% variables
ci=y(1);
ck=y(2);
e=y(3);
im=y(4);
is=y(5);
a=y(6);

% define the ODE
dy(1)=sigma*e;
dy(2)=rm*(1-pa)*(1-ps)*sigma*e+(1-pa)*ps*sigma*e;
dy(3)=beta*(1-(e+ci)./np)*(rhos*is+im+rhoa*a)-sigma*e;
dy(4)=(1-ps)*(1-pa)*sigma*e-gammam*im;
dy(5)=(1-pa)*ps*sigma*e-gammas*is;
dy(6)=pa*sigma*e-gammaa*a;

dy=dy(:);
end% make sure that we return a column vector