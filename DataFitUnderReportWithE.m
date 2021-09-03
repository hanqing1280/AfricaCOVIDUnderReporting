%----------------------------Introduction------------------------------%
% Countrywise Data Fit with EXPOSED STAGE
% Method: MCMC
% Data Set: Reported Cumulative Cases + Reported New Cases
% Input: Country Name (Mandantory)/ mindays (minimum days of data fitting)
%        / setdays (Selected;defaultly not used)
% Country Name are in:
% 'Algeria''Angola''Benin''Botswana''BurkinaFaso''Burundi''Cameroon'
% 'CapeVerde''CentralAfricanRepublic''Chad''Comoros''Congo''Cotedlvoire'
% 'DemocraticRepublicCongo''Djibouti''Egypt''EquatorialGuinea''Eritrea'
% 'Eswatini''Ethiopia''Gabon''Gambia''Ghana''Guinea''GuineaBissau'
% 'Kenya''Lesotho''Liberia''Libya''Madagascar''Malawi''Mali''Mauritania'
% 'Mauritius''Morocco''Mozambique''Namibia''Niger''Nigeria''Rwanda'
% 'SaoTomeAndPrincipe''Senegal''Seychelles''SierraLeone''Somalia'
% 'SouthAfrica''SouthSudan''Sudan''Tanzania''Togo''Tunisia''Uganda'
% 'Zamibia''Zimbabwe'
%----------------------------------------------------------------------%

clear model data params options
load('AfricaCovidData');
load('AfricaCountryData');

global country ck0 np days sigma gammas gammam gammaa nsimu paramstr nparam color mindays;

%--------------Input----------------%
country='CentralAfricanRepublic';% mandantory input
mindays=10;% defaulty 10 minimun days is used
setdays=0;% selected input;defaultly not used
%-----------------------------------%

%countrywise parameter values
index=find(strcmp(AfricaCountry,country));
name=AfricaCountryName{index};
np=AfricaPopulationSize(index);
ck0=CK0(index);

if (setdays==0)
    if (Days(index)>mindays)
        days=Days(index);
    else
        days=mindays;
        fprintf('Days are set to %d\n',mindays)
    end
else days=setdays;
     fprintf('Days are set to %d\n',setdays)
end    

%colors
color=[252 41 30;250 200 205;219 249 244;54 195 201;0 70 222;255 187 0;255 170 50]./255;

%data of cumulative cases and new cases
countrydata=[eval([country,'CumCase']) eval([country,'NewCase'])];
data.ydata=countrydata((1:days),:);
data.time=(0:length(data.ydata(:,1))-1)';
data.ydata=[data.time data.ydata];
sigma=1./5.1;
gammaa=1./5.88;
gammas=1./4.2;
gammam=1./4.2;
nsimu=5000;% number of MCMC simulations
paramstr={'R_0';'\rho_s';'\rho_a';'p_s';'p_a';'r_m';'E_0';'I_{m0}';'I_{s0}';'A_0';'R'};
nparam=length(paramstr);

theta00=[2.5,1,0.55,0.16,0.46,0.5,ck0,ck0,ck0,ck0,0]';
[theta0,ss0]=fminsearch(@urmodelss,theta00,[],data)
mse = ss0/(length(data.ydata)-nparam);

params = {
   %  name,                 init,       min,    max
    {'R_0',                theta00(1),   1,     3.9}
    {'\rho_s',             theta00(2),   1,     1.5}
    {'\rho_a',             theta00(3),   0,     1}
    {'p_s',                theta00(4),   0,     1}
    {'p_a',                theta00(5),   0,     1}
    {'r_m' ,               theta00(6),   0,     1}
    {'E_0' ,               theta00(7),  ck0,    np}
    {'I_{m0}' ,            theta00(8),  ck0,    np}
    {'I_{s0}' ,            theta00(9),  ck0,    np}
    {'A_0' ,               theta00(10), ck0,    np}
    {'R' ,                 theta00(11),  0,     np}
         };
         
model.ssfun = @urmodelss;
model.sigma2 = mse;

options.nsimu = nsimu;
options.updatesigma = 1; 
options.method= 'dram'; % adaptation method, 'mh', 'dr', 'am', or 'dram'
options.adaptint= 500;    % how often to adapt the proposal

[results,chain,s2chain] = mcmcrun(model,data,params,options);
r=mean(chain(:,6)).*(1-mean(chain(:,4))-mean(chain(:,5)))+mean(chain(:,4));% overall report rate


chainstats(chain,results)
fprintf('r mean=%f\n',r);


%-------------------------------Figures------------------------------------%
%--------Parameter Chain Plot--------%
clf 
close all
figure(1);
mcmcplot(chain,[],results,'chainpanel')
%mcmcplot(sqrt(s2chain),[],[],'dens',2)
%title('Error Std')
sgtitle(name);
%saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Latex/',[country,'ParaPanelWithE.png']]);

%-------Parameter Histogram Plot--------%
figure(2);
for i=1:nparam
    subplot(1,nparam,i)
    histogram(chain(:,i),nsimu);
    title(paramstr{i},'FontSize',20);
end
sgtitle(name,'FontSize',30);
set(gcf,'unit','centimeters','position',[3 5 55 10])
%saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Latex/',[country,'ParaDensityWithE.png']]);

%------------Fit Plot (NO CI)--------------%
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
hold on
set(gca, 'XTick', XTick(index,1:(length(XTick(XTick(index,:)>0))+1)),'XTickLabel',XTickLabel{index},'XTickLabelRotation',45,'FontSize',18);
xlabel('Date','FontSize',30);ylabel('Cases','FontSize',30);
legend('Data C_K','Modeled C_K','Data I_K','Modeled I_K','Modeled C_I','FontSize',20,'Location','northwest');
title(name,'FontSize',42);
set(gcf,'unit','centimeters','position',[10 15 20 15]);
%saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Latex/',[country,'FitWithE.png']]);

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
%plot(t,y(:,6),'-','Color',color(4,:),'LineWidth',2);
axis([0 inf 0 inf]);
set(gca, 'XTick', XTick(index,1:(length(XTick(XTick(index,:)>0))+1)),'XTickLabel',XTickLabel{index},'XTickLabelRotation',45,'FontSize',18);
xlabel('Date','FontSize',30);ylabel('Cases','FontSize',30);
legend('95% HDI','Modeled I_K','95% HDI','Modeled C_K','Data C_K','Data I_K','FontSize',20,'Location','northwest');
title(name,'FontSize',42);
set(gcf,'unit','centimeters','position',[10 15 20 15]);
%saveas(gcf,['/Users/Qing/Desktop/IDRC/Under-report/Latex/',[country,'FitWithE.png']]);

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
set(gca, 'XTick', XTick(index,1:(length(XTick(XTick(index,:)>0))+1)),'XTickLabel',XTickLabel{index},'XTickLabelRotation',45,'FontSize',18);
legend('95% HDI','Modeled C_I','Data C_K','Modeled C_K','FontSize',18,'Location','north west');
xlabel('Date','FontSize',30);ylabel('Cumulative cases','FontSize',30);
title(name,'FontSize',42);
set(gcf,'unit','centimeters','position',[10 15 20 15]);


%---------Parameter ACF Plot-------------%
figure(6);
mcmcplot(chain,[],results.names,'acf',100)
sgtitle(name);
set(gcf,'unit','centimeters','position',[3 5 30 40])

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

function imod=urmodeli(data,theta)
global sigma
time=data.time;
[t,y]=urmodel(time,theta);
imod=y(:,3)*sigma;
end

function [t,y]=urmodel(time,theta)
global ck0 gammas gammam gammaa np sigma;
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


function dy=urode(t,y,theta)
% known parameter values
global gammas gammam gammaa np sigma;

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