% %This script plots m-year recurrence values of the GPD model
% %foe electron fluxes based on EVT
close all
restoredefaultpath
rehash toolboxcache


%u=threshold, m=occurrence interval in years, N=sample size, n_u=sample >u
%d=observations per year
u=5e5;
m=1:1:100;
lm=ones(1,length(m));
N=118263;
n_u=26;
d=8447.357143; %Based on actual samnple divided by total measurement period

x_p=d*m*n_u/N;

%EVT parameter estimates for shape and scale, respectively
xi=1.229;
sgm=129228;

%m year-recurence level of x>u
x_mu=u+sgm/xi*((x_p).^xi-1);


%Partial derivatives
dxi_x_mu=sgm/xi^2*(x_p).^xi.*(xi*log(x_p)-1)+sgm/xi^2;
dsig_x_mu=x_p.^xi/xi-1/xi;

%s.d. of sigma and xi
var_sigma=2863386428;
var_xi=0.191093885;
cov_xi_sgm=11078.81585;

%Error of x_u
var_x_mu=dxi_x_mu.^2*var_xi+dsig_x_mu.^2*var_sigma+2*dxi_x_mu.*dsig_x_mu.*cov_xi_sgm;
sd_x_mu=sqrt(var_x_mu);

%Plot x_mu with 2 sigma confidence intervals
fig1=figure();
linew=2;
set(fig1,'DefaultLineLineWidth',linew);
ax=subplot(1,1,1);


plot(m,x_mu,'DisplayName','m-year recurrence flux')
hold on
set(ax,'YScale', 'log')
ax.set
ylabel('flux (cm-2s-1sr-1keV-1)')
xlabel('m (years)')
plot(m,x_mu+1.96*sd_x_mu,'b--','DisplayName','95% conf.int.')
plot(m,x_mu-1.96*sd_x_mu,'b--','HandleVisibility','off');
plot(m,lm*10^5,'r--','DisplayName','flux=10^5')
%Plot fluxes at 25, 50, 100 and 150 year marks
f25=x_mu(25);
f50=x_mu(50);
f100=x_mu(100);
%f150=x_mu(150);
%plot(m,lm*f25,'k','DisplayName','flux at m=25')
plot(m,lm*f50,'k:','DisplayName','flux at m=50')
plot(m,lm*f100,'k-.','DisplayName','flux at m=100')
%plot(m,lm*f150,'k-.','DisplayName','flux at m=150')

legend('location','east')
set(ax,'FontSize',14)
title('m-year recurrence values of the 130 keV flux')
ylim([0,max(x_mu+1.96*sd_x_mu)]);

fname=strcat('Figures/','m_year_rec_130keV_u=',num2str(u),'_xi=',num2str(xi),'_sgm=',num2str(sgm),'.png');
saveas(fig1,fname);

%Change y axis to linear scale and x axis to log scale and save the figure
hold off
ylim([-10^8,max(x_mu+1.96*sd_x_mu)]);
set(ax,'YScale', 'linear')
set(ax,'XScale', 'log')
legend('location','northwest')
fname=strcat('Figures/','m_year_rec_130keV_u=',num2str(u),'_xi=',num2str(xi),'_sgm=',num2str(sgm),'_xlog.png');
saveas(fig1,fname);