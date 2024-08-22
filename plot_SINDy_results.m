%This script plots vsw max (measured) vs. the e flux of 130 keV electrons
%predicted by the SINDy-PI model. The data is contained in a separate Excel
%sheet.
close all
restoredefaultpath
rehash toolboxcache

vsw_rising=[627.5824176,517.4847652,564.6291209,579.3986014,703.6433566,628.5856643,433.9160839,529.4410589,448.6785714];
vsw_stable=[554.79,561.60,489.65];
vsw_nan=[689.46,672.95,357.56,627.78];
flux_stable=[329306.2702,9204.495718,11574.44032];
flux_rising=[44963.52,27915.68,64112.42,34453.69,53312.38,56743.67,30385.33,49243.88,21635.31];
flux_nans=[60458.87,8637.74,100716.36,3218.10];

%Create tables
tbl_rising=table(vsw_rising',flux_rising');
tbl_rising_stable=table([vsw_rising,vsw_stable]',[flux_rising,flux_stable]');
tbl_rising_stable_nans=table([vsw_rising,vsw_stable,vsw_nan]',[flux_rising,flux_stable,flux_nans]');

% %Specify model and do fitting
% modelspec = 'Var2 ~ 1 + Var1';
% mdl_rising = fitlm(tbl_rising,modelspec,'RobustOpts','on');
% mdl_rising_stable = fitlm(tbl_rising_stable,modelspec,'RobustOpts','on');
% mdl_rising_stable_nans = fitlm(tbl_rising_stable_nans,modelspec,'RobustOpts','on');

% %Calculate 95% confidence intervals
% ci_rising=coefCI(mdl_rising);
% ci_rising_stable=coefCI(mdl_rising_stable);
% ci_rising_stable_nans=coefCI(mdl_rising_stable_nans);
% ci_level_flux='95%';

%Set KP limit for plotting
flux_KP_L5=2e4;
flux_KP_L6=1.3e4;
flux_KP_L67=9.3e4;
ci_KP_L5=[flux_KP_L5/3,flux_KP_L5*3];
ci_KP_L6=[flux_KP_L6/3,flux_KP_L6*3];
ci_KP_L67=[flux_KP_L67/3,flux_KP_L67*3];

%Set EVT results for plotting
flux_EVT=1.7e8;
sd_EVT=3.6e8;
ci_EVT=flux_EVT+2*[-sd_EVT,sd_EVT];
nsigma_EVT=2;

%Use the function plot_results for plotting
% plot_results('rising_flux',mdl_rising,flux_EVT,flux_KP_L5,flux_KP_L6,ci_EVT,ci_KP_L5,ci_KP_L6)
% plot_results('rising_stable_flux',mdl_rising_stable,flux_EVT,flux_KP_L5,flux_KP_L6,ci_EVT,ci_KP_L5,ci_KP_L6)
% plot_results('rising_stable_nans_flux',mdl_rising_stable_nans,flux_EVT,flux_KP_L5,flux_KP_L6,ci_EVT,ci_KP_L5,ci_KP_L6)
do_model('rising_flux',tbl_rising,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,false,false)
do_model('rising_stable_flux',tbl_rising_stable,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,false,false)
do_model('rising_stable_nans_flux',tbl_rising_stable_nans,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,false,false)
do_model('rising_flux',tbl_rising,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,true,false)
do_model('rising_stable_flux',tbl_rising_stable,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,true,false)
do_model('rising_stable_nans_flux',tbl_rising_stable_nans,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,true,false)
do_model('rising_flux',tbl_rising,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,false,true)
do_model('rising_stable_flux',tbl_rising_stable,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,false,true)
do_model('rising_stable_nans_flux',tbl_rising_stable_nans,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,false,true)
do_model('rising_flux',tbl_rising,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,true,true)
do_model('rising_stable_flux',tbl_rising_stable,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,true,true)
do_model('rising_stable_nans_flux',tbl_rising_stable_nans,flux_EVT,flux_KP_L5,flux_KP_L6,flux_KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,true,true)

function do_model(fluxes,tbl,EVT,KP_L5,KP_L6,KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,handle_outliers,robust_fitting)   
    %Specify model and do fitting
    modelspec = 'Var2 ~ 1 + Var1';
    mdl_fluxes=fitlm(tbl,modelspec);
%     mdl_rising = fitlm(tbl_rising,modelspec,'RobustOpts','on');
%     mdl_rising_stable = fitlm(tbl_rising_stable,modelspec,'RobustOpts','on');
%     mdl_rising_stable_nans = fitlm(tbl_rising_stable_nans,modelspec,'RobustOpts','on');

    %Remove outliers and do robust fitting, if desired
    if handle_outliers==false
        removed_outliers=false;
    else
        [outliers_leverage,outliers_Cook]=find_outliers(mdl_fluxes);
        removed_outliers=true;
        
        %Remove rows corresponding to outliers
        tbl([outliers_leverage],:)=[];
        tbl([outliers_Cook],:)=[];
%         tbl_rising_stable([outliers_leverage],:)=[];
%         tbl_rising_stable([outliers_Cook],:)=[];
%         tbl_rising_stable_nans([outliers_leverage],:)=[];
%         tbl_rising_stable_nans([outliers_Cook],:)=[];
    end

    if robust_fitting==true
        rob_ops='on';
    else
        rob_ops='off';
    end

    %Now do refitting based on removal of outliers and include robust
    %fitting, if selected
    ref_mdl_fluxes = fitlm(tbl,modelspec,'RobustOpts',rob_ops);
%     mdl_rising_stable = fitlm(tbl_rising_stable,modelspec,'RobustOpts',rob_ops);
%     mdl_rising_stable_nans = fitlm(tbl_rising_stable_nans,modelspec,'RobustOpts',rob_ops);

    plot_results(fluxes,ref_mdl_fluxes,EVT,KP_L5,KP_L6,KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,removed_outliers,robust_fitting,1)

    %Model coefficients and standard errors
%     c_vsw=mdl_fluxes.Coefficients.Estimate(2);
%     c_vsw_se=mdl_fluxes.Coefficients.SE(2);
%     icept=mdl_fluxes.Coefficients.Estimate(1);
%     icept_se=mdl_fluxes.Coefficients.SE(1);
    carr_vsw=1500;
%     carr_flux=c_vsw*carr_vsw+icept;
%     carr_flux_sd=c_vsw_se*carr_vsw+icept_se;

    %Print model name and R2 value and Carrington scale flux confidence
    %intervals
    fprintf(strcat('\n130keV_',fluxes,'_SINDy_vs_KP_vs_EVT_','removed_outliers=',num2str(removed_outliers),' robust fitting=',num2str(rob_ops)))
    fprintf('\n\nModel coefficients:\n')
    disp(mdl_fluxes.Coefficients)
    
    %Return predicted flux for Carrington scale vsw with confidence
    %intervals
    alpha=0.05;
    [carr_flux,carr_flux_sd]=predict(mdl_fluxes,carr_vsw,'Alpha',alpha);
    ci_level=100*(1-alpha);
    
    fprintf('\nFor a Carrington scale event with vsw=1500 km/s the electron flux, based on the linear model, is: %.4e',carr_flux)
    fprintf('\nWith a lower %.2i confidence interval bound: %.4e',ci_level,carr_flux_sd(1))
    fprintf('\nand the upper %.2i confidence interval bound: %.4e',ci_level,carr_flux_sd(2))
    fprintf('\nModel R^2=%.4f\n\n',mdl_fluxes.Rsquared.Ordinary)
        
end

%function plot_results(fluxes,mdl,EVT,KP_L5,KP_L6,ci_EVT,ci_KP_L5,ci_KP_L6,rem_outliers)
function plot_results(fluxes,mdl,EVT,KP_L5,KP_L6,KP_L67,ci_EVT,ci_KP_L5,ci_KP_L6,ci_KP_L67,removed_outliers,robust_fitting,p_SINDY)
    %Create x values for plotting
    xs_min=min(mdl.Variables.Var1);
    xs_max=max(mdl.Variables.Var1);
    xs_dif=xs_max-xs_min;
    xs=xs_min:xs_dif/99:xs_max;
    ones_xs=ones([length(xs),1]);
    lxs=length(xs);

    %Plot the results
    fig1=figure();
    ax=subplot(1,1,1);
    linew=2;

    set(fig1,'DefaultLineLineWidth',linew);
    
    %SINDy-PI results
    hold on
    if p_SINDY==0
        
    else
        plot(mdl,'DisplayName','SINDy-PI flux');
    end
    set(ax, 'YScale', 'log')
    set(ax, 'XScale', 'log')
    
    %EVT
    %Allow plotting negative data on log axis
    %warning off MATLAB:Axes:NegativeDataInLogAxis
    %semilogy(1500,EVT,'DisplayName','EVT flux','Color',"#0072BD");
    %fill([xs_min,xs_min,xs_max,xs_max],[ci_EVT,flip(ci_EVT)],"blue",'FaceAlpha',0.2,'HandleVisibility', 'off')
    %semilogy(1500,ci_EVT(1),'DisplayName','EVT CI','LineStyle',':','Color',"#0072BD")
    %semilogy(1500,ci_EVT(2),'DisplayName','none','LineStyle',':','Color',"#0072BD",'HandleVisibility', 'off')
    sd_EVT=(ci_EVT(2)-EVT)/2;
    %Enfroce the lower error bar to a flux value of 1e5.
    errb=errorbar(1500,EVT,min(EVT-1e5,2*sd_EVT),2*sd_EVT,'DisplayName','EVT flux','Color',"#0072BD",'marker','x','LineWidth',2);
    %plot(1500,EVT,'bo','DisplayName','none')
    %plot([1500,1500],[100,EVT+2*sd_EVT])

    %warning on MATLAB:Axes:NegativeDataInLogAxis
    
    %KP limits for L5, L6 and 6.7
%     semilogy(xs,ones_xs*KP_L5,'DisplayName','KP flux L5','Color','cyan');
%     fill([xs_min,xs_min,xs_max,xs_max],[ci_KP_L5,flip(ci_KP_L5)],"cyan",'FaceAlpha',0.2,'HandleVisibility', 'off')
%     semilogy(xs,ones_xs*ci_KP_L5(1),'DisplayName','KP CI L5','LineStyle',':','Color','cyan')
%     semilogy(xs,ones_xs*ci_KP_L5(2),'DisplayName','none','LineStyle',':','Color','cyan','HandleVisibility', 'off')
%     
%     semilogy(xs,ones_xs*KP_L6,'DisplayName','KP flux L6','Color',"#EDB120");
%     fill([xs_min,xs_min,xs_max,xs_max],[ci_KP_L6,flip(ci_KP_L6)],"red",'FaceAlpha',0.2,'HandleVisibility', 'off')
%     semilogy(xs,ones_xs*ci_KP_L6(1),'DisplayName','KP CI L6','LineStyle',':','Color',"#EDB120")
%     semilogy(xs,ones_xs*ci_KP_L6(2),'DisplayName','none','LineStyle',':','Color',"#EDB120",'HandleVisibility', 'off')

    semilogy(xs,ones_xs*KP_L67,'DisplayName','KP flux L=6.7','Color','cyan');
    fill([xs_min,xs_min,xs_max,xs_max],[ci_KP_L67,flip(ci_KP_L67)],"cyan",'FaceAlpha',0.2,'HandleVisibility', 'off')
    semilogy(xs,ones_xs*ci_KP_L67(1),'DisplayName','KP CI L=6.7','LineStyle',':','Color','cyan')
    semilogy(xs,ones_xs*ci_KP_L67(2),'DisplayName','none','LineStyle',':','Color','cyan','HandleVisibility', 'off')
   
    hold off
    %Set plot properties, title, legend and axis titles
    R2=mdl.Rsquared.Ordinary;
    title(strcat('130 keV fluxes fron SINDy-PI, EVT and KP limit, R2=',num2str(R2)))
    ylabel('flux (cm-2s-1sr-1keV-1)')
    xlabel('solar wind speed (km/s)')
    xlim([0,1600]);
    
    lgnd=legend('location','southeast');
    set(lgnd,'color','none');
               
    set(ax,'FontSize',14);
    %grid on
    set(fig1,'Position',[100 100 600 400]);
    set(fig1,'PaperPositionMode','auto');
    
    if removed_outliers==false
        rem_out="";
    else
        rem_out='removed_outliers_';
    end
    if robust_fitting==false
        rob_fit="";
    else
        rob_fit="robust_fit_";
    end

    if p_SINDY==1
        %Plot the SINDy-PI derived regression model
        plot_model(fluxes,mdl,rob_fit,rem_out)
    end

    %Save plot
    fprintf(fluxes)
    fname=strcat('Figures/','130keV_',fluxes,'_SINDy_vs_KP_vs_EVT_',rem_out,rob_fit,'.png');
    saveas(fig1,fname);
    
%     %Analyse diagnostics
%     t_leverage=2*mdl.NumCoefficients/mdl.NumObservations;
%     %Find observations that exceed t_leverage
%     outliers_leverage=find((mdl.Diagnostics.Leverage)>t_leverage);
%     %Find observation whose Cook's distance exceed 3 times the mean
%     %Cook's distance
%     outliers_Cook=find((mdl.Diagnostics.CooksDistance)>3*mean(mdl.Diagnostics.CooksDistance));
%     fprintf('Outliers based on leverage')
%     outliers_leverage
%     fprintf('Outliers based on Cook''s distance')
%     outliers_Cook
end

function [outliers_leverage,outliers_Cook]=find_outliers(mdl)
    %Analyse diagnostics
    t_leverage=2*mdl.NumCoefficients/mdl.NumObservations;
    %Find observations that exceed t_leverage
    outliers_leverage=find((mdl.Diagnostics.Leverage)>t_leverage);
    %Find observation whose Cook's distance exceed 3 times the mean
    %Cook's distance
    outliers_Cook=find((mdl.Diagnostics.CooksDistance)>3*mean(mdl.Diagnostics.CooksDistance));
%     fprintf('Outliers based on leverage')
%     outliers_leverage
%     fprintf('Outliers based on Cook''s distance')
%     outliers_Cook
end

function plot_model(fluxes,model,rob_fit,rem_out)
%Do a scatter plot of the data and the model
    fig15=figure();
    linew=2;
    set(fig15,'DefaultLineLineWidth',linew);
    

    ax1=subplot(1,1,1);
    set(ax1,'FontSize',14)
    
    plot(model)
    xlabel('Maximum of solar wind speed (kms/s)')
    ylabel('Maximum of electron flux (cm-2s-1sr-1keV-1)')
    title('Regression from SINDy-PI results')
    fname=strcat('Figures/','130keV_',fluxes,'_SINDy_model',rem_out,rob_fit,'.png');
    saveas(fig15,fname);
        
end