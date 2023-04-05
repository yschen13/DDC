function PSDBootstrap_example(cx,mdl_name,savename)
    %{
    INPUT:
        cx: time series: Number of time points x Number of nodes
        mdl_name: *.mat file where previously fit AR model was stored
    OUTPUT
        null connection mean and variance; saved as *.mat file
    %}
    load(mdl_name)
    TR = 0.72;
    boot_list = {}; boot_mean = {}; boot_var = {};
    for m = 1:11 % across methods
        boot_list{m} = [];
    end
    for j = 1:1000 % shuffling index
        SubMat = {};
        cx_PSDshuffle = cx;
        % generate null time series per node
        for i = 1:size(cx,2) % component index
            u = idinput(size(cx,1),'rgs');
            cx_PSDshuffle(:,i) = sim(mld_list{i},u);
        end
        % estimate FC from null time series
        V_obs = zscore(cx_PSDshuffle);
        [~,~,delta_p,delta_s,~] = FC_Est(V_obs,TR);
        [dCov1, dCov2,~,~] = dCov_numerical(V_obs,TR);
        [Cov,Precision,B,~] = estimators(V_obs,prctile(V_obs(:),50),TR);
        A = {Cov,Precision,inv(B),dCov1,dCov1*Precision,dCov1*inv(B),...
        dCov2,dCov2*Precision,dCov2*inv(B),delta_p,delta_s};
        for m = 1:11 % across methods
            boot_list{m} = cat(3,boot_list{m},A{m});
        end
    end
    for m = 1:11 % across methods
        boot_mean{m} = mean(boot_list{m},3);
        boot_var{m} = var(boot_list{m},0,3);
    end
    save(savename,'boot_mean','boot_var')
end
