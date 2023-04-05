function PSDBootstrap_Getmdl(cx, savename)
    %{
    INPUT:
        cx: time series: Number of time points x Number of nodes
    OUTPUT:
        saved as *.mat
    %}
    thres = -2; % threshold for rejecting higher BIC
    q_pool = 1:50; % the set of AR orders to test
    q_list = {}; mld_list = {};
    for i = 1:size(cx,2)
        x = cx(:,i);
        BIC = []; 
        for k = q_pool
            mdl = ar(x,k);
            BIC(k) = mdl.Report.Fit.BIC;
        end
        idx = find(diff(BIC)<thres == 0);
        q = q_pool(idx(1));
        q_list{i} = q;
        mld_list{i} = ar(x,q);
    end
    save(savename,'mld_list')
end
