function A_reg = dCov_linear_Reg(V,TR,lambda)
    % L2 Regularized version of deltaL
    % INPUT:
        % V: time series
        % TR: sampling interval
        % lambda: regularization strength
    % OUTPUT:
        % A_reg; linear DDC with ridge regularization
    [T, N] = size(V); % number of timepoints x number of nodes
    V_obs = zscore(V);
    [dCov1,dCov2,~,dCov_center] = dCov_numerical(V_obs,TR);
    [Cov,Precision,B,~] = estimators(V_obs,0,TR);
    C = dCov2;
    B = Cov;
    A_reg = zeros(size(C));
    % lamda = 1e-2;
    for i = 1:size(C,1)
        ci = C(i,:);
        % [coef,stats]=lasso(B,transpose(ci),'Lambda',lamda);
        coef = ridge(transpose(ci),B,lamda);
        A_reg(i,:) = transpose(coef);
    end
end