function [c_sens] = c_sensitivity_YC(FC,GT,flag)
%{
Estimate the c-sensitivity of the method (refer to Smith NeuroImage 2011)
Symmetrize the GT matrix and take the absolute value of the functional
connectivity matrix
INPUT:
    FC: functional connectivity matrix
    GT: GT matrix 
    flag: if symmetrize GT matrix or not
OUTPUT:
    c_sens: 
%}
	if nargin < 3
		flag = 1;
	end
	if flag == 1
    	GT = GT + transpose(GT);
    end
    TP_value = abs(FC(GT~=0));
    FP_value = abs(FC(GT==0));
    FP_95q = quantile(FP_value,0.95);
    c_sens = sum(TP_value>FP_95q)/length(TP_value);
end

% f = figure;
% subplot(221);imagesc(GT_sym);colorbar
% subplot(222);histogram(TP_value);title('TP values')
% subplot(223);histogram(FP_value);title('FP values')
% saveas(f,'Csens_test.png')