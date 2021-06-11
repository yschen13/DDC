function [FPR,TPR,AUC] = FC_GT_ROC(FC,W,cutoff)
%{
	Evaluate Connection Recovery Performance
	INPUT: 
		FC: estimated connection, weighted and directed , remove diagonal values
		W: sparse ground truth connection, binarized during the process
        cutoff: target sparsity (in percentile) of FC
		(s_count: sparsity of the FC matrix (by default: 0:20))
	OUTPUT:
		FPR: false postive rate: FP/Neg;
		TPR: true positive rate: TP/Pos
		AUC: area under the curve
	Written by Y.C. 09/30/2020
%}
	if nargin < 3
		s_count = 0:20; N = size(W,1);
        W_bi = double(W~=0);
        s_true = sum(W_bi,'all') / (N*N);
        cutoff = (1-s_true.*s_count)*100;
    end
	N = size(W,1);
    W_bi = double(W~=0);
	Pos = sum(W_bi(:)); Neg = N*N-Pos; 
	cutoff(cutoff<0) = 0; % make sure percentile is between 0 and 100
	mat = FC - diag(diag(FC));
	thres_list = prctile(abs(mat(:)),cutoff);
	FPR = zeros([length(thres_list) 1]); TPR = zeros([length(thres_list) 1]);
	for i = 1:length(thres_list)
		thres = thres_list(i);
		mat_bi = double(abs(mat)>thres);
		TP = sum(W_bi.*mat_bi,'all');
		FP = sum(double(mat_bi==1&W_bi==0),'all');
		FPR(i) = FP/Neg; TPR(i) = TP/Pos;
	end
	AUC = trapz(FPR,TPR);
end