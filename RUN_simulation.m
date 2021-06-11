%{
NOTE:
This is a script to simulate time traces for dCov_nonlinear analysis
%}

clear; close all
cd('/home/yuchen/Documents/dCov_nonlinear/')
addpath('/home/yuchen/Documents/fMRI_Real/SupportFunctions/FSLNets/')
addpath('/home/yuchen/Documents/fMRI_Real/SupportFunctions/L1precision/')
addpath('~/Documents/Balloon/Manuscript/SupportFunctions/');

%--------------------------------------------
% 3 NEURON SIMULATION
%--------------------------------------------

% savedir = '/nadata/cnl/data/yuchen/dCov_nonlinear/3neuron/';
savedir = '/nadata/cnl/data/yuchen/dCov_nonlinear/3neuron/forced/';
NTrial = 50;

%% Simulation parameters
% linear simulation
deltaT = 0.01;
Ttotal = 1000;
TR = deltaT;
params.deltaT = deltaT;
params.TR = TR;
params.NTrial = NTrial;
params.Ttotal = Ttotal;
T = round(Ttotal/deltaT); N = 3;
u = zeros(T,N); 
u(T/4:T/4*3,1) = 1;

S_list = -[0.1 0.2 0.5 0.8];
RandInput_list = [0.01 0.1 1 10 20 30 40 50 60 70 80 90 100];
V_pre_list = {}; 
parpool
parfor Input_idx = 1:length(RandInput_list)
	tmp = {};
	for S_idx = 1:length(S_list)
		V_pre_multi = zeros(1000/deltaT,3,NTrial);
		for Trial = 1:NTrial
			S = S_list(S_idx);
			RandInput = RandInput_list(Input_idx);
			G_confounder = -1*eye(3);
			G_confounder(2,1) = S; G_confounder(3,1) = S;
			G_chain = -1*eye(3);
			G_chain(2,1) = S; G_chain(3,2) = S;
			G = G_confounder;
			% G = G_chain;
			% V_pre = Linear_simulation(G,deltaT,RandInput);
			V_pre = Linear_simulation(G,deltaT,RandInput,Ttotal,u);
			V_pre_multi(:,:,Trial) = V_pre;
		end
		tmp{S_idx} = V_pre_multi;
	end
	V_pre_list(Input_idx,:) = tmp;
	disp(['Progress: ' num2str(Input_idx) '/' num2str(length(RandInput_list))])
end
% save([savedir 'Confounder_Ts.mat'],'V_pre_list','S_list','RandInput_list','params','-v7.3')
% save([savedir 'Chain_Ts.mat'],'V_pre_list','S_list','RandInput_list','params','-v7.3')
save([savedir 'Confounder_HalfON_Ts.mat'],'V_pre_list','S_list','RandInput_list','params','-v7.3')


%%  nonlinear simulation
deltaT = 0.0001; 
TR = 0.01; 
S = -0.5;
params.deltaT = deltaT;
params.TR = TR;
params.NTrial = NTrial;
params.S = S;
RandInput_list = [0.01 0.1 1 10 20 30 40 50 60 70 80 90 100];
alpha_list = [1 5 10 15 50]; 
V_pre_list = {}; 

% parpool
parpool('threads')
parfor Input_idx = 1:length(RandInput_list)
	tmp = {};
	for alpha_idx = 1:length(alpha_list)
		V_pre_multi = zeros(1000/TR,3,NTrial);
		for Trial = 1:NTrial
			Alpha = alpha_list(alpha_idx);
			RandInput = RandInput_list(Input_idx);
			G_confounder = -1*eye(3);
			G_confounder(2,1) = S; G_confounder(3,1) = S;
			G_chain = -1*eye(3);
			G_chain(2,1) = S; G_chain(3,2) = S;
			% G = G_confounder;
			G = G_chain;
			V_pre = Nonlinear_simulation(G,deltaT,RandInput,'nonlinearity','sigmoid_sym','parameter',Alpha);
			V = V_pre(1:TR/deltaT:size(V_pre,1),:);
			V_pre_multi(:,:,Trial) = V;
		end
		tmp{alpha_idx} = V_pre_multi;
		disp(['Progress: ' num2str(alpha_idx) '/' num2str(length(alpha_list))])
	end
	V_pre_list(Input_idx,:) = tmp;
end

% save([savedir 'SigmoidSymSim_Confounder_Ts.mat'],'V_pre_list','alpha_list','RandInput_list','params','-v7.3')
save([savedir 'SigmoidSymSim_Chain_Ts.mat'],'V_pre_list','alpha_list','RandInput_list','params','-v7.3')


%----------------------------------------------------------------
% 50 NEURON SIMULATION 
%----------------------------------------------------------------
savedir = '/nadata/cnl/data/yuchen/dCov_nonlinear/50neuron/';

G = load('GT_cxcx56789_uni_latent.mat');
G = -G.G(1:50,1:50);
NTrial = 50;

%% linear simulation
deltaT = 0.01; TR = 0.01; Ttotal = 5000; 
RandInput = 50;
params.deltaT = deltaT;
params.TR = TR;
params.RandInput = RandInput;
params.Ttotal = Ttotal;
V_pre = Linear_simulation(G,deltaT,RandInput,Ttotal);
V = zscore(V_pre);
save([savedir 'Linear_Ts.mat'],'V','G','params')

% linear simulation: multiple trial
deltaT = 0.01; TR = 0.01; Ttotal = 10000; RandInput = 50;
params.deltaT = deltaT;
params.TR = TR;
params.RandInput = RandInput;
params.Ttotal = Ttotal;
save([savedir 'Linear_RandInput50_MultiTrial/Paras.mat'],'params','G')
for Trial = 1:NTrial
	V_pre = Linear_simulation(G,deltaT,RandInput,Ttotal);
	Trial
	save([savedir 'Linear_RandInput50_MultiTrial/' num2str(Trial) '.mat'],'V_pre')
end

% linear simulation: multiple trial across RandInput
RandInput_list = [0.01 0.1 1 10 50 100];
deltaT = 0.01; TR = 0.01; Ttotal = 2000; 
params.deltaT = deltaT;
params.TR = TR;
params.Ttotal = Ttotal;
for i = 1:length(RandInput_list)
	RandInput = RandInput_list(i);
	params.RandInput = RandInput;
	V_pre_list = {};
	for Trial = 1:NTrial
		V_pre = Linear_simulation(G,deltaT,RandInput,Ttotal);
		V_pre_list{Trial} = V_pre;
	end
	save([savedir 'Linear_RandInput_MultiTrial/RandInput' num2str(RandInput) '.mat'],...
		'V_pre_list','params','G','-v7.3')
	['Progress: ' num2str(i/length(RandInput_list))]
end







%% nonlinear simulation

% simulate one situation; long duration; multiple trials
deltaT = 0.0001;TR = 0.01;
RandInput = 50;
Alpha = 1;
Ttotal = 10000;
params.deltaT = deltaT;
params.TR = TR;
params.RandInput = RandInput;
params.Alpha = Alpha;
params.Ttotal = Ttotal;
save([savedir 'SigmoidSymSim_RandInput50_alpha1_MultiTrial/Paras.mat'],'G','params')

NTrial = 50;
V_list = {};
for Trial = 1:NTrial
	V_pre = Nonlinear_simulation(G,deltaT,RandInput,'nonlinearity','sigmoid_sym','parameter',Alpha,'Ttotal',Ttotal);
	% V_pre = Nonlinear_simulation_v2018(G,deltaT,RandInput,'sigmoid_sym',Alpha,Ttotal);
	V = V_pre(1:TR/deltaT:size(V_pre,1),:);
	save([savedir 'SigmoidSymSim_RandInput50_alpha1_MultiTrial/' num2str(Trial) '.mat'],'V')
end


% simulate across randinput
deltaT = 0.0001;TR = 0.01;
% Alpha = 1;
Alpha = 50; 
Ttotal = 10000;
params.deltaT = deltaT;
params.TR = TR;
params.Alpha = Alpha;
params.Ttotal = Ttotal;
RandInput_list = [0.01 0.1 1 10 50 100];
NTrial = 50;
for Rand_idx = 1:length(RandInput_list)
	RandInput = RandInput_list(Rand_idx);
	params.RandInput = RandInput;
	folder = ['SigmoidSymSim_RandInput' num2str(RandInput) '_alpha' num2str(Alpha) '_MultiTrial/'];
	if ~exist([savedir folder], 'dir')
		mkdir([savedir folder])
	end
	save([savedir folder 'Paras.mat'],'G','params')
	V_list = {};
	for Trial = 1:NTrial
		% V_pre = Nonlinear_simulation(G,deltaT,RandInput,'nonlinearity','sigmoid_sym','parameter',Alpha,'Ttotal',Ttotal);
		V_pre = Nonlinear_simulation_v2018(G,deltaT,RandInput,'sigmoid_sym',Alpha,Ttotal);
		V = V_pre(1:TR/deltaT:size(V_pre,1),:);
		save([savedir folder num2str(Trial) '.mat'],'V')
	end
end



