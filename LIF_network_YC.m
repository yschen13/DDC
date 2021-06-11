 % Modifed based on  /cnl/chaos/ROBERT/wm_intrinsic_timescales/code/matlab/sim/LIF_network_fnc.m 
% Written by Y.C. 09/30/2020
% Add options 02/02/2021

function [REC, spk, rs, fr, params] = LIF_network_YC(W,ext_stim,BIAS,options)

%{
  INPUT
    W: weight matrix
    ext_stim: external stimulation: sampled at 200Hz (100xdt) (N x t)
    BIAS: adjust the excitability of the network
  OUTPUT
    REC: membrane voltage from all the units (N x t)
    spk: binary matrix indicating spikes (N x t)
    rs: firing rates (filtered spikes) from all the units (N x t)
    fr: overall firing rate (N x 1)
%}

	arguments
		W double
		ext_stim double
		BIAS double = -40
		options.tau (1,1) {mustBeNumeric} = 0.01
		options.td(1,1) {mustBeNumeric} = 0.125
	end


	%------------------------------------------------------
	% LIF network parameters
	%------------------------------------------------------
	N = size(W,1);
	dt = 0.00005;   % sampling rate 20000Hz
	T = (size(ext_stim, 2)-1)*dt*100; % external stimulus is at 200Hz
	nt = round(T/dt);
	tref = 0.002;   % refractory time constant (in sec)
	% tm = 0.010;      % membrane time constant (in sec)
	tm = options.tau;
	vreset = -65;   % voltage reset (in mV)
	vpeak = -40;    % voltage peak (in mV) for linear LIF

	% Synaptic parameters
	IPSC = zeros(N,1);      % post synaptic current storage variable
	h = zeros(N,1);         % storage variable for filtered firing rates
	r = zeros(N,1);         % second storage variable for filtered rates
	hr = zeros(N,1);        % third variable for filtered rates
	JD = 0*IPSC;            % storage variable required for each spike time
	tspike = zeros(4*nt,2); % storage variable for spike times
	ns = 0;                 % number of spikes, counts during simulation

	v = vreset + rand(N,1)*(30-vreset); % initialize voltage with random distribtuions
	%v = vreset;
	v_ = v;   % v_ is the voltage at previous time steps
	v0 = v;   % store the initial voltage values

	% Record REC (membrane voltage), Is (input currents), 
	% spk (spike raster), rs (firing rates) from all the units
	REC = zeros(nt,N);  % membrane voltage (in mV) values
	Is = zeros(N, nt);  % input currents from the ext_stim
	IPSCs = zeros(N, nt); % IPSC over time
	%IPSCs = zeros(N, 1); % IPSC over time
	spk = zeros(N, nt); % spikes
	rs = zeros(N, nt);  % firing rates
	hs = zeros(N, nt); % filtered firing rates

	% used to set the refractory times
	tlast = zeros(N,1); 

	% set the BIAS current, can help decrease/increase firing rates. 0 is fine.
	if nargin < 3
	  BIAS = vpeak; % for linear LIF
	end
	%BIAS = 0; % for quadratic LIF
	%BIAS = vreset;

	OMEGA = W;
	td = options.td; tr = 0.002; % parameters to synaptic filtering

	%------------------------------------------------------
	% Start the simulation
	%------------------------------------------------------
	for i = 1:nt
	    IPSCs(:, i) = IPSC; % record the IPSC over time (comment out if not used to save time)

	    I = IPSC + BIAS; % synaptic current

	    % Apply external input stim if there is any
	    I = I + ext_stim(:, round(i/100)+1);
	    Is(:, i) = ext_stim(:, round(i/100)+1);

	    % LIF voltage equation with refractory period
	    dv = (dt*i > tlast+tref).*(-v+I)/tm; % linear LIF
	    v = v + dt*(dv) + randn(N, 1)/10;

	    % find the neurons that have spiked
	    index = find(v>=vpeak);  

	    % store spike times, and get the weight matrix column sum of spikers
	    if length(index)>0
	        JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking
	        tspike(ns+1:ns+length(index),:) = [index, [0*index+dt*i]];
	        ns = ns + length(index);  % total number of spikes so far
	    end

	    % used to set the refractory period of LIF neurons
	    tlast = tlast + (dt*i -tlast).*(v>=vpeak);      

	    % if the rise time is 0, then use the single synaptic filter,
	    % otherwise (i.e. rise time is positive) use the double filter
	    if tr == 0
	        IPSC = IPSC.*exp(-dt./td)+JD*(length(index)>0)./(td);
	        r = r.*exp(-dt./td) + (v>=vpeak)./td;
	        rs(:, i) = r;
	    else
	        IPSC = IPSC.*exp(-dt./td) + h*dt;
	        h = h*exp(-dt/tr) + JD*(length(index)>0)./(tr*td);  %Integrate the current
	        
	        r = r.*exp(-dt./td) + hr*dt;
	        hr = hr*exp(-dt/tr) + (v>=vpeak)./(tr.*td);
	        rs(:, i) = r;
	    end

	    % record the spikes
	    spk(:, i) = v >= vpeak;
	    
	    % spike to 30mV
	    v = v + (30 - v).*(v>=vpeak);

	    % record the membrane voltage tracings from all the units
	    REC(i,:) = v(1:N); 

	    % reset with spike time interpolant implemented.
	    v = v + (vreset - v).*(v>=vpeak); 

	    if mod(i,round(nt/10)) == 0
	      disp(['Simulation progress ' num2str(i/nt)])
	    end
	end

	time = 1:1:nt;

	REC = REC';

	params = {};
	params.dt =  dt;
	params.T = T;
	params.nt = nt;
	params.td = td;
	params.IPSCs = IPSCs;

	fr = zeros(N,1);
	for i = 1:N
	  fr(i) = sum(spk(i,:))/T;
	end
end