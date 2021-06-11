function [V_pre] = Linear_simulation(G,deltaT,RandInput,Ttotal,ext)
%{
	INPUT: 	
		G: connectivity pattern, G(i,j): connections from j to i
		deltaT: integration step
		RandInput: magnitude of noisy input
		Ttotal: total simulation duration (sec), default to 1000 sec
		ext: (T x N) external input, default to zero
	OUTPUT:
		V_pre: time points x variables
	Y.C. 04/20/2021: 
		modified to simulate forced system
%}

	arguments
		G double
		deltaT double
		RandInput double
		Ttotal double = 1000
		ext double = NaN
	end
	N = size(G,1);
	T = Ttotal/deltaT; % by default, simulate 1000 seconds
	V_pre = zeros(T, N);
	I = zeros(T, N);
	if all(isnan(ext),'all')
		ext = zeros(T,N);
	end
	tic
	for t= 2:T
		W = randn(N,1)*RandInput; 
		u = W' + ext(t,:);
		I(t,:) = (G*V_pre(t-1, :)')'; % inhibitory coupling
		V_pre(t, :) = V_pre(t-1, :) + (I(t,:)+ u)*deltaT;
		if any(V_pre(t,:)>10000)
			disp('Simulation exploded')
			break
		end
		if mod(t,T/10)==0
	        toc
	        disp(['Simulation iteration: ' num2str(t/T)])
	        tic    
	    end
	end
end
