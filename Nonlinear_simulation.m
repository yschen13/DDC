function [V_pre] = Nonlinear_simulation(G,deltaT,RandInput,options)
%{
	By default, simulate 1000 seconds
	dx/dt = G*F(x) + u 
	INPUT: 	
		G: connectivity pattern, G(i,j): connections from j to i
		deltaT: integration step
		RandInput: magnitude of noisy input
		options.nonlinearity: relu or sigmoid or sigmoid_sym
		options.parameter: parameter controlling relu offset or sigmoid slope
	OUTPUT:
		V_pre: time points x variables
%}
	arguments
		G double
		deltaT double
		RandInput double
		options.nonlinearity (1,1) string = 'relu'
		options.parameter (1,1) double = 1
		options.Ttotal (1,1) double = 1000
	end

	if options.nonlinearity == 'relu'
		disp('ReLu nonlinearity')
		F = @relu;
	elseif options.nonlinearity == 'sigmoid'
		disp('sigmoid nonlinearity')
		F = @sigmoid;
	elseif options.nonlinearity == 'sigmoid_sym'
		disp('sigmoid sym nonlinearity')
		F = @sigmoid_sym;
	end
	N = size(G,1);
	T = options.Ttotal / deltaT; % by default, simulate 1000 seconds
	V_pre = zeros(T, N);
	I = zeros(T, N);
	tic
	for t= 2:T
		u = randn(1,N)*RandInput; ;
		I(t,:) = (G*F(V_pre(t-1,:))')'; 
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

	function Fx = relu(x) % x is a vector
		Fx = x - options.parameter; Fx(Fx<0) = 0; 
	end
	function Fx = sigmoid(x)
		Fx = 1./(1+exp(-options.parameter*x));
	end
	function Fx = sigmoid_sym(x)
		Fx = 1./(1+exp(-options.parameter*x)) - 0.5;
	end
end