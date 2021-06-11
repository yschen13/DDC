function [Cov,precision,B,dCov] = estimators(V_obs,thres,TR)
%{
	INPUT: 	
		V_obs: time points x variables
		thres: Relu offset
		TR: sampling interval (seconds)
	OUTPUT:
		B: <ReLu(x),x>
		dCov: <dx/dt,x>
    Version2: change the implementation of ReLu, get rid of the translation
%}


	[T,N] = size(V_obs);
	Cov = cov(V_obs);
	precision = inv(Cov);
	Fx = V_obs-thres; Fx(Fx<0) = 0;
 	% Fx = V_obs; Fx(Fx<thres) = 0;
	tmp = cov([Fx V_obs]);
	B = tmp(1:N,N+1:end);
	dV = (-1/2*V_obs(1:end-2,:) + 1/2*V_obs(3:end,:))/TR; % (z(t+1)-z(t-1))/2
	dV = [mean(dV);dV;mean(dV)]; % substitute the first and last row with mean
	tmp = cov([dV V_obs]);
	dCov = tmp(1:N,N+1:end);

end
