% generate nonstationary (time-variant covariance) data using mixed SDE
RandInput = 10;
N = size(G,1);
D1 = eye(3); % noise structure 
D2 = [1 0 1;1 1 0;0 1 1];
% numerically calculated as covariance
Q = eye(N)*RandInput^2*deltaT;
Cov1_v = -inv(kron(eye(N),G) + kron(G,eye(N)))*(kron(D1,D1))*reshape(eye(N),[N^2 1]);
Cov1 = reshape(Cov1_v,[N N]);
Cov2_v = -inv(kron(eye(N),G) + kron(G,eye(N)))*(kron(D2,D2))*reshape(eye(N),[N^2 1]);
Cov2 = reshape(Cov2_v,[N N]);


T = Ttotal/deltaT; % by default, simulate 1000 seconds
V_pre = zeros(T, N);
I = zeros(T, N);
Noise = zeros(T,N);
for t= 2:T
	if t < T/2
		D = D1;
	else
		D = D2;
	end
	W = D*randn(N,1)*RandInput;
	u = W';
	Noise(t-1,:) = u;
	I(t,:) = (G*V_pre(t-1, :)')'; % inhibitory coupling
	V_pre(t, :) = V_pre(t-1, :) + (I(t,:)+ u)*deltaT;
	if mod(t,T/10)==0
        toc
        disp(['Simulation iteration: ' num2str(t/T)])
        tic    
    end
end

