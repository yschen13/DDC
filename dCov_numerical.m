function [dCov1,dCov2,dCov5,dCov_center] = dCov_numerical(cx,h,dm)
%{
	Linear dCov computation
	1st order derivative computation: dv/dt = (v(t+1)-v(t))/dt
	INPUT:
		cx: T x N
		h: sampling interval
        dm: averaging window: parameter related to dCov_center
	OUTPUT:
		dCov1: 1st order
		dCov2: 2nd order 
		dCov5: five-point stencil approximation
        dCov_center: centered derivative from Taylor expansion
	NOTE: 
		dCov = <dv/dt,v>
		covariance is computed through cov()
%}

  	addpath('/home/claudia/TOOLS/')

    if nargin < 3
        dm = 4;
    end

    [T,N] = size(cx);
    diff_cx = (cx(2:end,:) - cx(1:end-1,:))/h; % (v(t+1)-v(t))/h
    diff_cx = [diff_cx; mean(diff_cx)]; % last row with mean
    Csample = cov([diff_cx cx]);
    dCov1 = Csample(1:N,N+1:N+N);

    diff_cx = (1/2*cx(3:end,:) - 1/2*cx(1:end-2,:))/h; % (v(t+1)-v(t))/h
    diff_cx = [mean(diff_cx); diff_cx; mean(diff_cx)]; % last row with mean
    Csample = cov([diff_cx cx]);
    dCov2 = Csample(1:N,N+1:N+N);

    diff_cx = (-cx(5:end,:) + 8*cx(4:end-1,:) - 8*cx(2:end-3,:) + cx(1:end-4,:))/(12*h);
    diff_cx = [mean(diff_cx); mean(diff_cx); diff_cx; mean(diff_cx); mean(diff_cx)];
    Csample = cov([diff_cx cx]);
    dCov5 = Csample(1:N,N+1:N+N);

    diff_cx = [];
    for i = 1:N
        [dx,~,~] = derivative_123(cx(:,i),dm,h);
        diff_cx = [diff_cx dx]; 
    end
    cx_trunc = cx(1+dm:T-dm,:);
    Csample = cov([diff_cx cx_trunc]);
    dCov_center = Csample(1:N,N+1:N+N);

end
