function [D1,D2,D3]=derivative_123(f,dm,dt)

%% Claudia  27.6.2012


t = 1+dm:length(f)-dm; 

D1=0;d=0;
for n1=1:dm,
    for n2=n1+1:dm,
        d=d+1;
        D1=D1 + (-((f(t-n2)*n1.^3-f(t+n2)*n1.^3-f(t-n1)*n2.^3+ ...
                    f(t+n1)*n2.^3)/(2*dt*n1.^3*n2-2*dt*n1*n2.^3)));
    end;
end;
D1=D1./d;  clear d;


D2=0;d=0;
for n1=1:dm,
    for n2=n1+1:dm,
        d=d+1;
        D2=D2 + ((f(t-n2)*n1.^4+f(t+n2)*n1.^4-f(t-n1)*n2.^4- ...
                 f(t+n1)*n2.^4-2*f(t)*(n1.^4-n2.^4))/ ...
                 (dt.^2*n2.^2*(n1.^4-n1.^2*n2.^2)));
    end;
end;
D2=D2./d;  clear d;

D3=0;d=0;
for n1=1:dm,
  for n2=n1+1:dm,
    for n3=n2+1:dm,
        d=d+1;
        D3=D3 + ((3*(f(t-n3)*n1*n2*(n1.^4-n2.^4)+ ...
                 f(t+n3)*(-(n1.^5*n2)+n1*n2.^5)+ ...
                 n3*((f(t-n1)-f(t+n1))*n2*(n2.^4-n3.^4)+ ...
                 f(t+n2)*(n1.^5-n1*n3.^4)+f(t-n2)*(-n1.^5+n1*n3.^4))))/ ...
                 (dt.^3*n1*(n1.^2-n2.^2)*n3*(n1.^2-n3.^2)*(n2.^3-n2*n3.^2)));
    end;
  end;
end;
D3=D3./d;  clear d;







