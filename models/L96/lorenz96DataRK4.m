function [time,data]=lorenz96DataRK4(tspan,y0,n,F)
%function [time,data]=lorenz96Data(odefun,tspan,y0,options,n,F);
%This rountine is developed to generate Lorenz 96 attractor
%%%%By Xiaodong Luo, June 2007%%%%

% if nargin<6,
%     F=8;
% end
% if nargin<5,
%     n=40;
% end
% if nargin<4,
%     options=odeset('RelTol',1e-4,'AbsTol',1e-8);
% end
% if nargin<3,
%     y0=ones(n,1)+rand(n,1);
% end
% if nargin<2,
%     tspan=0:0.01:40;
% end
% if nargin<1,
%     odefun=@lorenz96ODE;
% end


time=tspan;
data(:,1)=y0;
tStep=(tspan(end)-tspan(1))/(length(tspan)-1);

for i=1:length(tspan)-1,   % 4th order Runge-Kutta method 
    
    k1=lorenz96ODE(data(:,i),n,F);
    tmpSt=data(:,i)+0.5*tStep*k1;
    
    k2=lorenz96ODE(tmpSt,n,F);
    tmpSt=data(:,i)+0.5*tStep*k2;
    
    k3=lorenz96ODE(tmpSt,n,F);
    tmpSt=data(:,i)+tStep*k3;
    
    k4=lorenz96ODE(tmpSt,n,F);
    
    data(:,i+1) = data(:,i) + tStep*(k1+2*k2+2*k3+k4)/6;
end




