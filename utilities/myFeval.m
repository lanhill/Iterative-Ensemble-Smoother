function [time,data]=myFeval(modelInfo)
% function [time,data]=myFeval(modelInfo)
% script to evaluate a function handle given the input "modelInfo"
% essentially this is the forward numerical simulator
%
% Please note that script is made in a Linux enviroment.
% By Xiaodong Luo, Oct. 2019
% Copyright (c) 2019 Norwegian Research Centre (NORCE), All Rights Reserved.
%

fhandle=modelInfo.propagator; % fhandle: function handle
nPar=modelInfo.nPar; % number of model parameters
tspan=modelInfo.tspan; % integration time step
iniCon=modelInfo.iniCon; % initial condition for integration

for i=1:nPar
    varargin{i}=modelInfo.val{i}; %#ok<AGROW> value of model parameters
end

[time,data]=feval(fhandle,tspan,iniCon,varargin{:});

