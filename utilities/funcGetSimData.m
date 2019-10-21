function [simData varargout]= funcGetSimData(modelInfo,obvInfo,iniCon) %#ok<*NCOMMA>
%: function simData = funcGetSimData(modelInfo,obvInfo)
%: compute simulated observations
%
% Please note that script is made in a Linux enviroment. 
% By Xiaodong Luo, Oct. 2019
% Copyright (c) 2019 Norwegian Research Centre (NORCE), All Rights Reserved.
%
nstep = size(obvInfo.truth,2);
integrationStep = modelInfo.integrationStep;
assimilationTime = (nstep -1) * integrationStep; 
tWindow = 0:integrationStep:assimilationTime; % integration time window 

%--generating the truths--
tmpModelInfo = modelInfo;
tmpModelInfo.iniCon = iniCon;
tmpModelInfo.tspan = tWindow;
[~,simTrajectory]=myFeval(tmpModelInfo);

tmpData = obvInfo.observer(simTrajectory);
time_index = 1:obvInfo.skipStep:size(tmpData,2); % time instants at which observations are made
simData = tmpData(:,time_index);

if nargout > 1
   varargout{1} = simTrajectory; 
end