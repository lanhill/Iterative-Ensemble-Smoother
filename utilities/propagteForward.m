function backgroundEnsemble=propagteForward(modelInfo,analysisEnsemble)
% function backgroundEnsemble=propagteForward(modelInfo,analysisEnsemble)
%
% fields in modelInfo
% --nin: number of model parameter (excluding tspan, iniCon). Integer
% --tspan: time step of numerical integration (for ODEs only). Real
% --iniCon: initial conditions for numerical integration. One-dimensional array
% --val: values of model parameters. One-dimensional array
%
% Please note that script is made in a Linux enviroment.
% By Xiaodong Luo, Oct. 2019
% Copyright (c) 2019 Norwegian Research Centre (NORCE), All Rights Reserved.
%

n=length(analysisEnsemble(1,:)); %ensemble size

% propagating analysis ensemble forward
tmpModelInfo=modelInfo;
intgralStep=modelInfo.tspan(2)-modelInfo.tspan(1);
tmpModelInfo.tspan=0:intgralStep:intgralStep;
for j=1:n,
    tmpModelInfo.iniCon=analysisEnsemble(:,j);
    [tt,td]=myFeval(tmpModelInfo); %#ok<ASGLU>
    backgroundEnsemble(:,j)=td(:,end); %#ok<AGROW>
end


