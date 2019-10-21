function [modelInfo,obvInfo,methodInfo] = setupCase(obvSkip,skipStep,obv_nLevel,ensize,targeted_nstep,varargin)
%%
% function [modelInfo,obvInfo,methodInfo] = setupCase(obvSkip,skipStep,obv_nLevel,ensize,targeted_nstep)
%
% "setupCase.m" is used to specify the experiment settings in a data
% assimilation (DA) problem. Such a script is case depent. The current
% script is made for the Lorentz 96 (L96) model
%
% -- INPUTS--
% obvSkip:          parameter specifying the spatial density/frequency of observations at a given time instance, 
%                   in the form of $1:obvSkip:dynDim$, where $dynDim$ is the dimension of the dynamical system (L96 in this case)  
% skipStep:         Roughly speaking, this is the parameter telling the system to conduct DA every $skipStep$ step.
%                   specify the frequency of conducting temporal frequency of observations
% obv_nLevel:       parameter specifying the noise level in the observations, in the form of obv_nLevel * eye(obvDim),
%                   where $obvDim$ represents the dimension of observations at a given time instance
% ensize:           parameter specifying the ensemble size
% targeted_nstep:   parameter specifying the length of time window where DA is conducted, in the form of $integrationStep * targeted_nstep$, 
%                   where $integrationStep$ is the time-step size for ODE model integration.  
% varargin:         variable input reserved for future development
% 
% -- OUTPUTS--
% modelInfo,obvInfo,methodInfo: information with regard to dynamical model, observation system, and DA method, respectively.
%
%
% Please note that script is made in a Linux enviroment. 
% By Xiaodong Luo, Oct. 2019
% Copyright (c) 2019 Norwegian Research Centre (NORCE), All Rights Reserved.
%
%%
% set up path for MATLAB to locate correct scripts 
% please note that my ASSUMPTION is that this script ("setupCase.m")
% is run from its parent folder ("L96 in this case")
%
% Please make sure go to the folder "./run_example/L96" before running the 
% script "setupCase.m". Otherwise there may be a problem in using "cd ../.."
% and "addpath" afterwards. 
cd ../.. % move two levers up to get the correct path information
rootPath = pwd; %
addpath(genpath(rootPath)); % add the folder "IES" into MATLAB search path
cd([rootPath '/run_example/L96'])

% long term statistics (mean and cov) generated from a long free-run of the
% L96 model
data_dir = [rootPath '/models/L96/']; 
lt_mean = load([data_dir 'lt_mean.dat']); % long term mean
lt_cov = load([data_dir 'lt_cov.dat']); % long term mean

%% default input values
if nargin < 5,
   targeted_nstep = 20; % 5-day forcast (4 steps per day) 
end
if nargin < 4,
    ensize = 100;
end
if nargin < 3,
    obv_nLevel = 1;
end
if nargin < 2,
    skipStep = 1;
end
if nargin < 1,
    obvSkip = 1;
end

%% modelInfo
%
modelInfo.propagator=@lorenz96DataRK4; % dynamical system
modelInfo.nPar=2; % number of model parameters 

transitionTime=25; % transition time
integrationStep=0.05; % integration step
assimilationTime = integrationStep * targeted_nstep; % assimilation time window
modelInfo.integrationStep=integrationStep;
numAssimilationCycles=assimilationTime/integrationStep; % number of integration steps within the assimilation time window
modelInfo.tspan=0:integrationStep:integrationStep; % transition time window
modelInfo.tWindow=0:integrationStep:(transitionTime+assimilationTime); % integration time window 

nDOF=40;
F=8;
dynDim=nDOF; 
modelInfo.dynDim=dynDim; % dimension of the dynamical system
sigma=eye(dynDim); % cov of multivariate normal distribution
mu=zeros(dynDim,1); % mean of multivariate normal distribution
rng('default'); % random seed
modelInfo.iniCon=mvnrnd(mu,sigma,1)'; % initial condition for integration 
modelInfo.val={nDOF,F}; % values of model parameters
modelInfo.Q=zeros(nDOF,nDOF); % no stochastic model error 

%% obvInfo
%
obvInfo.targeted_nstep = targeted_nstep;
obvInfo.obvSkip = obvSkip;
obvInfo.isLinObs = 0; % 1 = linear observation; 0 = nonlinear observation

% observation operator
%state_index = [1:15 22 27 31 34 36]; % irregular obs network
state_index = 1:obvInfo.obvSkip:dynDim; % regular obs network   
if obvInfo.isLinObs
    obvInfo.linOpt = 'lin'; 
    obvInfo.observer=@(x) x(state_index,:); 
else
    obvInfo.linOpt = 'nln';
    obvInfo.observer=@(x) x(state_index,:).^3/5; % in case of nonlinear observation operator there is no matrix
end
%obvInfo.H = sparse(obvInfo.observer(eye(dynDim))); % get the observatio operator H
obvInfo.obs_site_index = state_index;

% observation dimension 
obvDim=length(obvInfo.observer(randn(dynDim,1)));
obvInfo.R = obv_nLevel * eye(obvDim); % covariance of the observation noise
%obvInfo.inv_R = (obvInfo.R)^(-1);
[U,S,V] = svd(obvInfo.R);
obvInfo.principal_sqrtR = U * sqrt(S) * (V');
obvInfo.principal_sqrtR_inv = U * sqrt(S^(-1)) * (V');

% generate the true trajectory
tmpModelInfo=modelInfo;
tmpModelInfo.tspan=modelInfo.tWindow;
[~,truth]=myFeval(tmpModelInfo);
nStep=numAssimilationCycles; % number of steps in assimilation
obvInfo.truth=truth(:,end-(nStep-1):end); % assimilating the last nStep time steps

% generate the initial background ensemble
obvInfo.enSize=ensize; %
obvInfo.iniBGEnsemble = mvnrnd(lt_mean,lt_cov,obvInfo.enSize)';

% generate observations
obvMu=zeros(obvDim,1); % mean of multivariate normal distribution
obv = obvInfo.observer(obvInfo.truth)+mvnrnd(obvMu,obvInfo.R,length(obvInfo.truth(1,:)))';
obvInfo.skipStep = skipStep; % conducting assimilation for every $skipStep$ steps
obvInfo.obv = obv(:,1:skipStep:end); %  observations used in assimilation

%% methodInfo
% localization setting
methodInfo.doLoc=0; % do localization? (1 = yes)
if methodInfo.doLoc
    % localization is not included here, but can be done if needed
    % see, for example, the recent paper
    % "Automatic and adaptive localization for ensemble-based history
    % matching", JPSE, 2019. https://doi.org/10.1016/j.petrol.2019.106559
    % and the references therein
end

% configuration of the iterative ensemble smoother, for more information, see
% "Iterative ensemble smoother as an approximate solution to a regularized 
% minimum-average-cost problem: theory and applications", SPE J. 2015.
% https://www.onepetro.org/journal-paper/SPE-176023-PA

methodInfo.beta = 0; % $beta$ determines the threshold value in one of the stopping criteria
methodInfo.maxOuterIter = 20; % maximum iteration number in the outer loop
methodInfo.maxInnerIter = 5; % maximum iteration number in the inner loop
methodInfo.init_lambda = 1; % initial lambda value
methodInfo.lambda_reduction_factor = 0.9; % reduction factor in case to reduce gamma
methodInfo.lambda_increment_factor = 2; % increment factor in case to increase gamma 
methodInfo.doTSVD = 1; % do a TSVD on the cov of simulated obs? (1 = yes)
if methodInfo.doTSVD
    methodInfo.tsvdCut = 0.99; % discard eigenvalues/eigenvectors if they are not among the truncated leading ones
end
methodInfo.min_RN_change = 1; % minimum residual norm (RN) change (in percentage); RN(k) - RN(k+1) > RN(k) * min_RN_change / 100

%
methodInfo.rootPath = rootPath;
methodInfo.casePath = [rootPath '/run_example/L96'];
methodInfo.results_dir = [methodInfo.casePath '/results'];
if ~exist(methodInfo.results_dir,'dir')
    mkdir(methodInfo.results_dir);
end

save([methodInfo.results_dir '/experiment_settings.mat'],'modelInfo','obvInfo','methodInfo');


