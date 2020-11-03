function output = iES(modelInfo,obvInfo,methodInfo)
% output = iES(modelInfo,obvInfo,methodInfo)
%
% We need to run a "setupCase.m" script first to get modelInfo,obvInfo,methodInfo
%  
% This scipt implements the iterative ensemble smoother (iES) in the
% paper "Iterative ensemble smoother as an approximate solution to a regularized 
% minimum-average-cost problem: theory and applications", SPE J. 2015.
% https://www.onepetro.org/journal-paper/SPE-176023-PA
%
% -- INPUTS --
% See a "setupCase.m" script for the meanings of modelInfo,obvInfo,methodInfo
%
% The scipt "iES.m" implements the iterative ensemble smoother (iES) in the
% paper "Iterative ensemble smoother as an approximate solution to a regularized 
% minimum-average-cost problem: theory and applications", SPE J. 2015.
% https://www.onepetro.org/journal-paper/SPE-176023-PA
%
% The implementation here follows the paper "Levenberg–Marquardt forms of the
% iterative ensemble smoother for efficient history matching and uncertainty 
% quantification" by Chen and Oliver, Computational Geosciences, 2013. 
% 
% Please note that script is made in a Linux enviroment. 
% By Xiaodong Luo, Oct. 2019
% Copyright (c) 2019 Norwegian Research Centre (NORCE), All Rights Reserved.
%

%%
% record time and activate diary
ticID = tic;
    
% specify the other dir 
cd(methodInfo.results_dir)

% "debug.txt" is a log file. We need to remove the previous log before
% starting a new experiment
if exist('./debug.txt','file')
   delete('./debug.txt');
end

obj_dir = './obj';
if ~exist(obj_dir,'dir')
    mkdir(obj_dir);
end

simData_dir = './simData';
if ~exist(simData_dir,'dir')
    mkdir(simData_dir);
end

ensemble_dir = './ensemble';
if ~exist(ensemble_dir,'dir')
    mkdir(ensemble_dir);
end

%--  initialization  --%
iter = 0;
lambda = methodInfo.init_lambda;
isDisplay = 1; % display iteration information or not? (1=yes)

% load information from obvInfo 
measurement = reshape(obvInfo.obv,numel(obvInfo.obv),1);% augment the obs at different times into a "super" vector
if  ~exist([simData_dir '/measurement.mat'],'file')
    save([simData_dir '/measurement.mat'],'measurement');
end
nd = length(measurement); % dim of obs

%
ensemble = obvInfo.iniBGEnsemble;
ne = size(ensemble,2); % ensemble size
ensemble = [ensemble mean(ensemble,2)]; % note that ensemble mean is appended to the ensemble at the end
if  ~exist([ensemble_dir '/ensemble','-iter',num2str(iter)],'file')
    save([ensemble_dir '/ensemble','-iter',num2str(iter)],'ensemble','lambda');
end

%--
Wbase = []; % augmented observation error std (in a vector)
for i = 1 : size(obvInfo.obv,2)
    Wbase = [Wbase;diag(obvInfo.principal_sqrtR)];
end
if  ~exist([simData_dir '/Wbase.mat'],'file')
    save([simData_dir '/Wbase.mat'],'Wbase');
end
measurement = measurement ./ Wbase; % normalized measurement
perturbedData=zeros(nd,ne); % peturbations to observations
W = ones(size(measurement));
for j=1 : ne
    perturbedData(:,j) = measurement + W.*randn(size(measurement)); % normalized perturbations 
end

% simulated observations of the ensemble members 
nm = size(ensemble,1); % dim of system state  
simData = [];
for i = 1 : size(ensemble,2)
    tmpData = funcGetSimData(modelInfo,obvInfo,ensemble(:,i));
    simData = [simData reshape(tmpData,numel(tmpData),1)]; %#ok<*AGROW>
    clear tmpData;
end

outer_simData_name = [simData_dir '/simData','-iter',num2str(iter),'.mat'];
if  ~exist(outer_simData_name,'file')
    save(outer_simData_name,'simData');
end

simData = simData ./ repmat(Wbase,1,size(simData,2));% normalized simulated obs;


% data mismatch w.r.t the ensemble 
[obj,objStd,objReal]=funcGetDataMismatch(simData(:,1:ne),W,measurement);
myDisplay(['   obj=',num2str(obj),', objStd=',num2str(objStd)],isDisplay);
outer_objReal_name = [obj_dir '/objReal', '-iter',num2str(iter),'.mat'];
if  ~exist(outer_objReal_name,'file')
    save(outer_objReal_name,'objReal');
end
init_obj = obj;
init_objStd = objStd;

% load information from methodInfo to configure iES
beta = methodInfo.beta;
objThreshold = beta^2 * nd;
myDisplay(['   objThreshold=',num2str(objThreshold)],isDisplay);
maxOuterIter = methodInfo.maxOuterIter;
maxInnerIter = methodInfo.maxInnerIter;
init_lambda = methodInfo.init_lambda;
lambda_reduction_factor = methodInfo.lambda_reduction_factor;
lambda_increment_factor = methodInfo.lambda_increment_factor;
if methodInfo.doTSVD
    tsvdCut = methodInfo.tsvdCut;
end
min_RN_change = methodInfo.min_RN_change;
isTooSmallRNChange = 0;
exitFlag = [0 0 0]; % flags of iES termination status; 1st => maxOuterIter; 2nd => objThreshold; 3rd => min_RN_change

%-------------%
%-- run iES --%
%-------------%

% outer iteration loops 
while (iter < maxOuterIter) && (obj > objThreshold) 
    
    %--
    myDisplay(['-- Outer iteration step: ' int2str(iter) ' --'],isDisplay);
    myDisplay(['    number of measurement elements is ',num2str(numel(measurement))],isDisplay);
    
    % deviation of the ensemble 
    deltaM = ensemble(:,1:ne) - ensemble(:,ne+1) * ones(1,ne); %#ok<*NASGU>
    
    % deviation of the simulated observations 
    deltaD = simData(:,1:ne) - simData(:,ne+1)*ones(1,ne);

    %-- if do TSVD --%
    if methodInfo.doTSVD
        [Ud,Wd,Vd]=svd(deltaD,'econ');
        
        val=diag(Wd);
        total=sum(val);
        for j=1:ne
            svdPd=j;
            if (sum(val(1:j))/total > tsvdCut)
                break
            end
        end
        
        myDisplay(['    svdPd=',num2str(svdPd)],isDisplay);
        Vd=Vd(:,1:svdPd);
        Ud=Ud(:,1:svdPd);
        Wd=val(1:svdPd); % a vector
        clear val;
    end
    
    
    %-- initialization of the inner loop --%
    iterLambda=1;
    
    %-- inner iteration loop --%
    while iterLambda < maxInnerIter
        
        myDisplay(['    -- Inner iteration step: ' int2str(iterLambda) '--'],isDisplay);
        
        ensembleOld = ensemble; % keep a copy of old status
        simDataOld = simData;
        
        if methodInfo.doTSVD
            alpha = lambda * sum(Wd.^2) / svdPd;
            %alpha = lambda * sum(Wd) / svdPd; % alternative rule 
            x1 = Vd * spdiags( Wd ./ (Wd.^2 + alpha), 0, svdPd, svdPd) ;
            KGain=deltaM*x1*Ud'; % Kalman gain
        else
            alpha = lambda * sum(sum(deltaD.^2)) / nd; % sum(sum(deltaD.^2)) = trace (deltaD * deltaD')
            KGain = deltaM * deltaD / (deltaD * deltaD' + alpha * eye(nd));
        end
        
        % localization can be implemented here
        % if methodInfo.doLoc
        %     taperMtx = localization(); %
        %     KGain = taperMtx .* KGain;
        % end
        
        iterated_ensemble = ensemble(:,1:ne) - KGain * (simData(:,1:ne) - perturbedData);
        ensemble = [iterated_ensemble mean(iterated_ensemble,2)];
        
        % check the change of ensemble mean 
        changeM = sqrt( sum( ( ensemble(:,ne+1) - ensembleOld(:,ne+1) ).^2 ) / nm );
        myDisplay(['        average change (in RMSE) of the ensemble mean = ',num2str(changeM)],isDisplay);
        
        % produce simData 
        simData = [];
        for i = 1 : size(ensemble,2)
            tmpData = funcGetSimData(modelInfo,obvInfo,ensemble(:,i));
            simData = [simData reshape(tmpData,numel(tmpData),1)]; %#ok<*AGROW>
            clear tmpData;
        end
        
        inner_simData_name = [simData_dir '/simData','-iter',num2str(iter), '-lamdaIter', num2str(iterLambda) '.mat'];
        if  ~exist(inner_simData_name,'file')
            save(inner_simData_name,'simData');
        end
        
        simData = simData ./ repmat(Wbase,1,size(simData,2));% normalized simulated obs;
        
        % check data mismatch 
        [objNew,objStdNew,objRealNew]=funcGetDataMismatch(simData(:,1:ne),W,measurement);
        myDisplay(['        objNew=',num2str(objNew),' objStdNew=',num2str(objStdNew)],isDisplay);
        inner_objReal_name = [obj_dir '/objReal', '-iter', num2str(iter) '-lamdaIter', num2str(iterLambda), '.mat'];
        if  ~exist(inner_objReal_name,'file')
              tmp_objReal = objReal;
              objReal = objRealNew;
              save(inner_objReal_name,'objReal');
              objReal = tmp_objReal;
              clear tmp_objReal;
        end
        
        %--
        if objNew>obj
            lambda = lambda * lambda_increment_factor;
            myDisplay(['         increasing Lambda to ',num2str(lambda)],isDisplay);
            iterLambda = iterLambda + 1;
            simData    = simDataOld;
            ensemble   = ensembleOld;
        else
            changeStd=(objStdNew-objStd)/objStd;
            myDisplay(['        changeStd=',num2str(changeStd)],isDisplay);
            
            lambda = lambda * lambda_reduction_factor;%
            myDisplay(['        reducing Lambda to ',num2str(lambda)],isDisplay);
            
            iter=iter+1;
            
            outer_simData_name = [simData_dir '/simData','-iter',num2str(iter),'.mat'];
            outer_objReal_name = [obj_dir '/objReal', '-iter',num2str(iter),'.mat'];
            
            % save files
            % 1) entire ensemble of model varialbes
            if  ~exist([ensemble_dir '/ensemble','-iter',num2str(iter)],'file')
                save([ensemble_dir '/ensemble','-iter',num2str(iter)],'ensemble','lambda');
            end
            
            % 2) simulated data
            sysCom = ['mv ',inner_simData_name, ' ', outer_simData_name];
            status=unix(sysCom);
            if status ~= 0
                error('error in mv files simulatedDataIter');
            end
            
            % 3) data mismatch at each iteration
            sysCom = ['mv ',inner_objReal_name, ' ', outer_objReal_name];
            status=unix(sysCom);
            if status ~= 0
                error('error in mv files objRealIter');
            end
            
            if abs(objNew-obj)/abs(obj)*100 < min_RN_change
                isTooSmallRNChange=1;
            end
            
            simDataOld=simData;
            ensembleOld=ensemble;
            objStd=objStdNew;
            obj=objNew;
            objReal = objRealNew;
            break % break the inner loop over lambda
        end
    end % end inner loop

    % if a better update not successfully found  
    
    if iterLambda >= maxInnerIter
        
        lambda = lambda * lambda_increment_factor;
        if lambda < init_lambda
            lambda = init_lambda;
        end
        
        if  ~exist([ensemble_dir '/ensemble','-iter',num2str(iter)],'file')
            save([ensemble_dir '/ensemble','-iter',num2str(iter)],'ensemble','lambda');
        end
        
        if  ~exist([ensemble_dir '/objReal','-iter',num2str(iter)],'file')
            save([obj_dir '/objReal','-iter',num2str(iter)],'objReal');
        end
        
        if  ~exist([ensemble_dir '/simData','-iter',num2str(iter)],'file')
            save([simData_dir '/simData','-iter',num2str(iter)],'simData');
        end
        
        iter=iter+1;

        %--
        tline = '       terminating inner iterations: iterLambda >= maxInnerIter';
        myDisplay(tline,isDisplay);
    end
    
    
    
    if isTooSmallRNChange 
        tline = ['  terminating outer iterations: reduction of objective function is less than ', num2str(min_RN_change),'%'];
        myDisplay(tline,isDisplay);
        exitFlag(3) = 1;
        break
    end
end

if iter >= maxOuterIter    
    tline = '   terminating outer iterations: iter >= maxOuterIter';
    myDisplay(tline,isDisplay);
    exitFlag(1) = 1;
end
if obj <= objThreshold
    tline = '   terminating outer iterations: obj <= objThreshold';
    myDisplay(tline,isDisplay);
    exitFlag(2) = 1;
end

%---------------------%
%-- post processing --%
%---------------------%

%-- mean trajectory and obs --%
[meanObs,meanTrajectory] = funcGetSimData(modelInfo,obvInfo,ensemble(:,ne+1));
rmse_dir = './rmse';
if  ~exist(rmse_dir,'dir')
    mkdir(rmse_dir);
end
save([rmse_dir,'/rmse.mat'],'meanObs','meanTrajectory');

%--
diff_trajectory = meanTrajectory - obvInfo.truth;
diff_obs = meanObs - obvInfo.obv;
save([rmse_dir,'/rmse.mat'],'diff_obs','diff_trajectory','-append');

rmse.state = zeros(1,size(diff_trajectory,2));
rmse.obs = zeros(1,size(diff_obs,2));
for j = 1 : size(diff_trajectory,2)
    rmse.state(j) = norm(diff_trajectory(:,j)) / sqrt(nm);
end
for j = 1 : size(diff_obs,2)
    rmse.obs(j) = norm(diff_obs(:,j)) / sqrt(nd);
end
save([rmse_dir,'/rmse.mat'],'rmse','-append');

%--
output.init_obj = init_obj;
output.init_objStd = init_objStd;
output.state_rmse = rmse.state;
output.obs_rmse = rmse.obs;
output.obj = obj;
output.objStd = objStd;
output.exitFlag = exitFlag;
output.iter = iter;

%
elapsedTime = toc(ticID);
output.elapsedTime = elapsedTime; % in 'second'
myDisplay(['Runtime in minutes: ' num2str(elapsedTime/60) 'm'],isDisplay);

%
cd(methodInfo.casePath);
end
