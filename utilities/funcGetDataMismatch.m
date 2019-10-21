function [obj,objStd,objReal]=funcGetDataMismatch(simData,W,measurment)
%: function [obj,objStd,objReal]=funcGetDataMismatch(simData,W,measurment)
%: compute the mismatch between simData and measurement;
%: both simData and measurement are pre-normalized by Wbase
%
% Please note that script is made in a Linux enviroment. 
% By Xiaodong Luo, Oct. 2019
% Copyright (c) 2019 Norwegian Research Centre (NORCE), All Rights Reserved.


ne=size(simData,2);
objReal=zeros(ne,1);
for j=1:ne
    objReal(j)=sum(((simData(:,j)-measurment).^2)./(W.^2)); % for diagonal weight matrices only
end
 
obj=mean(objReal);
objStd=std(objReal);