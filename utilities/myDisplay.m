function myDisplay(s,isDisplay)
% display a message "s" (optional), and write it into a log file "debug.txt"
%
% Please note that script is made in a Linux enviroment. 
% By Xiaodong Luo, Oct. 2019
% Copyright (c) 2019 Norwegian Research Centre (NORCE), All Rights Reserved.

if nargin < 2
    isDisplay = 1;
end
    
fid = fopen('./debug.txt','a+');
fprintf(fid,'%s \n',s);
fclose(fid);

if isDisplay
    disp(s);
end