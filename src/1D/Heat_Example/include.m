function include( path )
%INCLUDE 用来包含自己的函数库文件
%   此函数用来包含当前路径下名为path的文件夹
curdir = fileparts(mfilename('fullpath'));  %获取当前正在运行函数的路径

% fullfile返回组合文件的全部路径
% genpath递归的产生该文件夹下所有文件的路径
% addpath路径加入到函数或文件搜索的范围内
addpath(genpath(fullfile(curdir, path)));
end

