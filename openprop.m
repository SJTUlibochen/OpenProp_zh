clear global;
clear variables;
close all;
%清除命令行
clc;
%将当前目录下的SourceCode文件夹添加到“搜索路径”（主页-设置路径 可查看当前的搜索路径）中
% ./ 表示当前目录下；../ 表示父级目录下
addpath ./SourceCode;
addpath ../OpenProp_results;
run('./SourceCode/MainPage.m')