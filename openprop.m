clear global;
clear all;
close all;
%清除命令行
clc;

%disp用于显示括号内变量的值 
disp('% ----------------------------------------------------------------------- %')
disp('%                                                                         %')
disp('%                              0111000                                    %')
disp('%                           100 1 100 001                                 %')
disp('%                         10    1  1  1 00                                %')
disp('%                        01  1  1  1      0                               %')
disp('%                       0 1  1  1   1  1 1 0                              %')
disp('%                       0   1   1   1  1  1 0                             %')
disp('%                       0 1     1   1  1  1 0                             %')
disp('%                       0 1  1  1   1  0  1 0                             %')
disp('%                       0 1  1  1   0  1    0                             %')
disp('%                       01 1        1  1 1 0                              %')
disp('%                        0    0  1  0 1   0                               %')
disp('%                         0         1    0                                %')
disp('%                    10010 0 1101111110111                                %')
disp('%                  10 1 1  1111111111 11 11                               %')
disp('%                 0 1 1 1 11111111101011010111                            %')
disp('%                01 11    11111111 1  1    1 110                          %')
disp('%               011    1 1 111111110011  1 1 1 110                        %')
disp('%               0   11 1 1 1 111      0  1 1 1   10                       %')
disp('%               0 1   11  1  0         1 1 1 1 1 1 0                      %')
disp('%               1  11 1 1   11          0  1 1 1 1 11                     %')
disp('%                0     1 1  0           011  1 1 1 10                     %')
disp('%                10 1   1  0             0  1 1 1  11                     %')
disp('%                 10     01               01      10                      %')
disp('%                   10001                   001 100                       %')
disp('%                                             111                         %')
disp('%                                                                         %')
disp('%             ____                   _____                                %')
disp('%            / __ \                 |  __ \                               %')
disp('%           | |  | |_ __   ___ _ __ | |__) | __ ___  _ __                 %')
disp('%           | |  | | ''_ \ / _ \ ''_ \|  ___/ ''__/ _ \| ''_ \                %')
disp('%           | |__| | |_) |  __/ | | | |   | | | (_) | |_) |               %')
disp('%            \____/| .__/ \___|_| |_|_|   |_|  \___/| .__/                %')
disp('%                  | |                              | |                   %')
disp('%                  |_|                              |_|                   %')
disp('%                                                                         %')
disp('%             An integrated rotor design and analysis tool.               %')
disp('%                                                                         %')
disp('%                                                                         %')
disp('% Copyright (C) 2011, Brenden Epps.                                       %')
disp('%                                                                         %')
disp('% This program is free software; you can redistribute it and/or modify it %')
disp('% under the terms of the GNU General Public License as published by the   %')
disp('% Free Software Foundation.                                               %')
disp('%                                                                         %')
disp('% This program is distributed in the hope that it will be useful, but     %')
disp('% WITHOUT ANY WARRANTY; without even the implied warranty of              %')
disp('% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %')
disp('% See the GNU General Public License for more details.                    %')
disp('%                                                                         %')
disp('% ----------------------------------------------------------------------- %')

%将当前目录下的SourceCode文件夹添加到“搜索路径”（主页-设置路径 可查看当前的搜索路径）中
% ./ 表示当前目录下；../ 表示父级目录下
addpath  ./SourceCode

run('./SourceCode/OpenPropSingle.m')

