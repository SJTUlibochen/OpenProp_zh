%该函数用于将运行OpenProp或执行数据保存时可能的三种错误的文件位置移动至正确位置(/OpenPropDirectory/)
%对于其他的错误位置，则提示运行OpenProp的正确位置
function ChangeDir(rest,OpenPropDirectory,filename)
    %~等同于not，这里~isempty(rest)当rest为空是返回0，非空时返回1
    while ~isempty(rest)
        %strtok(rest,'\')将最先遇到的'\'之前的字符串赋值给CurrentDirectory，剩余的部分保留于rest中
        [CurrentDirectory,rest] = strtok(rest,'\');
        %strcmp比较字符串，相同则返回1
        if strcmp(CurrentDirectory,OpenPropDirectory)
            if isempty(rest)
                %you are in /OpenPropDirectory/
                addpath ./SourceCode
                cd('../OpenProp_results')
                mkdir(['./',filename])
                cd(['./',filename])
            %rest(2:end)是去除了首个字符后的rest
            elseif strcmp(rest(2:end),'SourceCode')
                %you are in /OpenPropDirectory/SourceCode/
                addpath ./
                cd('../../OpenProp_results')
                mkdir(['./',filename])
                cd(['./',filename])
                rest = [];
            else
                %you are in /OpenPropDirectory/wrongfolder
                %当前目录错误
                disp('路径错误！请从OpenProp_zh文件夹内启动OpenProp')
                return
            end
        elseif strcmp(CurrentDirectory,'OpenProp_results')
            if isempty(rest)
                %you are in /OpenProp_results/
                addpath(['../',OpenPropDirectory,'/SourcecCode'])
                mkdir(['./',filename])
                cd(['./',filename])
            elseif strcmp(rest(2:end),filename)
                %you are in /OpenProp_results/filename/
                addpath(['../../',OpenPropDirectory,'/SourcecCode'])
                rest = [];
            else
                %you are in /OpenProp_results/wrongfolder
                %当前目录错误
                disp('路径错误！请从OpenProp_zh文件夹内启动OpenProp')
                return
            end
        end
    end
end