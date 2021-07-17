%该函数用于将当前位置移动到/OpenProp_results/下
function LoadDir(rest,OpenPropDirectory)
    %~等同于not，这里~isempty(rest)当rest为空是返回0，非空时返回1
    while ~isempty(rest)
        %strtok(rest,'\')将最先遇到的'\'之前的字符串赋值给CurrentDirectory，剩余的部分保留于rest中
        [CurrentDirectory,rest] = strtok(rest,'\');
        %strcmp比较字符串，相同则返回1
        if strcmp(CurrentDirectory,OpenPropDirectory)
            if isempty(rest)
                %you are in /OpenPropDirectory/
                cd('../OpenProp_results')
            %rest(2:end)是去除了首个字符后的rest
            elseif strcmp(rest(2:end),'SourceCode')
                %you are in /OpenPropDirectory/SourceCode/
                cd('../../OpenProp_results')
                rest = [];
            else
                %you are in /OpenPropDirectory/wrongfolder
                %当前目录错误
                disp('路径错误！')
                return
            end
        elseif strcmp(CurrentDirectory,'OpenProp_results')
            if isempty(rest)
                %you are in /OpenProp_results/
                return
            else
                %you are in /OpenProp_results/otherfolder
                cd('../')
                rest = [];
            end
        end
    end
end