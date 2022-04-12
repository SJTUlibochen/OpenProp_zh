% 该函数用于在若干可能位置执行运行OpenProp和保存文件的功能：
% 包括的位置有：
% \UAV(21Summer)\
% \OpenProp_results\
% \OpenProp_results\filename\
% \OpenProp_results\anotherfile\或\OpenProp_results\filename\somefile\
% \OpenProp_zh\及所有的子文件内
% 对于其他位置，则提示运行OpenProp的正确位置
% 该函数运行完毕后，工作目录位于\OpenProp_results\filename\中
function existence = ChangeDirectory(OpenPropDirectory,filename)
    DircStrings = split(pwd,'\');
    DircLength = length(DircStrings);
    % [lia locb]数组中lia为逻辑值，表示是否存在；locb表示若存在则其索引值，不存在默认为0
    [lia1,locb1] = ismember({'OpenProp_results'},DircStrings);
    [lia2,locb2] = ismember({filename},DircStrings);
    [lia3,locb3] = ismember({OpenPropDirectory},DircStrings);
    if lia1
        if locb1 == DircLength
            % 位于\OpenProp_results\中
            [status,msg] = mkdir(['./',filename]);
            cd(['./',filename]);
        elseif lia2 && locb2 == DircLength
            % 位于\OpenProp_results\filename\中
            status = 1;
            msg = '目录已存在。';
        else
            % 位于\OpenProp_results\anotherfile\或\OpenProp_results\filename\somefile\中
            temp = [];
            for index = 1:locb1
                temp = strcat(temp,DircStrings{index},'\');
            end
            cd(temp);
            [status,msg] = mkdir(['./',filename]);
            cd(['./',filename]); 
        end
    elseif lia3
        % 位于\OpenPropDirectory\或其子文件中
        temp = [];
        for index = 1:locb3-1
            temp = strcat(temp,DircStrings{index},'\');
        end
        cd([temp,'OpenProp_results\']);
        [status,msg] = mkdir(['./',filename]);
        cd(['./',filename]);
    elseif exist('OpenProp_results','dir')
        % 位于\UAV(21Summer)\中
        cd('./OpenProp_results');
        [status,msg] = mkdir(['./',filename]);
        cd(['./',filename]);
    else
        % 位于其他错误位置中
        disp('路径错误，请从正确位置启动');
    end  
    if status && strcmp(msg,'目录已存在。') && exist([filename,'.mat'],'file')
        existence = 1;
    else
        existence = 0;
    end
end