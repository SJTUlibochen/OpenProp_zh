%�ú������ڽ���ǰλ���ƶ���/OpenProp_results/��
function LoadDir(rest,OpenPropDirectory)
    %~��ͬ��not������~isempty(rest)��restΪ���Ƿ���0���ǿ�ʱ����1
    while ~isempty(rest)
        %strtok(rest,'\')������������'\'֮ǰ���ַ�����ֵ��CurrentDirectory��ʣ��Ĳ��ֱ�����rest��
        [CurrentDirectory,rest] = strtok(rest,'\');
        %strcmp�Ƚ��ַ�������ͬ�򷵻�1
        if strcmp(CurrentDirectory,OpenPropDirectory)
            if isempty(rest)
                %you are in /OpenPropDirectory/
                cd('../OpenProp_results')
            %rest(2:end)��ȥ�����׸��ַ����rest
            elseif strcmp(rest(2:end),'SourceCode')
                %you are in /OpenPropDirectory/SourceCode/
                cd('../../OpenProp_results')
                rest = [];
            else
                %you are in /OpenPropDirectory/wrongfolder
                %��ǰĿ¼����
                disp('·������')
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