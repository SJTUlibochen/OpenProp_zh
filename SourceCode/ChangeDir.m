%�ú������ڽ�����OpenProp��ִ�����ݱ���ʱ���ܵ����ִ�����ļ�λ���ƶ�����ȷλ��(/OpenPropDirectory/)
%���������Ĵ���λ�ã�����ʾ����OpenProp����ȷλ��
function ChangeDir(rest,OpenPropDirectory,filename)
    %~��ͬ��not������~isempty(rest)��restΪ���Ƿ���0���ǿ�ʱ����1
    while ~isempty(rest)
        %strtok(rest,'\')������������'\'֮ǰ���ַ�����ֵ��CurrentDirectory��ʣ��Ĳ��ֱ�����rest��
        [CurrentDirectory,rest] = strtok(rest,'\');
        %strcmp�Ƚ��ַ�������ͬ�򷵻�1
        if strcmp(CurrentDirectory,OpenPropDirectory)
            if isempty(rest)
                %you are in /OpenPropDirectory/
                addpath ./SourceCode
                cd('../OpenProp_results')
                mkdir(['./',filename])
                cd(['./',filename])
            %rest(2:end)��ȥ�����׸��ַ����rest
            elseif strcmp(rest(2:end),'SourceCode')
                %you are in /OpenPropDirectory/SourceCode/
                addpath ./
                cd('../../OpenProp_results')
                mkdir(['./',filename])
                cd(['./',filename])
                rest = [];
            else
                %you are in /OpenPropDirectory/wrongfolder
                %��ǰĿ¼����
                disp('·���������OpenProp_zh�ļ���������OpenProp')
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
                %��ǰĿ¼����
                disp('·���������OpenProp_zh�ļ���������OpenProp')
                return
            end
        end
    end
end