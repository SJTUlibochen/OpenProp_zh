% �ú������ڽ���ǰλ���ƶ���/OpenProp_results/��
function LoadDirectory(rest,OpenPropDirectory)
    % ~��ͬ��not������~isempty(rest)��restΪ���Ƿ���0���ǿ�ʱ����1
    while ~isempty(rest)
        % strtok(rest,'\')������������'\'֮ǰ���ַ�����ֵ��CurrentDirectory��ʣ��Ĳ��ֱ�����rest��
        [CurrentDirectory,rest] = strtok(rest,'\');
        % strcmp�Ƚ��ַ�������ͬ�򷵻�1
        if strcmp(CurrentDirectory,OpenPropDirectory)
            if isempty(rest)
                % ��ǰλ��Ϊ/OpenPropDirectory/
                cd('../OpenProp_results')
            % rest(2:end)��ȥ�����׸��ַ����rest
            elseif strcmp(rest(2:end),'SourceCode')
                % ��ǰλ��Ϊ/OpenPropDirectory/SourceCode/
                cd('../../OpenProp_results')
                rest = [];
            else
                % ��ǰλ��Ϊ/OpenPropDirectory/wrongfolder
                disp('·������')
                return
            end
        elseif strcmp(CurrentDirectory,'OpenProp_results')
            if isempty(rest)
                % ��ǰλ��Ϊ/OpenProp_results/
                return
            else
                % ��ǰλ��Ϊ/OpenProp_results/otherfolder
                cd('../')
                rest = [];
            end
        end
    end
end