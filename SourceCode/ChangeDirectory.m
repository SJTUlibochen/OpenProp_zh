% �ú������������ɿ���λ��ִ������OpenProp�ͱ����ļ��Ĺ��ܣ�
% ������λ���У�
% \UAV(21Summer)\
% \OpenProp_results\
% \OpenProp_results\filename\
% \OpenProp_results\anotherfile\��\OpenProp_results\filename\somefile\
% \OpenProp_zh\�����е����ļ���
% ��������λ�ã�����ʾ����OpenProp����ȷλ��
% �ú���������Ϻ󣬹���Ŀ¼λ��\OpenProp_results\filename\��
function existence = ChangeDirectory(OpenPropDirectory,filename)
    DircStrings = split(pwd,'\');
    DircLength = length(DircStrings);
    % [lia locb]������liaΪ�߼�ֵ����ʾ�Ƿ���ڣ�locb��ʾ��������������ֵ��������Ĭ��Ϊ0
    [lia1,locb1] = ismember({'OpenProp_results'},DircStrings);
    [lia2,locb2] = ismember({filename},DircStrings);
    [lia3,locb3] = ismember({OpenPropDirectory},DircStrings);
    if lia1
        if locb1 == DircLength
            % λ��\OpenProp_results\��
            [status,msg] = mkdir(['./',filename]);
            cd(['./',filename]);
        elseif lia2 && locb2 == DircLength
            % λ��\OpenProp_results\filename\��
            status = 1;
            msg = 'Ŀ¼�Ѵ��ڡ�';
        else
            % λ��\OpenProp_results\anotherfile\��\OpenProp_results\filename\somefile\��
            temp = [];
            for index = 1:locb1
                temp = strcat(temp,DircStrings{index},'\');
            end
            cd(temp);
            [status,msg] = mkdir(['./',filename]);
            cd(['./',filename]); 
        end
    elseif lia3
        % λ��\OpenPropDirectory\�������ļ���
        temp = [];
        for index = 1:locb3-1
            temp = strcat(temp,DircStrings{index},'\');
        end
        cd([temp,'OpenProp_results\']);
        [status,msg] = mkdir(['./',filename]);
        cd(['./',filename]);
    elseif exist('OpenProp_results','dir')
        % λ��\UAV(21Summer)\��
        cd('./OpenProp_results');
        [status,msg] = mkdir(['./',filename]);
        cd(['./',filename]);
    else
        % λ����������λ����
        disp('·�����������ȷλ������');
    end  
    if status && strcmp(msg,'Ŀ¼�Ѵ��ڡ�') && exist([filename,'.mat'],'file')
        existence = 1;
    else
        existence = 0;
    end
end