%���ذ�ť�Ƿ񱻰��µ��߼�ֵ�������鼰�����°�ť������
function [pushedlist,pushednum] = PushedFind(hObject,ED)
    %���Ը�ѡ��������£�
    %ҶƬ������������(��������)-ҶƬ������������(������ƽ��)    1-9
    %ҶƬ��ȷֲ�����(��������)-ҶƬ��ȷֲ�����(������ƽ��)    2-10
    %�����ٶȷֲ�����-�յ��ٶȷֲ�����      3-7
    %������еİ�ť��0�л���1ʱ������ڵİ�ťenable����Ҫ�䶯�����ఴť��Ϊ������
    global PlotsValues;
    buttonstatus = zeros(1,length(PlotsValues));
    for index = 1 : length(PlotsValues)
        buttonstatus(index) = get(PlotsValues(index),'value');       
    end
    %���ر����°�ť������ֵ����
    pushedlist = find(buttonstatus);
    pushednum = length(pushedlist);
end