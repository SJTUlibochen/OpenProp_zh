% ���ߺ�������ֱ�Ӵ�GUI��ȡ��Ԫ������ת��Ϊ��ͨ����
function array = Cell2Array(cell)
    celllen = length(cell);
    array = zeros(celllen,1);
    for index = 1 : celllen
        array(index) = cell{index};
    end
end