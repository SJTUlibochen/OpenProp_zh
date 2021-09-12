% 工具函数，将直接从GUI读取的元胞数组转化为普通数组
function array = Cell2Array(cell)
    celllen = length(cell);
    array = zeros(celllen,1);
    for index = 1 : celllen
        array(index) = cell{index};
    end
end