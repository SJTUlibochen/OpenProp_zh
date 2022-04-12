%返回按钮是否被按下的逻辑值索引数组及被按下按钮的数量
function [pushedlist,pushednum] = PushedFind(hObject,ED)
    %可以复选的组合如下：
    %叶片伸张轮廓曲线(来自输入)-叶片伸张轮廓曲线(来自设计结果)    1-9
    %叶片厚度分布曲线(来自输入)-叶片厚度分布曲线(来自设计结果)    2-10
    %流入速度分布曲线-诱导速度分布曲线      3-7
    %当组合中的按钮由0切换至1时，组合内的按钮enable不需要变动，其余按钮变为不可用
    global PlotsValues;
    buttonstatus = zeros(1,length(PlotsValues));
    for index = 1 : length(PlotsValues)
        buttonstatus(index) = get(PlotsValues(index),'value');       
    end
    %返回被按下按钮的索引值数组
    pushedlist = find(buttonstatus);
    pushednum = length(pushedlist);
end