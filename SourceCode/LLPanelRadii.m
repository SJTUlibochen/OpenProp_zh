% Last Modified: 10/21/2011, Brenden Epps
%
% Lifting Line panel radii layout
%
% RC == radius of control points / propeller radius
% RV == radius of vortex  points / propeller radius
% DR == difference in vortex radii / propeller radius
%求解方式与Kerwin的课堂笔记P179并不相同
function [RC,RV,DR] = LLPanelRadii(Mp,Rhub_oR,Hub_flag,Duct_flag)
    if Duct_flag == 0 && Hub_flag == 0
       % Constant spacing -- 1/4 panel inset at hub and tip  (use with no hub or duct image)
        RV = Rhub_oR + (1-Rhub_oR) * ((0:Mp)+0.25)/(Mp+0.50);
        RC = Rhub_oR + (1-Rhub_oR) * ((1:Mp)-0.25)/(Mp+0.50);
    elseif Duct_flag == 1 && Hub_flag == 0 
       % Constant spacing -- 1/4 panel inset at hub only  (use with no hub image but yes duct image)
        RV = Rhub_oR + (1-Rhub_oR) * ((0:Mp)+0.25)/(Mp+0.25);
        RC = Rhub_oR + (1-Rhub_oR) * ((1:Mp)-0.25)/(Mp+0.25);
    elseif Duct_flag == 1 && Hub_flag == 1    
        % Constant spacing -- no inset  (use with both duct image and hub image)
        RV = Rhub_oR + (1-Rhub_oR) * ((0:Mp)     )/Mp;
        RC = Rhub_oR + (1-Rhub_oR) * ((1:Mp)-0.50)/Mp;
    elseif Duct_flag == 0 && Hub_flag == 1
        % Constant spacing -- 1/4 panel inset at tip only  (use with hub image but no duct image)
        RV = Rhub_oR + (1-Rhub_oR) * ((0:Mp)     )/(Mp+0.25);
        RC = Rhub_oR + (1-Rhub_oR) * ((1:Mp)-0.50)/(Mp+0.25);
    end
    %diff函数可以表示数组中相邻两项的插值，DR长度比RV小1
    DR = diff(RV);  
end

