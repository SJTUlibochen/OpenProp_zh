% Last Modified: 10/21/2011, Brenden Epps
%
% Lifting Line panel radii layout
%
% RC == radius of control points / propeller radius
% RV == radius of vortex  points / propeller radius
% DR == difference in vortex radii / propeller radius
% 
% ��ⷽʽ��Kerwin�Ŀ��ñʼ�P179������ͬ
% -----------------------------------------------------
function [rc,rv,drv] = LLPanelRadii(Mp,RhoRp,Hub_flag,Duct_flag)
    if Duct_flag == 0 && Hub_flag == 0
       % Constant spacing -- 1/4 panel inset at hub and tip  (use with no hub or duct image)
        rv = RhoRp + (1-RhoRp) * ((0:Mp)+0.25)/(Mp+0.50);
        rc = RhoRp + (1-RhoRp) * ((1:Mp)-0.25)/(Mp+0.50);
    elseif Duct_flag == 1 && Hub_flag == 0 
       % Constant spacing -- 1/4 panel inset at hub only  (use with no hub image but yes duct image)
        rv = RhoRp + (1-RhoRp) * ((0:Mp)+0.25)/(Mp+0.25);
        rc = RhoRp + (1-RhoRp) * ((1:Mp)-0.25)/(Mp+0.25);
    elseif Duct_flag == 1 && Hub_flag == 1    
        % Constant spacing -- no inset  (use with both duct image and hub image)
        rv = RhoRp + (1-RhoRp) * ((0:Mp)     )/Mp;
        rc = RhoRp + (1-RhoRp) * ((1:Mp)-0.50)/Mp;
    elseif Duct_flag == 0 && Hub_flag == 1
        % Constant spacing -- 1/4 panel inset at tip only  (use with hub image but no duct image)
        rv = RhoRp + (1-RhoRp) * ((0:Mp)     )/(Mp+0.25);
        rc = RhoRp + (1-RhoRp) * ((1:Mp)-0.50)/(Mp+0.25);
    end
    
    % diff�������Ա�ʾ��������������Ĳ�ֵ��DR���ȱ�RVС1
    drv = diff(rv); 
    
    % �������������Ϊ������
    rc = rc';
    rv = rv';
    drv = drv';
end

