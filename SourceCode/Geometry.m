% This function determines the geometry of the propeller.
% It outputs the geometry as a 2D image, 3D image, and stl file.
%
% 该函数决定螺旋桨的几何结构，其输出形式为2D图片、3D实时渲染或者快速成型文件
% 
% Reference: 
%   [1]J.S. Carlton, "Marine Propellers & Propulsion", ch. 3, 1994.
%   [2]Abbott, I. H., and Von Doenhoff, A. E.; Theory of Wing Sections. 
%      Dover, 1959. 
% 
% -------------------------------------------------------------------------
% Input Variables: 
% 输入变量：
%
%   filename            file name prefix for all output files
%                       文件名，该文件名将用于对所有输出文件的命名
%   date                time and date to print on reports
%                       输出日期
%   Geometry_flag       flag for whether to make 2D/3D geometry plot
%                       决定是否生成2/3D演示模型
%   Coordinate_flag     flag for whether to make points output file
%                       决定是否输出点坐标文档
%   Meanline_flag       flag for choice of meanline  form
%                       决定凸缘线类型是否沿径向同一
%   Thickness_flag      flag for choice of thickness form
%                       决定厚度类型是否沿径向同一
%   Meanline_x          types of meanline along rx
%                       沿径向分布的凸缘线类型
%   Thickness_x         types of thickness along rx
%                       沿径向分布的厚度类型
%   rx                  input radii / propeller radius
%                       径向坐标
%   f0oc0               input camber / chord at each radius
%   t0oc0               input thickness / chord at each radius
%   
%   rc                  control point radii / propeller radius
%                       控制点坐标
%   CL                  section lift coefficients
%   BetaC               Beta  at the control points
%   BetaIC              BetaI at the control points
%   alphaI              ideal angle of attack
%
%   Dp                  propeller diameter
%   Z                   number of blades
%   N                   propeller speed
%   Dhub                hub diameter
%   Rhub                hub radius
%
%   CpoDp               chord / diameter at each control point radius
%   Rp                  propeller radius
%   Mp                  number of radial 2D cross-sections
%   Np                  number of points in each 2D section
%   Js                  advance coefficient based on ship speed
%
% Output Variables:
%
%   The function has graphical and file outputs, in addition to the 
%   geometry data structure.
%
% -------------------------------------------------------------------------

function [geometry] = Geometry(pt,rg)
    global Fig_Main PlotsValues;
    if nargin == 1
        % 函数只有pt一个输入参数
        RadiiGiven_flag = 0;
    else 
        RadiiGiven_flag = 1; 
    end
    %% 
    input = pt.input;
    design = pt.design;
    
    filename = pt.filename;
    
    % 从pt.input中导出数据
    Z = input.Z;
    Dp = input.Dp;
    Hub_flag = input.Hub_flag;
    Dh = input.Dh;
    Duct_flag = input.Duct_flag;
    Nx = input.Nx;
    rx = input.rx;
    Meanline_x = input.Meanline_x;
    Thickness_x = input.Thickness_x;
    Geometry_flag = input.Geometry_flag;
    Coordinate_flag = input.Coordinate_flag;
    Printing_flag = input.Printing_flag;
    Mp = input.Mp;
    Np = input.Np;
    
    % 从pt.design中导出数据
    CpoDp_c = design.CpoDp_c;  
    T0oDp_c = design.T0oDp_c;
    if Duct_flag == 1 
        RdoRp = design.RdoRp;
        CdoRp = design.CdoRp;
        XdoRp = design.XdoRp;
        Gd = design.Gd;
        VARING = design.VARING;
        % 涵道切面凸缘线类型
        Meanline_d = 'NACA a=0.8 (modified)';
        % 涵道切面厚度类型
        Thickness_d = 'NACA 65A010';
    else
        RdoRp = 1;
        CdoRp = 1;
        XdoRp = 0;
        Gd = 0;
        VARING = 0;
    end
    rc = design.rc;
    CL_c = design.CL_c;
    tanBetaI_c = design.tanBetaI_c;
    BetaI_c = atand(tanBetaI_c);
    Js = design.Js;
    
    % 引申数据
    RhoRp = Dh/Dp;
    
    % 其他参数
    % 升力线位于中心则为0，位于1/4弦长处则为1
    QuarterChord_flag = 0;
    % 升力面几何修正
    LSGeoCorr = 'none';
    % 左手螺旋桨时为1
    LeftHand_flag = 0;
    
    %%
    % 确定径向坐标
    if RadiiGiven_flag == 0
        % 没有额外提供用于插值的坐标
        rg = RhoRp+(1-RhoRp)*(sin((0:Mp)*pi/(2*Mp)));
        % 将插值坐标调整为列向量
        rg = rg';
    else
        Mp = length(rg)-1;
        if (rg(1)<RhoRp)||(rg(end)>1)
            message = sprintf('错误：径向位置输入值超出范围');
            uialert(Fig_Main,message,'输入错误','icon','error');
            return 
        end
    end
    
    % 对参数进行插值
    CL_g = pchip(rc,CL_c,rg); 
    BetaI_g = pchip(rc,BetaI_c,rg); 
    tanBetaI_g = pchip(rc,tanBetaI_c,rg);
    if (Duct_flag == 0)||(RdoRp>1.001)
        % 叶梢弦长为0
        CpoDp_g = InterpolateChord(rc,CpoDp_c,rg);   
    else
        % 叶梢弦长为有限值
        CpoDp_g = pchip(rc,CpoDp_c,rg);
    end
    T0oDp_g = pchip(rc,T0oDp_c,rg);

    % 将归一化的几何参数复原
    R_g = rg*Dp/2;
    C_g = CpoDp_g*Dp;
    T0_g = T0oDp_g*Dp;
    
    % 计算伸张面积比
    EAR = (2*Z/pi)*trapz(linspace(RhoRp,1,100),...
          interp1(rg,CpoDp_g,linspace(RhoRp,1,100),'spline','extrap'));  
    % 计算叶芯厚径比
    BTF = interp1(rg,T0oDp_g,0,'linear','extrap');
    
    % 设置升力面坐标系
    % x0    x/C distance along mid-chord line to interpolate geometry data.
    % xy    xy  distance along mid-chord line to interpolate geometry data.
    x = zeros(1,Np);
    C_gx = zeros(Mp+1,Np);

    for j = 1 : Np
        % 坐标点沿弦长方向为余弦分布，方向由前缘线向后缘线
        x(j) = 0.5*(1-cos(pi*(j-1)/(Np-1)));
    end
    for i = 1 : Mp+1
        if QuarterChord_flag == 1
            C_gx(i,:) = C_g(i)/4-C_g(i)*x;
        else
            % 前缘线处为C_g(i)/2，后缘线处为-C_g(i)/2
            C_gx(i,:) = C_g(i)/2-C_g(i)*x;
        end
    end
    
    % 获取切面凸缘线和厚度相关数据
    F0oCpdtilde_x = zeros(Nx,1);
    CLItilde_x = zeros(Nx,1);
    AlphaItilde_x = zeros(Nx,1);
    FoF0_x = zeros(Nx,Np);
    dFoF0dx_x = zeros(Nx,Np);
    ToT0_x = zeros(Nx,Np);
    FoF0_g = zeros(Mp+1,Np);
    dFoF0dx_g = zeros(Mp+1,Np);
    ToT0_g = zeros(Mp+1,Np);

    % 沿径向方向写入几何信息
    for i = 1 : Nx
        [F0oCpdtilde_x(i),CLItilde_x(i),AlphaItilde_x(i),FoF0_x(i,:),...
         dFoF0dx_x(i,:),ToT0_x(i,:)] = GeometryFoil2D(Meanline_x{i},Thickness_x{i},x);
    end

    % 将径向坐标修改为rg
    F0oCptilde_g = pchip(rx,F0oCpdtilde_x,rg);
    CLItilde_g = pchip(rx,CLItilde_x,rg);
    AlphaItilde_g = pchip(rx,AlphaItilde_x,rg);
    for j = 1 : Np
        FoF0_g(:,j) = pchip(rx,FoF0_x(:,j),rg);
        dFoF0dx_g(:,j) = pchip(rx,dFoF0dx_x(:,j),rg);
        ToT0_g(:,j) = pchip(rx,ToT0_x(:,j),rg);
    end
    
    % 根据tilde值确定实际值
    AlphaI_g = AlphaItilde_g.*CL_g./CLItilde_g;
    F0oCp_g = F0oCptilde_g.*CL_g./CLItilde_g;
    
    % 基于三维升力面相关工程经验对攻角和拱度进行修正
    if strcmp(LSGeoCorr,'none') 
        % 没有选择修正方法，不进行修正
    elseif strcmp(LSGeoCorr,'Morgan1968')
        % 基于Morgan的方法进行修正
        [Kc,Ka,Kt] = Morgan1968(rg,tanBetaI_g,EAR,Z);
        for m = 1 : Mp+1
            F0oCp_g(m) = Kc(m) * F0oCp_g(m);
            AlphaI_g(m) = (Ka(m) * (pi/180)*AlphaI_g(m)+Kt(m)*BTF)*(180/pi);
        end        
        
    elseif strcmp(LSGeoCorr,'EckhardtMorgan1955')
        % 基于Eckhard的方法进行修正
        [K1K2] = EckhardtMorgan1955(EAR,rg,tanBetaI_g);
        for m = 1 : Mp+1
            F0oCp_g(m) = K1K2(m) * F0oCp_g(m);
            AlphaI_g(m) = K1K2(m) * AlphaI_g(m);
        end     
    end
    
    % 求解每个坐标点处实际的厚度、拱度和拱度的导数值
    T_gx = zeros(Mp+1,Np);
    F_gx = zeros(Mp+1,Np);
    dFdx_gx = zeros(Mp+1,Np);
    for i = 1 : Mp+1
        % 注意F0oCd_g(i)和C_g(i)均为具体数值，不存在矩阵相乘的问题
        F_gx(i,:) = FoF0_g(i,:)*F0oCp_g(i)*C_g(i);
        dFdx_gx(i,:) = dFoF0dx_g(i,:)*F0oCp_g(i);
        T_gx(i,:) = ToT0_g(i,:)*T0_g(i);
    end      
    
    % 螺距角为水动力螺旋角和攻角之和
    Phi_g = BetaI_g+AlphaI_g;
    % 螺距/叶片直径
    PoDp = tand(Phi_g).*pi.*rg;
    % 叶片之间的角度
    theta_Z  = 0:360/Z:360;
    T0oCp_g = T0_g./C_g;
    
    % 叶切面模型上表面的横坐标
    x2D_u = zeros(Mp+1,Np);
    % 叶切面模型下表面的横坐标
    x2D_l = zeros(Mp+1,Np);
    % 叶切面模型上表面的纵坐标
    y2D_u = zeros(Mp+1,Np);
    % 叶切面模型下表面的纵坐标
    y2D_l = zeros(Mp+1,Np);
    % 计算2D叶切面模型上的点坐标
    for i = 1 : Mp+1
        % 径向坐标方向：由叶根向叶梢
        for j = 1:Np
            % 弦向坐标方向：由前缘线向后缘线
            x2D_u(i,j) = C_gx(i,j)+(T_gx(i,j)/2)*sin(atan(dFdx_gx(i,j)));
            x2D_l(i,j) = C_gx(i,j)-(T_gx(i,j)/2)*sin(atan(dFdx_gx(i,j)));
            y2D_u(i,j) = F_gx(i,j)+(T_gx(i,j)/2)*cos(atan(dFdx_gx(i,j)));
            y2D_l(i,j) = F_gx(i,j)-(T_gx(i,j)/2)*cos(atan(dFdx_gx(i,j)));
        end
    end
    
    % 将上下表面的横坐标和纵坐标分别合成至同一数组中，顺序：
    % 径向(i:1 - Mp+1)：叶根 - 叶梢
    % 弦向(j:1 - Np)：后缘线 - 上表面（吸力面） - 前缘线 - 下表面（压力面）
    x2D(:,1:Np) = x2D_u(:,Np:-1:1);
    x2D(:,Np+1:2*Np) = x2D_l(:,1:Np);
    y2D(:,1:Np) = y2D_u(:,Np:-1:1);
    y2D(:,Np+1:2*Np) = y2D_l(:,1:Np);
    
    % 计算加入螺距角后的叶切面2D点坐标
    x2Dr = zeros(Mp+1,2*Np);
    y2Dr = zeros(Mp+1,2*Np);
    for i = 1 : Mp+1
        % 沿径向位置螺距角变化
        x2Dr(i,:) = x2D(i,:)*cosd(Phi_g(i))-y2D(i,:)*sind(Phi_g(i));
        y2Dr(i,:) = x2D(i,:)*sind(Phi_g(i))+y2D(i,:)*cosd(Phi_g(i));
    end

    % 计算三维点坐标
    % X3D [m], X position in 3D space (corresponds to y position in 2D space)
    % Y2D [m], Y position in 3D space
    % Z3D [m], Z position in 3D space
    x3D = zeros(Mp+1,2*Np,Z);
    y3D = zeros(Mp+1,2*Np,Z);
    z3D = zeros(Mp+1,2*Np,Z);
    for i = 1:Mp+1
        for j = 1:2*Np
            for k = 1:Z
                x3D(i,j,k) = y2Dr(i,j);
                y3D(i,j,k) = R_g(i)*sind(-(180/pi)*x2Dr(i,j)/R_g(i)-theta_Z(k));
                z3D(i,j,k) = R_g(i)*cosd(-(180/pi)*x2Dr(i,j)/R_g(i)-theta_Z(k));
            end
        end
    end
    
    % 左手螺旋桨为右手螺旋将关于xoz平面对称
    if LeftHand_flag == 1
        y3D = -y3D;
    end
    
    %% 生成2/3D演示模型
    if Geometry_flag
        set(PlotsValues(11),'valuechangedfcn',{@Geometry2Dfcn,...
                                               Mp,Np,rg,x2Dr,y2Dr});
        set(PlotsValues(12),'valuechangedfcn',{@Geometry3Dfcn,Z,...
                                               x3D,y3D,z3D});
    end    
    
%     Make_3D_Blade_Image(x3D,y3D,z3D,Duct_flag,RdoRp,CdoRp,XdoRp,Hub_flag,...
%                         RhoRp,Js,BetaI_g,theta,TANBIV,RV,Dp/2,Geometry_flag,Plots,PlotPanels);
%     if Duct_flag == 1
%         [xd,rd,Xd,Yd,Zd] = Duct_Plot_120329(rc,TANBIC,G,VARING,RV,Z,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR,Cduct_oR,Xduct_oR,Gd,Meanline_d,Thickness_d,x0,t0oc_duct,Rp,Mp,Np);
%         surf(Xd,Yd,Zd);
%     end
%     
%     if Coordinate_flag
%         ExportPointCoordinates(filename,Np,Mp,x3D,y3D,z3D);
%     end
    
    %% 保存数据
    % 将几何结构数据保存至pt.geometry中
    geometry.part1 = '设计规格';
    geometry.Z = Z;                             % 叶片个数
    geometry.Dp = Dp;                           % 叶片直径
    geometry.Dh = Dh;                           % 桨毂直径
    geometry.Duct_flag = Duct_flag;             % 设计中加入涵道
    geometry.EAR = EAR;                         % 伸张面积比
    geometry.BTF = BTF;                         % 叶芯厚径比
    if Duct_flag
        geometry.RdoRp = RdoRp;                 % 涵道直径
        geometry.CdoRp = CdoRp;                 % 涵道弦长
        geometry.XdoRp = XdoRp;                 % 涵道位置
        geometry.Meanline_d = Meanline_d;       % 涵道凸缘线
        geometry.Thickness_d = Thickness_d;     % 涵道厚度
        geometry.Gd = Gd;                       % 涵道环量
        geometry.VARING = VARING;               % 涵道中线轴向来流速度
    end
    
    geometry.part2 = '叶切面几何信息';
    geometry.Nx = Nx;                           % 切面数量
    geometry.rx = rx;                           % 径向位置
    geometry.Meanline_x = Meanline_x;           % 凸缘线类型
    geometry.Thickness_x = Thickness_x;         % 厚度类型
    geometry.rc = rc;                           % 径向位置
    geometry.CL_c = CL_c;                       % 升力系数
    geometry.BetaI_c = BetaI_c;                 % 水动力螺旋角
    geometry.rg = rg;                           % 径向位置
    geometry.CpoDp_g = CpoDp_g;                 % 叶片弦径比
    geometry.T0oDp_g = T0oDp_g;                 % 叶片厚径比
    geometry.T0oCp_g = T0oCp_g;                 % 叶片厚弦比
    geometry.F0oCp_g = F0oCp_g;                 % 叶片拱弦比
    geometry.BetaI_g = BetaI_g;                 % 水动力螺旋角
    geometry.AlphaI_g = AlphaI_g;               % 攻角
    geometry.Phi_g = Phi_g;                     % 桨距角
    geometry.PoDp = PoDp;                       % 叶片螺距
    
    geometry.part3 = '其他参数';
    geometry.Geometry_flag = Geometry_flag;     % 生成2/3D演示模型
    geometry.Coordinate_flag = Coordinate_flag; % 输出点坐标文档
    geometry.Printing_flag = Printing_flag;     % 生成快速成型文件
    
    geometry.part4 = '点坐标';
    geometry.x2Dr = x2Dr;                       % 2D演示模型横坐标
    geometry.y2Dr = y2Dr;                       % 2D演示模型纵坐标
    geometry.x3D = x3D;                         % 3D演示模型x坐标
    geometry.y3D = y3D;                         % 3D演示模型y坐标
    geometry.z3D = z3D;                         % 3D演示模型z坐标
    
    
end

% “叶切面二维轮廓”的回调函数
function Geometry2Dfcn(hObject,~,Mp,Np,rg,x2Dr,y2Dr)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    
    % 设置图像面板标题
    set(PlotsPanel,'title','叶切面二维轮廓');
    
    % 查找被按下的按钮数量
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % 对按钮11进行操作后，若一个按钮也没有被按下，则所有按钮均可用
        set(PlotsValues,'enable','on');
        if ~strcmp(pt.input.ChordMethod,'CLmax')
            set(PlotsValues(1),'enable','off');
            if ~strcmp(pt.input.ChordMethod,'ConeyPLL')
                set(PlotsValues(2),'enable','off');
            end    
        end    
        if pt.input.Analyze_flag == 0
            set(PlotsValues(10),'enable','off');
        end 
    else
        % 对按钮11进行操作后，若存在被按下的按钮，则除11外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(11),'enable','on');
    end
    
    % 绘制图像
    if get(hObject,'value')
        % 避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        color = {'r','g','b','m','k'};
        cla(coordinate);
        j = 1;
        for index = 1 : ceil(Mp/3) : Mp
            Contour_line(j) = plot(coordinate,x2Dr(index,:),y2Dr(index,:),...
                                  'color',color{j},'linewidth',2);
            hold(coordinate,'on');
            plot(coordinate,...
                 [0.5*(x2Dr(index,1)+x2Dr(index,2*Np)),...
                  0.5*(x2Dr(index,Np)+x2Dr(index,Np+1))],...
                 [0.5*(y2Dr(index,1)+y2Dr(index,2*Np)),...
                  0.5*(y2Dr(index,Np)+y2Dr(index,Np+1))],...
                 'color',color{j},'linewidth',1);
            hold(coordinate,'on');
            LegendStrings{j} = ['rg = ',num2str(rg(index))];
            j = j+1;
        end
        grid(coordinate,'on');
        box(coordinate,'on');
        legend(Contour_line,LegendStrings);
        xlabel(coordinate,'x(2D)','fontsize',16,'fontname','Times');
        ylabel(coordinate,'y(2D)','fontsize',16,'fontname','Times');
    else
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');
        legend(coordinate,'off');
        set(PlotsPanel,'title','点击左侧按钮显示对应图像');
    end
end

% “螺旋桨三维模型”的回调函数
function Geometry3Dfcn(hObject,~,Z,x3D,y3D,z3D)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    
    % 设置图像面板标题
    set(PlotsPanel,'title','螺旋桨三维模型');
    
    % 查找被按下的按钮数量
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % 对按钮12进行操作后，若一个按钮也没有被按下，则所有按钮均可用
        set(PlotsValues,'enable','on');
        if ~strcmp(pt.input.ChordMethod,'CLmax')
            set(PlotsValues(1),'enable','off');
            if ~strcmp(pt.input.ChordMethod,'ConeyPLL')
                set(PlotsValues(2),'enable','off');
            end    
        end    
        if pt.input.Analyze_flag == 0
            set(PlotsValues(10),'enable','off');
        end 
    else
        % 对按钮12进行操作后，若存在被按下的按钮，则除12外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(12),'enable','on');
    end
    
    % 绘制图像
    if get(hObject,'value')
        % 避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        color = {'r','g','b','m','k'};
        cla(coordinate);
        for index = 1:Z
            Propeller_surf = surf(coordinate,x3D(:,:,1),y3D(:,:,index),z3D(:,:,index));
            hold(coordinate,'on');
        end
        
    else
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');
        legend(coordinate,'off');
        set(PlotsPanel,'title','点击左侧按钮显示对应图像');
    end    

end