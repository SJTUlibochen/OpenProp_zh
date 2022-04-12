% Determines the "optimum" circulation, chord, and thickness distributions
% that satisfy the input operating conditions, using a variational 
% optimization algorithm for the propelller case (Coney, 1989), or a hybrid
% blade element momentum / vortex lattice method for the turbine case
% (Epps and Kimball, 2011).
%
% Returns: "design" data structure with performance specs such as 
%         thrust coefficient and efficiency, as well as the optimized
%         circulation distribution, chord distribution, etc.
%
% This implementation provides:
%   -- Propeller optimization via vortex-lattice theory, solved by...
%           ... linearized system of equations          (EppsOptimizer02.m) "LL-Linear"
%           ... Newton solver                           (EppsOptimizer23.m) "LL-Newton"
%           ... Newton solver with hub drag variation   (EppsOptimizer53.m) 
%   -- Turbine optimization via hybrid vortex lattice / momentum theory, solved by...
%           ... Newton solver                           (EppsOptimizer06.m) "Robust"
%   -- Chord length optimization via:
%           ... Specified maximum lift coefficient (default)
%           ... Cavitation mitigation method from (Epps, FAST2011)
%   -- Improved wake model for theoretical accuracy and numerical stability.
%   -- Improved duct model implementation for both propeller and turbine cases.
%
function [design] = EppsOptimizer(input)
    global Fig_Main;
    %% 从pt中导出输入值
    % part1 “螺旋桨规格”
    N = input.N;
    VS = input.VS;
    T = input.T;
    Z = input.Z;
    Dp = input.Dp;
    Hub_flag = input.Hub_flag;
    Dh = input.Dh;
    Duct_flag = input.Duct_flag;
    Dd = input.Dd;
    Cd = input.Cd;
    Ld = input.Ld;
    
    % part2 “叶片切面参数”
    Nx = input.Nx;
    ChordMethod = input.ChordMethod;
    Meanline_flag = input.Meanline_flag;
    Thickness_flag = input.Thickness_flag;
    rx = input.rx;
    CDp_x = input.CDp_x;
    Meanline_x = input.Meanline_x;
    Thickness_x = input.Thickness_x;   
    CL_x = input.CL_x;
    T0oDp_x = input.T0oDp_x;
    
    % part3 “外部环境参数”
    Ni = input.Ni;
    rho = input.rho;
    ri = input.ri;
    VA_i = input.VA_i;
    VT_i = input.VT_i;
    % 通过拟合使“流入速度/装置行进速度”数据更加光滑
    VA_i = RepairSpline(ri,VA_i);  
    VT_i = RepairSpline(ri,VT_i);
    
    % part4 “涵道相关参数”
    TdoT = input.TdoT;
    CDd = input.CDd;
    
    % part5 “其他参数与工具”
    Analyze_flag = input.Analyze_flag;
    Geometry_flag = input.Geometry_flag;
    Coordinate_flag = input.Coordinate_flag;
    Printing_flag = input.Printing_flag;
    Mp = input.Mp;
    Np = input.Np;
    Nd = input.Nd;
    ITER = input.ITER;
    
    % 其余参数
    HUF = 0;
    TUF = 0;
    Rhv = 0.5;
    H = 3.048;
    g = 9.81;
    Patm = 101325;
    Pv = 2500;
    rx_def = [0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;0.95;1];
    CpoDp_x_def = [0.1600;0.1812;0.2024;0.2196;0.2305;...
                   0.2311;0.2173;0.1807;0.1388;0.0010];
    % 输入中没有指定扭矩，若不指定升力而指定扭矩则为1
    TorqueSpec_flag = 0;
    % 指定展开面积比(EAR)
    EARspec = 0;
    
    %% 引申参数
    RhoRp = Dh/Dp;
    
    RdoRp = Dd/Dp;
    CdoRp = 2*Cd/Dp;
    % XdoRp反映涵道轴向中点处距离原点的距离
    XdoRp = 2*(Ld-Cd/2)/Dp;
    
    TpoT = 1-TdoT;
    
    SIGMAs = (Patm+rho*g*H-Pv)/(0.5*rho*VS^2);
    
    Js = VS/((N/60)*Dp);
    L = pi/Js;
    KT = T/(rho*(N/60)^2*Dp^4);
    CT = T/(0.5*rho*VS^2*pi*(Dp/2)^2);
    
    % CT_desired为所要求的总推力系数，CT_desired = CTp_desired+CTd_desired
    CT_desired = CT;
    % 叶片所提供的推力系数
    CTp_desired = CT_desired*TpoT;
    % 涵道所提供的推力系数
    CTd_desired = CT_desired*TdoT;
    
    % 获取切面拱度和理想升力系数的数据
    F0oCptilde_x = zeros(Nx,1);
    CLItilde_x = zeros(Nx,1);
    for index = 1 : Nx
        [F0oCptilde_x(index),CLItilde_x(index)] = GeometryFoil2D(Meanline_x{index},...
                                                                Thickness_x{index});
    end
    
    % 生成弦径比的插值数据
    CpoDp_x = InterpolateChord(rx_def,CpoDp_x_def,rx);
    
    %% 计算控制点和涡格点的径向位置并进行插值
    % 生成[RhoRp,1]之间100个等距点组成的向量
    ri_expanded = linspace(RhoRp,1,100);
    VA_i_expanded = pchip(ri,VA_i,ri_expanded);
    % 计算平均速度，可以参考Kerwin讲义p.138的eqn 163
    VMIV = 2/(1^2-RhoRp^2)*trapz(ri_expanded,ri_expanded.*VA_i_expanded);
    
    % 计算涡格点和控制点的径向位置 rc:[Mp*1];rv:[Mp+1*1]
    [rc,rv,drv] = LLPanelRadii(Mp,RhoRp,Hub_flag,Duct_flag);
    
    % 在相应径向位置进行插值
    VA_c = pchip(ri,VA_i,rc);
    VT_c = pchip(ri,VT_i,rc);
    CDp_c = pchip(rx,CDp_x,rc);
    T0oDp_c = pchip(rx,T0oDp_x,rc);
    CL_c = pchip(rx,CL_x,rc);
    F0oCdtilde_c = pchip(rx,F0oCptilde_x,rc);
    CLItilde_c = pchip(rx,CLItilde_x,rc);
    
    if (abs(rx(end)-1) < 1e-4) && (CpoDp_x(end) <= 0.01)
        % 叶梢处要求厚度和弦长为0时采用下列方法
        CpoDp_c = InterpolateChord(rx,CpoDp_x,rc);
    else
        CpoDp_c = pchip(rx,CpoDp_x,rc);
    end
    T0oCp_c = T0oDp_c./CpoDp_c;
    
    %% 设置迭代初始值
    % 设置归一化环量初始值 G = Gamma/2*pi*Rp*VS
    Gp = zeros(Mp,1);
    
    % 设置轴向和切向诱导速度初始值 Kerwin & Hadler (2010), eqn (4.26)
    UASTAR = 0*rc+0.5*(sqrt(1+CTp_desired)-1);
    UTSTAR = 0*rc;
    % 周围流场低流速情况 (Epps, 2013)
    if VMIV < 0.05
        UASTAR = 0*rc + 0.5*sqrt(CTp_desired);
    end
    
    % 根据诱导速度初始值计算进角和水动力螺旋角初始值的正切值
    tanBeta_c = VA_c./(L*rc+VT_c);
    tanBetaI_c = (VA_c+UASTAR)./(L*rc+VT_c+UTSTAR);
    
    % 计算总来流速度的初始值
    VSTAR = sqrt((VA_c+UASTAR).^2+(L*rc+VT_c+UTSTAR).^2);
    
    % 设置弦径比的初始值
    CpoDp_c = 0.1*ones(Mp,1);
    
    % 初始化叶片上的诱导速度
    [UAHIF,UTHIF] = Horseshoe(Mp,Z,tanBetaI_c,rc,rv,Hub_flag,...
                              RhoRp,Duct_flag,RdoRp);
    
    % 初始化涵道相关参数
    if Duct_flag == 1
        % 进行涵道对涵道内流场产生诱导速度的初始化
        [XdRING,GdRING,UADIF] = Duct_Influence(RdoRp,CdoRp,XdoRp,rc,Nd);
        
        % 基于pchip插值的方法确定涵道中线处的轴向来流速度
        VARING = pchip(ri,VA_i,RdoRp);  

        % 设置归一化涵道总环量的初始值
        Gd = 0;
        
        % 涵道内流场由涡环产生的总轴向诱导速度
        UADUCT = UADIF * Gd;
        
        % 进行叶片对涵道产生诱导速度的初始化
        [DAHIF_times_TANBIC,~,DRHIF_times_TANBIC] = Horseshoe_intr_110830(XdRING,RdoRp, rc,ones(size(tanBetaI_c)),rv,Z,Hub_flag,RhoRp,Duct_flag,RdoRp); 
        % DAHIF(n,m)单位强度环量下叶片上第m个涡格对涵道上第n个控制点的轴向诱导速度
        DAHIF = zeros(Nd,Mp);
        DRHIF = zeros(Nd,Mp);
        for m = 1 : Mp
            % (:,m)矩阵的第m列
            DAHIF(:,m) = DAHIF_times_TANBIC(:,m)/tanBetaI_c(m);
            DRHIF(:,m) = DRHIF_times_TANBIC(:,m)/tanBetaI_c(m);
        end  
    else
        Gd = 1;  % if set Gd==0, then Gd_res would be Inf  
        UADUCT = zeros(Mp,1);
    end
    
    % 平滑矩阵 X_smooth = Bsmooth*X;
    Bsmooth = RepairSplineMatrix(rc);
    
    %% 进行环量优化
    % 初始化迭代参数
    LM = -1;            % 拉格朗日乘子（用于螺旋桨）
    LM_last = LM;       % 拉格朗日乘子上一次的迭代值
    Gp_last = 0*Gp;     % 螺旋桨环量上一次的迭代值
    Gd_last = 0;        % 涵道环量上一次的迭代值
    G_iter = 1;         % 环量迭代次数
    Gp_res = 1;         % 螺旋桨环量两次迭代之间的差值
    Gd_res = 0;         % 涵道环量两次迭代之间的差值
    C_res = 0;          % 弦长两次迭代之间的差值  
    G_TOL = 1e-4;       % 环量达到收敛条件的最大残差值
    
    % 循环中方程求解直接所得：
    % Gp LM
    % 依据方程解和已知或上一代求解结果所得：
    % UASTAR UTSTAR tanBetaI_c DAHIF DRHIF UARING URRING Gd UADUCT
    % VSTAR UAHIF UTHIF
    while G_iter <= ITER && any([Gp_res;Gd_res] > G_TOL)
        % A为变量{Γ(1),Γ(2),...,Γ(Mp),λ}的系数，也就是方程组等号左侧的系数矩阵
        A = zeros(Mp+1,Mp+1);
        % B为方程组等号右侧的常数项向量
        B = zeros(Mp+1,1);
        for i = 1:Mp                           % for each equation for G(i)
            for m = 1:Mp                       % for each vortex panel, m        
                A(i,m) = UAHIF(m,i)*rc(m)*drv(m)+UAHIF(i,m)*rc(i)*drv(i)...
                         +LM_last*UTHIF(m,i)*drv(m)+LM_last*UTHIF(i,m)*drv(i);
            end  
            B(i) = -(VA_c(i)+UADUCT(i))*rc(i)*drv(i);
            A(i,Mp+1) = (L*rc(i)+VT_c(i))*drv(i);
        end

        if TorqueSpec_flag == 0 
            % 要求升力满足输入条件
            for m = 1:Mp                          
                A(Mp+1,m) = (L*rc(m)+VT_c(m)+UTSTAR(m))*drv(m);  
            end

            B(Mp+1) = CTp_desired/(4*Z)+(1/(2*pi))*...
                      sum(CDp_c.*VSTAR.*CpoDp_c.*(VA_c+UADUCT+UASTAR).*drv);
            if Hub_flag == 1
                % 加入桨毂阻力
                B(Mp+1) = B(Mp+1) + (Z/8)*(log(1/Rhv)+3)*(Gp_last(1)^2);
            end
        else
            % 要求扭矩满足条件
            for m = 1:Mp
                A(Mp+1,m) = (VA_c(m)+UADUCT(m)+UASTAR(m))*rc(m)*drv(m);
            end
            B(Mp+1) = CQdes/(4*Z)-(1/(2*pi))*...
                      sum(CDp_c.*VSTAR.*CpoDp_c.*(L*rc + VT_c + UTSTAR).*rc.*drv);
        end

        % 求解线性方程组
        GLM = linsolve(A,B);
        % 解的前Mp项为当前本次迭代求得的环量值
        Gp = GLM(1:Mp);
        % 解的第Mp+1项为本次迭代求得的拉格朗日乘子
        LM = GLM(Mp+1);

        % 更新诱导速度
        UASTAR = UAHIF*Gp;  
        UTSTAR = UTHIF*Gp;  

        % 更新水动力螺旋角的正切值
        tanBetaI_c = (VA_c+UADUCT+UASTAR)./(L*rc+VT_c+UTSTAR);
        % 优化水动力螺旋角
        tanBetaI_c_smooth = Bsmooth*tanBetaI_c;

        if any(isnan(GLM))||~isreal(GLM)||max(Gp-Gp_last) > 10
            % 解无效或非实数解或环量突变时进行警告
            Gp = zeros(Mp,1);
            UASTAR = zeros(Mp,1);
            UTSTAR = zeros(Mp,1);
            tanBetaI_c = tanBeta_c;
            % 弹出错误框
            error = sprintf('使用拉格朗日乘子法进行环量优化时，拉格朗日方程组无解或解无效');
            uialert(Fig_Main,error,'求解失败');
        end   

        if Duct_flag == 1
            % 更新叶片对涵道的诱导速度
            for m = 1:Mp
                DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / tanBetaI_c_smooth(m);
                DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / tanBetaI_c_smooth(m);
            end
            UARING = DAHIF*Gp;  
            URRING = DRHIF*Gp; 

            % 更新涵道上的环量
            [~,Gd] = Duct_Thrust(XdRING,RdoRp,VARING,UARING,URRING,...
                                 GdRING,Gd,CDd,CTd_desired);

            % 更新涵道在涵道内流场产生的轴向诱导速度
            UADUCT = UADIF*Gd;                        
        end
        
        % 更新总来流速度
        VSTAR = sqrt((VA_c+UADUCT+UASTAR).^2+(L*rc+VT_c+UTSTAR).^2);       
        
        % 更新叶片内部影响
        [UAHIF,UTHIF] = Horseshoe(Mp,Z,tanBetaI_c,rc,rv,Hub_flag,...
                                  RhoRp,Duct_flag,RdoRp);     
        
        % 更新弦长分布
        if strcmp(ChordMethod,'CLmax')
            % 直接由输入获得升力系数和厚径比，由此得到弦径比
            CpoDp_c = 2*pi*Gp./(VSTAR.*CL_c);
            
        elseif strcmp(ChordMethod,'ConeyPLL')
            % 直接由输入获取厚径比，根据计算获得升力系数
            SIGMA = SIGMAs./VSTAR.^2;    % local cavitation number
            F0oDp = (2*pi*Gp./VSTAR).*F0oCdtilde_c./CLItilde_c;
            CpoDp_c = (8.09*F0oDp+3.033*T0oDp_c)./(2*SIGMA)+...
                       sqrt((8.09*F0oDp+3.033*T0oDp_c).^2+...
                       4*SIGMA.*(26.67*F0oDp.^2+10*F0oDp.*T0oDp_c))./(2*SIGMA);
                   
        end

        % 如果伸张面积比有特殊要求，进行比例缩放
        if EARspec > 0
           EAR = (2*Z/pi)*trapz(linspace(RhoRp,1,100),...
                 interp1(rc,CpoDp_c,linspace(RhoRp,1,100), 'spline','extrap'));  
           CpoDp_c = (EARspec/EAR)*CpoDp_c; 
        end
        
        if all(CpoDp_c == 0) || any(isnan(CpoDp_c))
             Gp = 0*rc;
             UASTAR = 0*rc;
             UTSTAR = 0*rc;
             tanBetaI_c = tanBeta_c;
             CpoDp_c = 0*rc;
             % 弹出错误框
             error = sprintf(['使用',ChordMethod,...
                             '法求解厚径比失败 \n 请切换至其他方法']);
             uialert(Fig_Main,error,'求解失败');
             G_iter = 999;
        end
        
        % 更新迭代参数
        Gp_res = abs((Gp-Gp_last)./Gp);
        Gp_last = Gp;
        Gd_res = abs((Gd-Gd_last)/Gd);
        Gd_last = Gd;
        LM_last = LM;

        if G_iter < 10
            disp(['The max G_res for iteration ',num2str(G_iter),...
                 ' is:',num2str(max(Gp_res))]),  
        else
            disp(['The max G_res for iteration ',num2str(G_iter),...
                 'is: ',num2str(max(Gp_res))]),  
        end    
        
        % 迭代次数+1
        G_iter = G_iter+1;
    end
    
    % 进行优化结果说明
    if G_iter > ITER
        disp('警告：结果未收敛');
        % 弹出错误框
        error = sprintf('求解结果未收敛');
        uialert(Fig_Main,error,'求解失败');
        Converge_flag = 0;
    else
        disp('完成优化，结果收敛');
        Converge_flag = 1;
    end
    
    % If required, unload the hub and tip, then rescale the circulation
    % distribution to get the desired value of the thrust coefficient.
    if Hub_flag && (HUF > 0 || TUF > 0)
        if Duct_flag == 0
            DAHIF_times_TANBIC = 0;
            DRHIF_times_TANBIC = 0;
            XdRING = 0;
            RdoRp = 1;
            VARING = 0;
            GdRING = 0;
            Gd = 0;
            CDd = 0;
            CTd_desired = 0;
        end

        [Gp,UASTAR,UTSTAR,tanBetaI_c,UARING,URRING,Gd,UADUCT] = ...
            Unload_Blade(HUF,TUF,rc,RhoRp, Gp,  VA_c,VT_c, tanBetaI_c,rv,drv,L,Mp,Z, ...
                         Hub_flag,ITER,Bsmooth,CDp_c,CpoDp_c,Js,VMIV,Rhv,CTp_desired,...
                         Duct_flag,UADUCT,XdRING,RdoRp,VARING,GdRING,UADIF,Gd,CDd,CTd_desired,...
                         DAHIF_times_TANBIC,DRHIF_times_TANBIC,...
                         Plot_flag,Hgamma,HHgamma,Hvel,HHvel,Hbeta,HHbeta);    

    end
    
    % 计算无量纲参数和各类系数
    % 根据优化后的涵道环量计算涵道实际所占的推力系数
    if Duct_flag == 1
        [CTd_desired,~] = Duct_Thrust(XdRING,RdoRp,VARING,UARING,URRING,...
                                      GdRING,Gd,CDd,CTd_desired);
    end

    [CT,CQ,CP,KT,KQ,CTH,TpoT,Ja,~,VMWV,EFFYo,EFFY,ADEFFY,QF,QFo,QFw] = ...
        Forces(rc,drv,VA_c,VT_c,UASTAR,UTSTAR,UADUCT,CDp_c,CpoDp_c,Gp,Z,Js,...
               VMIV,Hub_flag,RhoRp,Rhv,CTd_desired);
    
    % 根据优化后的环量和速度计算实际的升力系数
    CL_c = 2*pi*Gp./(VSTAR.*CpoDp_c);

    % Expanded Area Ratio
    EAR = (2*Z/pi)*trapz(linspace(RhoRp,1,100),...
          interp1(rc,CpoDp_c,linspace(RhoRp,1,100),'spline','extrap'));  
    
    % 更新厚弦比
    T0oCp_c = T0oDp_c ./ CpoDp_c;    
    
    % 计算消耗的扭矩
    Q = CQ*0.5*rho*VS^2*pi*Dp^2/4*Dp/2;
    
    % 计算螺旋桨消耗的功率
    P = Q*2*pi*N/60;
    
    %% 打印和保存计算结果
    design.part1 = '叶片相关的计算结果';
    design.rc = rc;                         % [Mp*1] 控制点坐标
    design.rv = rv;                         % [Mp+1*1] 涡格点坐标
    design.drv = drv;                       % [Mp*1] 相邻两涡格点坐标差值
    design.Gp = Gp;                         % [Mp*1] 控制点/涡格环量
    design.VA_c = VA_c;                     % [Mp*1] 轴向来流速度 
    design.VT_c = VT_c;                     % [Mp*1] 切向来流速度 
    design.UASTAR = UASTAR;                 % [Mp*1] 控制点处总轴向诱导速度 
    design.UTSTAR = UTSTAR;                 % [Mp*1] 控制点处总切向诱导速度 
    design.VSTAR = VSTAR;                   % [Mp*1] 总相对来流速度  
    design.tanBeta_c = tanBeta_c;           % [Mp*1] 进角
    design.tanBetaI_c = tanBetaI_c;         % [Mp*1] 水动力螺旋角 
    design.CL_c = CL_c;                     % [Mp*1] 升力系数 
    design.CDp_c = CDp_c;                   % [Mp*1] 阻力系数 
    design.CpoDp_c = CpoDp_c;               % [Mp*1] 弦径比
    design.T0oCp_c = T0oCp_c;               % [Mp*1] 厚弦比 
    design.T0oDp_c = T0oDp_c;               % [Mp*1] 厚径比 

    design.part2 = '其余计算结果';
    design.Converge_flag = Converge_flag;   % 计算结果是否收敛
    design.iteration = G_iter;              % 优化总迭代次数
    design.RhoRp = RhoRp;                   % 桨毂半径/叶片半径
    design.EAR = EAR;                       % 伸张面积比
    design.LM = LM;                         % 拉格朗日乘子
    design.VMIV = VMIV;                     % 体平均速度
    design.VMWV = VMWV;                     % 
    design.SIGMAs = SIGMAs;                 % 
    design.Q = Q;                           % 螺旋桨扭矩
    design.P = P;                           % 螺旋桨功率
    
    if Duct_flag == 1  
        design.part3 = '涵道计算结果';
        design.RdoRp = RdoRp;               % 涵道半径/叶片半径
        design.CdoRp = CdoRp;               % 涵道弦长/叶片半径
        design.XdoRp = XdoRp;               % 涵道位置/叶片半径
        design.Gd = Gd;                     % 涵道涡环环量
        design.VARING = VARING;             % 涵道中线处轴向来流速度
        
        design.XdRING = XdRING;             % [Nd*1] 涡环轴向坐标
        design.UARING = UARING;             % [Nd*1] 涵道受到叶片轴向的诱导速度
        design.URRING = URRING;             % [Nd*1] 涵道受到叶片径向的诱导速度
        design.GdRING = GdRING;             % [Nd*1] 涡环环量沿轴向的占比

        design.DAHIFtT = DAHIF_times_TANBIC;% [Nd,Mp]
        design.DRHIFtT = DRHIF_times_TANBIC;% [Nd,Mp]
        
        design.UADIF = UADIF;               % [Mp*1] 涵道对涵道内流场的轴向诱导速度
        design.UADUCT = UADUCT;             % [Mp*1] 涵道对涵道内流场的轴向诱导速度
        design.CTp_desired = CTp_desired;   % 叶片的推力系数
        design.CTd_desired = CTd_desired;   % 涵道的推力系数

        design.TpoT = TpoT;                 % 叶片升力占比
        design.CDd = CDd;                   % 涵道阻力系数
        design.part4 = '性能参数';
    else
        design.part3 = '性能参数';
    end
    
    design.L = L;
    design.Js = Js;
    design.KT = KT;
    design.KQ = KQ;
    design.CT = CT;
    design.CQ = CQ;
    design.CP = CP;
    design.CTH = CTH;

    if abs(VMIV-1) > 1e-8   
        % i.e. if VMIV is not equal to 1 
        design.EFFYo = EFFYo;
        design.Ja = Ja;
    end
    design.EFFY = EFFY;
    design.ADEFFY = ADEFFY;
    design.QF = QF;

    if (VMIV < 0.05)   
        % assume design for bollard pull
        design.QFo = QFo;
        design.QFw = QFw;
    end
end