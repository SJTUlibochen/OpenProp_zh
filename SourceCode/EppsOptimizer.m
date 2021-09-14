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
    
    % part2 “叶片切面参数”
    Nx = input.Nx;
    ChordMethod = input.ChordMethod;
    Meanline_flag = input.Meanline_flag;
    Thickness_flag = input.Thickness_flag;
    rx = input.rx;
    CDp_x = input.CDp_x;
    if Meanline_flag
        Meanline_x = input.Meanline_x{1};
    else
        Meanline_x = input.Meanline_x;
    end    
    if Thickness_flag
        Thickness_x = input.Thickness_x{1};
    else
        Thickness_x = input.Thickness_x;
    end    
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
    XdoRp = 0;
    HUF = 0;
    TUF = 0;
    Rhv = 0.5;
    H = 3.048;
    g = 9.81;
    Patm = 101325;
    Pv = 2500;
    
    %% 引申参数
    RhoRp = Dh/Dp;
    
    RdoRp = Dd/Dp;
    CdoRp = 2*Cd/Dp;
    
    x1 = ones(Nx);
    x0 = zeros(Nx);
    i1 = ones(Ni);
    i0 = zeros(Ni);
    
    dVA_i = diff(VA_i);
    
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
    
    % 输入中没有指定扭矩，若不指定推力而指定扭矩则为1
    TorqueSpec_flag = 0;
    
    % 指定展开面积比 Specified Expanded Area Ratio (EAR)
    EARspec = 0;
    
    % 高速情况
    if strcmp(ChordMethod,'FAST2011dCTP')
        % Method: see (Epps et al., FAST'2011)
        if isfield(input,'VH')
            VH = input.VH;
        else
            VH = VS;
        end
        Jh = Js;  % initial guess for advance coefficient for high-speed state
        if isfield(input,'THRUSTh')
            CT_desired_high = input.THRUSTh/(0.5*rho*VH^2*pi*Rp^2);
        elseif isfield(input,'CT_desired_high')
            CT_desired_high = input.CT_desired_high; 
        else
            CT_desired_high = CT_desired;
        end
        SIGMAh = (Patm+rho*g*H-Pv)/(0.5*rho*VH^2);
    end

    % 获取切面拱度和理想升力系数的数据
    if ~Meanline_flag && ~Thickness_flag
        % 凸缘线和叶型均没有勾选同一性
        if length(Meanline_x) ~= Nx && length(Thickness_x) ~= Nx
            error = sprintf('输入数据中凸缘线或叶型的元胞数组大小与径向位置分段数量不同');
            uialert(Fig_Main,error,'计算失败');
            return
        else
            % 提前分配内存
            Xf0octilde = zeros(Nx,1);
            XCLItilde = zeros(Nx,1);
            for index = 1 : Nx
                [Xf0octilde(index),XCLItilde(index)] = GeometryFoil2D(Meanline_x{index},Thickness_x{index});
            end
        end
    elseif Meanline_flag && ~Thickness_flag
        % 凸缘线勾选同一性，叶型没有勾选同一性
        [f0octilde,CLItilde] = GeometryFoil2D(Meanline_x,Thickness_x);
        Xf0octilde = f0octilde * ones(size(rx));
        XCLItilde = CLItilde * ones(size(rx));
    end
    
    %% 计算控制点和涡格点的径向位置并进行插值
    %Compute the Volumetric Mean Inflow Velocity, eqn 163, p.138
    %生成[Rhub_oR,1]之间100个等距点组成的向量
    XRtemp = linspace(RhoRp,1,100);         % (temp) radius / propeller radius
    XVAtemp = pchip(ri,VA_i,XRtemp);            % (temp) axial inflow velocity /ship velocity
    VMIV = 2*trapz(XRtemp,XRtemp.*XVAtemp)/(1^2-RhoRp^2);  % [ ], VMIV/ship velocity
    %计算生成涡点和控制点的径向位置
    [RC,RV,DR] = LLPanelRadii(Mp,RhoRp,Hub_flag,Duct_flag);
    %在相应径向位置进行插值
    VAC = pchip(ri,VA_i,RC);   % axial inflow vel.  / VS at ctrl pts
    dVAC = pchip(ri,dVA_i,RC);   % delta( axial inflow vel.) / VS due to wake, for chord optimization
    VTC = pchip(ri,VT_i,RC);   % tangential inflow vel. / ship vel. at ctrl pts
    Cdp = pchip(rx,CDp_x,RC);   % section drag coefficient at ctrl pts
    t0oDp = pchip(rx,T0oDp_x,RC);   % section thickness / propeller dia. at ctrl pts
    CLmax = pchip(rx,CL_x,RC);   % maximum allowable lift coefficient at ctrl pts
    f0octilde = pchip(rx,Xf0octilde,RC);
    CLItilde = pchip(rx,XCLItilde,RC);
    %针对叶梢的情况选择不同的插值方法
    if (abs(rx(end)-1) < 1e-4) && (XCpoDp(end) <= 0.01)  % if XR == 1 and XCpoDp == 0
        CpoDp = InterpolateChord(rx,XCpoDp,RC);   % section chord / propeller diameter at ctrl pts
    else
        CpoDp = pchip(rx,XCpoDp,RC);   % section chord / propeller diameter at ctrl pts
    end
    t0oc = t0oDp./CpoDp; % input t0/c values
    %% 设置计算的初始值
    %基于桨盘理论初始化诱导速度和部分角度的估值
    if Propeller_flag == 1
        G = zeros(Mp,1);        %G = Gamma / 2*pi*Rp*VS
        UASTAR = 0*RC + 0.5*(sqrt(1+CTPdes)-1);     % Kerwin & Hadler (2010), eqn (4.26)
        UTSTAR = 0*RC;
        if (VMIV < 0.05)                       % assume bollard pull propeller design
            UASTAR = 0*RC + 0.5*sqrt(CTPdes);  % (Epps, 2013)
        end
    else                                                %  element is a turbine
        [CPBetz,BetzRC,BetzG,BetzUA,BetzUT,BetzTAN, BetzGRC] = Turbine_ADS_Theory(L,Z,CDpoCL,RC);
        G = BetzGRC';                              % G = Gamma / 2*pi*Rp*VS
        UASTAR = pchip(BetzRC,BetzUA,RC);
        UTSTAR = pchip(BetzRC,BetzUT,RC);
    end
    TANBC = VAC./(L*RC + VTC);
    TANBIC = (VAC + UASTAR)./(L*RC + VTC + UTSTAR);
    VSTAR = sqrt((VAC+UASTAR).^2 + (L*RC+VTC+UTSTAR).^2);        % V* / VS
    dVdG = zeros(Mp,Mp);					   	 % (2piR) * d(V*) /d(Gamma)
    
    %如果进行弦长优化，弦长迭代的初始值应设置较小
    CpoDp = 0.01*ones(size(RC));
    
    %初始化马蹄涡相关参数
    %Initialize vortex Horseshoe Influence Functions 
    [UAHIF,UTHIF] = Horseshoe(Mp,Z,TANBIC,RC,RV,Hub_flag,RhoRp,Duct_flag,RdoRp);
    %初始化涵道相关参数
    %Initialize duct variables
    if Duct_flag == 1
        %XdRING            [1,Nd],x/Rp location of each vortex ring downstream of propeller
        %VARING            [1,Nd],(axial free-stream velocity at duct)/VS
        %GdRING            [1,Nd],fraction of total non-dimensional duct circulation, sucht that sum(GdRING) = 1 
        %                         (i.e. non-dimensional circulation per unit Gd)
        %Gd                       total non-dimensional circulation about the duct, Gd == Gamma_d/(2*pi*Rp*VS),
        %                         such that the non-dimensional circulation of ring n is Gd*GdRING(n).
        %Gamma_d [m^2/s]          total dimensional circulation about the duct [m^2/s]
        %Cdd                      section drag coefficient for the duct
        %CTDdes                   desired duct thrust coefficient
        %CTD                      duct thrust coeff (viscous drag included) with total duct circulation of Gd 
        %
        %Influence of propeller on duct:
        %
        %(DAHIF ,DRHIF)     [Nd,Mp], (axial/radial horseshoe influence functions of prop on duct)*(2*pi*Rp)
        %(UARING,URRING)    [1,Nd],  (axial/radial velocity induced at duct (XdRING,Rduct_oR) by prop) / VS
        %
        %
        %Influence of duct on propeller:
        %
        %UARING      [1,Mp]  (axial velocity induced on PROP by duct)/VS
        %UADIF       [1,Mp]   axial velocity induced on PROP by duct per unit Gd, 
        %                    i.e. non-dimensional Duct Influence Function
        %                    == 2*pi*Rp * (duct influence function with unit dimensional duct circulation)   
        [XdRING,GdRING,UADIF] = Duct_Influence(RdoRp,CdoRp,Xduct_oR,RC); % Duct geometry, and influence of duct on propeller

        % ---------------------------- (axial inflow velocity / VS) at Rduct_oR
        VARING = pchip(ri,VA_i,RdoRp);  

        % ------------------------------------------------ Initial guess for Gd
        Gd     = 0;                  % total duct circulation / (2*pi*Rp*VS)
        UADUCT = UADIF * Gd;         % axial velocity induced at RC by the duct

        % ---------------------------------------------------------------------
        % Initialize Duct Horseshoe Influence Functions (influence of propeller on duct)
        %
        % DAHIF(n,m) = influence of m-th horseshoe vortex shed from propeller (Mp panels) 
        %                    on the n-th control point of the duct            (Nd rings) 
        %
        disp(' '), disp('Computing rotor-duct interaction...be patient...'), disp(' '),
        [DAHIF_times_TANBIC, DTHIF, DRHIF_times_TANBIC] = Horseshoe_intr_110830(XdRING,RdoRp, RC,ones(size(TANBIC)),RV,Z,Hub_flag,RhoRp,Duct_flag,RdoRp); 
        DAHIF = 0* DAHIF_times_TANBIC; % memory allocation
        DRHIF = 0* DRHIF_times_TANBIC; % memory allocation
        for m = 1:Mp                         
            DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBIC(m);
            DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBIC(m);
        end  
    else
        Gd = 1;  % if set Gd==0, then Gd_res would be Inf  
        UADUCT = 0*RC;
        CTD = 0;  
    end
    %Implement RepairSpline.m
    % To smooth data X(RC):  X_smooth = X*Bsmooth;
    Bsmooth = RepairSplineMatrix(RC);
    %% 进行环量优化（即螺旋桨优化）准备
    disp('开始环量优化');
    %螺旋桨优化方法
    %优化方法为LL-Linear
    if isfield(input,'EppsOptimizer02_flag')
        EppsOptimizer02_flag = input.EppsOptimizer02_flag;
    else
        EppsOptimizer02_flag = 1;
    end  % Propeller: LL-Linear
    %优化方法为LL-Newton (standard hub drag model)
    if isfield(input,'EppsOptimizer23_flag')
        EppsOptimizer23_flag = input.EppsOptimizer23_flag;
    else
        EppsOptimizer23_flag = 0; 
    end  % Propeller: LL-Newton (standard hub drag model)
    %优化方法为LL-Newton (variational hub drag model)
    if isfield(input,'EppsOptimizer53_flag')
        EppsOptimizer53_flag = input.EppsOptimizer53_flag; 
    else
        EppsOptimizer53_flag = 0; 
    end  % Propeller: LL-Newton (variational hub drag model)
    if (EppsOptimizer23_flag == 1 || EppsOptimizer53_flag == 1)
        EppsOptimizer02_flag = 0;
    end
    %初始化优化迭代所用参数
    LM = -1;        %拉格朗日乘子（用于螺旋桨）
    LM_last = LM;       %拉格朗日乘子上一次的迭代值
    G_last = 0*G;       %螺旋桨环量上一次的迭代值
    Gd_last = 0;        %涵道环量上一次的迭代值
    G_iter = 1;     %环量迭代次数
    G_res = 1;      %螺旋桨环量两次迭代之间的差值
    Gd_res = 0;     %涵道环量两次迭代之间的差值
    C_res = 0;      %弦长两次迭代之间的差值  
    relax = 0.9;        % Newton solver relaxation parameter
    G_TOL = 1e-4;       %环量达到收敛条件的最大残差值
    if (Chord_flag == 1) && ((strcmp(ChordMethod,'FAST2011dCTP') == 1) || (strcmp(ChordMethod,'FAST2011dVAC') == 1) || (strcmp(ChordMethod,'Brizzolara2007') == 1))
        ITER = ITER*3;  % allow three circulation iterations per chord iteration
    end
    %% 优化循环
    while G_iter <= ITER && any([G_res;Gd_res] > G_TOL)
        %循环的进入标准：没有达到最大迭代次数，且环量均没有收敛
        % UPDATE: G, UASTAR, UTSTAR, TANBIC, (duct stuff), LM
        if Propeller_flag == 1
            %% 迭代方法为LL-Linear
            %根据Coney的博士论文，拉格朗日乘子法必要条件，即偏导数为0所得的方程组中
            %前Mp个方程，也就是对Γ求偏导的化简结果中均包括所有的Mp+1个变量
            %而总共Mp+1个变量，共有Mp+1个方程，
            if EppsOptimizer02_flag == 1
                %A为变量{Γ(1),Γ(2),...,Γ(Mp),λ}的系数，也就是方程组等号左侧的系数矩阵
                A = zeros(Mp+1,Mp+1);
                %B为方程组等号右侧的常数项向量
                B = zeros(Mp+1,1);
                for i = 1:Mp                           % for each equation for G(i)
                    for m = 1:Mp                       % for each vortex panel, m        
                        A(i,m) = UAHIF(m,i)*RC(m)*DR(m)+UAHIF(i,m)*RC(i)*DR(i)...
                                 +LM_last*UTHIF(m,i)*DR(m)+LM_last*UTHIF(i,m)*DR(i);
                    end  
                    B(i) = -(VAC(i)+UADUCT(i))*RC(i)*DR(i);
                    A(i,Mp+1) = (L*RC(i)+VTC(i))*DR(i);
                end
                % The (Mp+1) equation is either the thrust constraint or torque constraint
                % Note: A(Mp+1,Mp+1) = 0            
                if TorqueSpec_flag == 0  % thrust is specified, CT_prop_inviscid/(4*Z) == CTPdes/(4*Z) + CT_prop_viscous/(4*Z) + CT_hub/(4*Z)
                    for m = 1:Mp                          
                        A(Mp+1,m) = (L*RC(m) + VTC(m) + UTSTAR(m))*DR(m);  
                    end
                    B(Mp+1) =  CTPdes/(4*Z)  +  (1/(2*pi))*sum(Cdp.*VSTAR.*CpoDp.*(VAC + UADUCT + UASTAR).*DR);
                    if Hub_flag == 1
                        B(Mp+1) = B(Mp+1) + (Z/8)*(log(1/Rhv)+3)*(G_last(1)^2);
                    end
                else % torque is specified, CQ_prop_inviscid/(4*Z) == CQdes/(4*Z) - CQ_prop_viscous/(4*Z)
                    for m = 1:Mp
                        A(Mp+1,m) = (VAC(m) + UADUCT(m) + UASTAR(m))*RC(m)*DR(m);
                    end
                    B(Mp+1) =  CQdes/(4*Z)  -  (1/(2*pi))*sum(Cdp.*VSTAR.*CpoDp.*(L*RC + VTC + UTSTAR).*RC.*DR);
                end
                % -------------------------------- Solve linear system of equations
                GLM = linsolve(A,B);             
                G = GLM(1:Mp);                     % G is the first Mp entries
                LM = GLM(Mp+1);                     % LM is the last entry
                % ----------------------------------------------------------------- 
                % Update induced velocities (influence of propeller on propeller)
                UASTAR = (UAHIF*G)';  
                UTSTAR = (UTHIF*G)';  
                TANBIC = (VAC + UADUCT + UASTAR)./(L*RC + VTC + UTSTAR);
                % Smooth the inflow angle for numerical stability:
                TANBICsmooth = TANBIC * Bsmooth;
                % -----------------------------------------------------------------

                if any(isnan(GLM))  || ~isreal(GLM) || max(G-G_last) > 10
                    G = 0*RC';
                    UASTAR = 0*RC;
                    UTSTAR = 0*RC;
                    TANBIC = TANBC;

                    disp(' ')
                    disp('<WARNING>')
                    disp('<WARNING> GLM == NaN or imaginary... crash avoided...')
                    disp('<WARNING>')
                    disp(' ')
                    disp('Switching numerical method to Newton solver...')
                    disp(' ')

                    EppsOptimizer02_flag = 0;
                    EppsOptimizer23_flag = 1;
                end   
                % -----------------------------------------------------------------      

                % ----------------------------------------------------------------- 
                if Duct_flag == 1
                    % Update the influence of propeller on duct
                    for m = 1:Mp
                        DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBICsmooth(m);
                        DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBICsmooth(m);
                    end

                    % Update induced velocities at the duct (influence of propeller on duct)
                    UARING = (DAHIF*G)';  
                    URRING = (DRHIF*G)'; 

                    % If propeller torque specified, then need to update desired duct thrust based on current prop thrust and specified thrust ratio, TAU
                    if TorqueSpec_flag == 1

                        CTP = 4*Z*sum(  G'.*(L*RC + VTC + UTSTAR).*DR  -  (1/(2*pi)).*VSTAR.*CpoDp.*Cdp.*(VAC + UADUCT + UASTAR).*DR  );

                        if Hub_flag == 1
                            CTP = CTP - 0.5*(log(1/Rhv)+3)*(Z*G(1))^2;   
                        end

                        CTDdes = ( (1-TpoT)/TpoT ) * CTP;
                    end

                    % Update duct circulation such that CTD == CTDdes            
                    [junk,Gd] = Duct_Thrust(XdRING,RdoRp,VARING,UARING,URRING,GdRING,Gd,CDd,CTDdes);

                    % Update the induced velocities at the propeller (influence of duct on propeller)
                    UADUCT =  UADIF*Gd;                        
                end
                % ----------------------------------------------------------------- 

            end
            %% 迭代方法为LL-Newton
            if (EppsOptimizer23_flag == 1 || EppsOptimizer53_flag == 1)
                % ------------------------------------------ Execute Newton solver   
                RNS = zeros(4*Mp+1,     1);     % Residual vector
                JNS = zeros(4*Mp+1,4*Mp+1);     % Jacobian

                % Evaluate induced velocities
                UASTARtemp = (UAHIF*G)';  
                UTSTARtemp = (UTHIF*G)';

                % ---------------------------------------------- Evaluate residuals        
                for i = 1:Mp

                    % 11/16/2011 BEPPS: This was the implementation in EppsOptimizer23.m, but the inclusion
                    %       of these drag terms has negligable effect ( O(10^-4) ) on the converged circulation distribution.
                    %            
                    % RNS(i) =       (VAC(i) + UADUCT(i) + UASTAR(i))     *RC(i)*DR(i)  ...
                    %          + sum( UAHIF(:,i)'.*G'.*RC  .*DR   ) ...
                    %          + sum( (1/(2*pi))*Cdp.*CpoDp .* dVdG(:,i)' .* (L*RC + VTC + UTSTAR) .* RC .* DR ) ...
                    %          + sum( (1/(2*pi))*Cdp.*CpoDp .*  VSTAR     .* (      UTHIF(:,i)'      ) .* RC .* DR ) ...
                    %          ...    
                    %          + LM * ( ...
                    %                        (L*RC(i) + VTC(i) + UTSTAR(i)) * DR(i)  ...
                    %                  + sum( UTHIF(:,i)' .* G'                .* DR   ) ...
                    %                  - sum( (1/(2*pi))*Cdp.*CpoDp .* dVdG(:,i)' .* (  VAC + UADUCT + UASTAR     ) .* DR ) ...
                    %                  - sum( (1/(2*pi))*Cdp.*CpoDp .*  VSTAR     .* (                 UAHIF(:,i)') .* DR ) ...
                    %                 );
                    %
                    % ------
                    RNS(i) =  (VAC(i) + UADUCT(i) + UASTAR(i)) *RC(i) *DR(i)  ...
                             +        sum( UAHIF(:,i)' .* G'  .*RC   .*DR   ) ...
                             ...
                             + LM * ( sum( UTHIF(:,i)' .* G'        .* DR   ) ...
                                     +(L*RC(i) + VTC(i) + UTSTAR(i)) * DR(i)  ...
                                    );            
                    % ------

                    RNS(i+  Mp) = UASTAR(i) - UASTARtemp(i);
                    RNS(i+2*Mp) = UTSTAR(i) - UTSTARtemp(i);
                    RNS(i+3*Mp) = TANBIC(i) - (VAC(i) + UADUCT(i) + UASTAR(i))/(L*RC(i) + VTC(i) + UTSTAR(i)); 
                end    

                % The (Mp+1) equation is either the thrust constraint or torque constraint
                %
                if TorqueSpec_flag == 0  % thrust is specified

                        RNS(1+4*Mp) = sum((L*RC + VTC + UTSTAR).*G'.*DR - (1/(2*pi))*Cdp.*CpoDp.*VSTAR.*(VAC + UADUCT + UASTAR).*DR)  - CTPdes/(4*Z);                

                    if     Hub_flag == 1 && EppsOptimizer23_flag == 1

                        RNS(1+4*Mp) = RNS(1+4*Mp) - (Z/8)*(log(1/Rhv)+3)*(G_last(1)^2);  % EppsOptimizer23.m hub drag treatment  (do NOT include hub drag in variational optimization)

                    elseif Hub_flag == 1 && EppsOptimizer53_flag == 1

                        RNS(1)      = RNS(1)      - (Z/4)*(log(1/Rhv)+3)*(G(1)*LM);      % EppsOptimizer53.m hub drag treatment  (include hub drag in variational optimization) 

                        RNS(1+4*Mp) = RNS(1+4*Mp) - (Z/8)*(log(1/Rhv)+3)*(G(1)^2);       % EppsOptimizer53.m hub drag treatment  (include hub drag in variational optimization)   
                    end

                else  % torque is specified

                        RNS(1+4*Mp) = sum( (VAC + UADUCT + UASTAR).*G'.*RC.*DR + (1/(2*pi))*Cdp.*CpoDp.*VSTAR.*(L*RC + VTC + UTSTAR).*RC.*DR )  - CQdes/(4*Z);
                end

                % ----------------------------------------------- Evaluate Jacobian        
                for i = 1:Mp 
                    JNS(i     ,(1:Mp)     ) = UAHIF(:,i)' .* RC .* DR  +  LM * UTHIF(:,i)' .* DR;

                    JNS(i     , i    +  Mp) = JNS(i,i+  Mp) + RC(i)*DR(i);
                    JNS(i     , i    +2*Mp) = JNS(i,i+2*Mp) + LM   *DR(i);

                    % 11/16/2011 BEPPS: This was the implementation in EppsOptimizer23.m, but the inclusion
                    %       of these drag terms has negligable effect ( O(10^-4) ) on the converged circulation distribution.
                    %
                    % JNS(i     ,(1:Mp)+  Mp) = - LM * (1/(2*pi))*Cdp.*CpoDp.*dVdG(:,i)'    .*DR;  
                    % JNS(i     ,(1:Mp)+2*Mp) =        (1/(2*pi))*Cdp.*CpoDp.*dVdG(:,i)'.*RC.*DR;  
                    %
                    % JNS(i ,1+4*Mp)      =       (L*RC(i) + VTC(i) + UTSTAR(i)) * DR(i)  ...
                    %                       + sum( UTHIF(:,i)' .* G'                .* DR   ) ...
                    %                       - sum( (1/(2*pi))*Cdp.*CpoDp .* dVdG(:,i)' .* (  VAC + UADUCT + UASTAR     ) .* DR ) ...
                    %                       - sum( (1/(2*pi))*Cdp.*CpoDp .*  VSTAR     .* (        UAHIF(:,i)') .* DR );
                    %
                    % ------
                    JNS(i ,1+4*Mp)      =       (L*RC(i) + VTC(i) + UTSTAR(i)) * DR(i)  ...
                                          + sum( UTHIF(:,i)' .* G'            .* DR   );
                    % ------


                    JNS(i+  Mp,1:Mp  ) = - UAHIF(i,1:Mp);
                    JNS(i+  Mp,i+  Mp) = 1;

                    JNS(i+2*Mp,1:Mp  ) = - UTHIF(i,1:Mp);
                    JNS(i+2*Mp,i+2*Mp) = 1;

                    JNS(i+3*Mp,i+  Mp) = -1                          /(L*RC(i) + VTC(i) + UTSTAR(i));
                    JNS(i+3*Mp,i+2*Mp) = (VAC(i)+UADUCT(i)+UASTAR(i))/(L*RC(i) + VTC(i) + UTSTAR(i))^2;
                    JNS(i+3*Mp,i+3*Mp) = 1;


                    % The (Mp+1) equation is either the thrust constraint or torque constraint
                    %
                    if TorqueSpec_flag == 0  % thrust is specified

                        JNS(1+4*Mp,i     ) = (L*RC(i) + VTC(i) + UTSTAR(i))*DR(i);
                        JNS(1+4*Mp,i+  Mp) = - (1/(2*pi))*Cdp(i)*CpoDp(i)*VSTAR(i)*DR(i);    
                        JNS(1+4*Mp,i+2*Mp) = G(i)*DR(i);


                        if Hub_flag == 1 && EppsOptimizer53_flag == 1     
                            JNS(1,1)      = JNS(1,1)      - (Z/4)*(log(1/Rhv)+3)*LM;   % EppsOptimizer53.m hub drag treatment     
                            JNS(1,1+4*Mp) = JNS(1,1+4*Mp) - (Z/4)*(log(1/Rhv)+3)*G(1); % EppsOptimizer53.m hub drag treatment 
                            JNS(1+4*Mp,1) = JNS(1+4*Mp,1) - (Z/4)*(log(1/Rhv)+3)*G(1); % EppsOptimizer53.m hub drag treatment
                        end            

                    else  % torque is specified

                        JNS(1+4*Mp,i     ) = (VAC(i) + UADUCT(i) + UASTAR(i))*RC(i)*DR(i);
                        JNS(1+4*Mp,i+  Mp) =                            G(i) *RC(i)*DR(i);
                        JNS(1+4*Mp,i+2*Mp) = (1/(2*pi))*Cdp(i)*CpoDp(i)*VSTAR(i)*RC(i)*DR(i);             
                    end

                end

                % ----------------------------- Update Newton solver vector of unknowns
                DX     = linsolve(JNS,-RNS);
                G      = G      + relax*DX( 1:Mp      ) ;
                UASTAR = UASTAR + relax*DX((1:Mp)+  Mp)';
                UTSTAR = UTSTAR + relax*DX((1:Mp)+2*Mp)';
                TANBIC = TANBIC + relax*DX((1:Mp)+3*Mp)';
                LM     = LM     + relax*DX(     1+4*Mp) ;

                % Smooth the inflow angle for numerical stability:
                TANBICsmooth = TANBIC * Bsmooth;
                % -----------------------------------------------------------------  

                if any(isnan(DX))  || ~isreal(DX)  || any(DX > 999)
                     G      = 0*RC';
                     UASTAR = 0*RC;
                     UTSTAR = 0*RC;
                     TANBIC = TANBC;

                     disp(' ')
                     disp('<WARNING>')
                     disp('<WARNING> DX == NaN or imaginary... crash avoided...')
                     disp('<WARNING>')
                     disp(' ')
                     G_iter = 999;
                end
                % ----------------------------------------------------------------- 

                % -----------------------------------------------------------------
                if Duct_flag == 1
                    % Update the influence of propeller on duct
                    for m = 1:Mp;
                        DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBICsmooth(m);
                        DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBICsmooth(m);
                    end

                    % Update induced velocities at the duct (influence of propeller on duct)
                    UARING = (DAHIF*G)';  
                    URRING = (DRHIF*G)'; 


                    % If propeller torque specified, then need to update desired duct thrust based on current prop thrust and specified thrust ratio, TAU
                    if TorqueSpec_flag == 1

                        CTP = 4*Z*sum(  G'.*(L*RC + VTC + UTSTAR).*DR  -  (1/(2*pi)).*VSTAR.*CpoDp.*Cdp.*(VAC + UADUCT + UASTAR).*DR  );

                        if Hub_flag == 1
                            CTP = CTP - 0.5*(log(1/Rhv)+3)*(Z*G(1))^2;   
                        end

                        CTDdes = ( (1-TpoT)/TpoT ) * CTP;
                    end

                    % Update duct circulation such that CTD == CTDdes            
                    [junk,Gd] = Duct_Thrust(XdRING,RdoRp,VARING,UARING,URRING,GdRING,Gd,CDd,CTDdes);

                    % Update the induced velocities at the propeller (influence of duct on propeller)
                    UADUCT =  UADIF*Gd;                        
                end
            end   
        else
            RNS = zeros(4*Mp,1);     % Residual vector
            JNS = zeros(4*Mp,4*Mp);  % Jacobian

            % Evaluate induced velocities
            UASTARtemp = (UAHIF*G)';  
            UTSTARtemp = (UTHIF*G)';

            % ---------------------------------------------- Evaluate residuals        
            for i = 1:Mp

                if (Chord_flag == 1) && (strcmp(ChordMethod,'CLmax') == 1)
                    % ==========
                    % 11/16/2011 BEPPS: This new drag treatment is consistent with actuator disk theory.
                    %                   The drag term here assumes d(CL)/d(G) == 0, which is true in chord optimization with CL == CLmax.
                    % 11/17/2011 BEPPS: However, in the no-chord-optimization case, d(CL)/d(G) is not zero, so the drag term is zero to the leading order.  
                    %                   Furthermore, this code as is crashes for the no-chord-optimization case, so use the formulation below instead.
                    RNS(i) =   ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * ( VAC(i) + UADUCT(i) +   UASTAR(i))             ...    
                             -                                       (L*RC(i) +    VTC(i) + 2*UTSTAR(i)) * UTSTAR(i) ...
                             + ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * (L*RC(i) +    VTC(i) + 2*UTSTAR(i)) * Cdp(i)/CLmax(i);
                    % ==========
                else
                    % 11/16/2011 BEPPS: This was the "Robust method" implementation in EppsOptimizer06.m, but the inclusion
                    %                   of these drag terms has negligable effect ( O(10^-4) ) on the converged circulation distribution.
                    %            
                    % Note, the first line is equivalent to:
                    %
                    %             (VAC(i) + UADUCT(i))^2 +  (3*(VAC(i) + UADUCT(i))+ 2*UASTAR(i))*UASTAR(i) ...       
                    %
                    % RNS(i) =   ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * (VAC(i) + UADUCT(i) + UASTAR(i))  ...    
                    %          - (L*RC(i) +    VTC(i) + 2*UTSTAR(i)) * UTSTAR(i) ... 
                    %          + ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * (1/(2*pi))*Cdp(i)*CpoDp(i)*dVdG(i,i)*(L*RC(i) + VTC(i) + UTSTAR(i)) ...
                    %          + ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * (1/(2*pi))*Cdp(i)*CpoDp(i)*VSTAR(i) * UTHIF(i,i);  
                    % ------  
                    % Inviscid terms:    
                    RNS(i) =   ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * (VAC(i) + UADUCT(i) + UASTAR(i))  ...    
                             - (L*RC(i) +    VTC(i) + 2*UTSTAR(i)) * UTSTAR(i);                    
                    % ------   
                end

                RNS(i+  Mp) = UASTAR(i) - UASTARtemp(i);
                RNS(i+2*Mp) = UTSTAR(i) - UTSTARtemp(i);
                RNS(i+3*Mp) = TANBIC(i) - (VAC(i) + UADUCT(i) + UASTAR(i))/(L*RC(i) + VTC(i) + UTSTAR(i));            

            end

            % ----------------------------------------------- Evaluate Jacobian        
            for i = 1:Mp 

                if (Chord_flag == 1) && (strcmp(ChordMethod,'CLmax') == 1)
                    % ==========
                    % 11/16/2011 BEPPS: New drag treatment...   
                    JNS(i     ,i+  Mp) =  3*(VAC(i) + UADUCT(i)) + 4*UASTAR(i)   + 2 * (L*RC(i) +    VTC(i) + 2*UTSTAR(i)) * Cdp(i)/CLmax(i);                              

                    JNS(i     ,i+2*Mp) = - (L*RC(i) +    VTC(i)  + 4*UTSTAR(i))  + 2 *  (VAC(i) + UADUCT(i) + 2*UASTAR(i)) * Cdp(i)/CLmax(i);               
                    % ==========
                else                  
                    % 11/16/2011 BEPPS: This was the implementation in EppsOptimizer06.m...
                    %                  
                    % JNS(i     ,i+  Mp) =    3*(VAC(i) + UADUCT(i)) + 4*UASTAR(i) ...
                    %                      + (1/pi)*Cdp(i)*CpoDp(i)*dVdG(i,i)*(L*RC(i) + VTC(i) + UTSTAR(i)) ...
                    %                      + (1/pi)*Cdp(i)*CpoDp(i)*VSTAR(i) * UTHIF(i,i); 
                    % 
                    % JNS(i     ,i+2*Mp) = - (L*RC(i) +    VTC(i) + 4*UTSTAR(i)) ... 
                    %                      + ( VAC(i) + UADUCT(i) + 2*UASTAR(i))*(1/(2*pi))*Cdp(i)*CpoDp(i)*dVdG(i,i);                 
                    % ------ 
                    % Inviscid terms:           
                    JNS(i     ,i+  Mp) =  3*(VAC(i) + UADUCT(i)) + 4*UASTAR(i);                              

                    JNS(i     ,i+2*Mp) = - (L*RC(i) +    VTC(i)  + 4*UTSTAR(i));               
                    % ------
                end

                JNS(i+  Mp,1:Mp  ) = - UAHIF(i,1:Mp);
                JNS(i+  Mp,i+  Mp) = 1;
                JNS(i+2*Mp,1:Mp  ) = - UTHIF(i,1:Mp);
                JNS(i+2*Mp,i+2*Mp) = 1;

                JNS(i+3*Mp,i+  Mp) = -1       /(L*RC(i) + VTC(i) + UTSTAR(i));
                JNS(i+3*Mp,i+2*Mp) = TANBIC(i)/(L*RC(i) + VTC(i) + UTSTAR(i));
                JNS(i+3*Mp,i+3*Mp) = 1;
            end

            % ----------------------------- Update Newton solver vector of unknowns
            DX     = linsolve(JNS,-RNS);

            G      = G      + relax*DX( 1:Mp      ) ;
            UASTAR = UASTAR + relax*DX((1:Mp)+  Mp)';
            UTSTAR = UTSTAR + relax*DX((1:Mp)+2*Mp)';
            TANBIC = TANBIC + relax*DX((1:Mp)+3*Mp)'; 

            % Smooth the inflow angle for numerical stability:
            TANBICsmooth = TANBIC * Bsmooth;
            % -----------------------------------------------------------------

            % -----------------------------------------------------------------
            if any(isnan(DX))  || ~isreal(DX)
                 G      = 0*RC';
                 UASTAR = 0*RC;
                 UTSTAR = 0*RC;
                 TANBIC = TANBC;

                 disp(' ')
                 disp('<WARNING>')
                 disp('<WARNING> DX == NaN or imaginary... crash avoided...')
                 disp('<WARNING>')
                 disp(' ')
                 G_iter = 999;
            end
            % ----------------------------------------------------------------- 

            % ----------------------------------------------------------------- 
            if Duct_flag == 1            
                % Update the influence of propeller on duct
                for m = 1:Mp;
                    DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBICsmooth(m);
                    DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBICsmooth(m);
                end

                % Update induced velocities at the duct (influence of propeller on duct)
                UARING = (DAHIF*G)';  
                URRING = (DRHIF*G)'; 

                % Update duct circulation such that CTD == CTDdes            
                [junk,GdNEW] = Duct_Thrust(XdRING,RdoRp,VARING,UARING,URRING,GdRING,Gd,CDd,CTDdes);

                Gd = 0.5*GdNEW + (1-0.5)*Gd;

                % Update the induced velocities at the propeller (influence of duct on propeller)
                UADUCT =  UADIF*Gd;                        
            end
            % ----------------------------------------------------------------- 

            % -----------------------------------------------------------------
            % END: "Robust" method (EppsOptimizer06.m)
            % -!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!- 


            % % -----------------------------------------------------------------        
            % % "SIMPLE (INCORRECT) OPTIMIZER"        
            % % -----------------------------------------------------------------
            % % ---------------------- Set up simultaneous equations for G and LM
            % A  = zeros(Mp,Mp);         % A matrix for linear system of equations
            % B  = zeros(Mp,1);          % B matrix for linear system of equations
            % 
            % for i = 1:Mp                           % for each equation for G(i)
            %     for m = 1:Mp                       % for each vortex panel, m        
            %         A(i,m) =  UAHIF(m,i)*RC(m)*DR(m)    ...                
            %                 + UAHIF(i,m)*RC(i)*DR(i); 
            %     end   
            % 
            %     B(i)  = -VAC(i)*RC(i)*DR(i) ...
            %             -(1/(2*pi))*sum(Cdp.*dVdG(:,i)'.*CpoDp.*(L*RC + VTC + UTSTAR).*RC.*DR) ...
            %             -(1/(2*pi))*sum(Cdp.*VSTAR.*CpoDp.*UTHIF(:,i)'.*RC.*DR);             
            % end
            % 
            % % -------------------------------- Solve linear system of equations
            % GLM = linsolve(A,B);             
            % G   = GLM(1:Mp);                     % G is the first Mp entries
            % 
            % % Update induced velocities (influence of propeller on propeller)
            % UASTAR = (UAHIF*G)';  UTSTAR = (UTHIF*G)';  
            % 
            % % -------------------------  Repair outlier {UASTAR, UTSTAR, URSTAR} values
            % UASTAR = RepairSpline(RC,UASTAR,'UASTAR',Plot_flag*Hvel);
            % UTSTAR = RepairSpline(RC,UTSTAR,'UTSTAR',Plot_flag*Hvel);
            %
            % % -------------------------------------------- Update TANBIC
            % TANBIC = (VAC + UASTAR)./(L*RC + VTC + UTSTAR);
            %          
            % % -----------------------------------------------------------------        
            % % END "SIMPLE (INCORRECT) OPTIMIZER"        
            % % -----------------------------------------------------------------           
        end 
        % ---------------------------------------------------------------------
        % END UPDATE: G, UASTAR, UTSTAR, TANBIC, (duct stuff), LM
        % =@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=
        %% 更新其他变量
        % ------------------------------------ Update VSTAR and its derivatives
        VSTAR  = sqrt((VAC+UADUCT+UASTAR).^2 + (L*RC+VTC+UTSTAR).^2);       
        % --------------------- Update the vortex Horseshoe Influence Functions
        % 11/17/2011 BEPPS: If you use TANBICsmooth here, then the circulation and 
        %                   induced velocities in the turbine case will be wavy!
        [UAHIF,UTHIF] = Horseshoe(Mp,Z,TANBIC,RC,RV,Hub_flag,RhoRp,Duct_flag,RdoRp);     
        % ---------------------------------------------------------------------    


        % ------------------------------------------- Update chord distribution
        if (Chord_flag == 1) && (G_iter <= ITER)

            if strcmp(ChordMethod,'CLmax') == 1
                % -----------------------------------------------------------------  
                % Method: Choose chord length based on CLmax
                % -----------------------------------------------------------------        
                CpoDp = 2*pi*G'./(VSTAR.*CLmax); % scale CpoDp to keep CL == CLmax
                % ----------------------------------------------------------------- 

            elseif strcmp(ChordMethod,'ConeyPLL') == 1
                % -----------------------------------------------------------------
                % Method: (Coney, 1989) cavitation method -- ASSUMES GIVEN THICKNESS DISTRIBUTION t0oDp      
                % -----------------------------------------------------------------   
                SIGMA = SIGMAs./VSTAR.^2;    % local cavitation number

                f0oD = (2*pi*G'./ VSTAR) .* f0octilde ./ CLItilde;

                CpoDp  = (8.09*f0oD+3.033*t0oDp)./(2*SIGMA) + sqrt((8.09*f0oD+3.033*t0oDp).^2 + 4* SIGMA .* (26.67*f0oD.^2 + 10*f0oD.*t0oDp) )./(2*SIGMA);
                % ----------------------------------------------------------------- 


            elseif strcmp(ChordMethod,'FAST2011dCTP') == 1
                % -----------------------------------------------------------------  
                % Method: see (Epps et al., FAST'2011)
                % -----------------------------------------------------------------
                % Every third iteration, update the chord and thickness.
                if mod(G_iter,3) == 0

                   [CpoDp, t0oc, C_res] = Chord_FAST2011_dCTP(SIGMAh,CT_desired_high, Jh,CpoDp,t0oc,...      
                                         Propeller_flag,Viscous_flag,Hub_flag,Duct_flag, ...
                                         Z,Mp,ITER,Rhv,RC,RV,DR,RhoRp,VMIV,CTD,...
                                         L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,Cdp,...
                                         Rp,rho,VS,N);            
                end

            elseif strcmp(ChordMethod,'FAST2011dVAC') == 1
                % -----------------------------------------------------------------  
                % Method: see (Epps et al., FAST'2011)
                % -----------------------------------------------------------------
                % Every third iteration, update the chord and thickness.
                if mod(G_iter,3) == 0

                   [CpoDp, t0oc, C_res] = Chord_FAST2011_dVAC(dVAC,SIGMAs,L,RC,VAC,VTC,UADUCT,UASTAR,UTSTAR,G,CpoDp,t0oc,...
                                               DR,Cdp,Z,Js,VMIV,Hub_flag,RhoRp,Rhv,CTD,Mp,Rp,rho,VS,N);
                end

            elseif strcmp(ChordMethod,'Brizzolara2007') == 1
                % -----------------------------------------------------------------  
                % Method: see (Brizzolara et al., 2007)
                % -----------------------------------------------------------------
                % Every third iteration, update the chord and thickness.
                if mod(G_iter,3) == 0

    %                 save temp
    %                 return

                   [CT,CQ,CP,KT,KQ, CTH,TpoT, Ja,Jw,VMWV, EFFYo, EFFY,ADEFFY,QF] = ...
                        Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,Cdp,CpoDp,G,Z,Js,VMIV,Hub_flag,RhoRp,Rhv,CTD);

                   [CpoDp, t0oc, C_res] = Chord_Brizzolara(rho,n,Dp,VS,H,Mp,Z,KT,KQ, RC,G,VSTAR,TANBIC,CpoDp,t0oc);
                end
            end

            % -----------------------------------------------------------------
            % Scale CpoDp to give specified Expanded Area Ratio (EAR == EARspec)
            %如果伸张面积比有特殊要求，进行比例缩放
            if EARspec > 0
               EAR = (2*Z/pi) * trapz(   linspace(RhoRp,1,100), interp1(RC,CpoDp, linspace(RhoRp,1,100), 'spline','extrap')   );  
               CpoDp = (EARspec/EAR) * CpoDp; 
            end
            % -----------------------------------------------------------------

            % -----------------------------------------------------------------
            if all(CpoDp == 0) || any(isnan(CpoDp))
                 G      = 0*RC';
                 UASTAR = 0*RC;
                 UTSTAR = 0*RC;
                 TANBIC = TANBC;
                 CpoDp    = 0*RC;

                 disp(' ')
                 disp('<WARNING>')
                 disp('<WARNING> CpoDp == NaN or zero... crash avoided...')
                 disp('<WARNING>')
                 disp(' ')
                 G_iter = 999;
            end
        end
        % ---------------------------------- Prepare for the next iteration
        G_res   = abs((G - G_last)./G);    % residual G
        G_last  = G;                       % the last value of G
        Gd_res  = abs((Gd - Gd_last)/Gd);  % residual Gd
        Gd_last = Gd;                      % the last value of Gd
        LM_last = LM;                      % last value of the Lagrange Multiplier    

        if G_iter < 10
            disp(['The max  G_res for iteration  ',num2str(G_iter),' is: ',num2str(max(G_res))]),  
        else
            disp(['The max  G_res for iteration ',num2str(G_iter),' is: ',num2str(max(G_res))]),  
        end    

    %     if Duct_flag == 1
    %         disp(['The     Gd_res for iteration ',num2str(G_iter),' is: ',num2str(Gd_res)]),        
    %     end
        % ---------------------------------------------------------------------

        % ---------------------------------------------------------------------
        if (Chord_flag == 1)  &&  ( (strcmp(ChordMethod,'FAST2011dCTP') == 1) || (strcmp(ChordMethod,'FAST2011dVAC') == 1) || (strcmp(ChordMethod,'Brizzolara2007') == 1) )
            if mod(G_iter,3) ~= 0
                G_res = 1;
            else
                G_res = max([G_res;C_res/10]);
            end
        end
        % ---------------------------------------------------------------------

        % ---------------------------------------------------------------------
        G_iter  = G_iter + 1;              % iteration in the G loop
        % ---------------------------------------------------------------------   
    end
    %% 进行优化结果说明
    if G_iter > ITER
        disp(' '),
        disp('警告：结果未收敛'),
        Converge_flag = 0;
    else
        disp(' '),
        disp('完成优化，结果收敛'),
        Converge_flag = 1;
    end
    %%

    % =========================================================================
    % 11/17/2011 BEPPS: Note that the OpenProp 3.2.0 version of Unload_Blade.m 
    %                   does not do a very good job of scaling the circulation
    %                   such that CT == CTdes.  This code should be improved.
    %
    % If required, unload the hub and tip, then rescale the circulation
    % distribution to get the desired value of the thrust coefficient.
    if Hub_flag && (HUF > 0 || TUF > 0)                       % (IF STATEMENT U1)

        if Duct_flag == 0
            DAHIF_times_TANBIC  = 0;
            DRHIF_times_TANBIC 	= 0;
            XdRING              = 0;
            RdoRp            = 1;
            VARING              = 0;
            GdRING              = 0;
            Gd                  = 0;
            CDd                 = 0;
            CTDdes              = 0;
        end

        [G,UASTAR,UTSTAR,TANBIC,UARING,URRING,Gd,UADUCT] = ...
                                        Unload_Blade(HUF,TUF,RC,RhoRp, G,  VAC,VTC, TANBIC,RV,DR,L,Mp,Z, ...
                                                     Hub_flag,ITER,Bsmooth,Cdp,CpoDp,Js,VMIV,Rhv,CTPdes,...
                                                     Duct_flag,UADUCT,XdRING,RdoRp,VARING,GdRING,UADIF,Gd,CDd,CTDdes,...
                                                     DAHIF_times_TANBIC,DRHIF_times_TANBIC,...
                                                     Plot_flag,Hgamma,HHgamma,Hvel,HHvel,Hbeta,HHbeta);    

    end
    % =========================================================================

    % if Rroot > Rhub 
    %     % --------------------------------------------------------------------- 
    %     % Develop code for turbine case with circular blade sections near the root.
    %     % ---------------------------------------------------------------------
    % end


    % --------------------- Compute thrust & torque coefficients and efficiency
    % Compute the actual CT for the duct with the current circulation Gd
    if Duct_flag == 1
        [CTD,junk] = Duct_Thrust(XdRING,RdoRp,VARING,UARING,URRING,GdRING,Gd,CDd,CTDdes);
    end

    [CT,CQ,CP,KT,KQ, CTH,TpoT, Ja,Jw,VMWV, EFFYo, EFFY,ADEFFY,QF,QFo,QFw] = ...
        Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,Cdp,CpoDp,G,Z,Js,VMIV,Hub_flag,RhoRp,Rhv,CTD);


    % -------------------------------------------------------------------------
    CL   = 2*pi*G'./(VSTAR.*CpoDp);     % lift coefficient

    % Expanded Area Ratio
    EAR = (2*Z/pi) * trapz(   linspace(RhoRp,1,100), interp1(RC,CpoDp, linspace(RhoRp,1,100), 'spline','extrap')   );  

    % -------------------------------------------------------------------------


    % ------------------------------------------- Update chord distribution
    if (Chord_flag == 1) && ( (strcmp(ChordMethod,'FAST2011dCTP') == 1) || (strcmp(ChordMethod,'FAST2011dVAC') == 1) || (strcmp(ChordMethod,'Brizzolara2007') == 1) )
        % Use with chord optimization methods that optimize t0oc
        t0oDp = t0oc.*CpoDp;
    else   
        % Use when no chord optimization or with chord optimization methods 
        % that do not optimize t0oc (e.g. CLmax or Coney use given t0oDp)
        t0oc = t0oDp ./ CpoDp;    
    end
    % -------------------------------------------------------------------------



    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
        disp(' '),
        disp('Forces after circulation optimization:')

    if Propeller_flag == 1

        disp(['    Js = ',num2str(Js)]),    
        disp(['    KT = ',num2str(KT)]),   
        disp(['    KQ = ',num2str(KQ)]),   
        disp(['    CT = ',num2str(CT)]),
        if abs(VMIV - 1) > 1e-8     % i.e. if VMIV is not equal to 1 
        disp(['    Ja = ',num2str(Ja)]),
        end
        disp(['  EFFY = ',num2str(EFFY)]),
        disp(['ADEFFY = ',num2str(ADEFFY)]),    
        disp(['    QF = ',num2str(QF)]),    

    else
        disp(['L      =  ',num2str(L)]),    
        disp(['CP     = ' ,num2str(CP)]),  
        disp(['CPBetz = ' ,num2str(CPBetz)]),  
        disp(['QF     =  ',num2str(CP/CPBetz)]),      
    end
        disp(' ')     
        disp(' ')     
    %% 保存计算结果
    design.part1      = '------ Section properties, size (1,Mp) ------';
    design.RC         = RC;                 % [1 x Mp] control point radii
    design.DR         = DR;                 % [1 x Mp] difference in vortex point radii
    design.G          = G';                 % [1 x Mp] circulation distribution
    design.VAC        = VAC;                % [1 x Mp] 
    design.VTC        = VTC;                % [1 x Mp] 
    design.UASTAR     = UASTAR;             % [1 x Mp] 
    design.UTSTAR     = UTSTAR;             % [1 x Mp] 
    design.VSTAR      = VSTAR;              % [1 x Mp]  
    design.TANBC      = TANBC;              % [1 x Mp] 
    design.TANBIC     = TANBIC;             % [1 x Mp] 
    design.CL         = CL;                 % [1 x Mp] 
    design.Cdp         = Cdp;                 % [1 x Mp] 
    design.CpoDp        = CpoDp;                % [1 x Mp]
    design.t0oc       = t0oc;               % [1 x Mp] 
    design.t0oDp       = t0oDp;               % [1 x Mp] 

    design.part2      = '------ Other properties  ------';
    design.converged  = Converge_flag;
    design.iteration  = G_iter;
    design.RV         = RV;                 % [1 x Mp+1] vortex point radii
    design.Rhub_oR    = RhoRp;            % [1 x 1]
    design.EAR        = EAR;
    design.LM         = LM;                 % [1 x 1]
    design.VMIV       = VMIV;               % [1 x 1]
    design.VMWV       = VMWV;               % [1 x 1]
    design.SIGMAs     = SIGMAs;


    if Duct_flag == 1  
        design.part3      = '------ Duct parameters ------';
        design.Rduct_oR   = RdoRp;           % [1 x 1]
        design.Cduct_oR   = CdoRp;           % [1 x 1]
        design.Xduct_oR   = Xduct_oR;           % [1 x 1]
        design.Gd         = Gd;                 % [1 x 1], duct circulation   
        design.VARING     = VARING;             % [1 x  1]

        design.XdRING     = XdRING;             % [1 x Nd], Nd=12 duct vortex rings 
        design.UARING     = UARING;             % [1 x Nd], Nd=12 duct vortex rings
        design.URRING     = URRING;             % [1 x Nd], Nd=12 duct vortex rings
        design.GdRING     = GdRING;             % [1 x Nd], Nd=12 duct vortex rings

        design.DAHIFtT    = DAHIF_times_TANBIC; % [Nd,Mp], Duct Horseshoe Influence Functions (influence of propeller on duct)
        design.DRHIFtT    = DRHIF_times_TANBIC; % [Nd,Mp], Duct Horseshoe Influence Functions (influence of propeller on duct)



        design.UADIF      = UADIF;              % [1 x Mp] 
        design.UADUCT     = UADUCT;             % [1 x Mp] 
        design.CTPdes     = CTPdes;             % [1 x 1], desired propeller CT
        design.CTDdes     = CTDdes;             % [1 x 1], desired duct      CT

        design.TAU        = TpoT;                  % [1 x 1]
        design.CTD        = CTD;                  % [1 x 1]
        design.part4      = '------ Performance metrics ------';
    else
        design.part3      = '------ Performance metrics ------';
    end


    if Propeller_flag == 1
        design.L      = L;
        design.Js     = Js;
        design.KT     = KT;                   % [1 x 1]
        design.KQ     = KQ;                   % [1 x 1]
        design.CT     = CT;                   % [1 x 1]
        design.CQ     = CQ;                   % [1 x 1]
        design.CP     = CP;                   % [1 x 1]
        design.CTH    = CTH;                  % [1 x 1]

        if abs(VMIV - 1) > 1e-8   % i.e. if VMIV is not equal to 1 
        design.EFFYo  = EFFYo;                 % [1 x 1]
        design.Ja     = Ja;
        end
        design.EFFY   = EFFY;                 % [1 x 1],  EFFY = EFFYa by convention (Kerwin)
        design.ADEFFY = ADEFFY;               % [1 x 1], prop.   actuator disk
        design.QF     = QF;                   % [1 x 1]

        if (VMIV < 0.05)   % assume design for bollard pull
            design.QFo    = QFo;                  % [1 x 1]
            design.QFw    = QFw;                  % [1 x 1]
        end
    else    
        design.L      = L;
        design.Js     = Js;
    %     design.KT     = KT;                   % [1 x 1]
    %     design.KQ     = KQ;                   % [1 x 1] 
        design.CT     = CT;                   % [1 x 1]
        design.CQ     = CQ;                   % [1 x 1]
        design.CP     = CP;                   % [1 x 1]
        design.CPBetz = CPBetz;               % [1 x 1], turbine actuator disk
        design.QF     = CP/CPBetz;
    end
end