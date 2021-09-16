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
    %% ��pt�е�������ֵ
    % part1 �����������
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
    
    % part2 ��ҶƬ���������
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
    
    % part3 ���ⲿ����������
    Ni = input.Ni;
    rho = input.rho;
    ri = input.ri;
    VA_i = input.VA_i;
    VT_i = input.VT_i;
    % ͨ�����ʹ�������ٶ�/װ���н��ٶȡ����ݸ��ӹ⻬
    VA_i = RepairSpline(ri,VA_i);  
    VT_i = RepairSpline(ri,VT_i);
    
    % part4 ��������ز�����
    TdoT = input.TdoT;
    CDd = input.CDd;
    
    % part5 �����������빤�ߡ�
    Analyze_flag = input.Analyze_flag;
    Geometry_flag = input.Geometry_flag;
    Coordinate_flag = input.Coordinate_flag;
    Printing_flag = input.Printing_flag;
    Mp = input.Mp;
    Np = input.Np;
    Nd = input.Nd;
    ITER = input.ITER;
    
    % �������
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
    % ������û��ָ��Ť�أ�����ָ��������ָ��Ť����Ϊ1
    TorqueSpec_flag = 0;
    % ָ��չ�������(EAR)
    EARspec = 0;
    
    %% �������
    RhoRp = Dh/Dp;
    
    RdoRp = Dd/Dp;
    CdoRp = 2*Cd/Dp;
    % XdoRp��ӳ���������е㴦����ԭ��ľ���
    XdoRp = 2*(Ld-Cd/2)/Dp;
    
    TpoT = 1-TdoT;
    
    SIGMAs = (Patm+rho*g*H-Pv)/(0.5*rho*VS^2);
    
    Js = VS/((N/60)*Dp);
    L = pi/Js;
    KT = T/(rho*(N/60)^2*Dp^4);
    CT = T/(0.5*rho*VS^2*pi*(Dp/2)^2);
    
    % CT_desiredΪ��Ҫ���������ϵ����CT_desired = CTp_desired+CTd_desired
    CT_desired = CT;
    % ҶƬ���ṩ������ϵ��
    CTp_desired = CT_desired*TpoT;
    % �������ṩ������ϵ��
    CTd_desired = CT_desired*TdoT;
    
    % ��ȡ���湰�Ⱥ���������ϵ��������
    F0oCptilde_x = zeros(Nx,1);
    CLItilde_x = zeros(Nx,1);
    for index = 1 : Nx
        [F0oCptilde_x(index),CLItilde_x(index)] = GeometryFoil2D(Meanline_x{index},...
                                                                Thickness_x{index});
    end
    
    % �����Ҿ��ȵĲ�ֵ����
    CpoDp_x = InterpolateChord(rx_def,CpoDp_x_def,rx);
    
    %% ������Ƶ���и��ľ���λ�ò����в�ֵ
    % ����[RhoRp,1]֮��100���Ⱦ����ɵ�����
    ri_expanded = linspace(RhoRp,1,100);
    VA_i_expanded = pchip(ri,VA_i,ri_expanded);
    % ����ƽ���ٶȣ����Բο�Kerwin����p.138��eqn 163
    VMIV = 2/(1^2-RhoRp^2)*trapz(ri_expanded,ri_expanded.*VA_i_expanded);
    
    % �����и��Ϳ��Ƶ�ľ���λ�� rc:[Mp*1];rv:[Mp+1*1]
    [rc,rv,drv] = LLPanelRadii(Mp,RhoRp,Hub_flag,Duct_flag);
    
    % ����Ӧ����λ�ý��в�ֵ
    VA_c = pchip(ri,VA_i,rc);
    VT_c = pchip(ri,VT_i,rc);
    CDp_c = pchip(rx,CDp_x,rc);
    T0oDp_c = pchip(rx,T0oDp_x,rc);
    CL_c = pchip(rx,CL_x,rc);
    F0oCdtilde_c = pchip(rx,F0oCptilde_x,rc);
    CLItilde_c = pchip(rx,CLItilde_x,rc);
    
    if (abs(rx(end)-1) < 1e-4) && (CpoDp_x(end) <= 0.01)
        % Ҷ�Ҵ�Ҫ���Ⱥ��ҳ�Ϊ0ʱ�������з���
        CpoDp_c = InterpolateChord(rx,CpoDp_x,rc);
    else
        CpoDp_c = pchip(rx,CpoDp_x,rc);
    end
    T0oCp_c = T0oDp_c./CpoDp_c;
    
    %% ���õ�����ʼֵ
    % ���ù�һ��������ʼֵ G = Gamma/2*pi*Rp*VS
    Gp = zeros(Mp,1);
    
    % ��������������յ��ٶȳ�ʼֵ Kerwin & Hadler (2010), eqn (4.26)
    UASTAR = 0*rc+0.5*(sqrt(1+CTp_desired)-1);
    UTSTAR = 0*rc;
    % ��Χ������������� (Epps, 2013)
    if VMIV < 0.05
        UASTAR = 0*rc + 0.5*sqrt(CTp_desired);
    end
    
    % �����յ��ٶȳ�ʼֵ������Ǻ�ˮ���������ǳ�ʼֵ������ֵ
    tanBeta_c = VA_c./(L*rc+VT_c);
    tanBetaI_c = (VA_c+UASTAR)./(L*rc+VT_c+UTSTAR);
    
    % �����������ٶȵĳ�ʼֵ
    VSTAR = sqrt((VA_c+UASTAR).^2+(L*rc+VT_c+UTSTAR).^2);
    
    % �����Ҿ��ȵĳ�ʼֵ
    CpoDp_c = 0.1*ones(Mp,1);
    
    % ��ʼ��ҶƬ�ϵ��յ��ٶ�
    [UAHIF,UTHIF] = Horseshoe(Mp,Z,tanBetaI_c,rc,rv,Hub_flag,...
                              RhoRp,Duct_flag,RdoRp);
    
    % ��ʼ��������ز���
    if Duct_flag == 1
        % ���к����Ժ��������������յ��ٶȵĳ�ʼ��
        [XdRING,GdRING,UADIF] = Duct_Influence(RdoRp,CdoRp,XdoRp,rc,Nd);
        
        % ����pchip��ֵ�ķ���ȷ���������ߴ������������ٶ�
        VARING = pchip(ri,VA_i,RdoRp);  

        % ���ù�һ�������ܻ����ĳ�ʼֵ
        Gd = 0;
        
        % �������������л��������������յ��ٶ�
        UADUCT = UADIF * Gd;
        
        % ����ҶƬ�Ժ��������յ��ٶȵĳ�ʼ��
        [DAHIF_times_TANBIC,~,DRHIF_times_TANBIC] = Horseshoe_intr_110830(XdRING,RdoRp, rc,ones(size(tanBetaI_c)),rv,Z,Hub_flag,RhoRp,Duct_flag,RdoRp); 
        % DAHIF(n,m)��λǿ�Ȼ�����ҶƬ�ϵ�m���и�Ժ����ϵ�n�����Ƶ�������յ��ٶ�
        DAHIF = zeros(Nd,Mp);
        DRHIF = zeros(Nd,Mp);
        for m = 1 : Mp
            % (:,m)����ĵ�m��
            DAHIF(:,m) = DAHIF_times_TANBIC(:,m)/tanBetaI_c(m);
            DRHIF(:,m) = DRHIF_times_TANBIC(:,m)/tanBetaI_c(m);
        end  
    else
        Gd = 1;  % if set Gd==0, then Gd_res would be Inf  
        UADUCT = zeros(Mp,1);
    end
    
    % ƽ������ X_smooth = Bsmooth*X;
    Bsmooth = RepairSplineMatrix(rc);
    
    %% ���л����Ż�
    % ��ʼ����������
    LM = -1;            % �������ճ��ӣ�������������
    LM_last = LM;       % �������ճ�����һ�εĵ���ֵ
    Gp_last = 0*Gp;     % ������������һ�εĵ���ֵ
    Gd_last = 0;        % ����������һ�εĵ���ֵ
    G_iter = 1;         % ������������
    Gp_res = 1;         % �������������ε���֮��Ĳ�ֵ
    Gd_res = 0;         % �����������ε���֮��Ĳ�ֵ
    C_res = 0;          % �ҳ����ε���֮��Ĳ�ֵ  
    G_TOL = 1e-4;       % �����ﵽ�������������в�ֵ
    
    % ѭ���з������ֱ�����ã�
    % Gp LM
    % ���ݷ��̽����֪����һ����������ã�
    % UASTAR UTSTAR tanBetaI_c DAHIF DRHIF UARING URRING Gd UADUCT
    % VSTAR UAHIF UTHIF
    while G_iter <= ITER && any([Gp_res;Gd_res] > G_TOL)
        % AΪ����{��(1),��(2),...,��(Mp),��}��ϵ����Ҳ���Ƿ�����Ⱥ�����ϵ������
        A = zeros(Mp+1,Mp+1);
        % BΪ������Ⱥ��Ҳ�ĳ���������
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
            % Ҫ������������������
            for m = 1:Mp                          
                A(Mp+1,m) = (L*rc(m)+VT_c(m)+UTSTAR(m))*drv(m);  
            end

            B(Mp+1) = CTp_desired/(4*Z)+(1/(2*pi))*...
                      sum(CDp_c.*VSTAR.*CpoDp_c.*(VA_c+UADUCT+UASTAR).*drv);
            if Hub_flag == 1
                % ���뽰�����
                B(Mp+1) = B(Mp+1) + (Z/8)*(log(1/Rhv)+3)*(Gp_last(1)^2);
            end
        else
            % Ҫ��Ť����������
            for m = 1:Mp
                A(Mp+1,m) = (VA_c(m)+UADUCT(m)+UASTAR(m))*rc(m)*drv(m);
            end
            B(Mp+1) = CQdes/(4*Z)-(1/(2*pi))*...
                      sum(CDp_c.*VSTAR.*CpoDp_c.*(L*rc + VT_c + UTSTAR).*rc.*drv);
        end

        % ������Է�����
        GLM = linsolve(A,B);
        % ���ǰMp��Ϊ��ǰ���ε�����õĻ���ֵ
        Gp = GLM(1:Mp);
        % ��ĵ�Mp+1��Ϊ���ε�����õ��������ճ���
        LM = GLM(Mp+1);

        % �����յ��ٶ�
        UASTAR = UAHIF*Gp;  
        UTSTAR = UTHIF*Gp;  

        % ����ˮ���������ǵ�����ֵ
        tanBetaI_c = (VA_c+UADUCT+UASTAR)./(L*rc+VT_c+UTSTAR);
        % �Ż�ˮ����������
        tanBetaI_c_smooth = Bsmooth*tanBetaI_c;

        if any(isnan(GLM))||~isreal(GLM)||max(Gp-Gp_last) > 10
            % ����Ч���ʵ�������ͻ��ʱ���о���
            Gp = zeros(Mp,1);
            UASTAR = zeros(Mp,1);
            UTSTAR = zeros(Mp,1);
            tanBetaI_c = tanBeta_c;
            % ���������
            error = sprintf('ʹ���������ճ��ӷ����л����Ż�ʱ���������շ������޽�����Ч');
            uialert(Fig_Main,error,'���ʧ��');
        end   

        if Duct_flag == 1
            % ����ҶƬ�Ժ������յ��ٶ�
            for m = 1:Mp
                DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / tanBetaI_c_smooth(m);
                DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / tanBetaI_c_smooth(m);
            end
            UARING = DAHIF*Gp;  
            URRING = DRHIF*Gp; 

            % ���º����ϵĻ���
            [~,Gd] = Duct_Thrust(XdRING,RdoRp,VARING,UARING,URRING,...
                                 GdRING,Gd,CDd,CTd_desired);

            % ���º����ں��������������������յ��ٶ�
            UADUCT = UADIF*Gd;                        
        end
        
        % �����������ٶ�
        VSTAR = sqrt((VA_c+UADUCT+UASTAR).^2+(L*rc+VT_c+UTSTAR).^2);       
        
        % ����ҶƬ�ڲ�Ӱ��
        [UAHIF,UTHIF] = Horseshoe(Mp,Z,tanBetaI_c,rc,rv,Hub_flag,...
                                  RhoRp,Duct_flag,RdoRp);     
        
        % �����ҳ��ֲ�
        if strcmp(ChordMethod,'CLmax')
            % ֱ��������������ϵ���ͺ񾶱ȣ��ɴ˵õ��Ҿ���
            CpoDp_c = 2*pi*Gp./(VSTAR.*CL_c);
            
        elseif strcmp(ChordMethod,'ConeyPLL')
            % ֱ���������ȡ�񾶱ȣ����ݼ���������ϵ��
            SIGMA = SIGMAs./VSTAR.^2;    % local cavitation number
            F0oDp = (2*pi*Gp./VSTAR).*F0oCdtilde_c./CLItilde_c;
            CpoDp_c = (8.09*F0oDp+3.033*T0oDp_c)./(2*SIGMA)+...
                       sqrt((8.09*F0oDp+3.033*T0oDp_c).^2+...
                       4*SIGMA.*(26.67*F0oDp.^2+10*F0oDp.*T0oDp_c))./(2*SIGMA);
                   
        end

        % ������������������Ҫ�󣬽��б�������
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
             % ���������
             error = sprintf(['ʹ��',ChordMethod,...
                             '�����񾶱�ʧ�� \n ���л�����������']);
             uialert(Fig_Main,error,'���ʧ��');
             G_iter = 999;
        end
        
        % ���µ�������
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
        
        % ��������+1
        G_iter = G_iter+1;
    end
    
    % �����Ż����˵��
    if G_iter > ITER
        disp('���棺���δ����');
        % ���������
        error = sprintf('�����δ����');
        uialert(Fig_Main,error,'���ʧ��');
        Converge_flag = 0;
    else
        disp('����Ż����������');
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
    
    % ���������ٲ����͸���ϵ��
    % �����Ż���ĺ����������㺭��ʵ����ռ������ϵ��
    if Duct_flag == 1
        [CTd_desired,~] = Duct_Thrust(XdRING,RdoRp,VARING,UARING,URRING,...
                                      GdRING,Gd,CDd,CTd_desired);
    end

    [CT,CQ,CP,KT,KQ,CTH,TpoT,Ja,~,VMWV,EFFYo,EFFY,ADEFFY,QF,QFo,QFw] = ...
        Forces(rc,drv,VA_c,VT_c,UASTAR,UTSTAR,UADUCT,CDp_c,CpoDp_c,Gp,Z,Js,...
               VMIV,Hub_flag,RhoRp,Rhv,CTd_desired);
    
    % �����Ż���Ļ������ٶȼ���ʵ�ʵ�����ϵ��
    CL_c = 2*pi*Gp./(VSTAR.*CpoDp_c);

    % Expanded Area Ratio
    EAR = (2*Z/pi)*trapz(linspace(RhoRp,1,100),...
          interp1(rc,CpoDp_c,linspace(RhoRp,1,100),'spline','extrap'));  
    
    % ���º��ұ�
    T0oCp_c = T0oDp_c ./ CpoDp_c;    
    
    % �������ĵ�Ť��
    Q = CQ*0.5*rho*VS^2*pi*Dp^2/4*Dp/2;
    
    % �������������ĵĹ���
    P = Q*2*pi*N/60;
    
    %% ��ӡ�ͱ��������
    design.part1 = 'ҶƬ��صļ�����';
    design.rc = rc;                         % [Mp*1] ���Ƶ�����
    design.rv = rv;                         % [Mp+1*1] �и������
    design.drv = drv;                       % [Mp*1] �������и�������ֵ
    design.Gp = Gp;                         % [Mp*1] ���Ƶ�/�и���
    design.VA_c = VA_c;                     % [Mp*1] ���������ٶ� 
    design.VT_c = VT_c;                     % [Mp*1] ���������ٶ� 
    design.UASTAR = UASTAR;                 % [Mp*1] ���Ƶ㴦�������յ��ٶ� 
    design.UTSTAR = UTSTAR;                 % [Mp*1] ���Ƶ㴦�������յ��ٶ� 
    design.VSTAR = VSTAR;                   % [Mp*1] ����������ٶ�  
    design.tanBeta_c = tanBeta_c;           % [Mp*1] ����
    design.tanBetaI_c = tanBetaI_c;         % [Mp*1] ˮ���������� 
    design.CL_c = CL_c;                     % [Mp*1] ����ϵ�� 
    design.CDp_c = CDp_c;                   % [Mp*1] ����ϵ�� 
    design.CpoDp_c = CpoDp_c;               % [Mp*1] �Ҿ���
    design.T0oCp_c = T0oCp_c;               % [Mp*1] ���ұ� 
    design.T0oDp_c = T0oDp_c;               % [Mp*1] �񾶱� 

    design.part2 = '���������';
    design.Converge_flag = Converge_flag;   % �������Ƿ�����
    design.iteration = G_iter;              % �Ż��ܵ�������
    design.RhoRp = RhoRp;                   % ��챰뾶/ҶƬ�뾶
    design.EAR = EAR;                       % ���������
    design.LM = LM;                         % �������ճ���
    design.VMIV = VMIV;                     % ��ƽ���ٶ�
    design.VMWV = VMWV;                     % 
    design.SIGMAs = SIGMAs;                 % 
    design.Q = Q;                           % ������Ť��
    design.P = P;                           % ����������
    
    if Duct_flag == 1  
        design.part3 = '����������';
        design.RdoRp = RdoRp;               % �����뾶/ҶƬ�뾶
        design.CdoRp = CdoRp;               % �����ҳ�/ҶƬ�뾶
        design.XdoRp = XdoRp;               % ����λ��/ҶƬ�뾶
        design.Gd = Gd;                     % �����л�����
        design.VARING = VARING;             % �������ߴ����������ٶ�
        
        design.XdRING = XdRING;             % [Nd*1] �л���������
        design.UARING = UARING;             % [Nd*1] �����ܵ�ҶƬ������յ��ٶ�
        design.URRING = URRING;             % [Nd*1] �����ܵ�ҶƬ������յ��ٶ�
        design.GdRING = GdRING;             % [Nd*1] �л������������ռ��

        design.DAHIFtT = DAHIF_times_TANBIC;% [Nd,Mp]
        design.DRHIFtT = DRHIF_times_TANBIC;% [Nd,Mp]
        
        design.UADIF = UADIF;               % [Mp*1] �����Ժ����������������յ��ٶ�
        design.UADUCT = UADUCT;             % [Mp*1] �����Ժ����������������յ��ٶ�
        design.CTp_desired = CTp_desired;   % ҶƬ������ϵ��
        design.CTd_desired = CTd_desired;   % ����������ϵ��

        design.TpoT = TpoT;                 % ҶƬ����ռ��
        design.CDd = CDd;                   % ��������ϵ��
        design.part4 = '���ܲ���';
    else
        design.part3 = '���ܲ���';
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