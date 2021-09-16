% This function determines the geometry of the propeller.
% It outputs the geometry as a 2D image, 3D image, and stl file.
%
% �ú��������������ļ��νṹ���������ʽΪ2DͼƬ��3Dʵʱ��Ⱦ���߿��ٳ����ļ�
% 
% Reference: 
%   [1]J.S. Carlton, "Marine Propellers & Propulsion", ch. 3, 1994.
%   [2]Abbott, I. H., and Von Doenhoff, A. E.; Theory of Wing Sections. 
%      Dover, 1959. 
% 
% -------------------------------------------------------------------------
% Input Variables: 
% ���������
%
%   filename            file name prefix for all output files
%                       �ļ��������ļ��������ڶ���������ļ�������
%   date                time and date to print on reports
%                       �������
%   Geometry_flag       flag for whether to make 2D/3D geometry plot
%                       �����Ƿ�����2/3D��ʾģ��
%   Coordinate_flag     flag for whether to make points output file
%                       �����Ƿ�����������ĵ�
%   Meanline_flag       flag for choice of meanline  form
%                       ����͹Ե�������Ƿ��ؾ���ͬһ
%   Thickness_flag      flag for choice of thickness form
%                       ������������Ƿ��ؾ���ͬһ
%   Meanline_x          types of meanline along rx
%                       �ؾ���ֲ���͹Ե������
%   Thickness_x         types of thickness along rx
%                       �ؾ���ֲ��ĺ������
%   rx                  input radii / propeller radius
%                       ��������
%   f0oc0               input camber / chord at each radius
%   t0oc0               input thickness / chord at each radius
%   
%   rc                  control point radii / propeller radius
%                       ���Ƶ�����
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
        % ����ֻ��ptһ���������
        RadiiGiven_flag = 0;
    else 
        RadiiGiven_flag = 1; 
    end
    %% 
    input = pt.input;
    design = pt.design;
    
    filename = pt.filename;
    
    % ��pt.input�е�������
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
    
    % ��pt.design�е�������
    CpoDp_c = design.CpoDp_c;  
    T0oDp_c = design.T0oDp_c;
    if Duct_flag == 1 
        RdoRp = design.RdoRp;
        CdoRp = design.CdoRp;
        XdoRp = design.XdoRp;
        Gd = design.Gd;
        VARING = design.VARING;
        % ��������͹Ե������
        Meanline_d = 'NACA a=0.8 (modified)';
        % ��������������
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
    
    % ��������
    RhoRp = Dh/Dp;
    
    % ��������
    % ������λ��������Ϊ0��λ��1/4�ҳ�����Ϊ1
    QuarterChord_flag = 0;
    % �����漸������
    LSGeoCorr = 'none';
    % ����������ʱΪ1
    LeftHand_flag = 0;
    
    %%
    % ȷ����������
    if RadiiGiven_flag == 0
        % û�ж����ṩ���ڲ�ֵ������
        rg = RhoRp+(1-RhoRp)*(sin((0:Mp)*pi/(2*Mp)));
        % ����ֵ�������Ϊ������
        rg = rg';
    else
        Mp = length(rg)-1;
        if (rg(1)<RhoRp)||(rg(end)>1)
            message = sprintf('���󣺾���λ������ֵ������Χ');
            uialert(Fig_Main,message,'�������','icon','error');
            return 
        end
    end
    
    % �Բ������в�ֵ
    CL_g = pchip(rc,CL_c,rg); 
    BetaI_g = pchip(rc,BetaI_c,rg); 
    tanBetaI_g = pchip(rc,tanBetaI_c,rg);
    if (Duct_flag == 0)||(RdoRp>1.001)
        % Ҷ���ҳ�Ϊ0
        CpoDp_g = InterpolateChord(rc,CpoDp_c,rg);   
    else
        % Ҷ���ҳ�Ϊ����ֵ
        CpoDp_g = pchip(rc,CpoDp_c,rg);
    end
    T0oDp_g = pchip(rc,T0oDp_c,rg);

    % ����һ���ļ��β�����ԭ
    R_g = rg*Dp/2;
    C_g = CpoDp_g*Dp;
    T0_g = T0oDp_g*Dp;
    
    % �������������
    EAR = (2*Z/pi)*trapz(linspace(RhoRp,1,100),...
          interp1(rg,CpoDp_g,linspace(RhoRp,1,100),'spline','extrap'));  
    % ����Ҷо�񾶱�
    BTF = interp1(rg,T0oDp_g,0,'linear','extrap');
    
    % ��������������ϵ
    % x0    x/C distance along mid-chord line to interpolate geometry data.
    % xy    xy  distance along mid-chord line to interpolate geometry data.
    x = zeros(1,Np);
    C_gx = zeros(Mp+1,Np);

    for j = 1 : Np
        % ��������ҳ�����Ϊ���ҷֲ���������ǰԵ�����Ե��
        x(j) = 0.5*(1-cos(pi*(j-1)/(Np-1)));
    end
    for i = 1 : Mp+1
        if QuarterChord_flag == 1
            C_gx(i,:) = C_g(i)/4-C_g(i)*x;
        else
            % ǰԵ�ߴ�ΪC_g(i)/2����Ե�ߴ�Ϊ-C_g(i)/2
            C_gx(i,:) = C_g(i)/2-C_g(i)*x;
        end
    end
    
    % ��ȡ����͹Ե�ߺͺ���������
    F0oCpdtilde_x = zeros(Nx,1);
    CLItilde_x = zeros(Nx,1);
    AlphaItilde_x = zeros(Nx,1);
    FoF0_x = zeros(Nx,Np);
    dFoF0dx_x = zeros(Nx,Np);
    ToT0_x = zeros(Nx,Np);
    FoF0_g = zeros(Mp+1,Np);
    dFoF0dx_g = zeros(Mp+1,Np);
    ToT0_g = zeros(Mp+1,Np);

    % �ؾ�����д�뼸����Ϣ
    for i = 1 : Nx
        [F0oCpdtilde_x(i),CLItilde_x(i),AlphaItilde_x(i),FoF0_x(i,:),...
         dFoF0dx_x(i,:),ToT0_x(i,:)] = GeometryFoil2D(Meanline_x{i},Thickness_x{i},x);
    end

    % �����������޸�Ϊrg
    F0oCptilde_g = pchip(rx,F0oCpdtilde_x,rg);
    CLItilde_g = pchip(rx,CLItilde_x,rg);
    AlphaItilde_g = pchip(rx,AlphaItilde_x,rg);
    for j = 1 : Np
        FoF0_g(:,j) = pchip(rx,FoF0_x(:,j),rg);
        dFoF0dx_g(:,j) = pchip(rx,dFoF0dx_x(:,j),rg);
        ToT0_g(:,j) = pchip(rx,ToT0_x(:,j),rg);
    end
    
    % ����tildeֵȷ��ʵ��ֵ
    AlphaI_g = AlphaItilde_g.*CL_g./CLItilde_g;
    F0oCp_g = F0oCptilde_g.*CL_g./CLItilde_g;
    
    % ������ά��������ع��̾���Թ��Ǻ͹��Ƚ�������
    if strcmp(LSGeoCorr,'none') 
        % û��ѡ����������������������
    elseif strcmp(LSGeoCorr,'Morgan1968')
        % ����Morgan�ķ�����������
        [Kc,Ka,Kt] = Morgan1968(rg,tanBetaI_g,EAR,Z);
        for m = 1 : Mp+1
            F0oCp_g(m) = Kc(m) * F0oCp_g(m);
            AlphaI_g(m) = (Ka(m) * (pi/180)*AlphaI_g(m)+Kt(m)*BTF)*(180/pi);
        end        
        
    elseif strcmp(LSGeoCorr,'EckhardtMorgan1955')
        % ����Eckhard�ķ�����������
        [K1K2] = EckhardtMorgan1955(EAR,rg,tanBetaI_g);
        for m = 1 : Mp+1
            F0oCp_g(m) = K1K2(m) * F0oCp_g(m);
            AlphaI_g(m) = K1K2(m) * AlphaI_g(m);
        end     
    end
    
    % ���ÿ������㴦ʵ�ʵĺ�ȡ����Ⱥ͹��ȵĵ���ֵ
    T_gx = zeros(Mp+1,Np);
    F_gx = zeros(Mp+1,Np);
    dFdx_gx = zeros(Mp+1,Np);
    for i = 1 : Mp+1
        % ע��F0oCd_g(i)��C_g(i)��Ϊ������ֵ�������ھ�����˵�����
        F_gx(i,:) = FoF0_g(i,:)*F0oCp_g(i)*C_g(i);
        dFdx_gx(i,:) = dFoF0dx_g(i,:)*F0oCp_g(i);
        T_gx(i,:) = ToT0_g(i,:)*T0_g(i);
    end      
    
    % �ݾ��Ϊˮ���������Ǻ͹���֮��
    Phi_g = BetaI_g+AlphaI_g;
    % �ݾ�/ҶƬֱ��
    PoDp = tand(Phi_g).*pi.*rg;
    % ҶƬ֮��ĽǶ�
    theta_Z  = 0:360/Z:360;
    T0oCp_g = T0_g./C_g;
    
    % Ҷ����ģ���ϱ���ĺ�����
    x2D_u = zeros(Mp+1,Np);
    % Ҷ����ģ���±���ĺ�����
    x2D_l = zeros(Mp+1,Np);
    % Ҷ����ģ���ϱ����������
    y2D_u = zeros(Mp+1,Np);
    % Ҷ����ģ���±����������
    y2D_l = zeros(Mp+1,Np);
    % ����2DҶ����ģ���ϵĵ�����
    for i = 1 : Mp+1
        % �������귽����Ҷ����Ҷ��
        for j = 1:Np
            % �������귽����ǰԵ�����Ե��
            x2D_u(i,j) = C_gx(i,j)+(T_gx(i,j)/2)*sin(atan(dFdx_gx(i,j)));
            x2D_l(i,j) = C_gx(i,j)-(T_gx(i,j)/2)*sin(atan(dFdx_gx(i,j)));
            y2D_u(i,j) = F_gx(i,j)+(T_gx(i,j)/2)*cos(atan(dFdx_gx(i,j)));
            y2D_l(i,j) = F_gx(i,j)-(T_gx(i,j)/2)*cos(atan(dFdx_gx(i,j)));
        end
    end
    
    % �����±���ĺ������������ֱ�ϳ���ͬһ�����У�˳��
    % ����(i:1 - Mp+1)��Ҷ�� - Ҷ��
    % ����(j:1 - Np)����Ե�� - �ϱ��棨�����棩 - ǰԵ�� - �±��棨ѹ���棩
    x2D(:,1:Np) = x2D_u(:,Np:-1:1);
    x2D(:,Np+1:2*Np) = x2D_l(:,1:Np);
    y2D(:,1:Np) = y2D_u(:,Np:-1:1);
    y2D(:,Np+1:2*Np) = y2D_l(:,1:Np);
    
    % ��������ݾ�Ǻ��Ҷ����2D������
    x2Dr = zeros(Mp+1,2*Np);
    y2Dr = zeros(Mp+1,2*Np);
    for i = 1 : Mp+1
        % �ؾ���λ���ݾ�Ǳ仯
        x2Dr(i,:) = x2D(i,:)*cosd(Phi_g(i))-y2D(i,:)*sind(Phi_g(i));
        y2Dr(i,:) = x2D(i,:)*sind(Phi_g(i))+y2D(i,:)*cosd(Phi_g(i));
    end

    % ������ά������
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
    
    % ����������Ϊ��������������xozƽ��Գ�
    if LeftHand_flag == 1
        y3D = -y3D;
    end
    
    %% ����2/3D��ʾģ��
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
    
    %% ��������
    % �����νṹ���ݱ�����pt.geometry��
    geometry.part1 = '��ƹ��';
    geometry.Z = Z;                             % ҶƬ����
    geometry.Dp = Dp;                           % ҶƬֱ��
    geometry.Dh = Dh;                           % ���ֱ��
    geometry.Duct_flag = Duct_flag;             % ����м��뺭��
    geometry.EAR = EAR;                         % ���������
    geometry.BTF = BTF;                         % Ҷо�񾶱�
    if Duct_flag
        geometry.RdoRp = RdoRp;                 % ����ֱ��
        geometry.CdoRp = CdoRp;                 % �����ҳ�
        geometry.XdoRp = XdoRp;                 % ����λ��
        geometry.Meanline_d = Meanline_d;       % ����͹Ե��
        geometry.Thickness_d = Thickness_d;     % �������
        geometry.Gd = Gd;                       % ��������
        geometry.VARING = VARING;               % �����������������ٶ�
    end
    
    geometry.part2 = 'Ҷ���漸����Ϣ';
    geometry.Nx = Nx;                           % ��������
    geometry.rx = rx;                           % ����λ��
    geometry.Meanline_x = Meanline_x;           % ͹Ե������
    geometry.Thickness_x = Thickness_x;         % �������
    geometry.rc = rc;                           % ����λ��
    geometry.CL_c = CL_c;                       % ����ϵ��
    geometry.BetaI_c = BetaI_c;                 % ˮ����������
    geometry.rg = rg;                           % ����λ��
    geometry.CpoDp_g = CpoDp_g;                 % ҶƬ�Ҿ���
    geometry.T0oDp_g = T0oDp_g;                 % ҶƬ�񾶱�
    geometry.T0oCp_g = T0oCp_g;                 % ҶƬ���ұ�
    geometry.F0oCp_g = F0oCp_g;                 % ҶƬ���ұ�
    geometry.BetaI_g = BetaI_g;                 % ˮ����������
    geometry.AlphaI_g = AlphaI_g;               % ����
    geometry.Phi_g = Phi_g;                     % �����
    geometry.PoDp = PoDp;                       % ҶƬ�ݾ�
    
    geometry.part3 = '��������';
    geometry.Geometry_flag = Geometry_flag;     % ����2/3D��ʾģ��
    geometry.Coordinate_flag = Coordinate_flag; % ����������ĵ�
    geometry.Printing_flag = Printing_flag;     % ���ɿ��ٳ����ļ�
    
    geometry.part4 = '������';
    geometry.x2Dr = x2Dr;                       % 2D��ʾģ�ͺ�����
    geometry.y2Dr = y2Dr;                       % 2D��ʾģ��������
    geometry.x3D = x3D;                         % 3D��ʾģ��x����
    geometry.y3D = y3D;                         % 3D��ʾģ��y����
    geometry.z3D = z3D;                         % 3D��ʾģ��z����
    
    
end

% ��Ҷ�����ά�������Ļص�����
function Geometry2Dfcn(hObject,~,Mp,Np,rg,x2Dr,y2Dr)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    
    % ����ͼ��������
    set(PlotsPanel,'title','Ҷ�����ά����');
    
    % ���ұ����µİ�ť����
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % �԰�ť11���в�������һ����ťҲû�б����£������а�ť������
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
        % �԰�ť11���в����������ڱ����µİ�ť�����11������а�ť��������
        set(PlotsValues,'enable','off');
        set(PlotsValues(11),'enable','on');
    end
    
    % ����ͼ��
    if get(hObject,'value')
        % �����ͼʱ����
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
        set(PlotsPanel,'title','�����ఴť��ʾ��Ӧͼ��');
    end
end

% ����������άģ�͡��Ļص�����
function Geometry3Dfcn(hObject,~,Z,x3D,y3D,z3D)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    
    % ����ͼ��������
    set(PlotsPanel,'title','��������άģ��');
    
    % ���ұ����µİ�ť����
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % �԰�ť12���в�������һ����ťҲû�б����£������а�ť������
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
        % �԰�ť12���в����������ڱ����µİ�ť�����12������а�ť��������
        set(PlotsValues,'enable','off');
        set(PlotsValues(12),'enable','on');
    end
    
    % ����ͼ��
    if get(hObject,'value')
        % �����ͼʱ����
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
        set(PlotsPanel,'title','�����ఴť��ʾ��Ӧͼ��');
    end    

end