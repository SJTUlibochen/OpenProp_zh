%This function determines the geometry of the propeller.
%It output the geometry as a 2D image, 3D image, and Rhino CAD file.
%
%Reference: 
%[1]J.S. Carlton, "Marine Propellers & Propulsion", ch. 3, 1994.
%[2]Abbott, I. H., and Von Doenhoff, A. E.; Theory of Wing Sections. 
%Dover, 1959. 
%
% -------------------------------------------------------------------------
% Input Variables:
%
%   filename            file name prefix for all output files
%   date                time and date to print on reports
%   Geometry_flag       flag for whether to make 2D/3D geometry plot
%   csv_flag            flag for whetehr to make points output file
%   Meanline            flag for choice of meanline  form
%   Thickness           flag for choice of thickness form
%
%   Xr          [ ],    input radii / propeller radius
%   f0oc0       [ ],    input camber    / chord at each radius
%   t0oc0       [ ],    input thickness / chord at each radius
%   skew0       [deg],  input skew              at each radius
%   rake0       [ ],    input rake / diameter   at each radius
%
%   RC          [ ],    control point radii / propeller radius
%   CL          [ ],    section lift coefficients
%   BetaC       [deg],  Beta  at the control points
%   BetaIC      [deg],  BetaI at the control points
%   alphaI      [deg],  ideal angle of attack
%
%   Dp           [m],    propeller diameter
%   Z           [ ],    number of blades
%   N           [RPM],  propeller speed
%   Dhub        [m],    hub diameter
%   Rhub        [m],    hub radius
%
%   CpoDp         [ ],    chord / diameter at each control point radius
%   Rp           [m],    propeller radius
%   Mp          [ ],    number of radial 2D cross-sections
%   Np          [ ],    number of points in each 2D section
%   Js          [ ],    advance coefficient based on ship speed
%
% Output Variables:
%
% The function has graphical and file outputs, in addition to the geometry 
% data structure.
%
%该函数用于生成叶片模型、几何结果文档以及csv点坐标文档
% -------------------------------------------------------------------------

function [geometry] = Geometry(pt,RG)
    global Fig_Main PlotsValues;
    %nargin表示当前执行函数的输入参数个数
    if nargin == 1
        %nargin = 1时即函数只接受了pt一个参数，没有RG
        %RG是用户自行设定的一套径向位置
        RadiiGiven_flag = 0;
    else 
        RadiiGiven_flag = 1; 
    end
    Date_string = pt.date;
    filename = pt.filename;
    Z = pt.input.Z;
    N = pt.input.N;
    Dp = pt.input.Dp;
    Rp = pt.input.Dp/2;
    Vs = pt.input.Vs;
    Dhub = pt.input.Dhub;
    Rhub = pt.input.Dhub/2;
    Rhub_oR = Rhub/Rp;
    Mp = pt.input.Mp;
    Np = pt.input.Np;
    Meanline = pt.input.Meanline;
    Thickness = pt.input.Thickness;
    Hub_flag = pt.input.Hub_flag;  
    Duct_flag = pt.input.Duct_flag;  
    Chord_flag = pt.input.Chord_flag;  
    Geometry_flag = pt.input.Geometry_flag;
    txt_flag = pt.input.txt_flag;
    csv_flag = pt.input.csv_flag;
    stl_flag = pt.input.stl_flag;
    if isfield(pt.input,'Js')
        Js = pt.input.Js;   
    elseif isfield(pt.input,'L' )
        Js = pi/pt.input.L;
    end
    %【删除请求】QuarterChord_flag不存在
    % 0 == lifting line at mid-chord, 1 == lifting line at quarter-chord (i.e. skew offset == 0.25*C/r)
    if isfield(pt.input,'QuarterChord_flag')
        QuarterChord_flag = pt.input.QuarterChord_flag;
    else
        QuarterChord_flag = 0;
    end
    %【删除请求】LSGeoCorr不存在
    % Lifting Surface Geometry Corrections
    if isfield(pt.input,'LSGeoCorr')
        LSGeoCorr = pt.input.LSGeoCorr;  
    else
        LSGeoCorr = 'none';
    end
    %
    if isfield(pt.input,'Xr')
        Xr = pt.input.Xr; 
        X1 = ones(size(Xr));
        X0 = zeros(size(Xr));
        if isfield(pt.input,'skew0')
            skew0 = pt.input.skew0;
        else
            skew0 = X0;
        end % [deg]
        if isfield(pt.input,'rake0')
            rake0 = pt.input.rake0;
        else
            rake0 = X0;
        end
        if isfield(pt.input,'XCpoDp') 
            XCpoDp = pt.input.XCpoDp;  
            if isfield(pt.input,'Xt0oDp')
                Xt0oDp = pt.input.Xt0oDp;          
            elseif isfield(pt.input,'t0oc0')
                Xt0oDp = pt.input.t0oc0 .* XCpoDp; 
            else
                Xt0oDp = interp1(pt.design.RC,pt.design.t0oDp,Xr,'spline','extrap');
            end 
        else
            XXR = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];
            XXCoD = [0.1600 0.1818 0.2024 0.2196 0.2305 ...
                     0.2311 0.2173 0.1806 0.1387 0.0100];
            Xt0oc0 = [0.2056 0.1551 0.1181 0.0902 0.0694 ...
                      0.0541 0.0419 0.0332 0.0324 0.0100];
            XCpoDp = interp1(XXR,XXCoD ,Xr,'pchip','extrap');
            t0oc0 = interp1(XXR,Xt0oc0,Xr,'pchip','extrap');
            if  isfield(pt.input,'Xt0oDp')
                Xt0oDp = pt.input.Xt0oDp;
            else
                Xt0oDp = t0oc0 .* XCpoDp;   
            end
        end 
    else
        Xr = RC;
        skew0 = 0*RC;
        rake0 = 0*RC;
        Xt0oDp = pt.design.t0oDp;
    end
    %
    if Duct_flag == 1 
        Rduct_oR = pt.design.Rduct_oR;
        Cduct_oR = pt.design.Cduct_oR;
        Xduct_oR = pt.design.Xduct_oR;
        Gd = pt.design.Gd;
        VARING = pt.design.VARING;
        %确认输入中是否设置了涵道的凸缘线和叶型
        %实际上输入中并不存在这一选项
        if isfield(pt.input,'Meanline_d')
            Meanline_d = pt.input.Meanline_d;
        else
            Meanline_d = 'NACA a=0.8 (modified)';
        end
        if isfield(pt.input,'Thickness_d')
            Thickness_d = pt.input.Thickness_d;
        else
            Thickness_d = 'NACA 65A010';
        end
        if isfield(pt.input,'t0oc_duct')
            t0oc_duct = pt.input.t0oc_duct;
        else
            t0oc_duct = 0.12;
        end
    else
        Rduct_oR = 1;
        Cduct_oR = 1;
        Xduct_oR = 0;
        Gd = 0;
        VARING = 0;
    end   
    RC = pt.design.RC;
    RV = pt.design.RV;
    G = pt.design.G';
    CL = pt.design.CL;
    TANBIC = pt.design.TANBIC;
    BetaIC = atand(pt.design.TANBIC);
    BetaC = atand(pt.design.TANBC );
    %将RC对应的βi在RV处进行插值
    TANBIV = pchip(RC,pt.design.TANBIC,RV);
    % Model Diameter
    % This is the diameter of the output SolidWorks/Rhino geometry
    if isfield(pt.input,'Dm')
        Dm = pt.input.Dm;
    else
        Dm = Dp;
    end
    Rm = Dm/2;
    %% 在选定的径向位置对输入的几何体进行插值
    if RadiiGiven_flag == 0
        %没有额外提供用于插值的径向位置
        RG = Rhub_oR + (1-Rhub_oR)*(sin((0:Mp)*pi/(2*Mp)));
    else
        %径向叶片分区数量 = 径向分区线数量-1 = 纬线数量-1
        Mp = length(RG)-1;
        if (RG(1)<Rhub_oR)||(RG(end)>1)
            message = sprintf('错误：径向位置输入值超出范围');
            uialert(Fig_Main,message,'输入错误','icon','error');
            return 
        end
    end
    CL = pchip(RC,CL,RG); 
    BetaIG = pchip(RC,BetaIC,RG); 
    TANBIG = pchip(RC,TANBIC,RG);
    %对弦长和叶片厚度进行插值
    if Chord_flag == 0  
        %没有进行弦长优化，弦长插值方法与EppsOptimizers.m中相同
        if (abs(Xr(end)-1)<1e-4) && (XCpoDp(end)<=0.01)  %Xr(end)=1且XCpoDp(end)=0
            CpoDp = InterpolateChord(Xr,XCpoDp,RG);   % section chord / propeller diameter at ctrl pts
        else
            CpoDp = pchip(Xr, XCpoDp,RG);   % section chord / propeller diameter at ctrl pts
        end
        t0oDp = pchip(Xr,Xt0oDp,RG);   % section thickness / propeller dia. at ctrl pts
    else
        %进行了弦长优化，因此此时弦长的插值方法应注意保护
        %自由螺旋桨叶尖的弦长为0以及涵道螺旋桨叶尖的弦长为有限值
        if (Duct_flag == 0)||(Rduct_oR>1.001)   
            CpoDp = InterpolateChord(RC,pt.design.CpoDp,RG);  % yields zero   chord at the tip       
        else                                                
            CpoDp = pchip(RC,pt.design.CpoDp,RG);  % yields finite chord at the tip   
        end
        t0oDp = pchip(RC,pt.design.t0oDp,RG);
    end
    % Thickness, rake, chord, and radius scaled for the model diameter
    r = RG*Rm;                          % radius of the RG sections [m]
    C = CpoDp*Dm;                          % section chord at the RG sections [m]
    t0 = t0oDp*Dm;
    skew = pchip(Xr,skew0,RG);       % [deg], angular translation along mid-chord helix
    rake = pchip(Xr,rake0,RG)*Dm;     % [m],   translation along propeller axis (3D X-axis)
    %计算伸张面积比 Compute Expanded Area Ratio
    EAR = (2*Z/pi)*trapz(linspace(Rhub_oR,1,100),interp1(RG,CpoDp,linspace(Rhub_oR,1,100),'spline','extrap'));  
    %计算BTF Compute Blade Thickness Fraction (BTF = t0oD at the prop centerline)
    BTF = interp1(RG,t0oDp,0,'linear','extrap');
    % ---------------------------------------- Lay out the 2D coordinate system
    % x0   [ ], x/C distance along mid-chord line to interpolate geometry data.
    % x1   [m], x   distance along mid-chord line to interpolate geometry data.
    %               By definition, x1 == C/2 - C*x0.

    %               At the Leading  Edge: x1 =  C/2, x0 = 0
    %               At the Trailing Edge: x1 = -C/2, x0 = 1
    %
    x0 = zeros(   1,Np);
    x1 = zeros(Mp+1,Np);

    for j = 1:Np                               % for each point along the chord
      % % Even spacing along the chord   
      % x0(j) =                 (j-1)/(Np-1);  % [0 : 1]

        % Cosine spacing along the chord
        x0(j) = 0.5*(1-cos(pi*(j-1)/(Np-1)));  % [0 : 1]
    end
    for i = 1:Mp+1                     % for each radial section along the span
        if QuarterChord_flag == 1
            x1(i,:) = C(i)/4 - C(i)*x0;    % lifting line at quarter-chord
        else        
            x1(i,:) = C(i)/2 - C(i)*x0;    % lifting line at mid-chord
        end
    end
    % ---------------------- Find normalized 2D foil geometry (at x0 positions)
    %   f0octilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
    %    CLItilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
    % alphaItilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
    if iscell(Meanline)  % Assume Meanline is given versus Xr

          f0octilde = zeros(length(Xr),1);  % can either be scalars or given versus RG(1:Mp+1)
           CLItilde = zeros(length(Xr),1);  % can either be scalars or given versus RG(1:Mp+1)
        alphaItilde = zeros(length(Xr),1);  % can either be scalars or given versus RG(1:Mp+1)
        fof0_temp        = zeros(length(Xr),Np); % eventually given versus [RG(1:Mp+1),x/C(1:Np)]
        dfof0dxoc_temp   = zeros(length(Xr),Np); % eventually given versus [RG(1:Mp+1),x/C(1:Np)]
        tot0_temp        = zeros(length(Xr),Np); % eventually given versus [RG(1:Mp+1),x/C(1:Np)]    

        fof0        = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/C(1:Np)]
        dfof0dxoc   = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/C(1:Np)]
        tot0        = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/C(1:Np)]    


        for i = 1:length(Xr)
            [f0octilde(i), CLItilde(i), alphaItilde(i), fof0_temp(i,:), dfof0dxoc_temp(i,:), tot0_temp(i,:)] = GeometryFoil2D(Meanline{i},Thickness{i},x0);
        end

          f0octilde = pchip(Xr,   f0octilde, RG);
           CLItilde = pchip(Xr,    CLItilde, RG);
        alphaItilde = pchip(Xr, alphaItilde, RG);


        for j = 1:Np
                 fof0(:,j) = pchip(Xr,      fof0_temp(:,j), RG');
            dfof0dxoc(:,j) = pchip(Xr, dfof0dxoc_temp(:,j), RG');
                 tot0(:,j) = pchip(Xr,      tot0_temp(:,j), RG');
        end

    else
        [f0octilde, CLItilde, alphaItilde, fof0_temp, dfof0dxoc_temp, tot0_temp] = GeometryFoil2D(Meanline,Thickness,x0);

        fof0        = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/C(1:Np)]
        dfof0dxoc   = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/C(1:Np)]
        tot0        = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/C(1:Np)]

        for i = 1:Mp+1
            %(i,:)矩阵的第i行
                 fof0(i,:) =      fof0_temp; 
            dfof0dxoc(i,:) = dfof0dxoc_temp;
                 tot0(i,:) =      tot0_temp;
        end

    end
    % -------------------------------------------------------------------------


    % ----- Scale camber ratio and ideal angle of attack by 2D lift coefficient
    alphaI = alphaItilde .* CL ./ CLItilde;        % [deg], ideal angle of attack
      f0oc =   f0octilde .* CL ./ CLItilde;        % f0/C, scaled for CL at RG    


    % -------------------------------------------------------------------------
    % Modify camber ratio and ideal angle of attack for 3D lifting surface effects
    if     strcmp(LSGeoCorr,'none') 
        % no geometry corrections

    elseif strcmp(LSGeoCorr,'Morgan1968')
        % Kc = 1; Ka = 1; Kt = 0;
        [Kc,Ka,Kt] = Morgan1968(RG,TANBIG,EAR,Z);

        for m = 1:Mp+1                 % for each radial section along the span
            f0oc(m)   = Kc(m) * f0oc(m);

            alphaI(m) = (Ka(m) * (pi/180)*alphaI(m) + Kt(m) * BTF)*(180/pi); % [deg], ideal angle of attack
        end        

    elseif strcmp(LSGeoCorr,'EckhardtMorgan1955')
        % K1K2 = 1;

        [K1K2] = EckhardtMorgan1955(EAR,RG,TANBIG);

        for m = 1:Mp+1                 % for each radial section along the span
            f0oc(m)   = K1K2(m) * f0oc(m);

            alphaI(m) = K1K2(m) * alphaI(m); % [deg], ideal angle of attack
        end     
    end
    % -------------------------------------------------------------------------


    % ------------------ Find meanline and thickness profiles (at x1 positions)
    % f      = camber               at x1 positions
    % dfdx   = slope of camber line at x1 positions
    % t      = thickness            at x1 positions
    t    = zeros(Mp+1,Np);
    f    = zeros(Mp+1,Np);
    dfdx = zeros(Mp+1,Np);

    for i = 1:Mp+1                 % for each radial section along the span
           f(i,:) =  fof0(i,:)    *f0oc(i)*C(i);
        dfdx(i,:) = dfof0dxoc(i,:)*f0oc(i);

           t(i,:) =  tot0(i,:) * t0(i);
    end      
    % -------------------------------------------------------------------------


    % -------------------------------------------------------------------------
    % ------------------------------------- Find 2D unroatated section profiles
    % x2D  [m], x   position in 2D space on upper (x2D_u) and lower (x2D_l) foil surfaces
    % y2D  [m], y   position in 2D space on upper (x2D_u) and lower (x2D_l) foil surfaces
    x2D_u = zeros(Mp+1,Np);     x2D_l = zeros(Mp+1,Np);
    y2D_u = zeros(Mp+1,Np);     y2D_l = zeros(Mp+1,Np);

    for i = 1:Mp+1                           % for each section along the span
        for j = 1:Np                         % for each point   along the chord
            x2D_u(i,j) = x1(i,j) + (t(i,j)/2)*sin(atan(dfdx(i,j))); % 2D upper surface x
            x2D_l(i,j) = x1(i,j) - (t(i,j)/2)*sin(atan(dfdx(i,j))); % 2D lower surface x
            y2D_u(i,j) =  f(i,j) + (t(i,j)/2)*cos(atan(dfdx(i,j))); % 2D upper surface y
            y2D_l(i,j) =  f(i,j) - (t(i,j)/2)*cos(atan(dfdx(i,j))); % 2D lower surface y
        end
    end



    % % -------------------------------------- Compute leading edge radius points
    % phiLEC = atan(dfdxLE);
    % NLE    = 3; % must be odd to capture leading edge point
    % phiLEs = 3*pi/8;
    % phiLE  = phiLEs:(pi-2*phiLEs)/(NLE-1):pi-phiLEs; 
    % xLEC  = x1(:,1)' - rLE.*cos(phiLEC);
    % yLEC  =            rLE.*sin(phiLEC);
    % 
    % xLE = zeros(Mp+1,NLE);
    % yLE = zeros(Mp+1,NLE);
    % 
    % for i = 1:Mp+1                           % for each section along the span
    %     xLE(i,:) = xLEC(i) + rLE(i)*sin(phiLE+phiLEC(i));
    %     yLE(i,:) = yLEC(i) - rLE(i)*cos(phiLE+phiLEC(i));
    % end


    % ----------------------------------------- Put all the numbers in one list
    % % Nose -> suctioin side -> tail -> pressure side -> nose
    % x2D(:,   1:Np   ) = x2D_u(:,1:Np);     % The first Np values are the upper surface (suction side),
    % x2D(:,1+Np:Np+Np) = x2D_l(:,Np:-1:1);  % and the second Np values are the lower surface (pressure side).
    % y2D(:,   1:Np   ) = y2D_u(:,1:Np);
    % y2D(:,1+Np:Np+Np) = y2D_l(:,Np:-1:1);

    % % j = 1          == tail
    % % j = 1:Np       == suction side
    % % j = Np         == nose
    % % j = Np + 1     == nose
    % % j = Np+ 1:2*Np == pressure side
    % % j = 2*Np       == tail
    % % Tail -> suctioin side -> nose, nose -> pressure side -> tail
    x2D(:,   1:Np   ) = x2D_u(:,Np:-1:1);   % The first Np values are the upper surface (suction side),
    x2D(:,Np+1:Np+Np) = x2D_l(:,1:Np);      % and the second Np values are the lower surface (pressure side).
    y2D(:,   1:Np   ) = y2D_u(:,Np:-1:1);
    y2D(:,Np+1:Np+Np) = y2D_l(:,1:Np);


    % % % Arrange points as follows:
    % % %     Tail -> suctioin side -> leading edge (with radius and nose) -> pressure side -> tail
    % % % j = 1                               == [1         point ] tail
    % % % j = 1              : Np-1           == [Np-1      points] suction side (tail to point aft of leading edge radius)
    % % % j = Np-1+1         : Np-1+(NLE-1)/2 == [(NLE-1)/2 points] suction side along leading edge radius
    % % % j = Np+(NLE-1)/2                    == [1         point ] nose
    % % % j = Np+(NLE-1)/2+1 :    Np-1+NLE    == [(NLE-1)/2 points] pressure side along leading edge radius
    % % % j = Np-1+NLE+1     : 2*(Np-1)+NLE   == [Np-1      points] pressure side (point aft of leading edge radius to tail)
    % % % j = 2*(Np-1)+NLE                    == [1         point ] tail
    % % %
    % % % j =            1   : Np+(NLE-1)/2   == suction  side (tail to nose)
    % % % j = Np+(NLE-1)/2   : 2*(Np-1)+NLE   == pressure side (nose to tail)
    % % 
    % % x2D(:,1              :    Np-1     ) = x2D_u(:,Np:-1:2); % [Np-1 points] suction  side (tail to point aft of leading edge radius)
    % % x2D(:,Np-1+1         :    Np-1 +NLE) =   xLE(:,NLE:-1:1);     % [NLE  points] leading edge radius  
    % % x2D(:,Np-1+NLE+1     : 2*(Np-1)+NLE) = x2D_l(:,2:Np);    % [Np-1 points] pressure side (point aft of leading edge radius to tail)
    % % 
    % % y2D(:,1              :    Np-1     ) = y2D_u(:,Np:-1:2); % [Np-1 points] suction  side (tail to point aft of leading edge radius)
    % % y2D(:,Np-1+1         :    Np-1 +NLE) =   yLE(:,NLE:-1:1);     % [NLE  points] leading edge radius  
    % % y2D(:,Np-1+NLE+1     : 2*(Np-1)+NLE) = y2D_l(:,2:Np);    % [Np-1 points] pressure side (point aft of leading edge radius to tail)



    % %--------------------------------------- plot unrotated blade
    %  Fig2_S = figure('units','normalized','position',[0.31 .06 .4 .3],'name',...
    %         'Blade Image','numbertitle','off');
    %     style=['r' 'g' 'b' 'm' 'k'];
    %     str_prefix = {'r/Rp = '};
    %     flag=1;
    %     for i = 1:ceil(Mp/5):Mp     % for five radial sections from root to tip
    %         plot(x2D(i,:),y2D(i,:),style(flag));
    %         
    %         for j = 1:Np
    %             plot([x2D_l(i,j),x2D_u(i,j)],[y2D_l(i,j),y2D_u(i,j)],style(flag));
    %         end
    %         
    %         str_legend(flag)=strcat(str_prefix,num2str(RC(i)));
    %         hold on;
    %         flag = flag+1;
    %     end
    %     legend(str_legend,'location','northwest');
    %     axis equal;     grid on;
    %     title('2D Blade Image');  xlabel('X (2D) [m]');  ylabel('Y (2D) [m]');
    % %---------------------------------------



    % ---------------------------------------------- Find pitch angle and pitch
    theta    = BetaIG + alphaI;               % Nose-tail pitch angle, [deg]
    PoD      = tand(theta).*pi.*RG;           % Pitch / propeller diameter, [ ]
    theta_Z  = 0:360/Z:360;                   % angle between blades [deg]



    % --------------------------------------- Find 2D roatated section profiles
    % x2Dr [m], x position in 2D space after rotation for pitch angle
    % y2Dr [m], y position in 2D space after rotation for pitch angle
    % x2Dr = zeros(Mp+1,2*(Np-1)+NLE);
    % y2Dr = zeros(Mp+1,2*(Np-1)+NLE);
    x2Dr = zeros(Mp+1,2*Np);
    y2Dr = zeros(Mp+1,2*Np);
    % for i = 1:Mp        % for each section along the span
    for i = 1:Mp+1        % for each section along the span
        x2Dr(i,:) = x2D(i,:)*cosd(theta(i)) - y2D(i,:)*sind(theta(i)); % rotated 2D upper and lower surface x
        y2Dr(i,:) = x2D(i,:)*sind(theta(i)) + y2D(i,:)*cosd(theta(i)); % rotated 2D upper and lower surface y
    end

    % --------------------------- Invoke skew and rake, and find 3D coordinates
    % X3D [m], X position in 3D space (corresponds to y position in 2D space)
    % Y2D [m], Y position in 3D space
    % Z3D [m], Z position in 3D space
    % X3D = zeros(Mp+1,2*(Np-1)+NLE,Z);
    % Y3D = zeros(Mp+1,2*(Np-1)+NLE,Z);
    % Z3D = zeros(Mp+1,2*(Np-1)+NLE,Z);
    X3D = zeros(Mp+1,2*Np,Z);
    Y3D = zeros(Mp+1,2*Np,Z);
    Z3D = zeros(Mp+1,2*Np,Z);
    % for i = 1:Mp        % for each section along the span
    for i = 1:Mp+1        % for each section along the span
    %     for j = 1:2*(Np-1)+NLE    % for each point   along the upper and lower surfaces
        for j = 1:2*Np    % for each point   along the upper and lower surfaces
            for k = 1:Z   % for each blade
                X3D(i,j,k) = - rake(i) - r(i)*(pi*skew(i)/180)*tand(theta(i)) + y2Dr(i,j);
                Y3D(i,j,k) = r(i)*sind(skew(i) - (180/pi)*x2Dr(i,j)/r(i) - theta_Z(k));
                Z3D(i,j,k) = r(i)*cosd(skew(i) - (180/pi)*x2Dr(i,j)/r(i) - theta_Z(k));
            end
        end
    end
    %%
    % % ---- Find axial length of blade
    % max(max(max(X3D))) - min(min(min(X3D)))



    % --------------------- If left-hand screw, then mirror the Y3D coordinates
    LeftHand_flag = 0;  % 1 == left-handed propeller, 0 == right-handed propeller

    if LeftHand_flag == 1
        Y3D = -Y3D;
    end
    % -------------------------------------------------------------------------

    % save geometry RC x2Dr y2Dr X3D Y3D Z3D

    %
    % =========================================================================
    % ============================ Pack up geometry data at the geometry points
    t0oc    = t0 ./ C;                   % [ ],   t0/C

    geometry.Meanline  = Meanline;
    geometry.Thickness = Thickness;

    geometry.Z         = Z;
    geometry.Dp         = Dp;                          % [m]
    geometry.Dhub      = Dhub;                       % [m]

    geometry.EAR       = EAR;
    geometry.BTF       = BTF;

    geometry.RG        = RG;                         % r/Rp
    geometry.CpoDp       = CpoDp;  
    geometry.t0oD      = t0oDp; 
    geometry.t0oc      = t0oc;           
    geometry.f0oc      = f0oc;      
    geometry.BetaI     = BetaIG;   
    geometry.alpha     = alphaI; 
    geometry.theta     = theta;
    geometry.PoD       = PoD; 
    geometry.skew      = skew;                       % [deg]
    geometry.rake      = rake/Dp;                     % 
    %% 绘制叶片模型(2D/3D)
    set(PlotsValues(13),'valuechangedfcn',{@Geometry2Dfcn,RG,x2Dr,y2Dr});
    set(PlotsValues(14),'valuechangedfcn',@Geometry3Dfcn);
%     Make_3D_Blade_Image(X3D,Y3D,Z3D,Duct_flag,Rduct_oR,Cduct_oR,Xduct_oR, Hub_flag,Rhub_oR,  Js,BetaIG,theta,TANBIV,RV,Rm, Geometry_flag,Plots,PlotPanels);
%     if Duct_flag == 1
%         [xd,rd,Xd,Yd,Zd] = Duct_Plot_120329(RC,TANBIC,G,VARING,RV,Z,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR,Cduct_oR,Xduct_oR,Gd,Meanline_d,Thickness_d,x0,t0oc_duct,Rp,Mp,Np);
%         surf(Xd,Yd,Zd);
%     end
    %% 导出几何文档
    if txt_flag
        MakeGeometryTxt(Mp,filename,Date_string,Dp,Z,N,Dhub,Meanline,Thickness,LSGeoCorr,...
                        RG,CpoDp,t0oc,f0oc,PoD,theta,skew,rake);
    end
    %% 导出点坐标csv文档
    if csv_flag
        ExportPointCoordinates(filename,Np,Mp,X3D,Y3D,Z3D);
    end
end
%生成几何文档
function MakeGeometryTxt(Mp,filename,Date_string,Dp,Z,N,Dhub,Meanline,Thickness,LSGeoCorr,...
                         RG,CpoDp,t0oc,f0oc,PoD,theta,skew,rake)
    filename_geometry = strcat(filename,'_Geometry.txt');
    fid = fopen(filename_geometry,'wt');
    fprintf(fid,'%s \n',filename_geometry);
    fprintf(fid,'螺旋桨几何文档\n');
    fprintf(fid,'日期与时间：%s \n',Date_string);
    fprintf(fid,'螺旋桨直径 = %.4f m\n',Dp);
    fprintf(fid,'叶片数量 = %.0f\n',Z);
    fprintf(fid,'装置行进速度 = %.0f RPM\n',N);
    fprintf(fid,'桨毂直径 = %.4f m\n',Dhub);
    if iscell(Meanline)
        fprintf(fid,'Meanline Type: ');
        for j = 1:length(Meanline)
            fprintf(fid,[Meanline{j},', ']);
        end
        fprintf(fid,'\n');
        fprintf(fid,'Thickness Type: ');
        for j = 1:length(Thickness)
            fprintf(fid,[Thickness{j},', ']);
        end
        fprintf(fid,'\n');   
    else
        fprintf(fid,['凸缘线类型：',Meanline,'\n']);
        fprintf(fid,['叶型类型：',Thickness,'\n']);
    end
    fprintf(fid,['Lifting Surface Geometry Corrections: ',LSGeoCorr,'\n']);
    fprintf(fid,' \n');
    fprintf(fid,' \n');
    fprintf(fid,' r/Rp\t  C/Dp\t  t0/C\t   f0/C\t P/Dp\t pitch\t skew \t rake/Dp\n');
    fprintf(fid,'    \t     \t      \t       \t    \t (deg)\t (deg)\t       \n');
    for i = 1:Mp
        fprintf(fid, '%5.4f  %5.4f  %5.4f  %5.4f  %5.4f  %7.4f  %5.4f  %5.4f\n',...
                RG(i),CpoDp(i),t0oc(i),f0oc(i),PoD(i),theta(i),skew(i),rake(i)/Dp);
    end
    fprintf(fid,'\n\n');
    fprintf(fid,'r/Rp  \t [ ], radial position / propeller radius \n');
    fprintf(fid,'Cp/Dp  \t [ ], chord length    / propeller diameter \n');
    fprintf(fid,'P/Dp  \t [ ], pitch           / propeller diameter \n');
    fprintf(fid,'fo/C \t [ ], max camber      / chord length \n');
    fprintf(fid,'to/C \t [ ], max thickness   / chord length \n');
    fclose(fid);
end
%按下“叶切面二维图像”后的回调函数
function Geometry2Dfcn(hObject,ED,RG,x2Dr,y2Dr)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate;
    %% 设置图像面板标题
    set(PlotsPanels,'title',PlotsStrings{13});
    %% 设置相应按钮的可用状态
    %查找被按下的按钮数量
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %对按钮13进行操作后，若一个按钮也没有被按下
        %则除4、5外的所有按钮均可用
        set(PlotsValues,'enable','on');
        if pt.input.Chord_flag ~= 0
            set(PlotsValues(1),'enable','off');
        end
        if pt.input.Analyze_flag == 0
            set(PlotsValues(12),'enable','off');
        end 
        set(PlotsValues(4),'enable','off');
        set(PlotsValues(5),'enable','off');
    else
        %对按钮13进行操作后，若存在被按下的按钮
        %则除13外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(13),'enable','on');
    end
    %% 绘制曲线
    %按钮按下后才绘制图像
    if get(PlotsValues(13),'value')
        %避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        Mp = size(x2Dr,1)-1;
        Np = size(x2Dr,2)/2;
        color = {'r','g','b','m','k'};
        str_prefix = {'r/R = '};  
        cla(coordinate);
        j = 1;
        for index = 1:ceil(Mp/5):Mp
            Contour_line(j) = plot(coordinate,x2Dr(index,:),y2Dr(index,:),...
                                  'color',color{j},'linewidth',2);
            hold(coordinate,'on');
            plot(coordinate,...
                 [0.5*(x2Dr(index,1)+x2Dr(index,2*Np)),0.5*(x2Dr(index,Np)+x2Dr(index,Np+1))],...
                 [0.5*(y2Dr(index,1)+y2Dr(index,2*Np)),0.5*(y2Dr(index,Np)+y2Dr(index,Np+1))],...
                 'color',color{j},'linewidth',1);
            hold(coordinate,'on');
            LegendStrings(j) = strcat(str_prefix,num2str(RG(index)));
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
        set(PlotsPanels,'title','点击左侧按钮显示对应图像');
    end
end
%按下“叶片三维图形”后的回调函数
function Geometry3Dfcn(hObject,ED)
    
end