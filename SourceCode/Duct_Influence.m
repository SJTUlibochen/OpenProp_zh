% This function:
%   1) Forms a discrete vortex ring representation of a duct with 
%           -- chord length CdoRp
%           -- axial offset XdoRp
%           -- radius RdoRp
%           -- a NACA a=0.8 meanline 
%
%   2) Calculates axial velocity duct influence function (UADIF)
%      on propeller lifting line control points
%
% -------------------------------------------------------------------------
% Model parameters:
%
%     RdoRp                  duct radius / propeller radius
%     CdoRp                  duct chord  / propeller radius
%     XdoRp                  location of duct mid-chord downstream of propeller (Xduct = 0 by default)
%
%     XdRING            [1,Nd], x/R location of each vortex ring downstream of propeller
%     VARING            [1,Nd], (axial free-stream velocity at duct)/Vs
%
%     GdRING            [1,Nd],  fraction of total non-dimensional duct circulation, sucht that sum(GdRING) = 1 
%                                (i.e. non-dimensional circulation per unit Gd)
%     Gd                         total non-dimensional circulation about the duct, Gd == Gamma_d/(2*pi*R*Vs),
%                                such that the non-dimensional circulation of ring n is Gd*GdRING(n).
%     Gamma_d [m^2/s]            total     dimensional circulation about the duct [m^2/s]
%
%     CDd                       section drag coefficient for the duct
%     CTDDES                    desired duct thrust coefficient
%     CTD                       duct thrust coeff (viscous drag included) with total duct circulation of Gd 
%
% Influence of propeller on duct:
%
%   (DAHIF ,DRHIF)     [Nd,Mp], (axial/radial horseshoe influence functions of prop on duct)*(2*pi*R)
%   (UARING,URRING)    [1,Nd],  (axial/radial velocity induced at duct (XdRING,RdoRp) by prop) / Vs
%
%
% Influence of duct on propeller:
%
%   UADUCT      [1,Mp]  (axial velocity induced on PROP by duct)/Vs
%   UADIF       [1,Mp]   axial velocity induced on PROP by duct per unit Gd, 
%                        i.e. non-dimensional Duct Influence Function
%                         == 2*pi*R * (duct influence function with unit dimensional duct circulation)
%
% -------------------------------------------------------------------------
%
% Inputs:
%     RdoRp,CdoRp,XdoRp
%     RC                        radius of control points on lifting line / R
%
% Outputs:
%     XdRING,GdRING,UADIF
%
% =========================================================================

function [XdRING,GdRING,UADIF] = Duct_Influence(RdoRp,CdoRp,XdoRp,rc,Nd)
    Mp = length(rc);
    % Setup vortex ring axial spacing
    % Note that the code is numerically unstable for Nd > 20 and generally 
    % slow for Nd > 10. Nd must be even.
    dS  = 1/Nd;          % (even spacing between vortex rings) / CdoRp
    hdS = 0.5*dS;        % half of dS

    % -------------------------------------------------------------------------
    % Compute the circulation on the vortex rings which each represent 
    % a section of length dS located at position XdRING of a NACA a=0.8 
    % mean line (L.E.=0.0, T.E.=1.0) that has unit circulation 1 [m^2/s].
    %
    XdRING  = zeros(Nd,1);
    GdRING  = zeros(Nd,1);

    % Note that since Gamma_d*GdRING(n) is the circulation of ring (n), then
    %    the circulation distribution of the rings, GdRING, is actually unitless!
    %    If Gamma_d is scaled or normalized, then GdRING remains unitless and
    %    always gives the correct circulation distribution.
    for index = 1 : Nd
        % 涡环的坐标
        XdRING(index) = (index-1)*dS+hdS;
        % 涡环的两个端点
        x2 = XdRING(index)+hdS;
        x1 = XdRING(index)-hdS;
        if x2 <= 0.8
            GdRING(index) = dS/0.9;
        elseif x1 >= 0.8
            y1 = 1.0-(x1-0.8)/0.2;
            y2 = 1.0-(x2-0.8)/0.2;
            GdRING(index) = dS*0.5*(y1 + y2)/0.9;
        else
            y2 = 1.0-(x2-0.8)/0.2;
            front = 0.8-x1;
            back = 0.5*(1.0+y2)*(x2-0.8);
            GdRING(index) = (front+back)/0.9;
        end
    end
    
    % 默认叶片位于涵道轴向中心，因此前缘边在原点上游LED处
    LED = -(Nd/2)*dS;
    % 将涵道沿轴向的坐标转化为以叶片为原点的坐标
    XdRING = XdRING+LED;
    % 统一使用Rp进行归一化
    XdRING = XdRING*CdoRp-XdoRp;
    
    % ------------------------------------------------------------------------- 
    % Calculate duct influence function (UADIF) on the UASTAR
    % axial velocity at propeller lifting line control points
        % Note: No tangential influence
        % Note: Radial influence does not create a force on radial lifting line
    
    % 单位环量强度下，涵道中线上轴向坐标XdRING的点对控制点产生的轴向诱导速度
    UAD = zeros(Mp,Nd);
    % 单位环量强度下，涵道对控制点产生的总轴向诱导速度
    UADIF = zeros(Mp,1);
    
    for m=1:Mp                      % for each control point on lifting line    
        for n=1:Nd                  % cycle thru all vortex rings on duct
            UAD(m,n) = vRing(XdRING(n),RdoRp,0,rc(m),GdRING(n));
            UADIF(m) = UADIF(m) + UAD(m,n);
        end
    end

    % The discrete vortex ring formulation breaks down near the duct, so 
    % extrapolate UADIF for RC > 0.9*Rduct/R
    indices = find(rc <= 0.9*RdoRp);
    UADIF = interp1(rc(indices),UADIF(indices),rc,'linear','extrap');


    % Since GdRING is non-dimensionalized with respect to Gd, then
    % UADIF represents the axial flow velocity induced per unit Gd, which
    % is the correct influence function.

end  

% =========================================================================
% ===================================================== vRing Function
% Created: J.M. Stubblefield, 2008
% Edited:  B.P. Epps, 2010
%
% Returns axial and radial velocity at field point (XF,RF) induced by 
% vortex ring of NON-DIMENSIONAL circulation Gring = Gamma/(2*pi*R*Vs) 
% at x-axis location (Xring) and radius (Rring).
%
% NOTE: The inputs MUST be non-dimensionalized as follows.  For dimensional
%       inputs, the output velocity must be multiplied by 2*pi*Vs
%
% Example: for R = 1 m, Vs = 1 m/s, Gamma = 1 m^2/s
%              Xring = XF = 0, Rring = 1, RF = 0 (i.e. center of ring)
%          then Gring = 1/(2*pi)
%          so  [UA,UR] = vRing(0,1,0,0,1/(2*pi))
%          check dimensional output: UA*Vs = Gamma / (2*Rring*R)
%
% Variables:
%     R  [m]     propeller radius
%     Vs [m/s]   ship speed
%     Gring      vortex ring circulation / (2*pi*R*Vs)
%     Xring      x/R, axial location of vortex ring / R
%     Rring      r/R, radius         of vortex ring / R
%     XF         x/R, field point at which velocity is induced
%     RF         r/R, field point at which velocity is induced
%
% Returns: 
%     UA         axial  velocity / Vs     at field point
%     UR         radial velocity / Vs     at field point
%
% Ref: Kuchemann and Weber, Aerodynamics of Propulsion p 305.
% =========================================================================

function [UA,UR] = vRing(Xring,Rring,XF,RF,Gring)
    if Rring == 0               %stops function if Rring = 0
        UA = 0;
        UR = 0;
        return
    end
    if XF == Xring && RF == Rring          % stop if field point on vortex ring
        UA = 0;
        UR = 0;
        return
    end
    % --------------------------------------- Non-dimensional coordinates (x,r)
    x = (XF-Xring)/Rring;                       %x/r' from Kuchemann
    r =  RF       /Rring;                       %r/r' from Kuchemann

    %Elliptic integral method (Kuchemann p. 305)
    %uses parameter k where k^2 = m for elliptic integrals

    k     = sqrt(4*r/(x^2+(r+1)^2));
    % [K,E] = ellipke(k^2);
    [K,E] = elliptic12(pi/2,k^2);

    UA = Gring/(Rring)/sqrt(x^2+(r+1)^2)*(K-(1+2*(r-1)/(x^2+(r-1)^2))*E);

    if r==0
        UR = 0;
    else
        UR = Gring/(Rring)*(-x)/r/sqrt(x^2+(r+1)^2)*(K-(1+2*r/(x^2+(r-1)^2))*E);
    end

end