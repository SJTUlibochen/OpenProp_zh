% -------------------------------------------------------------------------
% OpenProp_v3 duct model:
%
% The duct is modeled as a constant-diameter cylinder, composed of Nd 
% vortex rings at axial locations, XdRING(1:Nd).  These vortex rings have 
% non-dimensional strength Gd*GdRING(1:Nd), such that the total 
% non-dimensional circulation of the duct is Gd.
%
% The influence of the propeller on the duct (DAHIF,DRHIF) is computed
% using Horseshoe_intr_110830.m, which implements the circumferential
% average induced velocities of (Hough and Ordway, 1965).
%
% The influence of the duct on the propeller (UADIF) is computed using 
% elliptic integrals in Duct_Influence.m.
%
% Model parameters:
%
%     Rduct_oR                  duct radius / propeller radius
%     Cduct_oR                  duct chord  / propeller radius
%     Xduct_oR                  location of duct mid-chord downstream of propeller (Xduct = 0 by default)
%
%     XdRING            [1,Nd], x/R location of each vortex ring downstream of propeller
%     VARING            [1,1],  (axial free-stream velocity at duct)/Vs
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
% Thrust allocation between propeller and duct:
%
%     CTPdes = CTdes*   TAU;                       % CT desired for the propeller
%     CTDdes = CTdes*(1-TAU);                      % CT desired for the duct
%      CTdes = CTPdes + CTDdes;                    % total thrust coefficient for propeller + duct
%       TAU  = CTPdes / CTdes;                     % ratio of propeller thrust to total thrust
%
% Influence of propeller on duct:
%
%   (DAHIF ,DRHIF)     [Nd,Mp], (axial/radial horseshoe influence functions of prop on duct)*(2*pi*R)
%   (UARING,URRING)    [1,Nd],  (axial/radial velocity induced at duct (XdRING,Rduct_oR) by prop) / Vs
%
%
% Influence of duct on propeller:
%
%   UADUCT      [1,Mp]  (axial velocity induced on PROP by duct)/Vs
%   UADIF       [1,Mp]   axial velocity induced on PROP by duct per unit Gd, 
%                        i.e. non-dimensional Duct Influence Function
%                         == 2*pi*R * (duct influence function with unit dimensional duct circulation)
%
% Relations:
%
%   sum(GdRING) = 1
%
%   UADUCT =  UADIF*Gd;     % influence of duct      on propeller 
%
%   UASTAR = (UAHIF*G)';    % influence of propeller on propeller
%   UTSTAR = (UTHIF*G)';
%
%   UARING = (DAHIF * G)';  % influence of propeller on duct
%   URRING = (DRHIF * G)';
%
% Where:
%   Since, the axial and radial influence functions are proportional to 1/TANBIC, one can call...
%
%       [DAHIF_times_TANBIC, DTHIF, DRHIF_times_TANBIC] = Horseshoe_intr_110830(...,TANBIC == 1,...)
%
% ...and then complete the calculation by:
%
%     for m = 1:Mp;                               
%         DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBIC(m);
%         DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBIC(m);
%     end   
%
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Code implementation changes from OpenProp v2 to v3:
%
% OpenProp v2        >>  OpenProp v3
% -----------            -----------
%   VAD              >>  VARING               (nomenglature change)
% (UADUCT,URDUCT)    >>  (UARING,URRING)      (nomenglature change)
%  UASDI             >>  UADIF                (nomenglature change)
%  Xduct_oR          >>  XdRING               (nomenglature change)
%  UASTARd           >>  UADUCT               (nomenglature change)
%
% ductVelocity.m     >>  Horseshoe_intr_110830.m  + update DAHIF,DRHIF per above
%
% Induced_Velocity.m >>  UASTAR = (UAHIF*G)';  UTSTAR = (UTHIF*G)';  % influence of propeller on propeller  
%                        UARING = (DAHIF*G)';  URRING = (DRHIF*G)';  % influence of propeller on duct
%                        UADUCT =  UADIF*Gd;                         % influence of duct      on propeller
%                        
% 
% ductThrust.m       >>  Duct_Thrust.m    + nomenglature changes above
% ductInfluence.m    >>  Duct_Influence.m + nomenglature changes above
%
% % -------------------------------------------------------------------------
% % DAHIF(n,m) = influence of m-th horseshoe vortex shed from propeller (Mp panels) 
% %                    on the n-th control point of the duct            (Nd rings) 
% %
% [DAHIF_times_TANBIC, DTHIF, DRHIF_times_TANBIC] = Horseshoe_intr_110830(XdRING,Rduct_oR, RC,ones(size(TANBIC)),RV,Z,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR); 
% 
% for m = 1:Mp;                               
%     DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBIC(m);
%     DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBIC(m);
% end  
% % -------------------------------------------------------------------------

% =========================================================================





