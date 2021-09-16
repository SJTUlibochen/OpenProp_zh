% Created: 6/28/2011, Brenden Epps
%
% This function computes the vortex Horseshoe Influence Functions given
% in Kerwin, p.179.
%
% UAHIF = 2*pi*R*(HIF in Kerwin)
% UTHIF = 2*pi*R*(HIF in Kerwin)
%
% UAHIF(n,m) = influence of mth horseshoe vortex on nth control point
%
% -------------------------------------------------------------------------
% Locally-constant-pitch wake assumption is made:
%
%             |
%             |            
%             |
%             *---------------       
%             |            
%             |  <--- RC(n)
%             |
%             *---------------
%             |
%             |
%             |        
%             *--------------- r1 = RV(m+1);  TANBIV1 = TANBIC(m)*RC(m)/RV(m+1) 
%             |
%             |  <--- RC(m), TANBIC(m), G(m)
%             |
%             *--------------- r2 = RV(m  );  TANBIV2 = TANBIC(m)*RC(m)/RV(m)
%             |
%             |
%             |
%
% -------------------------------------------------------------------------
% 该函数的最终输出为：
% 单位环量强度下，由每一个涡格对每一个控制点轴向/切向的诱导速度构成的Mp*Mp矩阵
% 其中第n行第m列的元素为单位环量强度下第m个涡格对第n个控制点的轴向/切向诱导速度
% -------------------------------------------------------------------------
% Sign convention:
%
% The circulation is defined positive when it is directed AWAY from 
% the lifting line.
%
% UASTAR is defined positive in the downstream direction
% UTSTAR is defined positive in the direction of the apparent inflow, omega*r
% -------------------------------------------------------------------------

function [UAHIF,UTHIF] = Horseshoe(Mp,Z,TANBIC,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR)
    UAHIF = zeros(Mp,Mp);
    UTHIF = zeros(Mp,Mp);
    for n = 1:Mp                 % for each control point, n     (FOR LOOP MF2)
        for m = 1:Mp             % for each vortex  panel, m     (FOR LOOP MF3)
            % rc = RC(n);
            % r1 = RV(m+1);  TANBIV1 = TANBIC(m)*RC(m)/RV(m+1)
            % r2 = RV(m  );  TANBIV2 = TANBIC(m)*RC(m)/RV(m)
            [UAW1,UTW1] = Wrench(Z,TANBIC(m)*RC(m)/RV(m+1),RC(n),RV(m+1));  % Velocity induced at RC(n) by a unit vortex shed at RV(m+1), (Wrench returns 2*pi*R*u_bar)
            [UAW2,UTW2] = Wrench(Z,TANBIC(m)*RC(m)/RV(m)  ,RC(n),RV(m)  );  % Velocity induced at RC(n) by a unit vortex shed at RV(m)  , (Wrench returns 2*pi*R*u_bar)
            % ---------------------------- Find hub-image effects, Kerwin p.181
            if Hub_flag == 1
                [UAWh1,UTWh1] = Wrench(Z,TANBIC(m)*RC(m)/(Rhub_oR^2/RV(m+1)) ,RC(n), Rhub_oR^2/RV(m+1) ); 
                [UAWh2,UTWh2] = Wrench(Z,TANBIC(m)*RC(m)/(Rhub_oR^2/RV(m)) ,RC(n), Rhub_oR^2/RV(m  ) ); 
                UAW1 = UAW1 - UAWh1; 
                UAW2 = UAW2 - UAWh2;
                UTW1 = UTW1 - UTWh1;
                UTW2 = UTW2 - UTWh2;
            end
            % ----------------------------------------- Find duct-image effects
            if Duct_flag == 1
                [UAWd1,UTWd1] = Wrench(Z, TANBIC(m)*RC(m)/(Rduct_oR^2/RV(m+1)) ,RC(n), Rduct_oR^2/RV(m+1) );
                [UAWd2,UTWd2] = Wrench(Z, TANBIC(m)*RC(m)/(Rduct_oR^2/RV(m  )) ,RC(n), Rduct_oR^2/RV(m  ) );
                UAW1 = UAW1 - UAWd1;
                UAW2 = UAW2 - UAWd2;
                UTW1 = UTW1 - UTWd1;
                UTW2 = UTW2 - UTWd2;
            end        
        % -------------------------- Determine the Horseshoe Influence Function
        % The Horseshoe Influence Function for vortex panel m is the
        % effect of the induction by a helical trailing vortex at:
        % vortex point m   with circulation -Gamma(m) and another at
        % vortex point m+1 with circulation +Gamma(m).
        % UAHIF(n,m) = u_barA horseshoe influence function in eqn 254.
        % UAW(m)     = u_barA Wrench velocity given in eqn 202-203.
        % Note that the Wrench velocity assumes that a positive circulation is
        % directed AWAY from the lifting line.    
            UAHIF(n,m) = UAW1 - UAW2;           % 2*pi*R*(HIF)
            UTHIF(n,m) = UTW1 - UTW2;           % 2*pi*R*(HIF)    
        end                                                % (END FOR LOOP MF3)
    end                                                    % (END FOR LOOP MF2)
end
% Any escape might help to smooth the unattractive truth                  %
% But the suburbs have no charms to soothe the restless dreams of youth   %
% -------------------------------------------------------------------------


