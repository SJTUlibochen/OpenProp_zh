% Find meanline and thickness profiles at given x0 == x/c positions
%
% Inputs:
%   Meanline              meanline  type (string)
%   Thickness             thickness type (string)
%   x0                    x/c distance along chord line to interpolate NACA foil table data.
%
% Variables: 
%         f(xoc)              camber    at x/c
%         f0              max camber
%         t(xoc)              thickness at x/c
%         t0              max thickness
%         c               chord
%       xoc               x/c coordinate in 2D NACA foil tables
%                               At the Leading  Edge: xoc = x0 = 0
%                               At the Trailing Edge: xoc = x0 = 1
%
%       foc               camber / chord ratio (NACA data at xoc positions)
%       dfdxoc            slope of camber line (NACA data at xoc positions)
% 
% Outputs:
%     fof0                f/f0         at x0 == x/c positions 
%     dfof0dxoc         d(f/f0)/d(x/c) at x0 == x/c positions 
%
%     tot0                t/t0         at x0 == x/c positions 
%          
%     alphaItilde         [deg] NACA data ideal angle of attack
%        CLItilde         [ ]   NACA data ideal lift coefficient
%       f0octilde         f0/c  NACA data for CLI == CLItilde   
%
% -------------------------------------------------------------------------

function [f0octilde, CLItilde, alphaItilde, fof0, dfof0dxoc, tot0, As0] = GeometryFoil2D(Meanline,Thickness,x0)
    %% 输出OpenProp支持的叶型和凸缘线类型
    if nargin == 0  
        disp(' ')
        disp('Available meanline and thickness profiles:')
        disp(' ')
        disp('Available meanline profiles:')
        disp('  NACA a=0.8')
        disp('  NACA a=0.8 (modified)')
        disp('  parabolic')
        disp('  flat')
        disp(' ')
        disp('Available thickness profiles:')
        disp('  NACA 65A010')
        disp('  NACA 65A010 (Epps modified)') 
        disp('  NACA 66 (DTRC modified)') 
        disp('  NACA 00xx')
        disp('  elliptic') 
        disp('  parabolic') 
        return
    elseif nargin == 2
        x0 = 0:0.1:1;
    end
    %% 检查GUI输入是否在支持范围内
    Print_list = 0;
    if strcmp(Meanline,'NACA a=0.8 (modified)')
        Meanline = 'NACA a=0.8 (modified)';
    elseif strcmp(Meanline,'NACA a=0.8')
        Meanline = 'NACA a=0.8';
    elseif strcmp(Meanline,'parabolic') || strcmp(Meanline,'Parabolic')
        Meanline = 'parabolic';
    elseif strcmp(Meanline,'flat')
        Meanline = 'flat';
    else
        Print_list = 1;
    end
    if strcmp(Thickness,'NACA 65A010')
        Thickness = 'NACA 65A010';
    elseif strcmp(Thickness,'elliptic') || strcmp(Thickness,'elliptical')
        Thickness = 'elliptic';
    elseif strcmp(Thickness,'parabolic')
        Thickness = 'parabolic';
    elseif strcmp(Thickness,'NACA 65A010 (Epps modified)') || strcmp(Thickness,'NACA 65A010 (modified)')
        Thickness = 'NACA 65A010 (Epps modified)';
    elseif strcmp(Thickness,'NACA 66 (DTRC modified)') || strcmp(Thickness,'NACA66 (DTRC modified)') || strcmp(Thickness,'NACA66 (DTRC Modified)')
        Thickness = 'NACA 66 (DTRC modified)';
    elseif strcmp(Thickness,'NACA 00xx')
        Thickness = 'NACA 00xx';
    else
        Print_list = 1;
    end
    if Print_list == 1
        disp(' ')
        disp('Available meanline and thickness profiles:')
        disp(' ')
        disp('Available meanline profiles:')
        disp('  NACA a=0.8')
        disp('  NACA a=0.8 (modified)')
        disp('  parabolic')
        disp('  flat')
        disp(' ')
        disp('Available thickness profiles:')
        disp('  NACA 65A010')
        disp('  NACA 65A010 (Epps modified)') 
        disp('  NACA 66 (DTRC modified)')
        disp('  elliptic') 
        disp('  parabolic') 
        return
    end
    %% 凸缘线类型相关参数
    %NACA a=0.8 (modified)
    if strcmp(Meanline,'NACA a=0.8 (modified)')  
        xoc = [0 .5 .75 1.25 2.5 5 7.5 10 15 20 25 30 35 40 45 50 ...
              55 60 65 70 75 80 85 90 95 100]./100;
        foc = [0 0.281 0.396 0.603 1.055 1.803 2.432 2.981 3.903 4.651 5.257 ...
                 5.742 6.120 6.394 6.571 6.651 6.631 6.508 6.274 5.913 5.401 ...
                 4.673 3.607 2.452 1.226 0  ]./100;
        dfdxoc = [  0.53699  0.47539  0.44004  0.39531  0.33404  0.27149  0.23378  0.20618 0.16546 ...
                    0.13452  0.10873  0.08595  0.06498  0.04507  0.02559  0.00607 ...
                   -0.01404 -0.03537 -0.05887 -0.08610 -0.12058 -0.18034 -0.23430 ...
                   -0.24521 -0.24521 -0.24521];
        alphaItilde  = 1.40;                     % [deg] NACA data ideal angle of attack
           CLItilde  = 1.00;                     % [ ]   NACA data ideal lift coefficient
          f0octilde  = max(foc);                 % f0/c  NACA data for CLI == CLItilde  
          
        fof0         =    foc/f0octilde;         %   f/f0         at xoc = x/c positions 
        dfof0dxoc    = dfdxoc/f0octilde;         % d(f/f0)/d(x/c) at xoc = x/c positions 

        fof0         = pchip(xoc, fof0    ,x0);  %   f/f0         at  x0 = x/c positions 
        dfof0dxoc    = pchip(xoc,dfof0dxoc,x0);  % d(f/f0)/d(x/c) at  x0 = x/c positions 
        
    %NACA a=0.8 meanline
    elseif strcmp(Meanline,'NACA a=0.8')
        xcn = x0;               % NOTE: indexed leading edge to trailing edge
        a = 0.8;    % For NACA a=0.8
        g = -1/(1-a) * (a^2*(log(a)/2-1/4)+1/4);
        h =  1/(1-a) * ((1-a)^2*log(1-a)/2 - (1-a)^2/4) + g;
        C1 = max(1 - xcn,1e-6);
        CA = a - xcn;
        N = length(xcn);
        for i = 1:N
            if abs(CA(i))<1e-6
               CA(i) = CA(i)+1e-5;
            end
        end
        P = 1/2*CA.^2.*log(abs(CA))-1/2*C1.^2.*log(C1)+1/4*(C1.^2-CA.^2);
        % Camber distribution -- excluding ideal angle of attack
        F = (P/(1-a)-xcn.*log(xcn)+g-h*xcn)/(2*pi*(a+1));    
        % Slope of meanline -- excluding ideal angle of attack
        B = ((-CA.*log(abs(CA))+C1.*log(C1))/(1-a) - log(xcn) -1-h)/(2*pi*(a+1)); 
        % == F at x/c = 0.50, which is assumed to be the max f/c in the table above!
        maxF = 0.067943423879;
        % camber of meanline / max camber   (   f/f0         ), indexed leading edge to trailing edge
        fof0 = F / maxF;
        % slope  of meanline / max camber   ( d(f/f0)/d(x/c) ), indexed leading edge to trailing edge
        dfdx = B / maxF;  
        if x0(1) == 0
            fof0(1) = 0;
            dfdx(1) = interp1(x0(2:end),dfdx(2:end),x0(1),'spline','extrap');
        end
        alphaItilde  = 1.54;                     % [deg] NACA data ideal angle of attack        
           CLItilde  = 1.00;                     % [ ]   NACA data ideal lift coefficient
          f0octilde  = maxF;                     % f0/c  NACA data for CLI == CLItilde      

          dfof0dxoc  = dfdx;
    %NACA a=0.8 meanline
    elseif strcmp(Meanline,'NACA a=0.8 (v1)')               
        xoc = [0 0.5   0.75  1.25  2.5   5     7.5   10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]./100;        
        foc = [0 0.287 0.404 0.616 1.077 1.841 2.483 3.043 3.985 4.748 5.367 5.863 6.248 6.528 6.709 6.790 6.770 6.644 6.405 6.037 5.514 4.771 3.683 2.435 1.163 0]./100;
        dfdxoc = [0.54823 .48535 .44925 .40359 .34104 .27718 .23868 .21050 .16892 .13734 .11101 .08775 .06634 .04601 .02613 .00620 -.01433 -.03611 -.06010 -.08790 -.12311   -.18412 -.23921 -.25583 -.24904 -.20385];
        alphaItilde  = 1.54;                     % [deg] NACA data ideal angle of attack        
           CLItilde  = 1.00;                     % [ ]   NACA data ideal lift coefficient
          f0octilde  = max(foc);                 % f0/c  NACA data for CLI == CLItilde      
        
        fof0         =    foc/f0octilde;         %   f/f0         at xoc = x/c positions 
        dfof0dxoc    = dfdxoc/f0octilde;         % d(f/f0)/d(x/c) at xoc = x/c positions 

        fof0         = pchip(xoc, fof0    ,x0);  %   f/f0         at  x0 = x/c positions 
        dfof0dxoc    = pchip(xoc,dfof0dxoc,x0);  % d(f/f0)/d(x/c) at  x0 = x/c positions 
    %NACA a=0.8 meanline
    elseif strcmp(Meanline,'NACA a=0.8 (v2)')
        alphaItilde  = 1.54;                     % [deg] NACA data ideal angle of attack        
           CLItilde  = 1.00;                     % [ ]   NACA data ideal lift coefficient
          f0octilde  = 0.0679;                 % f0/c  NACA data for CLI == CLItilde  

        xoc       = [0.0000 0.0125 0.0250 0.0500 0.0750 0.1000 0.1500 0.2000 0.3000 0.4000 0.4500 0.5000 0.6000 0.7000 0.8000 0.9000 0.9500 1.0000];
        fof0      = [0.0000 0.0907 0.1586 0.2712 0.3657 0.4482 0.5869 0.6993 0.8635 0.9615 0.9881 1.0000 0.9786 0.8892 0.7027 0.3586 0.1713 0.0000];
        dfof0dxoc = pchipslopes(xoc,fof0); % == d(f/f0)/d(x/c)  

        fof0         = pchip(xoc, fof0    ,x0);  %   f/f0         at  x0 = x/c positions 
        dfof0dxoc    = pchip(xoc,dfof0dxoc,x0);  % d(f/f0)/d(x/c) at  x0 = x/c positions 
    %parabolic meanline
    elseif strcmp(Meanline,'parabolic')
        % For parabolic meanline: alphaItilde == 0, CLItilde == 4*pi*f0oc 
        % Therefore, set CLItilde and f0octilde such that f0oc == f0octilde * CL / CLItilde
        alphaItilde = 0;
        CLItilde = 1; 
        f0octilde = 1  / (4*pi);
         fof0       = (1-(2*(x0-0.5)).^2);    
        dfof0dxoc   = -8   *(x0-0.5);    
    %flat meanline
    elseif strcmp(Meanline,'flat')
        % For flat plate meanline: alphaItilde == f0octilde == CLItilde == 0 
        % However, set CLItilde == 1 such that f0oc == f0octilde * CL / CLItilde does not crash
        alphaItilde = 0;
        CLItilde    = 1; 
        f0octilde   = 0;
         fof0       = 0*x0;  
        dfof0dxoc   = 0*x0;
    %其他凸缘线类型
    else
        disp('错误：不支持的凸缘线类型')
        disp('Available meanline profiles:')
        disp('  NACA a=0.8')
        disp('  NACA a=0.8 (modified)')
        disp('  parabolic')
        disp(' ')
    end
    %% 叶型类型
    %NACA 65A010 thickness form
    if strcmp(Thickness,'NACA 65A010')   
        xoc = [0 .5 .75 1.25 2.5 5 7.5 10 15 20 25 30 35 40 45 50 ...
              55 60 65 70 75 80 85 90 95 100]./100;
        toc = [0 .765 .928 1.183 1.623 2.182 2.65 3.04 3.658 4.127 ...
                  4.483 4.742 4.912 4.995 4.983 4.863 4.632 4.304     ...
                  3.899 3.432 2.912 2.352 1.771 1.188 .604 .021]./100;
        tot0 = toc / max(toc);          %   t/t0         at xoc = x/c positions
        tot0 = pchip(xoc,tot0,x0);      %   t/t0         at  x0 = x/c positions 
        As0  = 0.6771;  % == trapz(x0,tot0), with x0 = 0:0.00001:1;  cross-sectional area for c == t0 == 1
    %NACA 65A010 (modified) thickness form
    elseif strcmp(Thickness,'NACA 65A010 (Epps modified)') || strcmp(Thickness,'NACA 65A010 (modified)')   
        xoc = [0    0.005000000000000   0.007500000000000   0.012500000000000 ...
                    0.025000000000000   0.050000000000000   0.075000000000000   0.100000000000000 ...
                    0.150000000000000   0.200000000000000   0.250000000000000   0.300000000000000 ...
                    0.350000000000000   0.400000000000000   0.471204188481675   0.523560209424084 ...
                    0.575916230366492   0.628272251308901   0.680628272251309   0.732984293193717 ...
                    0.785340314136126   0.837696335078534   0.890052356020942   0.942408376963351 ...
                    0.968586387434555   0.981675392670157   0.989528795811518   0.994764397905759 ...
                    0.997382198952880   1.000000000000000];
        toc = [0    0.007650000000000   0.009280000000000   0.011830000000000 ...
                    0.016230000000000   0.021820000000000   0.026500000000000   0.030400000000000 ...
                    0.036580000000000   0.041270000000000   0.044830000000000   0.047420000000000 ...
                    0.049120000000000   0.049950000000000   0.049830000000000   0.048630000000000 ...
                    0.046320000000000   0.043040000000000   0.038990000000000   0.034320000000000 ...
                    0.029120000000000   0.023520000000000   0.017710000000000   0.011880000000000 ...
                    0.008960000000000   0.007499530848329   0.006623639691517   0.006040000000000 ...
                    0.004049015364794   0.000210000000000];        
        tot0 = toc / max(toc);          %   t/t0         at xoc = x/c positions
        tot0 = pchip(xoc,tot0,x0);      %   t/t0         at  x0 = x/c positions     
        As0  = 0.7107;  % == trapz(x0,tot0), with x0 = 0:0.00001:1;  cross-sectional area for c == t0 == 1
    %NACA 66 (DTRC Modified) thickness form
    elseif strcmp(Thickness,'NACA 66 (DTRC modified)') || strcmp(Thickness,'NACA66 (DTRC Modified)')
        xoc  = [0.0000 0.0125 0.0250 0.0500 0.0750 0.1000 0.1500 0.2000 0.3000 0.4000 0.4500 0.5000 0.6000 0.7000 0.8000 0.9000 0.9500 1.0000];
        tot0 = [0.0000 0.2088 0.2932 0.4132 0.5050 0.5814 0.7042 0.8000 0.9274 0.9904 1.0000 0.9924 0.9306 0.8070 0.6220 0.3754 0.2286 0.0666];
        tot0 = pchip(xoc,tot0,x0);      %   t/t0         at  x0 = x/c positions  
        As0  = 0.7207;  % == trapz(x0,tot0), with x0 = 0:0.00001:1;  cross-sectional area for c == t0 == 1
    %NACA 00xx (four-digit airfoil)
    elseif strcmp(Thickness,'NACA 00xx') || strcmp(Thickness,'NACA00xx')
        yot0 = (   0.29690*sqrt(x0)  ...    % This NACA 00xx formula gives the y ordinate  
                 - 0.12600*x0        ...    % of the upper and lower surface (i.e. the half-thickness).      
                 - 0.35160*x0.^2     ...    % Here, we assume t0/c = 1 and thus ignore it in the formula, 
                 + 0.28430*x0.^3     ...    % so this is in fact:  y/t0  at the  x0 = x/c positions.
                 - 0.10150*x0.^4) / 0.20;   % The thickess of the section is twice the y values:   (t/t0) = 2*(y/t0)
        tot0 = 2 * yot0;   % t/t0  at   x0 = x/c positions    
        As0  = 0.3425;  % == trapz(x0,tot0), with x0 = 0:0.00001:1;  cross-sectional area for c == t0 == 1
    %Elliptical thickness form
    elseif strcmp(Thickness,'elliptic')  
        tot0 = sqrt(1-(2*(x0-0.5)).^2);
        As0  = pi/4;
    %Parabolic thickness form
    elseif strcmp(Thickness,'parabolic') 
        tot0 = (1-(2*(x0-0.5)).^2);    
        As0  = 2/3;
    %其他叶型类型
    else
        disp('错误：不支持的叶型类型')
        disp('Available thickness profiles:')
        disp('  NACA 65A010')
        disp('  NACA 65A010 (Epps modified)') 
        disp('  NACA 66 (DTRC modified)')
        disp('  elliptic') 
        disp('  parabolic')   
    end
end