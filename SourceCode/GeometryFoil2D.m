% Find meanline and thickness profiles at given x0 == x/C positions
% 沿弦长方向的厚度、拱度等几何数据的分布
% 
% Inputs:
%   Meanline              meanline  type (string)
%   Thickness             thickness type (string)
%   x                     x/C distance along chord line to interpolate NACA foil table data.
%
% Variables: 
%   F                     camber at x/C
%   F0                    max camber
%   T                     thickness at x/C
%   T0                    max thickness
%   C                     chord
%   x                     x/C coordinate in 2D NACA foil tables
%                               At the Leading Edge: x = 0
%                               At the Trailing Edge: x = 1
%
%   FoC                   camber / chord ratio (NACA data at x positions)
%   dfdxoc                slope of camber line (NACA data at x positions)
% 
% Outputs:
%   fof0                  F/F0 at x0 == x/C positions 
%   dfof0dxoc             d(F/F0)/d(x/C) at x0 == x/C positions 
%
%   tot0                  T/T0 at x0 == x/C positions 
%          
%   alphaItilde           [deg] NACA data ideal angle of attack
%   CLItilde              NACA data ideal lift coefficient
%   f0octilde             F0/C  NACA data for CLI == CLItilde   
%
% 该函数目前还不支持Meanline_x和Thickness_x为元胞数组的情况，需要修改
% -------------------------------------------------------------------------

function [F0oCtilde,CLItilde,AlphaItilde,FoF0,dFoF0dx,ToT0,As0] = GeometryFoil2D(Meanline_x,Thickness_x,x0)
    %% 凸缘线类型相关参数
    if strcmp(Meanline_x,'NACA a=0.8 (modified)')  
        % NACA a=0.8 (modified)
        x = [0,0.5,0.75,1.25,2.5,5,7.5,10,15,20,25,30,35,40,45,50,55,60,...
             65,70,75,80,85,90,95,100]./100;
        
        FoC = [0,0.281,0.396,0.603,1.055,1.803,2.432,2.981,3.903,4.651,...
               5.257,5.742,6.120,6.394,6.571,6.651,6.631,6.508,6.274,5.913,...
               5.401,4.673,3.607,2.452,1.226,0]./100;
        
        dFdx = [0.53699,0.47539,0.44004,0.39531,0.33404,0.27149,0.23378,...
                0.20618,0.16546,0.13452,0.10873,0.08595,0.06498,0.04507,...
                0.02559,0.00607,-0.01404,-0.03537,-0.05887,-0.08610,...
                -0.12058,-0.18034,-0.23430,-0.24521,-0.24521,-0.24521];
            
        AlphaItilde = 1.40;             % [deg] NACA data ideal angle of attack
        
        CLItilde = 1.00;                % NACA data ideal lift coefficient
        
        F0oCtilde = max(FoC);           % F0/C NACA data for CLI == CLItilde  
          
        FoF0 = FoC/F0oCtilde;           % F/F0 at x = x/C positions
        
        dFoF0dx = dFdx/F0oCtilde;       % d(F/F0)/d(x/C) at x = x/C positions 

        FoF0 = pchip(x,FoF0,x0);        % F/F0 at x0 = x/C positions
        
        dFoF0dx = pchip(x,dFoF0dx,x0);  % d(F/F0)/d(x/C) at  x0 = x/C positions 
        
    elseif strcmp(Meanline_x,'NACA a=0.8')
        % NACA a=0.8
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
        % == F at x/C = 0.50, which is assumed to be the max F/C in the table above!
        maxF = 0.067943423879;
        % camber of meanline / max camber   (   F/F0         ), indexed leading edge to trailing edge
        FoF0 = F / maxF;
        % slope  of meanline / max camber   ( d(F/F0)/d(x/C) ), indexed leading edge to trailing edge
        dfdx = B / maxF;  
        if x0(1) == 0
            FoF0(1) = 0;
            dfdx(1) = interp1(x0(2:end),dfdx(2:end),x0(1),'spline','extrap');
        end
        AlphaItilde  = 1.54;                     % [deg] NACA data ideal angle of attack        
           CLItilde  = 1.00;                     % [ ]   NACA data ideal lift coefficient
          F0oCtilde  = maxF;                     % F0/C  NACA data for CLI == CLItilde      

          dFoF0dx  = dfdx;
          
    elseif strcmp(Meanline_x,'NACA a=0.8 (v1)')               
        % NACA a=0.8
        x = [0 0.5   0.75  1.25  2.5   5     7.5   10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]./100;        
        FoC = [0 0.287 0.404 0.616 1.077 1.841 2.483 3.043 3.985 4.748 5.367 5.863 6.248 6.528 6.709 6.790 6.770 6.644 6.405 6.037 5.514 4.771 3.683 2.435 1.163 0]./100;
        dFdx = [0.54823 .48535 .44925 .40359 .34104 .27718 .23868 .21050 .16892 .13734 .11101 .08775 .06634 .04601 .02613 .00620 -.01433 -.03611 -.06010 -.08790 -.12311   -.18412 -.23921 -.25583 -.24904 -.20385];
        AlphaItilde  = 1.54;                     % [deg] NACA data ideal angle of attack        
           CLItilde  = 1.00;                     % [ ]   NACA data ideal lift coefficient
          F0oCtilde  = max(FoC);                 % F0/C  NACA data for CLI == CLItilde      
        
        FoF0         =    FoC/F0oCtilde;         %   F/F0         at x = x/C positions 
        dFoF0dx    = dFdx/F0oCtilde;         % d(F/F0)/d(x/C) at x = x/C positions 

        FoF0         = pchip(x, FoF0    ,x0);  %   F/F0         at  x0 = x/C positions 
        dFoF0dx    = pchip(x,dFoF0dx,x0);  % d(F/F0)/d(x/C) at  x0 = x/C positions 
    
    elseif strcmp(Meanline_x,'NACA a=0.8 (v2)')
        % NACA a=0.8
        AlphaItilde  = 1.54;                     % [deg] NACA data ideal angle of attack        
           CLItilde  = 1.00;                     % [ ]   NACA data ideal lift coefficient
          F0oCtilde  = 0.0679;                 % F0/C  NACA data for CLI == CLItilde  

        x       = [0.0000 0.0125 0.0250 0.0500 0.0750 0.1000 0.1500 0.2000 0.3000 0.4000 0.4500 0.5000 0.6000 0.7000 0.8000 0.9000 0.9500 1.0000];
        FoF0      = [0.0000 0.0907 0.1586 0.2712 0.3657 0.4482 0.5869 0.6993 0.8635 0.9615 0.9881 1.0000 0.9786 0.8892 0.7027 0.3586 0.1713 0.0000];
        dFoF0dx = pchipslopes(x,FoF0); % == d(F/F0)/d(x/C)  

        FoF0         = pchip(x, FoF0    ,x0);  %   F/F0         at  x0 = x/C positions 
        dFoF0dx    = pchip(x,dFoF0dx,x0);  % d(F/F0)/d(x/C) at  x0 = x/C positions 
    
    elseif strcmp(Meanline_x,'parabolic')
        % parabolic
        % For parabolic meanline: alphaItilde == 0, CLItilde == 4*pi*f0oc 
        % Therefore, set CLItilde and f0octilde such that f0oc == f0octilde * CL / CLItilde
        AlphaItilde = 0;
        CLItilde = 1; 
        F0oCtilde = 1  / (4*pi);
         FoF0       = (1-(2*(x0-0.5)).^2);    
        dFoF0dx   = -8   *(x0-0.5);    
    
    elseif strcmp(Meanline_x,'flat')
        % flat
        % For flat plate meanline: alphaItilde == f0octilde == CLItilde == 0 
        % However, set CLItilde == 1 such that f0oc == f0octilde * CL / CLItilde does not crash
        AlphaItilde = 0;
        CLItilde    = 1; 
        F0oCtilde   = 0;
         FoF0       = 0*x0;  
        dFoF0dx   = 0*x0;
    
    else
        % 其他凸缘线类型
        disp('错误：不支持的凸缘线类型')
        disp('Available meanline profiles:')
        disp('  NACA a=0.8')
        disp('  NACA a=0.8 (modified)')
        disp('  parabolic')
        disp(' ')
    end
    %% 叶型类型
    % NACA 65A010
    if strcmp(Thickness_x,'NACA 65A010')   
        x = [0 .5 .75 1.25 2.5 5 7.5 10 15 20 25 30 35 40 45 50 ...
              55 60 65 70 75 80 85 90 95 100]./100;
        toc = [0 .765 .928 1.183 1.623 2.182 2.65 3.04 3.658 4.127 ...
                  4.483 4.742 4.912 4.995 4.983 4.863 4.632 4.304     ...
                  3.899 3.432 2.912 2.352 1.771 1.188 .604 .021]./100;
        ToT0 = toc / max(toc);          %   T/T0         at x = x/C positions
        ToT0 = pchip(x,ToT0,x0);      %   T/T0         at  x0 = x/C positions 
        As0  = 0.6771;  % == trapz(x0,tot0), with x0 = 0:0.00001:1;  cross-sectional area for C == T0 == 1
    % NACA 65A010 (modified)
    elseif strcmp(Thickness_x,'NACA 65A010 (Epps modified)') || strcmp(Thickness_x,'NACA 65A010 (modified)')   
        x = [0    0.005000000000000   0.007500000000000   0.012500000000000 ...
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
        ToT0 = toc / max(toc);          %   T/T0         at x = x/C positions
        ToT0 = pchip(x,ToT0,x0);      %   T/T0         at  x0 = x/C positions     
        As0  = 0.7107;  % == trapz(x0,tot0), with x0 = 0:0.00001:1;  cross-sectional area for C == T0 == 1
    % NACA 66 (DTRC Modified) thickness form
    elseif strcmp(Thickness_x,'NACA 66 (DTRC modified)') || strcmp(Thickness_x,'NACA66 (DTRC Modified)')
        x  = [0.0000 0.0125 0.0250 0.0500 0.0750 0.1000 0.1500 0.2000 0.3000 0.4000 0.4500 0.5000 0.6000 0.7000 0.8000 0.9000 0.9500 1.0000];
        ToT0 = [0.0000 0.2088 0.2932 0.4132 0.5050 0.5814 0.7042 0.8000 0.9274 0.9904 1.0000 0.9924 0.9306 0.8070 0.6220 0.3754 0.2286 0.0666];
        ToT0 = pchip(x,ToT0,x0);      %   T/T0         at  x0 = x/C positions  
        As0  = 0.7207;  % == trapz(x0,tot0), with x0 = 0:0.00001:1;  cross-sectional area for C == T0 == 1
    % NACA 00xx (four-digit airfoil)
    elseif strcmp(Thickness_x,'NACA 00xx') || strcmp(Thickness_x,'NACA00xx')
        yot0 = (   0.29690*sqrt(x0)  ...    % This NACA 00xx formula gives the y ordinate  
                 - 0.12600*x0        ...    % of the upper and lower surface (i.e. the half-thickness).      
                 - 0.35160*x0.^2     ...    % Here, we assume T0/C = 1 and thus ignore it in the formula, 
                 + 0.28430*x0.^3     ...    % so this is in fact:  y/T0  at the  x0 = x/C positions.
                 - 0.10150*x0.^4) / 0.20;   % The thickess of the section is twice the y values:   (T/T0) = 2*(y/T0)
        ToT0 = 2 * yot0;   % T/T0  at   x0 = x/C positions    
        As0  = 0.3425;  % == trapz(x0,tot0), with x0 = 0:0.00001:1;  cross-sectional area for C == T0 == 1
    % Elliptical thickness form
    elseif strcmp(Thickness_x,'elliptic')  
        ToT0 = sqrt(1-(2*(x0-0.5)).^2);
        As0  = pi/4;
    % Parabolic thickness form
    elseif strcmp(Thickness_x,'parabolic') 
        ToT0 = (1-(2*(x0-0.5)).^2);    
        As0  = 2/3;
    % 其他叶型类型
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