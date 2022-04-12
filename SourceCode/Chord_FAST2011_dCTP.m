% -------------------------------------------------------------------------
% 11/18/2011 Brenden Epps
%
% This function optimizes blade chord and thickness using the method of
% (Epps et al., FAST 2011), given the design point "endurance speed" 
% and a "high speed" state.
%
% This implementation enforces the ABS thickness rule.
%
% -------------------------------------------------------------------------


function [CoD, t0oc, C_res] = Chord_FAST2011_dCTP(SIGMAh,CTDESh,Jh,CoD,t0oc,...      
                                                  Propeller_flag,Viscous_flag,...
                                                  Hub_flag,Duct_flag,Z,Mp,...
                                                  ITER,Rhv,RC,RV,DR,Rhub_oR,...
                                                  VMIV,CTD,L,Gp,VAC,VTC,UASTAR,...
                                                  UTSTAR,VSTAR,TANBIC,CD,R,...
                                                  rho,Vs,N)
                     
    if Duct_flag == 1
        disp('Chord_FAST2011_dCTP.m 不支持考虑涵道的工况');
        return
    else 
        UADUCT = zeros(Mp,1);
    end
    
    % last value of CoD
    CoD_last = CoD;    
    
    disp(' ')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('Performing FAST2011 chord optimization...')
    disp(' ')

    % -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
    % Find off-design state (Gh, VSTARh, Jh) - h == high-speed state
    % -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
    %
    % Compute design state variables
    ALPHA = zeros(Mp,1);    
    CL = 2*pi*Gp./(VSTAR.*CoD);
    Js = pi/L;

    TANBIC0 = TANBIC;
    CL0     = CL;
    CD0     = CD;


    % Off-design analysis parameters
    ALPHAstall  = 8*pi/180 * ones(Mp,1);

    % -------------------------------------------------------------------------
    % Propeller Aspect Ratio (PAR)
    PAR = (1-Rhub_oR)^2 / trapz( linspace(Rhub_oR,1,100) , InterpolateChord(RC,CoD,linspace(Rhub_oR,1,100))  );

    % Lift curve slope:
    dCLdALPHA = (2*pi / (1 + 2/PAR)) * ones(Mp,1);
    % -------------------------------------------------------------------------

    % Compute off-design state                     
    [Gh, VSTARh, Jh, UASTARh,UTSTARh,CDh] = Chord_Find_High_Speed_State(CTDESh, Jh,...
                             ...       
                             Propeller_flag,Viscous_flag,Hub_flag,0,              ...
                             Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV,...
                             TANBIC0,CL0,CD0,                                             ...
                             L,Gp,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL);                     
    
    % If Find_High_Speed_State crashed...     
    if Jh == 0
        CoD   = 0*CoD;
        t0oc  = 0*t0oc;
        C_res = 0;    

        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        disp('Aborting FAST2011 chord optimization... (c/D = 0, t0/c = 0)'),
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        disp(' ')
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        disp('Continuing circulation optimization:')
        disp(' ')

        return
    end
    % -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
    % -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

    % -!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
    % UPDATE CHORD and THICKNESS DISTRIBUTIONS -- CoD and t0oc
    % -!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!- 

    % -------------------------------------------------------------------------
    % Parameters:
    NACAa       = 0.8;               % 'NACA a=0.8' meanline
    RHOLE       = 0.640279919559199; % == 0.00639/(2*0.04995)^2 == leading edge radius / chord, for t0oc==1 (RLE = RHOLE*t0oc^2)
    CAVmargin   = 1.0;               % allowable SIGMAs = CAVmargin*SIGMAs
    CoDmax      = 0.65;              % maximum allowable c/D
    % -------------------------------------------------------------------------

    SIGMA     = SIGMAh./VSTARh.^2;      % local cavitation number

    % -------------------------------------------------------------------------
    % CASE 2: Find optimum using Newton solver
    g1   = 4/pi;
    g2   = pi*Gp./((1+NACAa)*VSTAR);
    g3   =              (Gh./VSTARh - Gp./VSTAR);
    g3sq = (2/RHOLE) .* (Gh./VSTARh - Gp./VSTAR).^2;  % not exactly g3^2

    CoD_2 = zeros(Mp,1);
    for i = 1:Mp
        [CoD_2(i),t0oc(i)] = Chord_ctSolver(CAVmargin*SIGMA(i), g1, g2(i), g3sq(i) , CoD(i), t0oc(i));
    end


    % -----------------------------------------------------------------
    % CASE 1:
    a1   = 4/pi;
    a2   = pi*Gp./((1+NACAa)*VSTAR);
    a3   = (Gh./VSTARh - Gp./VSTAR);

    % Solve for CoD required by CASE 1 for given t0oc calculated using CASE 2 or 3.
    %          (1+CAVmargin*SIGMA-(1+a1*t0oc).^2)
    %         -(2*(1+a1*t0oc).*(a2+a3))
    %         -(a2.^2+2*a2.*a3)
    %
    CoD_1 = ( (2*(1+a1*t0oc).*(a2+a3)) + sqrt( (2*(1+a1*t0oc).*(a2+a3)).^2 + 4*(1+CAVmargin*SIGMA-(1+a1*t0oc).^2).*(a2.^2+2*a2.*a3) ) )./( 2*(1+CAVmargin*SIGMA-(1+a1*t0oc).^2) );


    % -----------------------------------------------------------------

    % Set CoD to maximum of CASE 1 or CASE 2
    CoD = max(CoD_1,CoD_2);

    % Enforce maximum allowable c/D 
    CoD = min(CoD,CoDmax);

    t0oD = t0oc.*CoD;
    % -----------------------------------------------------------------
    % -----------------------------------------------------------------
    % -----------------------------------------------------------------

    % save temp
    % disp('paused...saved temp')
    % pause,

    % -----------------------------------------------------------------
    % Augment CoD and t0oD using ABS_Blade method (where "rated speed" is the endurance speed, Vs)
    % -----------------------------------------------------------------
    % alphaItilde = 1.54; % [deg]
    %    CLItilde = 1;

    CL       = 2*pi*Gp./(VSTAR.*CoD);           % lift coefficient

    alphaI   = 1.54*CL;                         % == alphaItilde.*CL/CLItilde;        

    TANtheta = tand( atand(TANBIC) + alphaI );  % tan(pitch angle)


    % CP == power coefficient
    [junk1,junk2,CP] = Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,Gp,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD);


    % fig(1234); plot(RC,CoD,'k')
    % fig(1235); plot(RC,t0oD,'k')

    [CoD, t0oD] = Chord_ABS_Blade(RC,t0oD,t0oc,TANtheta,  CP,R,rho,Vs,N,Z);  

    disp('t0oD = ');
    disp(t0oD);
    
    CoD_2 = (g1.*t0oD + g2) ./ ( sqrt(1 + CAVmargin*SIGMA - g3sq./t0oD.^2) - 1 );

    CoD_1 = (g1.*t0oD + g2 + g3)./(CAVmargin*SIGMA) + sqrt( ((g1.*t0oD + g2 + g3)./(CAVmargin*SIGMA)).^2 + (2*g3.*(g1.*t0oD + g2) + (g1.*t0oD + g2).^2)./(CAVmargin*SIGMA) );

    t0oc_1 = t0oD./CoD_1;
    t0oc_2 = t0oD./CoD_2;


    % Set CoD to maximum of CASE 1 or CASE 2
    CoD = max(CoD_1,CoD_2);

    % Enforce maximum allowable c/D 
    CoD = min(CoD,CoDmax);

    t0oc = t0oD./CoD;

    C_res = max( abs(CoD-CoD_last) ); % residual CoD
    % -------------------------------------------------------------------------

    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['The max C_res is: ',num2str(C_res)]),
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(' ')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('Continuing circulation optimization:')
    disp(' ')
end


function [c,tau] = Chord_ctSolver(SIGMA, g1, g2, g3sq, c, tau)
    if	c < 0.2
        c = 0.2;
    end
    if  tau < 0.02
        tau = 0.02;
    elseif tau > 0.1
        tau = 0.1;
    end

    relax   = 0.9;

    J = zeros(2,2);
    R =  ones(2,1);

    while abs(R(1)) > 1e-6 || abs(R(2)) > 1e-6
        % given: SIGMA, g1, g2, g3sq, 

        % Residual equations:
        R1      =      (1 + g1*tau + g2/c)^2 +   g3sq/(c^2*tau^2) - 1 - SIGMA; % => 0

        R2      = 2*g1*(1 + g1*tau + g2/c)   - 2*g3sq/(c^2*tau^3);             % => 0

        R       = [R1;R2];

        % J = Jacobian = [dR1dc dR1dt; dR2dc dR2dt]

        dR1dc   = 2 * (1 + g1 * tau + g2 / c) * (-g2 / c^2) - 2 * g3sq / (c^3 * tau^2);

        dR1dt   = 2 * (1 + g1 * tau + g2 / c) * ( g1      ) - 2 * g3sq / (c^2 * tau^3);

        dR2dc   = -2 * g1 * g2 / (c^2) + 4 * g3sq / (c^3 * tau^3);

        dR2dt   =  2 * g1^2            + 6 * g3sq / (c^2 * tau^4);



        J       = [dR1dc dR1dt; ...
                   dR2dc dR2dt];

        % dx = delta[c; tau]

        dx      = linsolve(J, -R);

        % enforce a maximum step size
        dx(1) = sign(dx(1)) * min([abs(dx(1)),0.05]);
        dx(2) = sign(dx(2)) * min([abs(dx(2)),0.02]);

        % update the state variables
        if       (c  + relax * dx(1)) > 0
            c  =  c  + relax * dx(1);
        else
            c  =  c/2;
        end

        if         (tau  + relax * dx(2)) > 0
            tau  =  tau  + relax * dx(2);
        else
            tau  = tau/2;
        end
    end
end
