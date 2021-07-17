% Make graphical and text reports for a propeller design
%
% -------------------------------------------------------------------------

function [] = Make_Reports(pt)
global Plots PlotPanels;
%���˻�[0.988235 0.631372 0.023529]    ����ʯ��[0.341176 0.584313 0.447058]
%����[0.576470 0.709803 0.811764]     ɽ����[0.380392 0.392156 0.623529]
lc1 = [0.576470 0.709803 0.811764];
lc2 = [0.341176 0.584313 0.447058];
lc3 = [0.380392 0.392156 0.623529];
lc4 = [0.988235 0.631372 0.023529];
input = pt.input;
if isfield(input,'GUI_flag')
    GUI_flag = input.GUI_flag;
else
    GUI_flag = 0;
end
if isfield(pt,'date')
    Date_string = pt.date;
else
    Date_string = date;
end
if isfield(pt,'filename')
    filename = pt.filename;
elseif isfield(pt,'name')
    filename = pt.name;
elseif isfield(pt.input,'filename')
    filename = pt.input.filename; 
else
    filename = 'OpenProp';
end
Z = input.Z;
if isfield(input,'Np')
    Np = input.Np;
else
    Np = 20;
end
if isfield(input,'Vs')
    Vs = input.Vs;
else
    Vs = 1;
end
if isfield(input,'R')
    R = input.R;
else
    R = 1;
end
if isfield(input,'Rhub')
    Rhub = input.Rhub;
else
    Rhub = 0.2;
end
if isfield(input,'rho')
    rho = input.rho;
else
    rho = 1000;
end
%���û�и���ҶƬ�ļ�����ƻ������ٶȣ��ͼٶ�����������Ϊ4118
%4118ֻ�о����������룬�����������룬������ϵ��Ϊ0.008
%ע�⣬4118��Ҷ�Ҵ�����ϵ��Ϊ0.001��Ҷ���Ϊ0
if isfield(input,'XR'), 
    XR  = input.XR;
    X1  =  ones(size(XR));
    X0  = zeros(size(XR));   
    if isfield(input,'XCoD'),  XCoD  = input.XCoD;  else  XCoD = X0; end
    if isfield(input,'t0oc0'), t0oc0 = input.t0oc0; else t0oc0 = X0; end  
    if isfield(input,'XVA'  ),   XVA = input.XVA;   else   XVA = X1; end
    if isfield(input,'XVT'  ),   XVT = input.XVT;   else   XVT = X0; end
else 
    XR    = [];
    X1    = 1;
    X0    = 0;
    XCoD  = [];
    t0oc0 = [];
    XVA   = [];
    XVT   = [];
end
if isfield(input,'XCD')
    XCD = input.XCD; 
    if length(XCD) == 1
        XCD = XCD*X1;
    end
else
    XCD = 0.008*X1;
end
Propeller_flag = input.Propeller_flag;
Viscous_flag   = input.Viscous_flag;
if isfield(input,'Hub_flag')
    Hub_flag = input.Hub_flag;  
else
    Hub_flag = 1;
end
if isfield(input,'Duct_flag')
    Duct_flag = input.Duct_flag;  
else
    Duct_flag = 0;
end
if isfield(input,'Chord_flag')
    Chord_flag = input.Chord_flag;  
else
    Chord_flag = 0;
end
if isfield(input,'Plot_flag')
    Plot_flag = input.Plot_flag;  
else
    Plot_flag = 0;
end
if Propeller_flag == 1
    Js = input.Js;
    L = pi/Js;   
    %CTdes == desired total thrust coefficient == CTPdes + CTDdes
    %������ϵ�� = ����������ϵ�� + ��������ϵ��
    if isfield(input,'THRUST')
        CTdes = input.THRUST / (0.5*rho*Vs^2*pi*R^2);             
    elseif isfield(input,'CTDES')
        CTdes = input.CTDES;             
    elseif isfield(input,'CT')
        CTdes = input.CT;             
    elseif isfield(input,'KTDES')
        CTdes = input.KTDES * (8/pi)/Js^2;           
    elseif isfield(input,'KT')
        CTdes = input.KT * (8/pi)/Js^2;
    else
        CTdes = 0;
    end   
else
    L = input.L;
    Js = pi/L;
    CTdes = 0;     
end

if Viscous_flag == 0
    XCD = X0;
    CDoCL = 0;
end

if Duct_flag == 1
    if isfield(input,'TAU')
        TAU = input.TAU;        %TAU == thrust ratio == propeller thrust / total thrust
    else
        TAU = 1;
    end  
    if isfield(input,'Rduct_oR')
        Rduct_oR = input.Rduct_oR;      % duct radius
    elseif isfield(input,'Rduct')
        Rduct_oR = input.Rduct /R; 
    else
        Rduct_oR = 1;
    end
    if isfield(input,'Cduct_oR')
        Cduct_oR = input.Cduct_oR;      % duct chord length
    elseif isfield(input,'Cduct')
        Cduct_oR = input.Cduct /R; 
    else
        Cduct_oR = 1;
    end
    if isfield(input,'Xduct_oR')
        Xduct_oR = input.Xduct_oR;
    elseif isfield(input,'Xduct')
        Xduct_oR = input.Xduct /R;      % duct axial position downstream
    else
        Xduct_oR = 0;       % i.e. duct mid-chord centered at propeller
    end
    if isfield(input,'CDd')
        CDd = input.CDd;
    else
        CDd = 0.008;
    end    
else
    TAU = 1; 
    Rduct = 0;
    Cduct = 0;
    CDd = 0;
end

if isfield(input,'HUF' ), HUF  = input.HUF;   else  HUF  = 0;    end
if isfield(input,'TUF' ), TUF  = input.TUF;   else  TUF  = 0;    end
if isfield(input,'Rhv' ), Rhv  = input.Rhv;   else  Rhv  = 0.5;  end

if isfield(input,'dVs'),  dVs = input.dVs;  else dVs  = 0.30;  end % m/s
if isfield(input,'H'  ),  H   = input.H;    else H    = 3.048; end % m
if isfield(input,'g'  ),  g   = input.g;    else g    = 9.81;  end % m/s^2
if isfield(input,'Patm'), Patm= input.Patm; else Patm = 101325;end % Pa
if isfield(input,'Pv'),   Pv  = input.Pv;   else Pv   = 2500;  end % Pa

CT = pt.design.CT;
CQ = pt.design.CQ;
CP = pt.design.CP;
VMIV = pt.design.VMIV;
if Propeller_flag == 1
    KT = pt.design.KT;
    KQ = pt.design.KQ;
    EFFY = pt.design.EFFY;
end
RC = pt.design.RC;
Mp = length(RC);
RV = pt.design.RV;
G = pt.design.G';
VAC = pt.design.VAC;
VTC = pt.design.VTC;
UASTAR = pt.design.UASTAR;
UTSTAR = pt.design.UTSTAR;
TANBC = pt.design.TANBC;
TANBIC = pt.design.TANBIC;
CoD = pt.design.CoD;
CL = pt.design.CL;
CD = pt.design.CD;
if Duct_flag == 1    
    TAU = pt.design.TAU;
    XdRING = pt.design.XdRING;
    GdRING = pt.design.GdRING;
    UADIF = pt.design.UADIF;
    Gd = pt.design.Gd;
    UARING = pt.design.UARING;
    URRING = pt.design.URRING;
else
    TAU = 1;
    XdRING = 1;
    GdRING = 0;
    UADIF = 0*RC;
    Gd = 0;
    UARING = 0;
    URRING = 0;
    Rduct_oR = 1;
end
Beta_c = atand(pt.design.TANBC);
BetaI_c = atand(pt.design.TANBIC);

D = 2*R;
Dhub = 2*Rhub;
Rhub_oR = Rhub/R;
N = 60*Vs/(Js*D);
n = N/60;
%% ����'OpenPropSingle.m'��û�л��Ƶ�����
if GUI_flag
    %���ơ������ֲ����ߡ�
    set(0,'CurrentFigure',Plots);
    h = axes('parent',PlotPanels(7));
    axes(h);
    hold on;
    plot(RC,G,'color',lc1,'linestyle','-','LineWidth',2,'marker','.','MarkerSize',16);
    xlabel('r/R','FontSize',16,'FontName','Times');
    ylabel('��/2\piRVs','FontSize',16,'FontName','Times');
    set(gca,'FontSize',14,'FontName','Times');
    grid on; box on,
    ylimits = get(gca,'Ylim');
    if abs(ylimits(2)) > abs(ylimits(1))       
        set(gca,'Ylim',[0 ylimits(2)]);
    else
        set(gca,'Ylim',[ylimits(1) 0]);
    end
    %���ơ��յ��ٶȷֲ����ߡ�   
    set(0,'CurrentFigure',Plots);
    h = axes('parent',PlotPanels(8));
    axes(h);
    hold on;   
    plot(RC,VAC,'color',lc1,'linestyle','-','LineWidth',2,'marker','.','MarkerSize',16)
    plot(RC,VTC,'color',lc2,'linestyle','-','LineWidth',2,'marker','.','MarkerSize',16)
    plot(RC,UASTAR,'color',lc3,'linestyle','-','LineWidth',2,'marker','.','MarkerSize',16)
    plot(RC,UTSTAR,'color',lc4,'linestyle','-','LineWidth',2,'marker','.','MarkerSize',16);   
    legend('Va/Vs','Vt/Vs','Ua*/Vs','Ut*/Vs');   
    xlabel('r/R','FontSize',16,'FontName','Times');
    ylabel('velocity','FontSize',16,'FontName','Times');
    set(gca,'FontSize',14,'FontName','Times');
    grid on; box on,
    %���ơ���ؽǶȵķֲ����ߡ�  
    set(0,'CurrentFigure',Plots);
    h = axes('parent',PlotPanels(9));
    axes(h);
    hold on;
    plot(RC,Beta_c,'color',lc1,'linestyle','-','LineWidth',2,'marker','.','MarkerSize',16)
    plot(RC,BetaI_c,'color',lc2,'linestyle','-','LineWidth',2,'marker','.','MarkerSize',16)
    legend('��','��i'); 
    xlabel('r/R','FontSize',16,'FontName','Times');
    ylabel('Angle','FontSize',16,'FontName','Times');
    set(gca,'FontSize',14,'FontName','Times');
    grid on; box on,    
    ylimits = get(gca,'Ylim');
    set(gca,'Ylim',[0 ylimits(2)]);   
    %���ơ�ҶƬ������������(������ƽ��)��
    set(0,'CurrentFigure',Plots);
    h = axes('parent',PlotPanels(10));
    axes(h);
    hold on;
    XXRC = Rhub_oR + (1-Rhub_oR)*(sin((0:60)*pi/(2*60)));
    XXCoD = InterpolateChord(RC, CoD,XXRC);
    plot(XXRC, XXCoD,'color',lc4,'LineWidth',2);
    plot(XXRC,-XXCoD,'color',lc4,'LineWidth',2);
    plot(RC, CoD,'marker','.','markeredgecolor',lc4,'linestyle','none','MarkerSize',16);
    plot(RC,-CoD,'marker','.','markeredgecolor',lc4,'linestyle','none','MarkerSize',16); 
    xlabel('r/R','FontSize',16,'FontName','Times');
    ylabel('c/R','FontSize',16,'FontName','Times');
    set(gca,'FontSize',14,'FontName','Times');
    grid on; box on,   
    %���ơ�ҶƬ��ȷֲ�����(������ƽ��)��
    set(0,'CurrentFigure',Plots);
    h = axes('parent',PlotPanels(11));
    axes(h);
    hold on;
    %interp1Ϊһά���Բ�ֵ�����������'extrap'��ʾʹ��'spline'�����η��������������
    %�����ָ������������������ֵ�����ڵ�һ������֮��Ĳ���
    %�����޸ġ������RhubӦ��ΪRhub_oR
    plot([Rhub_oR,pt.design.RC,1],...
         interp1(pt.design.RC, pt.design.t0oD,[Rhub_oR,pt.design.RC,1],'spline','extrap'),...
         'color',lc1,'LineWidth',2);
    plot([Rhub_oR,pt.design.RC,1],...
         interp1(pt.design.RC,-pt.design.t0oD,[Rhub_oR,pt.design.RC,1],'spline','extrap'),...
         'color',lc1,'LineWidth',2);
    plot(pt.design.RC, pt.design.t0oD,'marker','.','markeredgecolor',lc1,'linestyle','none','MarkerSize',16);
    plot(pt.design.RC,-pt.design.t0oD,'marker','.','markeredgecolor',lc1,'linestyle','none','MarkerSize',16);
    xlabel('r/R','FontSize',16,'FontName','Times');   
    ylabel('t0/D','FontSize',16,'FontName','Times');            
    set(gca,'FontSize',14,'FontName','Times');
    grid on; box on,
    %���ơ�����ϵ���ֲ�����(������ƽ��)��
    set(0,'CurrentFigure',Plots);
    h = axes('parent',PlotPanels(12));
    axes(h);
    hold on;
    plot([Rhub_oR,pt.design.RC,1],...
         interp1(pt.design.RC, pt.design.CL,[Rhub_oR,pt.design.RC,1],'spline','extrap'),...
         'color',lc3,'LineWidth',2);
    plot(pt.design.RC,pt.design.CL,'marker','.','markeredgecolor',lc3,'linestyle','none','MarkerSize',16);
    xlabel('r/R','FontSize',16,'FontName','Times');   
    ylabel('CL','FontSize',16,'FontName','Times');            
    set(gca,'FontSize',14,'FontName','Times');
    grid on; box on,
    %������ϵ��ͼ������������Ϊ0
    ylimits = get(gca,'Ylim');
    if abs(ylimits(2)) > abs(ylimits(1))
        set(gca,'Ylim',[0 ylimits(2)]);
    else
        set(gca,'Ylim',[ylimits(1) 0]);
    end  
else
    
    Fig1_S = figure('units','normalized','position',[.01 .06 .4 .3],...
        'name','Graphical Report','numbertitle','off');
    
    subplot(2,2,1);
    plot(RC,G);
    xlabel('r/R');     ylabel('Gamma / 2piRVs');      grid on;
    
    if Propeller_flag == 1
        TitleString=strcat('J='  ,num2str(Js,        '%10.3f'),...
                           '; Ct='  ,num2str(CT,  '%10.3f'),...
                           '; Cq='  ,num2str(CQ,  '%10.3f'),...
                           '; Kt='  ,num2str(KT,  '%10.3f'),...
                           '; Kq='  ,num2str(KQ,  '%10.3f'),...
                           '; \eta=',num2str(EFFY,'%10.3f'),...
                           '; \tau=',num2str(TAU, '%10.3f'));
    else
        TitleString=strcat('J='     ,num2str(Js,  '%10.3f'),...
                           '; Ct='  ,num2str(CT,  '%10.3f'),...
                           '; Cq='  ,num2str(CQ,  '%10.3f'),...
                           '; \tau=',num2str(TAU, '%10.3f'));        
    end
    title(TitleString);
    
    subplot(2,2,2);
    plot(RC,VAC,'-b',RC,VTC,'--b',RC,UASTAR,'-.r',RC,UTSTAR,':r');
    xlabel('r/R');   grid on;    legend('Va/Vs','Vt/Vs','Ua*/Vs','Ut*/Vs');
    
    subplot(2,2,3);
    plot(RC,Beta_c,'--b',RC,BetaI_c,'-r');
    xlabel('r/R');   ylabel('Degrees');   grid on;  legend('Beta','BetaI');
    
    subplot(2,2,4);
    plot(RC,CoD);
    xlabel('r/R');   ylabel('c/D');       grid on;    
end
end