%�ó�����ʵ��������͡�������ҶƬ��ơ�����
%�������ݽṹ��ͳһ��ͳһʹ��char����Ҫ����ʱʹ��cell����ʹ��string
%���ں���������ͳһ���ؼ�������ĸ��д�����ֺ�������
%% �����沿��
function OpenPropSingle
    clear variables;
    clear global;
    %% ����ȫ�ֱ���
    global Fig_Main Tabs SingleDesign ParametricStudy TabPageStrings TabPageMenu Menu;
    global OpenPropDirectory Filename SpecificationsValues DuctValues FlagsValues FoilValues...
           N_R0 Xr_in Col_Label XCpoDp_in XCdp_in VAI_in VTI_in ri_in Xt0oDp_in skew0_in rake0_in...
           Meanline_cell Thickness_cell CalculatorsValues;
    global XCpoDp_values XCLmax_values XCdp_values;
    global zhfontsize numfontsize subcolumnfontsize buttonfontsize
    
    Meanline_cell = {'NACA a=0.8' 'NACA a=0.8 (modified)' 'Parabolic'};
    Thickness_cell = {'NACA 65A010' 'NACA 65A010 (modified)' 'Elliptical' 'Parabolic' 'NACA 66 (DTRC modified)'};
    
    %% ��ҶƬ��Ʋ������͡������ٶ�/װ���н��ٶȡ�����ȱʡֵ�趨
    %����Ҷ����뾶����Բ�뾶��ֵ����ľ���λ��
    Xr_def = [0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;0.95;1];
    N_R0 = length(Xr_def);      %����λ�õĸ���
    %���徶��λ�ô����Ҿ��ȣ����߳�������Բֱ���ı�ֵ
    XCpoDp_def = [0.1600;0.1812;0.2024;0.2196;0.2305;0.2311;0.2173;0.1807;0.1388;0.0010];
    %���徶��λ�ô�������ϵ��(Drag Coefficient)��ȱʡֵ��Ϊ������Ϊ0.008
    XCdp_def = ones(N_R0,1).*0.008;
    %���徶��λ�ô���Ҷ��Ȧģ�ҶƬ�ĺ�������߳��ȱ�
    %t0oc0_def = [0.2055;0.1553;0.1180;0.09016;0.06960;0.05418;0.04206;0.03321;0.03228;0.03160];
    %���徶��λ�ô��ĺ񾶱ȣ�ҶƬ�ĺ������԰ֱ���ı�ֵ Xt0oDp_def = t0oc0_def .* XCpoDp_def;
    Xt0oDp_def = [0.0329;0.0281;0.0239;0.0198;0.0160;0.0125;0.0091;0.0060;0.0045;0];
    %���徶��λ�ô�����б��
    skew0_def = zeros(N_R0,1);
    %���徶��λ�ô��Ĳ�б��
    %ʵ���ϸñ�������Ǹ�ֵ��б����
    rake0_def = zeros(N_R0,1);
    %�������ϵ���ֲ�
    %����ļ���ʽ�ǽ�Ҷ�ҵ�Ҷ�����������ϵ��������0.2-0.5֮�䣬�ؾ��������Էֲ�
    XCLmax_def = 0.5 + (0.2-0.5)*(Xr_def-Xr_def(1))/(Xr_def(end)-Xr_def(1));
    %Ϊ�仯�ص�����Ĭ��ֵ
    XCpoDp_values = XCpoDp_def;
    XCLmax_values = XCLmax_def;
    XCdp_values = XCdp_def;
    %���徶��λ�ô��������ֵ
    XVA_def = ones(N_R0,1);
    %���徶��λ�ô��������ֵ
    XVT_def = zeros(N_R0,1);   
    %% ��ҶƬ��񡱲���ȱʡֵ�趨
    Z_def = 2;      %ҶƬ����
    N_def = 3500;        %������ת��(RPM)
    Dp_def = 0.4;      %������ֱ��(m)
    Thrust_def = 1.85;      %װ��������Ҫ��(thrust,N)
    Vs_def = 5;     %װ���н��ٶ�(m/s)
    Dhub_def = 0.02;     %���ֱ��(m)
    rho_def = 1.29;     %�����ܶ�(kg/m^3)
    Mp_def = 20;        %����ҶƬ��������(Number of vortex panels over the radius)
    Np_def = 20;        %����ҶƬ��������(Number of points over the chord)
    %% �������ٲ������������ȱʡֵ�趨
    n_def = N_def/60;       %������ת��(r/s)
    lambda_def = n_def*2*pi*(Dp_def/2)/Vs_def;           %�������ٶ�=Ҷ��Բ���ٶ�/װ���н��ٶ�(Non-dimensional)
    Js_def = Vs_def/(n_def*Dp_def);      %����ϵ��
    KT_def = Thrust_def/(rho_def*n_def^2*Dp_def^4);       %����ϵ��(Thrust Coeffcient)����Ҷ���ٶ�Ϊ�����ٶ�
    CT_def = Thrust_def/(0.5*rho_def*Vs_def^2*pi*(Dp_def/2)^2);       %����ϵ��(Thrust Coefficient),��װ���н��ٶ�Ϊ�����ٶ�    
    %% ������������ز���Ducted Propeller variables
    TpoT_def = 1;        %����������ռ�ȣ����������ṩ������ռ�������ı���(Thrust ratio)
    Cdd_def = 0.008;        %�������������ϵ��(Duct section drag coefficient)
    DdoDp_def = 1;        %����ֱ��/������ֱ��(duct D / prop D)
    filename = 'DefaultPropeller';      %�趨�ļ���ǰ׺��ȱʡֵ
    %% GUI��ʹ�õ�Ԫ�سߴ�ĳ���ֵ
    %���������С
    zhfontsize = 15;
    %�༭�������ֺ͵�λ��С
    numfontsize = 13;
    %ҶƬ��� & ҶƬ��Ʋ��� & �����ٶ�/װ���н��ٶ� & ���ѡ�� & ���������� & �����ٲ��� & ����
    %�������������С
    subcolumnfontsize = 18;
    %��ť�����С  Load & Save & Run OpenProp
    buttonfontsize = 16;
    %���ڵ��ܸ߶Ⱥ��ܿ��
    Windowht = 780;
    Window = 1530;
    %�ļ�Ŀ¼
    OpenPropDirectory = 'OpenProp_zh';
    %% ���������    
    close all;
    Fig_Main = uifigure('position',[5 40 Window Windowht],...
                        'resize','on',...
                        'autoresizechildren','on',...
                        'numbertitle','off',...
                        'name','OpenProp');
    %% ����������
    %���ɵ�һ��������
    MenuStrings = {'�ļ�' '��ǩҳ'};
    MenuTips = {'�����ļ���ز���' '�����Ӧ��ǩҳ���ɹر�'};
    Menu = zeros(1,length(MenuStrings));
    for index = 1 : length(MenuStrings)
        Menu(index) = uimenu('parent',Fig_Main,...
                             'text',MenuStrings{index},...
                             'tooltip',MenuTips{index});
    end    
    %�����ӹ�����
    FileStrings = {'����' '���Ϊ' '����' '����'};
    FileMenu = zeros(1,length(FileStrings));
    for index = 1 : length(FileStrings)
        FileMenu(index) = uimenu('parent',Menu(1),'text',FileStrings{index});
    end    
    set(FileMenu(1),'menuselectedfcn',@SaveData);
    set(FileMenu(2),'menuselectedfcn',@SaveDataCopy);
    set(FileMenu(3),'menuselectedfcn',@LoadData);
    set(FileMenu(4),'menuselectedfcn',@Execute);
    TabPageStrings = {'������ҶƬ���' '�о�������Ӱ��'};
    for index = 1 : length(TabPageStrings)
        TabPageMenu(index) = uimenu('parent',Menu(2),...
                                    'text',TabPageStrings{index},...
                                    'enable','off');
    end    
    %% �����������ϵı�ǩҳ
    %�����������еı�ǩ��uitabgroup
    Tabs = uitabgroup('parent',Fig_Main,...
                      'position',[0 0 Window Windowht],...
                      'selectionchangedfcn',@ToggleMode);
    %����Tabs�ĵ�һ����ǩҳ��������ҶƬ��ơ�
    SingleDesign = uitab('parent',Tabs,...
                         'title','������ҶƬ���');
    %����Tabs�ĵڶ�����ǩҳ���о�������Ӱ�족
    ParametricStudy = uitab('parent',Tabs,...
                            'title','�о�������Ӱ��');
    
    %% ����SingleDesign�е�����ռ估�ռ��еĸ���Ŀ
    %����SingleDesign��2*4������ռ�
    SingleDesignGrid = uigridlayout(SingleDesign,[2 4],...
                                    'columnwidth',{'2x','4x','2x','1.05x'},...
                                    'rowheight',{'3x','1x'});
    
    %���ɡ�ҶƬ�������ԭ��'Specifications'
    Specifications = uipanel('parent',SingleDesignGrid,...
                             'title','ҶƬ���',...
                             'titleposition','centertop',...
                             'fontsize',subcolumnfontsize,...
                             'fontweight','bold',...
                             'scrollable','off');
    
    %���ɡ�ҶƬ��Ʋ���������ԭ��'Blade Design Values'
    BladeDesign = uipanel('parent',SingleDesignGrid,...
                          'title','ҶƬ��Ʋ���',...
                          'titleposition','centertop',...
                          'fontsize',subcolumnfontsize,...
                          'fontweight','bold',...
                          'scrollable','on');
    
    %���ɡ������ٶ�/װ���н��ٶȡ�����ԭ��'Inflow Profile Values'
    Inflow = uipanel('parent',SingleDesignGrid,...
                     'title','�����ٶ�/װ���н��ٶ�',...
                     'titleposition','centertop',...
                     'fontsize',subcolumnfontsize,...
                     'fontweight','bold',...
                     'scrollable','on');
    
    %���ɡ����ѡ�����ԭ��'Options'
    Flags = uipanel('parent',SingleDesignGrid,...
                    'title','���ѡ��',...
                    'titleposition','centertop',...
                    'fontsize',subcolumnfontsize,...
                    'fontweight','bold',...
                    'scrollable','on');
                         
    %���ɡ�����������������ԭ��'Ducted Propeller'
    Duct = uipanel('parent',SingleDesignGrid,...
                   'title','����������',...
                   'titleposition','centertop',...
                   'fontsize',subcolumnfontsize,...
                   'fontweight','bold',...
                   'scrollable','off');

    %���ɡ������ٲ���������ԭ��'Non-dimensional Parameters'
    Calculators = uipanel('parent',SingleDesignGrid,...
                          'title','�����ٲ���',...
                          'titleposition','centertop',...
                          'fontsize',subcolumnfontsize,...
                          'fontweight','bold',...
                          'scrollable','on');

    %���ɡ����ߡ�����ԭ��'Tools'
    Tools = uipanel('parent',SingleDesignGrid,...
                    'title','����',...
                    'titleposition','centertop',...
                    'fontsize',subcolumnfontsize,...
                    'fontweight','bold',...
                    'scrollable','off');
    Tools.Layout.Column = [3 4];        %�����ߡ�һ�����ڶ��е�3��4��
       
    %% ��ҶƬ��������е�����
    %����Specifications��11*2������ռ�
    SpecificationsGrid = uigridlayout(Specifications,[11 2],...
                                      'scrollable','off');    
    %Specifications�е�ǰ9�в���   
    SpecificationsStrings = {'ҶƬ������'...     %ԭ��'Number of blades'
                             '������ת�٣�'...        %ԭ��'Rotation speed'
                             '������ֱ����'...        %ԭ��'Rotor diameter'
                             'װ��������Ҫ��'...      %ԭ��'Required thrust'
                             'װ���н��ٶȣ�'...       %ԭ��'Ship speed'
                             '���ֱ����'...     %ԭ��'Hub diameter'
                             '�����ܶȣ�'...     %ԭ��'Fluid density'
                             '����ҶƬ����������'...     %ԭ��'radial panels'
                             '����ҶƬ����������'};      %ԭ��'chordwise panels'
    SpecificationsTips = {'Z' 'N' 'Dp' 'Thrust' 'Vs' 'Dhub' '��' 'Mp' 'Np'};
    SpecificationsValues_def = {Z_def N_def Dp_def Thrust_def Vs_def Dhub_def rho_def Mp_def Np_def};
    %����ѭ������ǰ9��
    for index = 1 : length(SpecificationsValues_def)
        SpecificationsTexts(index) = uilabel('parent',SpecificationsGrid,...
                                             'text',SpecificationsStrings{index},...
                                             'horizontalalignment','left',...
                                             'verticalalignment','center',...
                                             'tooltip',SpecificationsTips{index},...
                                             'fontsize',zhfontsize);

        SpecificationsValues(index) = uieditfield(SpecificationsGrid,'numeric',...
                                                  'fontsize',numfontsize,...
                                                  'horizontalalignment','center',...
                                                  'value',SpecificationsValues_def{index},...
                                                  'limit',[0 Inf],...
                                                  'lowerlimitinclusive','off',...
                                                  'valuechangedfcn',@UpdateSpecifications);
    end
    %����ǰ9�б༭���е���ʾ��ʽ�͵�λ
    set(SpecificationsValues(1),'valuedisplayformat','%.0f');
    set(SpecificationsValues(2),'valuedisplayformat','%.1f RPM');
    set(SpecificationsValues(3),'valuedisplayformat','%.2f m');
    set(SpecificationsValues(4),'valuedisplayformat','%.3f N');
    set(SpecificationsValues(5),'valuedisplayformat','%.2f m/s');
    set(SpecificationsValues(6),'valuedisplayformat','%.3f m');
    set(SpecificationsValues(7),'valuedisplayformat','%.3f kg/m^3');
    set(SpecificationsValues(8),'valuedisplayformat','%.0f');
    set(SpecificationsValues(9),'valuedisplayformat','%.0f');
    
    %Specifications�еĺ�2�в���
    FoilStrings = {'͹Ե�����ͣ�'...      %ԭ��'Meanline type'
                   'Ҷ�����ͣ�'};        %ԭ��'Thickness type' 
    %����ѭ�����ƺ�2��
    for index = 1 : 2
        FoilTexts(index) = uilabel('parent',SpecificationsGrid,...
                                   'text',FoilStrings{index},...
                                   'horizontalalignment','left',...
                                   'verticalalignment','center',...
                                   'fontsize',zhfontsize);
        FoilValues(index) = uidropdown('parent',SpecificationsGrid,...
                                       'fontsize',numfontsize);
    end
    %���ú�2��ѡ������
    set(FoilValues(1),'items',Meanline_cell);
    set(FoilValues(2),'items',Thickness_cell);
    set(FoilValues(1),'itemsdata',[1 2 3]);
    set(FoilValues(2),'itemsdata',[1 2 3 4 5]);
    %% �������������������е�����
    %����Duct��3*2������ռ�
    DuctGrid = uigridlayout(Duct,[3 2],...
                            'columnwidth',{'1x','1x'},...
                            'rowheight',{'1x','1x','1x'},...
                            'scrollable','off');
    %Duct�еĲ���
    DuctStrings = {'����������ռ�ȣ�'...        %ԭ��'Thrust Ratio'
                   '������������ϵ����'...     %ԭ��'Duct section drag (Cd)'
                   '����ֱ��/������ֱ����'};      %ԭ��'duct D / prop D'
    DuctTips = {'Tp/Thrust' 'Cdd' 'Dd/Dp'};
    DuctValues_def = {TpoT_def Cdd_def DdoDp_def};
    %����ѭ������3��
    for index = 1 : 3
        DuctTexts(index) = uilabel('parent',DuctGrid,...
                                   'text',DuctStrings{index},...
                                   'horizontalalignment','left',...
                                   'verticalalignment','center',...
                                   'tooltip',DuctTips{index},...
                                   'fontsize',zhfontsize);        
        DuctValues(index) = uieditfield(DuctGrid,'numeric',...
                                        'fontsize',numfontsize,...
                                        'horizontalalignment','center',...
                                        'value',DuctValues_def{index});
    end
    %���ò�����Χ
    set(DuctValues(1),'limit',[0 1]);
    set(DuctValues(1),'lowerlimitinclusive','off');
    set(DuctValues(2),'limit',[0 Inf]);
    set(DuctValues(3),'limit',[1 Inf]);
    %% ��ҶƬ��Ʋ����������е�����
    %����BladeDesign��11*6������ռ�
    BladeDesignGrid = uigridlayout(BladeDesign,[N_R0+1 6],...
                                   'columnspacing',0,...
                                   'rowspacing',0,...
                                   'scrollable','on');   
    %BladeDesign�еĲ���
    BladeDesignStrings = {'����λ��'...        %ԭ��'r/R'
                          '�Ҿ���'...     %ԭ��'c/D'
                          '����ϵ��'...        %ԭ��'Cd'
                          '�񾶱�'...     %ԭ��'t0/D'
                          '�����'...     %ԭ��'Skew'
                          'б����'};      %ԭ��'Xs/D'��Ӧ����Ϊ'Zr/D'
    BladeDesignTips = {'r/Rp' 'Cp/Dp' 'Cdp' 't0/Dp' 'skew' 'Zr/Dp'};
    %����ѭ�����Ʊ�����Ŀ
    for index = 1 : length(BladeDesignStrings)
        Col_Label(index) = uilabel('parent',BladeDesignGrid,...
                                   'text',BladeDesignStrings{index},...
                                   'horizontalalignment','center',...
                                   'verticalalignment','center',...
                                   'tooltip',BladeDesignTips{index},...
                                   'fontsize',zhfontsize);
    end
    %����ѭ�����Ʊ�������
    for index = 1 : N_R0
        Xr_in(index) = uieditfield(BladeDesignGrid,'numeric',...
                                   'fontsize',numfontsize,...
                                   'horizontalalignment','center',...
                                   'value',Xr_def(index),...
                                   'limit',[0 1]);

        XCpoDp_in(index) = uieditfield(BladeDesignGrid,'numeric',...
                                     'fontsize',numfontsize,...
                                     'horizontalalignment','center',...
                                     'value',XCpoDp_def(index),...
                                     'limit',[0 Inf]);

        XCdp_in(index) = uieditfield(BladeDesignGrid,'numeric',...
                                    'fontsize',numfontsize,...
                                    'horizontalalignment','center',...
                                    'value',XCdp_def(index),...
                                    'limit',[0 Inf]);

        Xt0oDp_in(index)	= uieditfield(BladeDesignGrid,'numeric',...
                                      'fontsize',numfontsize,...
                                      'horizontalalignment','center',...
                                      'value',Xt0oDp_def(index),...
                                      'limit',[0 Inf]);
        
        skew0_in(index)	= uieditfield(BladeDesignGrid,'numeric',...
                                      'fontsize',numfontsize,...
                                      'horizontalalignment','center',...
                                      'value',skew0_def(index),...
                                      'limit',[0 Inf]);

        rake0_in(index)	= uieditfield(BladeDesignGrid,'numeric',...
                                      'fontsize',numfontsize,...
                                      'horizontalalignment','center',...
                                      'value',rake0_def(index),...
                                      'limit',[0 Inf]);
    end   
    %% �������ٶ�/װ���н��ٶȡ������е�����
    %����Inflow��11*3������ռ�
    InflowGrid = uigridlayout(Inflow,[N_R0+1 3],...
                              'columnspacing',0,...
                              'rowspacing',0,...
                              'scrollable','on');  
    %Inflow�еĲ���
    InflowStrings = {'����λ��'...       %ԭ��'r'���Ʋ�Ϊ'r/R'
                     '�����ֵ'...       %ԭ��'Va/Vs'
                     '�����ֵ'};        %ԭ��'Vt/Vs'
    InflowTips = {'r/Rp' 'Va/Vs' 'Vt/Vs'};  
    %����ѭ�����Ʊ��
    for index = 1 : length(InflowStrings)
        Col_Label2(index) = uilabel('parent',InflowGrid,...
                                    'text',InflowStrings(index),...
                                    'horizontalalignment','center',...
                                    'verticalalignment','center',...
                                    'tooltip',InflowTips(index),...
                                    'fontsize',zhfontsize);
    end
    for index = 1 : N_R0
        ri_in(index) = uieditfield(InflowGrid,'numeric',...
                                   'fontsize',numfontsize,...
                                   'horizontalalignment','center',...
                                   'value',Xr_def(index),...
                                   'limit',[0 1]);

        VAI_in(index) = uieditfield(InflowGrid,'numeric',...
                                    'fontsize',numfontsize,...
                                    'horizontalalignment','center',...
                                    'value',XVA_def(index),...
                                    'limit',[0 Inf]);
                                
        VTI_in(index) = uieditfield(InflowGrid,'numeric',...
                                    'fontsize',numfontsize,...
                                    'horizontalalignment','center',...
                                    'value',XVT_def(index),...
                                    'limit',[0 Inf]);
    end  
    %% �������ٲ����������е�����
    %����Calculators��2*2������ռ�
    CalculatorsGrid = uigridlayout(Calculators,[3 4],...
                                  'scrollable','on');
    %Calculators�еĲ���
    CalculatorsStrings = {'����ϵ����'...     %ԭ��'Js = V/nD ='
                          'Ҷ���ٶ�/�н��ٶȣ�'...     %ԭ��'L = omega*R/V ='
                          '����ϵ��(Ҷ���ٶ�)��'...       %ԭ��'KT = T/(rho*n^2*D^4) ='
                          '����ϵ��(�н��ٶ�)��'};      %ԭ��'CT = T/(1/2*rho*V^2*pi*R^2) ='
    CalculatorsTips = {'Js = V/n*Dp' 'L = ��*Rp/V' 'KT = T/(��*n^2*Dp^4)' 'CT = T/(0.5�Ц�*V^2*Rp^2)'};
    CalculatorsValues_def = {Js_def lambda_def KT_def CT_def};  
    %����ѭ������4�������ٲ���
    for index = 1 : length(CalculatorsStrings)
        CalculatorsTexts(index) = uilabel('parent',CalculatorsGrid,...
                                          'text',CalculatorsStrings{index},...
                                          'horizontalalignment','left',...
                                          'verticalalignment','center',...
                                          'tooltip',CalculatorsTips{index},...
                                          'fontsize',zhfontsize);
        CalculatorsValues(index) = uieditfield(CalculatorsGrid,'numeric',...
                                               'fontsize',numfontsize,...
                                               'horizontalalignment','center',...
                                               'editable','off',...
                                               'enable','on',...
                                               'value',CalculatorsValues_def{index});
    end  
    %% �����ѡ������е�����
    %����Flags��10*1������ռ�
    FlagsGrid = uigridlayout(Flags,[10 1],...
                             'rowspacing',15,...
                             'scrollable','on');
    %Flags�еĲ���
    SelectStrings = {'������'...       %ԭ��'Propeller'
                     '����ҶƬ'};       %ԭ��'Turbine'
    FlagsStrings = {'�Ƿ�������'...     %ԭ��'Hub'���ø�ѡ��ť�����ж����ͼ�����Ƿ�������
                    '�Ƿ�Ϊ����������'...       %ԭ��'Ducted'���ø�ѡ��ť�����ж��Ƿ����ɵ���������
                    '�Ƿ�����ҳ��Ż�'...       %ԭ��'Chord Optimization'���ø�ѡ��ť�����ж��Ƿ�����ҳ��Ż�
                    '�Ƿ���ճ������'...       %ԭ��'Viscous forces'���ø�ѡ��ť�����ж��Ƿ���ճ����
                    '�Ƿ������������'...       %ԭ��'Performance curve'���ø�ѡ��ť�����ж��Ƿ������������
                    '�Ƿ����ҶƬģ��'...       %ԭ��'Geometry plots'���ø�ѡ��ť�����ж��Ƿ����2D��3Dͼ��
                    '�Ƿ񵼳������ĵ�'...       %������ѡ��ť���ø�ѡ��ť�����ж��Ƿ񵼳�����ҶƬ������Ϣ���ĵ�(txt)
                    '�Ƿ񵼳�������'...        %������ѡ��ť���ø�ѡ��ť�����ж��Ƿ񵼳��������ĵ�(csv)
                    '�Ƿ񵼳�stl�ļ�'};        %������ѡ��ť���ø�ѡ��ť�����ж��Ƿ񵼳�stl��ʽ�Ŀɴ�ӡ�ļ�
    FlagsTips = {'Hub_flag����ѡ�����ɽ��' 'Duct_flag����ѡ�����ɺ���' 'Chord_flag����ѡ������ҳ��Ż�'...
                 'Viscous_flag����ѡ����ճ������' 'Analyze_flag����ѡ�������������' 'Geometry_flag����ѡ����ʾҶƬģ��'...
                 'txt_flag����ѡ�����txt��ʽҶƬ������Ϣ' 'csv_flag����ѡ�����csv��ʽ����' 'stl_flag����ѡ�����stl��ʽģ��'};
    FlagsValues_def = {1 1 0 1 0 1 1 1 0};    
    %����ѡ��ʹ��'itemsdata'���Խ�ѡ���ı��ֱ���1��0��Ӧ��ѡ�����������򷵻�ֵΪ1
    FlagsValues(1) = uidropdown('parent',FlagsGrid,...
                                'items',SelectStrings,...
                                'itemsdata',[1 0],...
                                'fontsize',zhfontsize,...
                                'valuechangedfcn',@PropTurb);
    %����ѭ������9�и�ѡ��
    for index = 1 : length(FlagsStrings)
        FlagsValues(index+1) = uicheckbox('parent',FlagsGrid,...
                                          'text',FlagsStrings{index},...
                                          'value',FlagsValues_def{index},...
                                          'fontsize',zhfontsize,...
                                          'tooltip',FlagsTips{index});
    end
    %���ø�ѡ��Ļص�����
    set(FlagsValues(3),'valuechangedfcn',@ChangeDuct);
    set(FlagsValues(4),'valuechangedfcn',@ChangeChord);
    set(FlagsValues(5),'valuechangedfcn',@ChangeViscous);
    %% �����ߡ������е�����
    %����Tools��2*4������ռ�
    ToolsGrid = uigridlayout(Tools,[2 4],...
                             'columnwidth',{'1x','1x','1x','1x'},...
                             'rowheight',{'2x','3x'},...
                             'columnspacing',25,...
                             'rowspacing',25,...
                             'scrollable','off');
    %�����ļ����ı���
    ToolsText = uilabel('parent',ToolsGrid,...
                        'text','��ǰ��Ŀ���ƣ�',...
                        'fontsize',zhfontsize,...
                        'horizontalalignment','center',...
                        'verticalalignment','center');
    %���Ʊ༭��
    Filename = uieditfield(ToolsGrid,...
                           'fontsize',numfontsize,...
                           'horizontalalignment','center',...
                           'value',filename);
    Filename.Layout.Column = [2 4];
    %���Ƽ��ذ�ť
    LoadButton = uibutton('parent',ToolsGrid,...
                          'text','',...
                          'tooltip','��.mat�ļ�����GUI',...
                          'icon','load.png',...
                          'fontsize',buttonfontsize,...
                          'buttonpushedfcn',@LoadData);
    %���Ʊ��水ť
    SaveButton = uibutton('parent',ToolsGrid,...
                          'text','',...
                          'tooltip','��GUI������.mat�ļ�',...
                          'icon','save.png',...
                          'fontsize',buttonfontsize,...
                          'buttonpushedfcn',@SaveData);
    %�������Ϊ��ť
    SaveCopyButton = uibutton('parent',ToolsGrid,...
                              'text','',...
                              'tooltip','��GUI���Ϊ.mat�ļ�',...
                              'icon','savecopy.png',...
                              'fontsize',buttonfontsize,...
                              'buttonpushedfcn',@SaveDataCopy);
    %�������а�ť
    RunButton = uibutton('parent',ToolsGrid,...
                         'text','',...
                         'tooltip','����OpenProp',...
                         'icon','run.png',...
                         'fontsize',buttonfontsize,...
                         'buttonpushedfcn',@Execute);
end
%% �л���������ҶƬ��ơ��͡��о�������Ӱ�족��ص��ĺ���
function ToggleMode(hObject,ED)
    global Tabs SingleDesign ParametricStudy;
    global Filename;
    global SpecificationsValues DuctValues FlagsValues FoilValues CalculatorsValues...
           Xr_in XCpoDp_in XCdp_in VAI_in VTI_in ri_in Xt0oDp_in skew0_in rake0_in...
           Meanline_cell Thickness_cell XCpoDp_values XCLmax_values;
    if get(Tabs,'selectedtab') == SingleDesign
        %��ǰѡ���ģʽΪ��������ҶƬ��ơ���������в���
    elseif get(Tabs,'selectedtab') == ParametricStudy
        %��ǰѡ���ģʽΪ���о�������Ӱ�족����Ҫ����������ҶƬ��ơ�����������ݴ�����ʱ�ļ���
        filename = get(Filename,'value');
        Z = get(SpecificationsValues(1),'value');       %ҶƬ����
        N = get(SpecificationsValues(2),'value');       %������ת��(RPM)
        Dp = get(SpecificationsValues(3),'value');      %������ֱ��(m)
        Thrust = get(SpecificationsValues(4),'value');      %װ��������Ҫ��(N)
        Vs = get(SpecificationsValues(5),'value');      %װ���н��ٶ�(m/s)
        Dhub = get(SpecificationsValues(6),'value');        %���ֱ��(m)
        rho = get(SpecificationsValues(7),'value');     %�����ܶ�(kg/m^3)
        Mp = get(SpecificationsValues(8),'value');      %����ҶƬ��������
        Np = get(SpecificationsValues(9),'value');      %����ҶƬ��������
        Meanline_index = get(FoilValues(1),'value');
        Meanline = char(Meanline_cell(Meanline_index));     %͹Ե������
        Thickness_index	= get(FoilValues(2),'value');
        Thickness = char(Thickness_cell(Thickness_index));      %Ҷ������
        TpoT = get(DuctValues(1),'value');      %����������ռ��
        Cdd = get(DuctValues(2),'value');       %������������ϵ��
        Rduct_oR = get(DuctValues(3),'value');     %����ֱ��/������ֱ��
        Js = get(CalculatorsValues(1),'value');      %����ϵ��
        L = get(CalculatorsValues(2),'value');       %Ҷ���ٶ�/�н��ٶ�
        KT = get(CalculatorsValues(3),'value');      %����ϵ��(Ҷ���ٶ�)
        CT = get(CalculatorsValues(4),'value');      %����ϵ��(�н��ٶ�)
        Propeller_flag = get(FlagsValues(1),'value');       %������/����ҶƬ
        Hub_flag = get(FlagsValues(2),'value');     %�Ƿ�������
        Duct_flag = get(FlagsValues(3),'value');        %�Ƿ�Ϊ����������
        Chord_flag = get(FlagsValues(4),'value');       %�Ƿ�����ҳ��Ż�
        Viscous_flag = get(FlagsValues(5),'value');     %�Ƿ���ճ������
        Analyze_flag = get(FlagsValues(6),'value');     %�Ƿ������������
        Geometry_flag = get(FlagsValues(7),'value');        %�Ƿ����ҶƬģ��
        txt_flag = get(FlagsValues(8),'value');     %�Ƿ񵼳������ĵ�
        csv_flag = get(FlagsValues(9),'value');     %�Ƿ񵼳�������
        stl_flag = get(FlagsValues(10),'value');     %�Ƿ񵼳�stl�ļ�   
        Xr = get(Xr_in,'value');        %����λ��
        XCpoDp = get(XCpoDp_in,'value');        %�Ҿ���
        XCLmax = get(XCpoDp_in,'value');      %�������ϵ��
        XCdp = get(XCdp_in,'value');      %����ϵ��
        Xt0oDp = get(Xt0oDp_in,'value');      %�񾶱�
        skew0 = get(skew0_in,'value');      %�����
        rake0 = get(rake0_in,'value');      %б����
        ri = get(ri_in,'value');        %����λ��
        VAI = get(VAI_in,'value');      %�����ֵ
        VTI = get(VTI_in,'value');      %�����ֵ
        ri = ri (~isnan(ri));
        VAI = VAI(~isnan(VAI));
        VTI = VTI(~isnan(VTI));
        save('OpenPropTempFile','Z','N','Dp','Thrust','Vs','Dhub','rho','Mp','Np',...
             'Meanline','Meanline_index','Thickness','Thickness_index',...
             'TpoT','Cdd','Rduct_oR','Js','L','KT','CT',...
             'Propeller_flag','Hub_flag','Duct_flag','Chord_flag',...
             'Viscous_flag','Analyze_flag','Geometry_flag',...
             'txt_flag','csv_flag','stl_flag',...
             'filename','Xr','XCpoDp','XCdp','XCLmax','Xt0oDp','skew0','rake0',...
             'ri','VAI','VTI','XCpoDp_values','XCLmax_values');
        OpenPropParam;
    end
end
%% ��ҶƬ����в����ı����õĺ���
function UpdateSpecifications(hObject,ED)
    %���¡������ٲ������еı���ֵ���롰ҶƬ����еĲ�������һ��
    global SpecificationsValues CalculatorsValues Fig_Main;
    Z = round(get(SpecificationsValues(1),'value'));        %ҶƬ����
    set(SpecificationsValues(1),'value',Z);
    N = get(SpecificationsValues(2),'value');      %������ת��(RPM)
    Dp = get(SpecificationsValues(3),'value');      %������ֱ��(m)
    Thrust = get(SpecificationsValues(4),'value');     %װ��������Ҫ��(N)
    Vs = get(SpecificationsValues(5),'value');     %װ���н��ٶ�(m/s)
    Dhub = get(SpecificationsValues(6),'value');       %���ֱ��(m)
    rho = get(SpecificationsValues(7),'value');        %�����ܶ�(kg/m^3)
    Mp = round(get(SpecificationsValues(8),'value'));      %����ҶƬ��������
    set(SpecificationsValues(8),'value',Mp);
    Np = round(get(SpecificationsValues(9),'value'));      %����ҶƬ��������
    set(SpecificationsValues(9),'value',Np);
    %��齰�ֱ����������ֱ��֮��Ĺ�ϵ����ǰ�ߴ��ں��ߣ��򵯳����ѿ�
    if Dhub >= Dp
        message = sprintf('��ǰ���õĽ��ֱ������������ֱ�� \n ���������ã�');
        uialert(Fig_Main,message,'�������ô���','closefcn',@ResetDiameters);
    end
    n = N/60;       %������ת��(r/s)
    lambda = n*2*pi*(Dp/2)/Vs;       %�������ٶ�=Ҷ��Բ���ٶ�/װ���н��ٶ�(Non-dimensional)
    Js = Vs/(n*Dp);      %����ϵ��
    KT = Thrust/(rho*n^2*Dp^4);      %����ϵ��(Thrust Coeffcient)����Ҷ���ٶ�Ϊ�����ٶ�
    CT = Thrust/(0.5*rho*Vs^2*pi*(Dp/2)^2);      %����ϵ��(Thrust Coefficient)����װ���н��ٶ�Ϊ�����ٶ�
    set(CalculatorsValues(1),'value',Js);
    set(CalculatorsValues(2),'value',lambda);
    set(CalculatorsValues(3),'value',KT);
    set(CalculatorsValues(4),'value',CT);
end
%����ֱ�����ô������õĺ���
function ResetDiameters(hObject,ED)
    global SpecificationsValues;
    set(SpecificationsValues(3),'value',0.4);
    set(SpecificationsValues(6),'value',0.02);
end
%% �޸ġ���������/������ҶƬ����ť��ص��ĺ���
function PropTurb(hObject,ED)
    global FlagsValues SpecificationsValues;
    if get(FlagsValues(1),'value')
        set(SpecificationsValues(4),'editable','on')
    else
        %ѡ������ҶƬ����ť�󣬡�ҶƬ��������µġ�װ��������Ҫ�󡱱�Ϊ���ɱ༭
        set(SpecificationsValues(4),'editable','off')
    end
end
%% �޸ġ��Ƿ�Ϊ��������������ť��ص��ĺ���
function ChangeDuct(hObject,ED)
    global DuctValues FlagsValues;
    if get(FlagsValues(3),'value')
        set(DuctValues,'editable','on');
    else
        set(DuctValues,'editable','off');
    end
end
%% �޸ġ��Ƿ�����ҳ��Ż�����ť��ص��ĺ���
function ChangeChord(hObject,ED)
    global FlagsValues XCpoDp_in Col_Label;
    global XCpoDp_values XCLmax_values;
    if get(FlagsValues(4),'value')
        set(Col_Label(2),'text','�������ϵ��');
        for index = 1:length(XCpoDp_in)
            XCpoDp_values(index) = get(XCpoDp_in(index),'value');
            set(XCpoDp_in(index),'value',XCLmax_values(index));
        end
    else
        set(Col_Label(2),'text','�Ҿ���');
        for index = 1:length(XCpoDp_in)
            XCLmax_values(index) = get(XCpoDp_in(index),'value');
            set(XCpoDp_in(index),'value',XCpoDp_values(index));
        end        
    end
end
%% �޸ġ��Ƿ���ճ����������ť��ص��ĺ���
function ChangeViscous(hObject,ED)
    global FlagsValues XCdp_in XCdp_values;
    if get(FlagsValues(5),'value')
        set(XCdp_in,'editable','on');
        for index = 1 : length(XCdp_values)
            set(XCdp_in(index),'value',XCdp_values(index));
        end
    else
        for index = 1 : length(XCdp_in)
            XCdp_values(index) = get(XCdp_in(index),'value');
        end
        set(XCdp_in,'editable','off','value',0);
    end
end
%% ��������ء���ص��ĺ���
function LoadData(hObject,ED)
    global N_R0 OpenPropDirectory SpecificationsValues DuctValues FlagsValues FoilValues Filename...
           Xr_in XCpoDp_in XCdp_in VAI_in VTI_in ri_in Xt0oDp_in skew0_in rake0_in;
    global pt
    rest = pwd;
    LoadDir(rest,OpenPropDirectory);
    %uiload�������ڴ��ļ�ѡ��Ŀ¼
    uiload;
    set(SpecificationsValues(1),'value',pt.input.Z);
    set(SpecificationsValues(2),'value',pt.input.N);
    set(SpecificationsValues(3),'value',pt.input.Dp);
    set(SpecificationsValues(4),'value',pt.input.Thrust);
    set(SpecificationsValues(5),'value',pt.input.Vs);
    set(SpecificationsValues(6),'value',pt.input.Dhub);
    set(SpecificationsValues(7),'value',pt.input.rho);
    set(SpecificationsValues(8),'value',pt.input.Mp);
    set(SpecificationsValues(9),'value',pt.input.Np);
    set(DuctValues(1),'value',pt.input.TpoT);
    set(DuctValues(2),'value',pt.input.Cdd);
    set(DuctValues(3),'value',pt.input.Rduct_oR);
    set(FlagsValues(1),'value',pt.input.Propeller_flag);
    set(FlagsValues(2),'value',pt.input.Hub_flag);
    set(FlagsValues(3),'value',pt.input.Duct_flag);
    set(FlagsValues(4),'value',pt.input.Chord_flag);
    set(FlagsValues(5),'value',pt.input.Viscous_flag);
    set(FlagsValues(6),'value',pt.input.Analyze_flag);
    set(FlagsValues(7),'value',pt.input.Geometry_flag);
    set(FlagsValues(8),'value',pt.input.txt_flag);
    set(FlagsValues(9),'value',pt.input.csv_flag);
    set(FlagsValues(10),'value',pt.input.stl_flag);
    set(FoilValues(1),'value',pt.input.Meanline_index);
    set(FoilValues(2),'value',pt.input.Thickness_index);
    set(Filename,'value',pt.filename);
    for index = 1 : N_R0
        set(Xr_in(index),'value',pt.input.Xr(index));
        set(XCpoDp_in(index),'value',pt.input.XCpoDp(index));
        set(XCdp_in(index),'value',pt.input.XCdp(index));
        set(Xt0oDp_in(index),'value',pt.input.Xt0oDp(index));
        set(skew0_in(index),'value',pt.input.skew0(index));
        set(rake0_in(index),'value',pt.input.rake0(index));
    end    
    for index = 1 : N_R0
        if isnan(pt.input.ri(index))
            set(ri_in(index),'value','');
        else
            set(ri_in(index),'value',pt.input.ri(index));
        end
        if isnan(pt.input.VAI(index))
            set(VAI_in(index),'value',1);
        else
            set(VAI_in(index),'value',pt.input.VAI(index));
        end
        if isnan(pt.input.VTI(index))
            set(VTI_in(index),'value',0);
        else
            set(VTI_in(index),'value',pt.input.VTI(index));
        end
    end
end
%% ��������桱��ص��ĺ���
function SaveData(hObject,ED)
    global OpenPropDirectory SpecificationsValues DuctValues CalculatorsValues FlagsValues FoilValues Filename...
           Xr_in XCpoDp_in XCdp_in VAI_in VTI_in ri_in Xt0oDp_in skew0_in rake0_in...
           Meanline_cell Thickness_cell;
    global Fig_Main pt;
    %% ��GUI�����ж�ȡ����ֵ
    filename = get(Filename,'value');       %��ǰ��Ŀ����
    %existence��ʾ�Ƿ����ͬ����Ŀ
    existence = ChangeDirectory(OpenPropDirectory,filename);
    if existence
        message = sprintf('��ǰ��Ŀ���Ѿ����� \n ����������иò�����������ԭ���ļ� \n ȷ��������');
        SaveSelection = uiconfirm(Fig_Main,message,'���Ǿ���',...
                                  'options',{'����ԭ���ļ�','���Ϊ������','ȡ������'},...
                                  'defaultoption',2,...
                                  'canceloption',3,...
                                  'icon','warning');
        if strcmp(SaveSelection,'����ԭ���ļ�')
            %ѡ�񡰸���ԭʼ�ļ���,����Ҫ���в���
        elseif strcmp(SaveSelection,'���Ϊ������')
            %ѡ�����Ϊ�����ơ�������savecopy����
            SaveDataCopy;
        else
            %ѡ��ȡ�����������˳�
            return
        end    
    end    
    Z = get(SpecificationsValues(1),'value');       %ҶƬ����
    N = get(SpecificationsValues(2),'value');       %������ת��(RPM)
    Dp = get(SpecificationsValues(3),'value');      %������ֱ��(m)
    Thrust = get(SpecificationsValues(4),'value');      %װ��������Ҫ��(N)
    Vs = get(SpecificationsValues(5),'value');      %װ���н��ٶ�(m/s)
    Dhub = get(SpecificationsValues(6),'value');        %���ֱ��(m)
    rho = get(SpecificationsValues(7),'value');     %�����ܶ�(kg/m^3)
    Mp = get(SpecificationsValues(8),'value');      %����ҶƬ��������
    Np = get(SpecificationsValues(9),'value');      %����ҶƬ��������
    Meanline_index = get(FoilValues(1),'value');
    Meanline = char(Meanline_cell(Meanline_index));     %͹Ե������
    Thickness_index	= get(FoilValues(2),'value');
    Thickness = char(Thickness_cell(Thickness_index));      %Ҷ������
    TpoT = get(DuctValues(1),'value');      %����������ռ��
    Cdd = get(DuctValues(2),'value');       %������������ϵ��
    Rduct_oR = get(DuctValues(3),'value');     %����ֱ��/������ֱ��
    Js = get(CalculatorsValues(1),'value');      %����ϵ��
    L = get(CalculatorsValues(2),'value');       %Ҷ���ٶ�/�н��ٶ�
    KT = get(CalculatorsValues(3),'value');      %����ϵ��(Ҷ���ٶ�)
    CT = get(CalculatorsValues(4),'value');      %����ϵ��(�н��ٶ�)
    Propeller_flag = get(FlagsValues(1),'value');       %������/����ҶƬ
    Hub_flag = get(FlagsValues(2),'value');     %�Ƿ�������
    Duct_flag = get(FlagsValues(3),'value');        %�Ƿ�Ϊ����������
    Chord_flag = get(FlagsValues(4),'value');       %�Ƿ�����ҳ��Ż�
    Viscous_flag = get(FlagsValues(5),'value');     %�Ƿ���ճ������
    Analyze_flag = get(FlagsValues(6),'value');     %�Ƿ������������
    Geometry_flag = get(FlagsValues(7),'value');        %�Ƿ����ҶƬģ��
    txt_flag = get(FlagsValues(8),'value');     %�Ƿ񵼳������ĵ�
    csv_flag = get(FlagsValues(9),'value');     %�Ƿ񵼳�������
    stl_flag = get(FlagsValues(10),'value');     %�Ƿ񵼳�stl�ļ�   
    Xr = cell2array(get(Xr_in,'value'));        %����λ��
    XCpoDp = cell2array(get(XCpoDp_in,'value'));        %�Ҿ���
    XCLmax = cell2array(get(XCpoDp_in,'value'));      %�������ϵ��
    XCdp = cell2array(get(XCdp_in,'value'));      %����ϵ��
    Xt0oDp = cell2array(get(Xt0oDp_in,'value'));      %�񾶱�
    skew0 = cell2array(get(skew0_in,'value'));      %�����
    rake0 = cell2array(get(rake0_in,'value'));      %б����
    ri = cell2array(get(ri_in,'value'));        %����λ��
    VAI = cell2array(get(VAI_in,'value'));      %�����ֵ
    VTI = cell2array(get(VTI_in,'value'));      %�����ֵ
    ITER = 40;      %��������
    Rhv = 0.5;      %hub vortex radius / hub radius
    %% ��GUI���������ֵ�������ṹ������input
    input.part1 = 'ҶƬ���';
    input.Z = Z;
    input.N = N;
    input.Dp = Dp;
    input.Thrust = Thrust;
    input.Vs = Vs;
    input.Dhub = Dhub;
    input.rho = rho;
    input.Mp = Mp;
    input.Np = Np;
    input.Meanline = Meanline;
    input.Meanline_index = Meanline_index;
    input.Thickness = Thickness;
    input.Thickness_index = Thickness_index;
    input.part2 = 'ҶƬ��Ʋ���';
    input.Xr = Xr;
    input.XCdp = XCdp;
    input.XCLmax = XCLmax;
    input.XCpoDp = XCpoDp;
    input.Xt0oDp = Xt0oDp;
    input.skew0 = skew0;
    input.rake0 = rake0;
    input.part3 = '�����ٶ�/װ���н��ٶ�';
    input.ri = ri;
    input.VAI = VAI;
    input.VTI = VTI;
    input.part4 = '���ѡ��';
    input.Propeller_flag = Propeller_flag;
    input.Hub_flag = Hub_flag;
    input.Duct_flag = Duct_flag;
    input.Chord_flag = Chord_flag;
    input.Viscous_flag = Viscous_flag;
    input.Analyze_flag = Analyze_flag;
    input.Geometry_flag = Geometry_flag;
    input.txt_flag = txt_flag;
    input.csv_flag = csv_flag;
    input.stl_flag = stl_flag;
    input.part5 = '����������';
    input.TpoT = TpoT;
    input.Cdd = Cdd;
    input.Rduct_oR = Rduct_oR;
    input.part6 = '�����ٲ���';
    input.Js = Js;
    input.L = L;
    input.KT = KT;
    input.CT = CT;
    input.part7 = '�������������';
    input.ITER = ITER;
    input.Rhv = Rhv;
    %% ��������Ŀ��Ϣ��װ���ṹ������pt��
    pt.filename = filename;
    pt.date = date;
    pt.input = input;
    pt.design = [];     %�ṹ�����飺design conditions
    pt.geometry = [];       %�ṹ�����飺design geometry
    pt.states = [];     %�ṹ�����飺off-design state analysis
    %% ����
    save(filename,'pt');
end
%% ��������Ϊ����ص��ĺ���
function SaveDataCopy(hObject,ED)
    global Fig_Main SpecificationsValues DuctValues CalculatorsValues FlagsValues FoilValues Filename...
           Xr_in XCpoDp_in XCdp_in VAI_in VTI_in ri_in Xt0oDp_in skew0_in rake0_in...
           Meanline_cell Thickness_cell;
    global pt;
    %% ��GUI�����ж�ȡ����ֵ
    filename = get(Filename,'value');       %��ǰ��Ŀ����
    Z = get(SpecificationsValues(1),'value');       %ҶƬ����
    N = get(SpecificationsValues(2),'value');       %������ת��(RPM)
    Dp = get(SpecificationsValues(3),'value');      %������ֱ��(m)
    Thrust = get(SpecificationsValues(4),'value');      %װ��������Ҫ��(N)
    Vs = get(SpecificationsValues(5),'value');      %װ���н��ٶ�(m/s)
    Dhub = get(SpecificationsValues(6),'value');        %���ֱ��(m)
    rho = get(SpecificationsValues(7),'value');     %�����ܶ�(kg/m^3)
    Mp = get(SpecificationsValues(8),'value');      %����ҶƬ��������
    Np = get(SpecificationsValues(9),'value');      %����ҶƬ��������
    Meanline_index = get(FoilValues(1),'value');
    Meanline = char(Meanline_cell(Meanline_index));     %͹Ե������
    Thickness_index	= get(FoilValues(2),'value');
    Thickness = char(Thickness_cell(Thickness_index));      %Ҷ������
    TpoT = get(DuctValues(1),'value');      %����������ռ��
    Cdd = get(DuctValues(2),'value');       %������������ϵ��
    Rduct_oR = get(DuctValues(3),'value');     %����ֱ��/������ֱ��
    Js = get(CalculatorsValues(1),'value');      %����ϵ��
    L = get(CalculatorsValues(2),'value');       %Ҷ���ٶ�/�н��ٶ�
    KT = get(CalculatorsValues(3),'value');      %����ϵ��(Ҷ���ٶ�)
    CT = get(CalculatorsValues(4),'value');      %����ϵ��(�н��ٶ�)
    Propeller_flag = get(FlagsValues(1),'value');       %������/����ҶƬ
    Hub_flag = get(FlagsValues(2),'value');     %�Ƿ�������
    Duct_flag = get(FlagsValues(3),'value');        %�Ƿ�Ϊ����������
    Chord_flag = get(FlagsValues(4),'value');       %�Ƿ�����ҳ��Ż�
    Viscous_flag = get(FlagsValues(5),'value');     %�Ƿ���ճ������
    Analyze_flag = get(FlagsValues(6),'value');     %�Ƿ������������
    Geometry_flag = get(FlagsValues(7),'value');        %�Ƿ����ҶƬģ��
    txt_flag = get(FlagsValues(8),'value');        %�Ƿ񵼳������ĵ�
    csv_flag = get(FlagsValues(9),'value');        %�Ƿ񵼳�������
    stl_flag = get(FlagsValues(10),'value');     %�Ƿ񵼳�stl�ļ�   
    Xr = cell2array(get(Xr_in,'value'));        %����λ��
    XCpoDp = cell2array(get(XCpoDp_in,'value'));        %�Ҿ���
    XCLmax = cell2array(get(XCpoDp_in,'value'));      %�������ϵ��
    XCdp = cell2array(get(XCdp_in,'value'));      %����ϵ��
    Xt0oDp = cell2array(get(Xt0oDp_in,'value'));      %�񾶱�
    skew0 = cell2array(get(skew0_in,'value'));      %�����
    rake0 = cell2array(get(rake0_in,'value'));      %б����
    ri = cell2array(get(ri_in,'value'));        %����λ��
    VAI = cell2array(get(VAI_in,'value'));      %�����ֵ
    VTI = cell2array(get(VTI_in,'value'));      %�����ֵ
    ITER = 40;      %��������
    Rhv = 0.5;      %hub vortex radius / hub radius
    %% ��GUI���������ֵ�������ṹ������input
    input.part1 = 'ҶƬ���';
    input.Z = Z;
    input.N = N;
    input.Dp = Dp;
    input.Thrust = Thrust;
    input.Vs = Vs;
    input.Dhub = Dhub;
    input.rho = rho;
    input.Mp = Mp;
    input.Np = Np;
    input.Meanline = Meanline;
    input.Meanline_index = Meanline_index;
    input.Thickness = Thickness;
    input.Thickness_index = Thickness_index;
    input.part2 = 'ҶƬ��Ʋ���';
    input.Xr = Xr;
    input.XCdp = XCdp;
    input.XCLmax = XCLmax;
    input.XCpoDp = XCpoDp;
    input.Xt0oDp = Xt0oDp;
    input.skew0 = skew0;
    input.rake0 = rake0;
    input.part3 = '�����ٶ�/װ���н��ٶ�';
    input.ri = ri;
    input.VAI = VAI;
    input.VTI = VTI;
    input.part4 = '���ѡ��';
    input.Propeller_flag = Propeller_flag;
    input.Hub_flag = Hub_flag;
    input.Duct_flag = Duct_flag;
    input.Chord_flag = Chord_flag;
    input.Viscous_flag = Viscous_flag;
    input.Analyze_flag = Analyze_flag;
    input.Geometry_flag = Geometry_flag;
    input.txt_flag = txt_flag;
    input.csv_flag = csv_flag;
    input.stl_flag = stl_flag;
    input.part5 = '����������';
    input.TpoT = TpoT;
    input.Cdd = Cdd;
    input.Rduct_oR = Rduct_oR;
    input.part6 = '�����ٲ���';
    input.Js = Js;
    input.L = L;
    input.KT = KT;
    input.CT = CT;
    input.part7 = '�������������';
    input.ITER = ITER;
    input.Rhv = Rhv;
    %% ��������Ŀ��Ϣ��װ���ṹ������pt��
    pt.filename = [filename,'-copy'];
    pt.date = date;
    pt.input = input;
    pt.design = [];     %�ṹ�����飺design conditions
    pt.geometry = [];       %�ṹ�����飺design geometry
    pt.states = [];     %�ṹ�����飺off-design state analysis
    %Ŀǰ��������Ϊ��ɺ���Ļ�����ص�matlab�����������Ҫ���
    uisave('pt',[filename,'-copy']);
end
%��ֱ�Ӵ�GUI��ȡ��cellת��Ϊdouble����Ĺ��ߺ���
function array = cell2array(cell)
    celllen = length(cell);
    array = zeros(celllen,1);
    for index = 1 : celllen
        array(index) = cell{index};
    end
end
%% ִ�в���
function Execute(hObject,ED)
    global PlotsPanels NumbersValues;
    global pt
    %% �����������
    SaveData;
    %% �����½���
    newPlots;
    %% ��������ת��
    Xr = pt.input.Xr;
    XCdp = pt.input.XCdp;
    XCLmax = pt.input.XCLmax;
    XCpoDp = pt.input.XCpoDp;
    Xt0oDp = pt.input.Xt0oDp;
    skew0 = pt.input.skew0;
    rake0 = pt.input.rake0;
    ri = pt.input.ri;
    VAI = pt.input.VAI;
    VTI = pt.input.VTI;
    Propeller_flag = pt.input.Propeller_flag;
    Analyze_flag = pt.input.Analyze_flag;
    %% ͨ��EppsOptimizer�����õ�pt.design��ֵ
    pt.design = EppsOptimizer(pt.input);
    %������������Ť�أ�����CQΪŤ��ϵ��
    pt.design.Q = pt.design.CQ * 0.5 * pt.input.rho * pt.input.Vs^2 * pi*pt.input.Dp^2/4 * pt.input.Dp/2; % [Nm]  torque
    %������������ת��
    omega = 2*pi*pt.input.N/60; % [rad/s]
    %�������������ĵĹ���
    pt.design.P = pt.design.Q * omega;
    if Propeller_flag == 1
        set(NumbersValues(7),'value',pt.design.Js);
        set(NumbersValues(8),'value',pt.design.KT);
        set(NumbersValues(10),'value',pt.design.KQ);
        set(NumbersValues(13),'value',pt.design.EFFY);
        set(NumbersValues(14),'value',pt.design.ADEFFY);
    else
        set(NumbersValues(7),'value',pi/pt.design.L);
        set(NumbersValues(8),'value',' ');
        set(NumbersValues(10),'value',' ');
        set(NumbersValues(13),'value',' ');
        set(NumbersValues(14),'value',' ');
    end
    set(NumbersValues(1),'value',pt.input.Vs);
    set(NumbersValues(2),'value',pt.input.N);
    set(NumbersValues(3),'value',pt.input.Dp);
    set(NumbersValues(4),'value',pt.input.Thrust);
    set(NumbersValues(5),'value',pt.design.Q);
    set(NumbersValues(6),'value',pt.design.P);
    set(NumbersValues(9),'value',pt.design.CT);
    set(NumbersValues(11),'value',pt.design.CQ);
    set(NumbersValues(12),'value',pt.design.CP);
    %% ���ƽ������
    Make_Reports;
    %% ����ͼ��
    pt.geometry = Geometry(pt);
    %% �����Ƿ������������
    if Analyze_flag
        pt.states = AnalyzeAuto(pt);
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotsPanels(15));
        axes(h);
        hold on;
        if Propeller_flag == 1
            VMIV = pt.design.VMIV;
            Js_curve = linspace(min(pt.states.Js),max(pt.states.Js),100);
            EFFY_curve = (Js_curve/(2*pi)) * VMIV .* pchip(pt.states.Js,pt.states.KT,Js_curve)./pchip(pt.states.Js,pt.states.KQ,Js_curve); 
            % Efficiency (green squares)
            plot(Js_curve,EFFY_curve,'-','LineWidth',2,'Color',[0 0.8 0])
            Heffy = plot(pt.states.Js,pt.states.EFFY,'sk','MarkerSize',5,'LineWidth',1,'MarkerFaceColor',[0 0.8 0]); 
            % Thrust coefficient (blue diamonds)
            plot(pt.states.Js,pt.states.KT,'b-','LineWidth',2)
            Hkt   = plot(pt.states.Js,pt.states.KT,'dk','MarkerSize',5,'LineWidth',1,'MarkerFaceColor','b');
            % Torque coefficient (red circles)
            plot(pt.states.Js,10*pt.states.KQ,'r-','LineWidth',2)    
            Hkq = plot(pt.states.Js,10*pt.states.KQ,'ok','MarkerSize',5,'LineWidth',1,'MarkerFaceColor','r');
            % Design point
            plot(pt.design.Js*[1 1],[0 2],'k--','LineWidth',1);
            xlabel('Js','FontSize',16,'FontName','Times'), 
            ylabel('KT, 10*KQ, EFFY','FontSize',16,'FontName','Times')
            axis([min(pt.states.Js) max(pt.states.Js) 0 0.9])
            set(gca,'FontSize',14,'FontName','Times');
            box on, grid on,
            ylimits = get(gca,'Ylim');
            set(gca,'Ylim', [0  max(1,ylimits(2)) ] );
        else
            % Power coefficient (blue dots)
            plot(pt.states.L,-pt.states.CP,'b.-','LineWidth',2,'MarkerSize',12)
            % Design point
            plot(pt.design.L*[1 1],[0 0.6],'k--','LineWidth',1);
            % % Betz limit
            % plot([0 ceil(max(pt.states.L))],(16/27)*[1 1],'k--','LineWidth',1);
            xlabel('L','FontSize',16,'FontName','Times'), 
            ylabel('CP','FontSize',16,'FontName','Times'),
            set(gca,'Ytick',[0:0.1:0.6])
            axis([0 ceil(max(pt.states.L))  0 0.6])
            set(gca,'FontSize',14,'FontName','Times');
            box on, grid on,
        end
    end

    %% �����ǹ��ڲ����о��Ĳ���
    if exist([filename '.mat'],'file')        
                CLR = [     1       0       0;      ... % (1) Red
                            0       0.9     0;      ... % (2) Green
                            0       0       1;      ... % (3) Blue
                            0.75    0       0.75;   ... % (4) Purple
                            1       0.5     0;      ... % (5) Orange
                            0       1       1;      ... % (6) Cyan
                            1       0       1;      ... % (7) Magenta
                            0.75    0.5     0.25;   ... % (8) Brown
                            0.25    0.25    0.75;   ... % (9) Navy blue
                            0.25    0.5     0.75;   ... % (10) Steel blue
                            0.75    0.75    0];         % (11) Burnt Yellow        

        temp    = pt;
        load([filename '.mat']);
        if isfield(pt,'paroutput')
            paroutput   = pt.paroutput;
            parinput    = pt.parinput;
            % --- For EFFY vs N ---
            set(0,'CurrentFigure',Plots);
            h = axes('parent',PlotsPanels(4));
            axes(h);
            hold on, box on, grid on,
            ylim([0 1]),
            tempEFFY = paroutput.EFFY(1,:,1);
            plot(paroutput.N,tempEFFY(:),'-','color',CLR( (mod(0,11)+1) ,:));
            str_legend ={[num2str(parinput.Dp),' m  ']};
            legend(str_legend,'location','southwest');
            xlabel('Rotation Speed (RPM)  ');
            ylabel('Efficiency');
            title(['Number of Blades: ',num2str(Z),'  '])
            % --- For EFFY vs D ---
            set(0,'CurrentFigure',Plots);
            h = axes('parent',PlotsPanels(5));
            axes(h);
            hold on, box on, grid on,
            ylim([0 1]),
            tempEFFY = paroutput.EFFY(1,1,:);
            plot(paroutput.Dp,tempEFFY(:),'-','color',CLR( (mod(0,11)+1) ,:));
            str_legendD ={[num2str(parinput.N),' RPM  ']};
            legend(str_legendD,'location','southwest');
            xlabel('Propeller Diameter (m)  ');
            ylabel('Efficiency');
            title(['Number of Blades: ',num2str(Z),'  '])
        else
            set(Toggle(4),'enable','off');
            set(Toggle(5),'enable','off');
        end
        pt = temp;
    else
        set(Toggle(4),'enable','off');
        set(Toggle(5),'enable','off');
    end
    % pt overwrite sequence:
    if exist([filename '.mat'],'file')
        disp(['Found original file:  ',filename,'.mat']);
        temp  = pt;
        load(filename)
        disp(['Overwriting file:  ',filename,'.mat']);
        pt.filename = temp.filename;
        pt.date     = temp.date;
        pt.input	= temp.input;
        pt.design   = temp.design;
        pt.geometry = temp.geometry;
        pt.states   = temp.states;
    end
    pt
    save(filename,'pt');
    save([filename,'_GUIinput'],'Z','N','Dp','Thrust','Vs','Dhub','rho','Mp','Np',...
        'TpoT','Cdd','Propeller_flag','Hub_flag','Duct_flag',...
        'Chord_flag','Viscous_flag','Geometry_flag','Analyze_flag',...
        'Meanline','Meanline_index','Thickness',...
        'Thickness_index','filename','Xr','XCpoDp','XCdp','ri','VAI','VTI','Xt0oDp','skew0',...
        'rake0');
end
%% ����OpenProp���½��漸��Ԫ�ص����ɴ���
function newPlots
    global pt;
    global Tabs TabPageStrings TabPageMenu Menu;
    global Filename NumbersValues PlotsValues PlotsStrings PlotsPanels coordinate;
    global zhfontsize subcolumnfontsize numfontsize;
    global lc1_1 lc1_2 lc2_1 lc2_2 lc3_1 lc3_2 lc4_1 lc4_2;
    %����[0.933333 0.721568 0.764705]�������[0.709803 0.596078 0.631372]
    %���˻�[0.988235 0.631372 0.023529]������ʯ��[0.341176 0.584313 0.447058]
    %����[0.576470 0.709803 0.811764]��ɽ����[0.380392 0.392156 0.623529]
    lc1_1 = [0.047059 0.925490 0.866667];       %#0CECDD
    lc1_2 = [1.000000 0.952941 0.219608];       %#FFF338
    lc2_1 = [1.000000 0.403921 0.905882];       %#FF67E7
    lc2_2 = [0.768627 0.000000 1.000000];       %#C400FF
    lc3_1 = [0.627451 0.235294 0.470588];       %#A03C78
    lc3_2 = [0.929412 0.556863 0.486275];       %#ED8E7C
    lc4_1 = [0.803921 0.952941 0.635294];       %#CDF3A2
    lc4_2 = [0.576471 0.850980 0.639216];       %#93D9A3
    %% �����������ϵ�һ���±�ǩҳ
    filename = get(Filename,'value');
    DesignResult = uitab('parent',Tabs,...
                         'title',[filename,'����ƽ��']);
    %�л��������ɵĽ����ǩҳ
    set(Tabs,'selectedtab',DesignResult);
    %�������ɵĽ����ǩҳ������������ġ���ǩҳ��������
    TabPageStrings{length(TabPageStrings)+1} = [filename,'����ƽ��'];
    TabPageMenu(length(TabPageMenu)+1) = uimenu('parent',Menu(2),...
                                                'text',TabPageStrings{end});
    %% ����DesignResult�е�����ռ估�ռ��еĸ���Ŀ
    %����DesignResult��1*2������ռ�
    DesignResultGrid = uigridlayout(DesignResult,[1 2],...
                                    'columnwidth',{'2x','7x'});
    %���ɡ���ֵ�������
    Numbers = uipanel('parent',DesignResultGrid,...
                      'title','��ֵ���',...
                      'titleposition','centertop',...
                      'fontsize',subcolumnfontsize,...
                      'fontweight','bold');
    %���ɡ�ͼ��������
    Plots = uipanel('parent',DesignResultGrid,...
                    'title','ͼ����',...
                    'titleposition','centertop',...
                    'fontsize',subcolumnfontsize,...
                    'fontweight','bold');
    %% ����ֵ����������е�����
    %����Numbers��14*2������ռ�
    NumbersGrid = uigridlayout(Numbers,[14 2]);
    %Numbers�еĲ���
    NumbersStrings = {'װ���н��ٶȣ�' '������ת�٣�' '������ֱ����' 'װ����������' 'װ����Ť�أ�' 'װ���ܹ��ʣ�'...
                      '����ϵ����' '����ϵ��(Ҷ���ٶ�)��' '����ϵ��(�н��ٶ�)��' 'Ť��ϵ����' 'Ť��ϵ����' '����ϵ����'...
                      'EFFY��' 'ADEFFY��'};
    NumbersTips = {'Vs' 'N' 'Dp' 'Thrust' 'Torque' 'Power'...
                   'Js' 'KT' 'CT' 'KQ' 'CQ' 'CP'...
                   'EFFY' 'ADEFFY'};
    %����ѭ������14��
    for index = 1 : length(NumbersStrings)
        NumbersTexts(index) = uilabel('parent',NumbersGrid,...
                                      'text',NumbersStrings{index},...
                                      'horizontalalignment','left',...
                                      'verticalalignment','center',...
                                      'tooltip',NumbersTips{index},...
                                      'fontsize',zhfontsize);
        NumbersValues(index) = uieditfield(NumbersGrid,'numeric',...
                                           'fontsize',numfontsize,...
                                           'horizontalalignment','center',...
                                           'editable','off',...
                                           'value',0);
    end
    %����14�еı༭���е���ʾ��ʽ�͵�λ
    set(NumbersValues(1),'valuedisplayformat','%.2f m/s');
    set(NumbersValues(2),'valuedisplayformat','%.1f RPM');
    set(NumbersValues(3),'valuedisplayformat','%.2f m');
    set(NumbersValues(4),'valuedisplayformat','%.2f N');
    set(NumbersValues(5),'valuedisplayformat','%.2f N��m');
    set(NumbersValues(6),'valuedisplayformat','%.2f W');
    %% ��ͼ�����������е�����
    %����Plots��1*2������ռ�
    PlotsGrid = uigridlayout(Plots,[1 2],...
                             'columnwidth',{'1x','3x'});
    %ѡ���ͼ������Ҫʹ�õĲ���
    PlotsStrings = {'ҶƬ������������(��������)'...     %ԭ��'Expanded Blade (input blade)'
                    'ҶƬ��ȷֲ�����(��������)'...     %ԭ��'Blade Thickness (input blade)'
                    '�����ٶȷֲ�����'...       %ԭ��'Inflow Profile'
                    'Ч��-������ֱ��'...       %ԭ��'Efficiency vs Diameter'
                    'Ч��-������ת��'...       %ԭ��'Efficiency vs Rotation Speed'
                    '�����ֲ�����'...     %ԭ��'Circulation Distribution'
                    '�յ��ٶȷֲ�����'...       %ԭ��'Induced Velocity'
                    '��ؽǶȵķֲ�����'...      %ԭ��'Inflow Angle'
                    'ҶƬ������������(������ƽ��)'...       %ԭ��'Expanded Blade (as designed)'
                    'ҶƬ��ȷֲ�����(������ƽ��)'...       %ԭ��'Blade Thickness (as designed)'
                    '����ϵ���ֲ�����'...       %ԭ��'Lift Coefficient'
                    '��������'...      %ԭ��'Performance Curves'
                    'Ҷ�����ά����'...        %ԭ��'2D Geometry'
                    'ҶƬ��άģ��'};      %ԭ��'3D Geometry'
    %����14�е�ѡ��ť
    PlotsButtonBox = uibuttongroup('parent',PlotsGrid,...
                                   'title','ѡ����ʾ��ͼ��',...
                                   'titleposition','centertop',...
                                   'fontsize',zhfontsize,...
                                   'fontweight','bold',...
                                   'scrollable','on');
    ButtonBoxGrid = uigridlayout(PlotsButtonBox,[14 1]);
    for index = 1 : length(PlotsStrings)
        PlotsValues(index) = uibutton(ButtonBoxGrid,'state',...
                                      'text',PlotsStrings{index},...
                                      'fontsize',zhfontsize);
    end
    %���ø���ť�ص�����
    set(PlotsValues(1),'valuechangedfcn',@ExpandedBladeInfcn);
    set(PlotsValues(2),'valuechangedfcn',@BladeThicknessInfcn);
    set(PlotsValues(3),'valuechangedfcn',@InflowProfilefcn);
    %��Ч��-������ֱ�����͡�Ч��-������ת�١�������Ĭ�Ϲر�
    set(PlotsValues(4),'enable','off');
    set(PlotsValues(5),'enable','off');
    if pt.input.Chord_flag ~= 0
        %ѡ���ҳ��Ż�ʱ��ҶƬ������������(��������)������
        set(PlotsValues(1),'enable','off');
    end
    if pt.input.Analyze_flag == 0
        %û��ѡ�������������ʱ���������߲�����
        set(PlotsValues(12),'enable','off');
    end 
    if pt.input.Geometry_flag == 0
        %û��ѡ���������ͼ��ʱ����ά��������άģ�Ͳ�����
        set(PlotsValues(13),'enable','off');
        set(PlotsValues(14),'enable','off');
    end       
    %����ͼ�����
    PlotsPanels = uipanel('parent',PlotsGrid,...
                          'title','�����ఴť��ʾ��Ӧͼ��',...
                          'titleposition','centertop',...
                          'fontsize',zhfontsize,...
                          'fontweight','bold');
    PanelGrid = uigridlayout(PlotsPanels,[1 1]);
    %��������ϵ
    coordinate = uiaxes(PanelGrid);
end
%% ���°�ť��Ļص�����
%���¡�ҶƬ������������(��������)����Ļص�����
function ExpandedBladeInfcn(hObject,ED)
    global PlotsPanels PlotsStrings
    global PlotsValues
    global pt;
    global Fig_Main coordinate Xr_XCpoDp_In_line Xr_XCpoDp_De_line Xr_XCpoDp_In_point;
    global lc1_1 lc1_2;
    %% ����ͼ��������
    set(PlotsPanels,'title',PlotsStrings{1});
    %% ������Ӧ��ť�Ŀ���״̬
    %���ұ����µİ�ť����
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %�԰�ť1���в�������һ����ťҲû�б�����
        %���4��5������а�ť������
        set(PlotsValues,'enable','on');
        if pt.input.Analyze_flag == 0
            set(PlotsValues(12),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(13),'enable','off');
            set(PlotsValues(14),'enable','off');
        end
        set(PlotsValues(4),'enable','off');
        set(PlotsValues(5),'enable','off');
    else
        %�԰�ť1���в����������ڱ����µİ�ť
        %���1��9������а�ť��������
        set(PlotsValues,'enable','off');
        set(PlotsValues(1),'enable','on');
        set(PlotsValues(9),'enable','on');
    end
    %% ��������
    %���¡�ҶƬ������������(��������)����Ļص�����
    if get(PlotsValues(1),'value')
        %�����ͼʱ����
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.input.Xr;
        XCpoDp = pt.input.XCpoDp;
        %XXr��Xr��һ����ֵ���У��������Һ�����չXr����
        XXr = Rhub_oR+(1-Rhub_oR)*(sin((0:60)*pi/(2*60)));
        %XXCpoDpΪͨ����ֵ�õ���XXr����λ�ô����Ҿ��ȣ��ò�ֵ������matlab�Դ�
        XXCpoDp = InterpolateChord(Xr,XCpoDp,XXr);
        %�����ʱ�����µİ�ť��Ϊ1����˵��9Ҳ������
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        Xr_XCpoDp_In_line = plot(coordinate,XXr,XXCpoDp,...
                                 'linewidth',2,...
                                 'color',lc1_1);
        hold(coordinate,'on');
        grid(coordinate,'on');
        box(coordinate,'on');
        Xr_XCpoDp_In_point = scatter(coordinate,Xr,XCpoDp,...
                                     'marker','o',...
                                     'markeredgecolor',lc1_1,...
                                     'markerfacecolor',lc1_2,...
                                     'sizedata',25);
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        if pushednum == 1
            legend(Xr_XCpoDp_In_line,'Cp/Dp');
        else
            legend([Xr_XCpoDp_In_line,Xr_XCpoDp_De_line],'Input:Cp/Dp','Design:Cp/Dp');
            set(PlotsPanels,'title','ҶƬ������������(��������/������ƽ��)');
        end    
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'Cp/Dp','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanels,'title','�����ఴť��ʾ��Ӧͼ��');
    else
        delete(Xr_XCpoDp_In_line);
        delete(Xr_XCpoDp_In_point);
        legend(Xr_XCpoDp_De_line,'Cp/Dp');
        set(PlotsPanels,'title',PlotsStrings{9});
    end
end
%���¡�ҶƬ��ȷֲ�����(��������)����Ļص�����
function BladeThicknessInfcn(hObject,ED)
    global PlotsPanels PlotsStrings
    global PlotsValues
    global pt;
    global Fig_Main coordinate Xr_Xt0oDp_In_line Xr_Xt0oDp_De_line Xr_Xt0oDp_In_point;
    global lc2_1 lc2_2;
    %% ����ͼ��������
    set(PlotsPanels,'title',PlotsStrings{2});
    %% ������Ӧ��ť�Ŀ���״̬
    %���ұ����µİ�ť����
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %�԰�ť2���в�������һ����ťҲû�б�����
        %���4��5������а�ť������
        set(PlotsValues,'enable','on');
        if pt.input.Chord_flag ~= 0
            set(PlotsValues(1),'enable','off');
        end    
        if pt.input.Analyze_flag == 0
            set(PlotsValues(12),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(13),'enable','off');
            set(PlotsValues(14),'enable','off');
        end
        set(PlotsValues(4),'enable','off');
        set(PlotsValues(5),'enable','off');
    else
        %�԰�ť2���в����������ڱ����µİ�ť
        %���2��10������а�ť��������
        set(PlotsValues,'enable','off');
        set(PlotsValues(2),'enable','on');
        set(PlotsValues(10),'enable','on');
    end
    %% ��������
    if get(PlotsValues(2),'value')
        %�����ͼʱ����
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.input.Xr;
        Xt0oDp = pt.input.Xt0oDp;
        %%XXr��Xr��һ����ֵ���У��������Һ�����չXr����
        XXr = Rhub_oR+(1-Rhub_oR)*(sin((0:60)*pi/(2*60)));
        %XXt0oDpΪͨ����ֵ�õ���XXr����λ�ô��ĺ񾶱ȣ��ò�ֵ����Ϊpchip
        XXt0oDp = pchip(Xr,Xt0oDp,XXr);
        %�����ʱ�����µİ�ť��Ϊ1����˵��10Ҳ������
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        Xr_Xt0oDp_In_line = plot(coordinate,XXr,XXt0oDp,...
                                 'linewidth',2,...
                                 'color',lc2_1);
        hold(coordinate,'on');
        grid(coordinate,'on');
        box(coordinate,'on');
        Xr_Xt0oDp_In_point = scatter(coordinate,Xr,Xt0oDp,...
                                     'marker','o',...
                                     'markeredgecolor',lc2_1,...
                                     'markerfacecolor',lc2_2,...
                                     'sizedata',25);
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        if pushednum == 1
            legend(Xr_Xt0oDp_In_line,'t0/Dp');
        else
            legend([Xr_Xt0oDp_In_line,Xr_Xt0oDp_De_line],'Input:t0/Dp','Design:t0/Dp');
            set(PlotsPanels,'title','ҶƬ��ȷֲ�����(��������/������ƽ��)');
        end    
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'t0/Dp','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanels,'title','�����ఴť��ʾ��Ӧͼ��');
    else
        delete(Xr_Xt0oDp_In_line);
        delete(Xr_Xt0oDp_In_point);
        legend(Xr_Xt0oDp_De_line,'t0/Dp');
        set(PlotsPanels,'title',PlotsStrings{10});
    end
end
%���¡������ٶȷֲ����ߡ���Ļص�����
function InflowProfilefcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate Xr_VAI_line Xr_VTI_line Xr_UASTAR_line Xr_UTSTAR_line;
    global lc1_1 lc2_1;
    %% ����ͼ��������
    set(PlotsPanels,'title',PlotsStrings{3});
    %% ������Ӧ��ť�Ŀ���״̬
    %���ұ����µİ�ť����
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %�԰�ť3���в�������һ����ťҲû�б�����
        %���4��5������а�ť������
        set(PlotsValues,'enable','on');
        if pt.input.Chord_flag ~= 0
            set(PlotsValues(1),'enable','off');
        end
        if pt.input.Analyze_flag == 0
            set(PlotsValues(12),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(13),'enable','off');
            set(PlotsValues(14),'enable','off');
        end
        set(PlotsValues(4),'enable','off');
        set(PlotsValues(5),'enable','off');
    else
        %�԰�ť3���в����������ڱ����µİ�ť
        %���3��7������а�ť��������
        set(PlotsValues,'enable','off');
        set(PlotsValues(3),'enable','on');
        set(PlotsValues(7),'enable','on');
    end
    %% ��������
    if get(PlotsValues(3),'value')
        %�����ͼʱ����
        set(0,'CurrentFigure',Fig_Main);
        Xr = pt.input.ri;
        VAI = pt.input.VAI;
        VTI = pt.input.VTI;
        %�����ʱ�����µİ�ť��Ϊ1����˵��7Ҳ������
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        Xr_VAI_line = plot(coordinate,Xr,VAI,...
                           'linewidth',2,...
                           'color',lc1_1);
        hold(coordinate,'on');
        grid(coordinate,'on');
        box(coordinate,'on');
        Xr_VTI_line = plot(coordinate,Xr,VTI,...
                           'linewidth',2,...
                           'color',lc2_1);
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        if pushednum == 1
            legend([Xr_VAI_line,Xr_VTI_line],'VA/Vs','VT/Vs');
        else
            legend([Xr_VAI_line,Xr_VTI_line,Xr_UASTAR_line,Xr_UTSTAR_line],...
                   'VA/Vs','VT/Vs','Ua*/Vs','Ut*/Vs');
            set(PlotsPanels,'title','�ٶȷֲ�����');
        end    
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'V/Vs','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanels,'title','�����ఴť��ʾ��Ӧͼ��');
    else
        delete(Xr_VAI_line);
        delete(Xr_VTI_line);
        legend([Xr_UASTAR_line,Xr_UTSTAR_line],'Ua*/Vs','Ut*/Vs');
        set(PlotsPanels,'title',PlotsStrings{7});
    end    
end

