% �ó�����OpenProp���İ�������GUI������
%% ������
function NewGUI
    %% ǰ����
    clear variables;
    clear global;
    global titlefontsize zhfontsize subcolumnfontsize numfontsize;
    global Fig_Main;
    global SpecificationsValues FlagsValues Property1Values Table1PanelGrid...
           rx_in CDp_x_in Meanline_x_in Thickness_x_in CL_x_in T0oDp_x_in...
           Property2Values Table2PanelGrid ri_in VA_i_in VT_i_in DuctValues...
           ToolsValues FilenameValues
    %% �趨��������ȱʡֵ
    % ���������    
    N_def = 7000;
    VS_def = 5;
    T_def = 0.2;
    Z_def = 3;
    Dp_def = 0.2;
    Hub_flag_def = 1;
    Dh_def = 0.02;
    Duct_flag_def = 1;
    Dd_def = 0.2;
    Cd_def = 0.2;
    % ҶƬ�������
    Nx_def = 10;
    ChordMethod_def = 'CLmax';
    Meanline_flag_def = 1;
    Thickness_flag_def = 1;
    rx_def = [0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;0.95;1];
    CDp_x_def = ones(Nx_def,1)*0.008;
    Meanline_x_items = {'NACA a=0.8','NACA a=0.8 (modified)','Parabolic'};
    Thickness_x_items = {'NACA 65A010','NACA 65A010 (modified)',...
                         'Elliptical','Parabolic','NACA 66 (DTRC modified)'};
    CL_x_def = 0.5+(0.2-0.5)*(rx_def-rx_def(1))/(rx_def(end)-rx_def(1));
    T0oDp_x_def = [0.0329;0.0281;0.0239;0.0198;0.0160;0.0125;0.0091;0.0060;0.0045;0];
    % �ⲿ��������
    Ni_def = 10;
    rho_def = 1.29;
    ri_def = [0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;0.95;1];
    VA_i_def = ones(Ni_def,1);
    VT_i_def = zeros(Ni_def,1);
    % ������ز���
    TpoT_def = 1;
    TdoT_def = 1-TpoT_def;
    CDd_def = 0.008;
    % ���������빤��
    Analyze_flag_def = 1;
    Geometry_flag_def = 1;
    Coordinate_flag_def = 1;
    Printing_flag_def = 1;
    Mp_def = 20;
    Np_def = 20;
    Nd_def = 12;
    ITER_def = 50;
    filename_def = 'DefaultProject';
    %% GUIԪ�صĲ���
    % �����С
    titlefontsize = 60;
    zhfontsize = 15;
    subcolumnfontsize = 17;
    numfontsize = 13;
    % �������ʼ�߶ȺͿ��
    Windowht = 750;
    Window = 1530;
    % �趨�ļ�Ŀ¼
    OpenPropDirectory = 'OpenProp_zh';
    %% �������е�GUIԪ��
    % ����������
    Fig_Main = uifigure('position',[5 55 Window Windowht],...
                        'numbertitle','off',...
                        'name','OpenProp');
    % ���ɶ����˵���
    MenuStrings = {'�ļ�' '����' '��ǩҳ'};
    Menu = zeros(1,length(MenuStrings));
    for index = 1 : length(MenuStrings)
        Menu(index) = uimenu('parent',Fig_Main,...
                             'text',MenuStrings{index});
    end    
    % �����ӹ�����
    FileStrings = {'����' '���Ϊ' '����' '����'};
    FileMenu = zeros(1,length(FileStrings));
    for index = 1 : length(FileStrings)
        FileMenu(index) = uimenu('parent',Menu(1),...
                                 'text',FileStrings{index});
    end
    
    SettingStrings = {'�л�����' 'ҹ��ģʽ'};
    SettingMenu = zeros(1,length(SettingStrings));
    for index = 1 : length(SettingStrings)
        SettingMenu(index) = uimenu('parent',Menu(2),...
                                    'text',SettingStrings{index});
    end    
    
    TabPageStrings = {'������ҶƬ���' '�о�������Ӱ��'};
    for index = 1 : length(TabPageStrings)
        TabPageMenu(index) = uimenu('parent',Menu(3),...
                                    'text',TabPageStrings{index},...
                                    'enable','off');
    end    
    % �����������еı�ǩ��uitabgroup
    Tabs = uitabgroup('parent',Fig_Main,...
                      'position',[0 0 Window Windowht]);
    % ����Tabs�ĵ�һ����ǩҳ��������ҶƬ��ơ�
    SingleDesign = uitab('parent',Tabs,...
                         'title','������ҶƬ���');
    % ����Tabs�ĵڶ�����ǩҳ���о�������Ӱ�족
    ParametricStudy = uitab('parent',Tabs,...
                            'title','�о�������Ӱ��');
    %% ��ǩҳ��������ҶƬ��ơ��е�GUIԪ��
    % ����SingleDesign�е�����ռ�
    SingleDesignGrid = uigridlayout(SingleDesign,[3 3],...
                                    'columnwidth',{'1x','3x','1.5x'},...
                                    'rowheight',{'4x','1x','1.1x'},...
                                    'padding',[10 10 10 10]);
                                
    % ��������� Propeller Specification
    Specifications = uipanel('parent',SingleDesignGrid,...
                             'title','���������',...
                             'titleposition','centertop',...
                             'fontsize',subcolumnfontsize,...
                             'fontweight','bold',...
                             'scrollable','off');
                         
    % ҶƬ������� Blade Design Values
    BladeDesign = uipanel('parent',SingleDesignGrid,...
                          'title','ҶƬ�������',...
                          'titleposition','centertop',...
                          'fontsize',subcolumnfontsize,...
                          'fontweight','bold',...
                          'scrollable','on');
    BladeDesign.Layout.Row = [1 2];
    
    % �ⲿ�������� Inflow Profile Values
    Inflow = uipanel('parent',SingleDesignGrid,...
                     'title','�ⲿ��������',...
                     'titleposition','centertop',...
                     'fontsize',subcolumnfontsize,...
                     'fontweight','bold',...
                     'scrollable','on');
    Inflow.Layout.Row = [1 2];
    
    % ������ز��� Duct Design Values
    Duct = uipanel('parent',SingleDesignGrid,...
                   'title','������ز���',...
                   'titleposition','centertop',...
                   'fontsize',subcolumnfontsize,...
                   'fontweight','bold',...
                   'scrollable','off');
    Duct.Layout.Row = 2;
    Duct.Layout.Column = 1;
    
    % ���� Title
    Title = uilabel(SingleDesignGrid,...
                    'text','OpenProp',...
                    'horizontalalignment','center',...
                    'fontname','Bauhaus 93',...
                    'fontsize',titlefontsize,...
                    'fontweight','normal');
     
    % ���������빤�� Other Desgin Values & Tools
    Tools = uipanel('parent',SingleDesignGrid,...
                    'title','���������빤��',...
                    'titleposition','centertop',...
                    'fontsize',subcolumnfontsize,...
                    'fontweight','bold',...
                    'scrollable','off');
    Tools.Layout.Column = [2 3];
    
    %% ��ǩҳ���о�������Ӱ�족�е�GUIԪ��
    % ����ParametricStudy�е�����
    ParametricStudyGrid = uigridlayout(ParametricStudy,[1 1]);
    Announcement = uilabel(ParametricStudyGrid,...
                           'text','�ù�����δ��ɣ������ڴ�',...
                           'horizontalalignment','center',...
                           'verticalalignment','center',...
                           'fontsize',zhfontsize);
                           
    %% ����������������е�GUIԪ��
    % ����Specifications�е�����
    SpecificationsGrid = uigridlayout(Specifications,[4 1],...
                                      'rowheight',{'3x','2x','2x','3x'},...
                                      'columnspacing',5,...
                                      'rowspacing',5,...
                                      'padding',[5 5 5 5],...
                                      'scrollable','off');   
                                  
    % GUIԪ��
    SpecificationsStrings = {'���ת�٣�','�н��ٶȣ�','��������','ҶƬ������',...
                             'ҶƬֱ����','���ֱ����','����ֱ����','�����ҳ���'};
    SpecificationsTips = {'N','VS','T','Z','Dp','Dh','Dd','Cd'};
    SpecificationsValues_def = {N_def,VS_def,T_def,Z_def,Dp_def,Dh_def,Dd_def,Cd_def};
    
    % װ�����
    DevicePanel = uipanel(SpecificationsGrid);
    DevicePanelGrid = uigridlayout(DevicePanel,[3 2],...
                                   'columnwidth',{'2x','3.5x'},...
                                   'columnspacing',0,...
                                   'rowspacing',10,...
                                   'padding',[10 10 10 10]);    
    for index = 1 : 3
        SpecificationsTexts(index) = uilabel('parent',DevicePanelGrid,...
                                             'text',SpecificationsStrings{index},...
                                             'horizontalalignment','left',...
                                             'verticalalignment','center',...
                                             'tooltip',SpecificationsTips{index},...
                                             'fontsize',zhfontsize);

        SpecificationsValues(index) = uieditfield(DevicePanelGrid,'numeric',...
                                                  'fontsize',numfontsize,...
                                                  'horizontalalignment','center',...
                                                  'value',SpecificationsValues_def{index},...
                                                  'limit',[0 Inf],...
                                                  'lowerlimitinclusive','off');
    end 
    
    % ҶƬ���
    BladePanel = uipanel(SpecificationsGrid);
    BladePanelGrid = uigridlayout(BladePanel,[2 2],...
                                  'columnwidth',{'2x','3.5x'},...
                                  'columnspacing',0,...
                                  'rowspacing',10,...
                                  'padding',[10 10 10 10]);  
    for index = 4 : 5
        SpecificationsTexts(index) = uilabel('parent',BladePanelGrid,...
                                             'text',SpecificationsStrings{index},...
                                             'horizontalalignment','left',...
                                             'verticalalignment','center',...
                                             'tooltip',SpecificationsTips{index},...
                                             'fontsize',zhfontsize);

        SpecificationsValues(index) = uieditfield(BladePanelGrid,'numeric',...
                                                  'fontsize',numfontsize,...
                                                  'horizontalalignment','center',...
                                                  'value',SpecificationsValues_def{index},...
                                                  'limit',[0 Inf],...
                                                  'lowerlimitinclusive','off');
    end    
    
    % ������
    HubPanel = uipanel(SpecificationsGrid);
    HubPanelGrid = uigridlayout(HubPanel,[2 2],...
                                'columnwidth',{'2x','3.5x'},...
                                'columnspacing',0,...
                                'rowspacing',10,...
                                'padding',[10 10 10 10]); 
                            
    HubPanelUpSubGrid = uigridlayout(HubPanelGrid,[1 1],...
                                     'padding',[0 0 0 0]);
    HubPanelUpSubGrid.Layout.Column = [1 2]; 
    FlagsValues(1) = uicheckbox('parent',HubPanelUpSubGrid,...
                                'text','������Ƿ���뽰챣�',...
                                'value',Hub_flag_def,...
                                'tooltip','Hub_flag',...
                                'fontsize',zhfontsize);
    
    SpecificationsTexts(6) = uilabel('parent',HubPanelGrid,...
                                     'text',SpecificationsStrings{6},...
                                     'horizontalalignment','left',...
                                     'verticalalignment','center',...
                                     'tooltip',SpecificationsTips{6},...
                                     'fontsize',zhfontsize);
    SpecificationsValues(6) = uieditfield(HubPanelGrid,'numeric',...
                                          'fontsize',numfontsize,...
                                          'horizontalalignment','center',...
                                          'value',SpecificationsValues_def{6},...
                                          'limit',[0 Inf],...
                                          'lowerlimitinclusive','on');
                                      
    % �������
    DuctPanel = uipanel(SpecificationsGrid);
    DuctPanelGrid = uigridlayout(DuctPanel,[3 2],...
                                 'columnwidth',{'2x','3.5x'},...
                                 'columnspacing',0,...
                                 'rowspacing',10,...
                                 'padding',[10 10 10 10]);  
                             
    DuctPanelUpSubGrid = uigridlayout(DuctPanelGrid,[1 1],...
                                     'padding',[0 0 0 0]);
    DuctPanelUpSubGrid.Layout.Column = [1 2]; 
    FlagsValues(2) = uicheckbox('parent',DuctPanelUpSubGrid,...
                                'text','������Ƿ���뺭����',...
                                'value',Duct_flag_def,...
                                'tooltip','Duct_flag',...
                                'fontsize',zhfontsize);
    
    for index = 7 : 8
        SpecificationsTexts(index) = uilabel('parent',DuctPanelGrid,...
                                             'text',SpecificationsStrings{index},...
                                             'horizontalalignment','left',...
                                             'verticalalignment','center',...
                                             'tooltip',SpecificationsTips{index},...
                                             'fontsize',zhfontsize);

        SpecificationsValues(index) = uieditfield(DuctPanelGrid,'numeric',...
                                                  'fontsize',numfontsize,...
                                                  'horizontalalignment','center',...
                                                  'value',SpecificationsValues_def{index},...
                                                  'limit',[0 Inf],...
                                                  'lowerlimitinclusive','on');
    end    
    
    % ����Specifications�༭���е���ʾ��ʽ�͵�λ
    set(SpecificationsValues(1),'valuedisplayformat','%.0f RPM');
    set(SpecificationsValues(2),'valuedisplayformat','%.2f m/s');
    set(SpecificationsValues(3),'valuedisplayformat','%.3f N');
    set(SpecificationsValues(4),'valuedisplayformat','%.0f');
    set(SpecificationsValues(5),'valuedisplayformat','%.3f m');
    set(SpecificationsValues(6),'valuedisplayformat','%.3f m');
    set(SpecificationsValues(7),'valuedisplayformat','%.3f m');
    set(SpecificationsValues(8),'valuedisplayformat','%.3f m');
    
    %% ������ҶƬ����������е�GUIԪ��
    % ����BladeDesign�е�����
    BladeDesignGrid = uigridlayout(BladeDesign,[2 1],...
                                   'rowheight',{'2x','21x'},...
                                   'columnspacing',5,...
                                   'rowspacing',5,...
                                   'padding',[5 5 5 5],...
                                   'scrollable','off');
                               
    % ����1���
    Property1Panel = uipanel(BladeDesignGrid);
    Property1PanelGrid = uigridlayout(Property1Panel,[1 4],...
                                     'columnwidth',{'1x','1x','0.8x','0.8x'},...
                                     'columnspacing',10,...
                                     'rowspacing',10,...
                                     'padding',[10 10 10 10]);  
    Property1PanelSubGrid(1) = uigridlayout(Property1PanelGrid,[1 2],...
                                            'columnwidth',{'2x','3.5x'},...
                                            'columnspacing',0,...
                                            'rowspacing',0,...
                                            'padding',[0 0 0 0]);  
    Property1PanelSubGrid(2) = uigridlayout(Property1PanelGrid,[1 2],...
                                            'columnwidth',{'2x','3.5x'},...
                                            'columnspacing',0,...
                                            'rowspacing',0,...
                                            'padding',[0 0 0 0]);  
    Property1Texts(1) = uilabel('parent',Property1PanelSubGrid(1),...
                                'text','����������',...
                                'horizontalalignment','left',...
                                'verticalalignment','center',...
                                'tooltip','Nx�������������������',...
                                'fontsize',zhfontsize);
    Property1Values(1) = uieditfield(Property1PanelSubGrid(1),'numeric',...
                                     'fontsize',numfontsize,...
                                     'horizontalalignment','center',...
                                     'value',Nx_def,...
                                     'limit',[0 Inf],...
                                     'lowerlimitinclusive','off');
    Property1Texts(2) = uilabel('parent',Property1PanelSubGrid(2),...
                                'text','��Ʒ�����',...
                                'horizontalalignment','left',...
                                'verticalalignment','center',...
                                'tooltip','ChordMethod',...
                                'fontsize',zhfontsize);
    Property1Values(2) = uidropdown('parent',Property1PanelSubGrid(2),...
                                    'value',ChordMethod_def,...
                                    'items',{'CLmax','ConeyPLL',...
                                             'FAST2011dCTP','FAST2011dVAC',...
                                             'Brizzolara2007'},...
                                    'fontsize',numfontsize);
    FlagsValues(3) = uicheckbox('parent',Property1PanelGrid,...
                                'text','͹Ե���Ƿ�ͬһ��',...
                                'value',Meanline_flag_def,...
                                'tooltip','Meanline_flag',...
                                'fontsize',zhfontsize);
    FlagsValues(4) = uicheckbox('parent',Property1PanelGrid,...
                                'text','Ҷ���Ƿ�ͬһ��',...
                                'value',Thickness_flag_def,...
                                'tooltip','Thickness_flag',...
                                'fontsize',zhfontsize);
                                
    % ���1���
    Table1Panel = uipanel(BladeDesignGrid);
    Table1PanelGrid = uigridlayout(Table1Panel,[Nx_def+1 6],...
                                   'columnspacing',2,...
                                   'rowspacing',2,...
                                   'padding',[10 10 10 0]);
    Table1Strings = {'����λ��','����ϵ��','͹Ե��','Ҷ��','����ϵ��','�񾶱�'};
    Table1Tips = {'rx','CDp_x','Meanline','Thickness',...
                  'CL_x','T0oDp_x��Ҷ����ҶƬֱ��֮��'};
    for index = 1 : length(Table1Strings)
        Col_Label(index) = uilabel('parent',Table1PanelGrid,...
                                   'text',Table1Strings{index},...
                                   'horizontalalignment','center',...
                                   'verticalalignment','center',...
                                   'tooltip',Table1Tips{index},...
                                   'fontsize',zhfontsize);
    end
    for index = 1 : Nx_def
        rx_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                   'fontsize',numfontsize,...
                                   'horizontalalignment','center',...
                                   'value',rx_def(index),...
                                   'limit',[0 1]);
                               
        CDp_x_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                      'fontsize',numfontsize,...
                                      'horizontalalignment','center',...
                                      'value',CDp_x_def(index),...
                                      'limit',[0 Inf]);
                                  
        Meanline_x_in(index) = uidropdown('parent',Table1PanelGrid,...
                                          'items',Meanline_x_items,...
                                          'fontsize',numfontsize,...
                                          'enable','off');
                                      
        Thickness_x_in(index) = uidropdown('parent',Table1PanelGrid,...
                                           'items',Thickness_x_items,...
                                           'fontsize',numfontsize,...
                                           'enable','off');
                                       
        CL_x_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                     'fontsize',numfontsize,...
                                     'horizontalalignment','center',...
                                     'value',CL_x_def(index),...
                                     'limit',[0 Inf],...
                                     'editable','on',...
                                     'enable','on');
                                 
        T0oDp_x_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                       'fontsize',numfontsize,...
                                       'horizontalalignment','center',...
                                       'value',T0oDp_x_def(index),...
                                       'limit',[0 Inf],...
                                       'editable','on',...
                                       'enable','on');
    end   
    % ȱʡ״̬��Meanline��Thickness������ֻ�е�һ�п���
    set(Meanline_x_in(1),'enable','on');
    set(Thickness_x_in(1),'enable','on');
    %% �������ⲿ�����������е�GUIԪ��
    % ����Inflow�е�����
    InflowGrid = uigridlayout(Inflow,[2 1],...
                              'rowheight',{'2x','21x'},...
                              'columnspacing',5,...
                              'rowspacing',5,...
                              'padding',[5 5 5 5],...
                              'scrollable','off');
                          
    % ����2���
    Property2Panel = uipanel(InflowGrid);
    Property2PanelGrid = uigridlayout(Property2Panel,[1 2],...
                                     'columnwidth',{'1x','1x'},...
                                     'columnspacing',10,...
                                     'rowspacing',10,...
                                     'padding',[10 10 10 10]);  
    Property2PanelSubGrid(1) = uigridlayout(Property2PanelGrid,[1 2],...
                                            'columnwidth',{'2x','3x'},...
                                            'columnspacing',0,...
                                            'rowspacing',0,...
                                            'padding',[0 0 0 0]);  
    Property2PanelSubGrid(2) = uigridlayout(Property2PanelGrid,[1 2],...
                                            'columnwidth',{'2x','3x'},...
                                            'columnspacing',0,...
                                            'rowspacing',0,...
                                            'padding',[0 0 0 0]);  
    Property2Texts(1) = uilabel('parent',Property2PanelSubGrid(1),...
                                'text','����������',...
                                'horizontalalignment','left',...
                                'verticalalignment','center',...
                                'tooltip','Ni�������������������',...
                                'fontsize',zhfontsize);
    Property2Values(1) = uieditfield(Property2PanelSubGrid(1),'numeric',...
                                     'fontsize',numfontsize,...
                                     'horizontalalignment','center',...
                                     'value',Ni_def,...
                                     'limit',[0 Inf],...
                                     'lowerlimitinclusive','off');
    Property2Texts(2) = uilabel('parent',Property2PanelSubGrid(2),...
                                'text','�����ܶȣ�',...
                                'horizontalalignment','left',...
                                'verticalalignment','center',...
                                'tooltip','rho',...
                                'fontsize',zhfontsize);
    Property2Values(2) = uieditfield(Property2PanelSubGrid(2),'numeric',...
                                     'fontsize',numfontsize,...
                                     'horizontalalignment','center',...
                                     'value',rho_def,...
                                     'limit',[0 Inf],...
                                     'lowerlimitinclusive','off');
    set(Property2Values(2),'valuedisplayformat','%.2f kg/m^3');
                                
    % ���2���
    Table2Panel = uipanel(InflowGrid);
    Table2PanelGrid = uigridlayout(Table2Panel,[Ni_def+1 3],...
                                   'columnspacing',2,...
                                   'rowspacing',2,...
                                   'padding',[10 10 10 0]);
    Table2Strings = {'����λ��','���������ٶ�','���������ٶ�'};
    Table2Tips = {'ri','VA_i�����������ٶ����н��ٶȵı�ֵ',...
                  'VT_i�����������ٶ����н��ٶȵı�ֵ'};
    for index = 1 : length(Table2Strings)
        Col_Label(index) = uilabel('parent',Table2PanelGrid,...
                                   'text',Table2Strings{index},...
                                   'horizontalalignment','center',...
                                   'verticalalignment','center',...
                                   'tooltip',Table2Tips{index},...
                                   'fontsize',zhfontsize);
    end
    for index = 1 : Ni_def
        ri_in(index) = uieditfield(Table2PanelGrid,'numeric',...
                                   'fontsize',numfontsize,...
                                   'horizontalalignment','center',...
                                   'value',ri_def(index),...
                                   'limit',[0 1]);
                               
        VA_i_in(index) = uieditfield(Table2PanelGrid,'numeric',...
                                     'fontsize',numfontsize,...
                                     'horizontalalignment','center',...
                                     'value',VA_i_def(index),...
                                     'limit',[0 Inf]);
                                
        VT_i_in(index) = uieditfield(Table2PanelGrid,'numeric',...
                                     'fontsize',numfontsize,...
                                     'horizontalalignment','center',...
                                     'value',VT_i_def(index),...
                                     'limit',[0 Inf]);
    end    
    %% ������������ز������е�GUIԪ��
    % ����Duct�е�����
    DuctGrid = uigridlayout(Duct,[2 2],...
                            'columnwidth',{'2x','3.5x'},...
                            'rowheight',{'1x','1x'},...
                            'columnspacing',10,...
                            'rowspacing',10,...
                            'padding',[15 10 15 10],...
                            'scrollable','off');
    DuctStrings = {'�������أ�','����ϵ����'};
    DuctTips = {'TdoT���������ṩ����ռ�������ı�ֵ','CDd'};
    DuctValues_def = {TdoT_def,CDd_def};
    for index = 1 : 2
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
                        
    %% ������������Ʋ����빤�ߡ��е�GUIԪ��
    % ����Tools�е�����
    ToolsGrid = uigridlayout(Tools,[1,3],...
                             'columnwidth',{'1.2x','2x','1x'},...
                             'columnspacing',5,...
                             'rowspacing',5,...
                             'padding',[5 5 5 5],...
                             'scrollable','off');
    
    % ��ѡ�����
    FlagsPanel = uipanel(ToolsGrid);
    FlagsPanelGrid = uigridlayout(FlagsPanel,[2,2]);
    FlagsStrings = {'�Ƿ�����������ߣ�','�Ƿ������ʾģ�ͣ�',...
                    '�Ƿ���������ꣿ','�Ƿ����stl�ļ���'};
    FlagsTips = {'Analyze_flag','Geometry_flag',...
                 'Coordinate_flag','Printing_flag'};
    FlagsValues_def = {Analyze_flag_def,Geometry_flag_def,...
                       Coordinate_flag_def,Printing_flag_def};
    for index = 1 : 4
        FlagsValues(index+4) = uicheckbox('parent',FlagsPanelGrid,...
                                          'text',FlagsStrings{index},...
                                          'value',FlagsValues_def{index},...
                                          'tooltip',FlagsTips{index},...
                                          'fontsize',zhfontsize);
    end   
    % ���Ƶ����
    ControlPointPanel = uipanel(ToolsGrid);
    ControlPointPanelGrid = uigridlayout(ControlPointPanel,[2 2]);
    ControlPointPanelSubGrid(1) = uigridlayout(ControlPointPanelGrid,[1 2],...
                                               'columnwidth',{'2x','2.5x'},...
                                               'columnspacing',0,...
                                               'padding',[0 0 0 0]);
    ControlPointPanelSubGrid(2) = uigridlayout(ControlPointPanelGrid,[1 2],...
                                               'columnwidth',{'2x','2.5x'},...
                                               'columnspacing',0,...
                                               'padding',[0 0 0 0]);
    ControlPointPanelSubGrid(3) = uigridlayout(ControlPointPanelGrid,[1 2],...
                                               'columnwidth',{'2x','2.5x'},...
                                               'columnspacing',0,...
                                               'padding',[0 0 0 0]);
    ControlPointPanelSubGrid(4) = uigridlayout(ControlPointPanelGrid,[1 2],...
                                               'columnwidth',{'2x','2.5x'},...
                                               'columnspacing',0,...
                                               'padding',[0 0 0 0]);
    ToolsTexts(1) = uilabel('parent',ControlPointPanelSubGrid(1),...
                            'text','ҶƬ���Ƶ�������',...
                            'horizontalalignment','left',...
                            'verticalalignment','center',...
                            'tooltip','Mp',...
                            'fontsize',zhfontsize);
    ToolsValues(1) = uieditfield(ControlPointPanelSubGrid(1),'numeric',...
                                 'fontsize',numfontsize,...
                                 'horizontalalignment','center',...
                                 'value',Mp_def,...
                                 'limit',[0 Inf],...
                                 'lowerlimitinclusive','off');
    ToolsTexts(2) = uilabel('parent',ControlPointPanelSubGrid(2),...
                            'text','ҶƬ�����������',...
                            'horizontalalignment','left',...
                            'verticalalignment','center',...
                            'tooltip','Np',...
                            'fontsize',zhfontsize);
    ToolsValues(2) = uieditfield(ControlPointPanelSubGrid(2),'numeric',...
                                 'fontsize',numfontsize,...
                                 'horizontalalignment','center',...
                                 'value',Np_def,...
                                 'limit',[0 Inf],...
                                 'lowerlimitinclusive','off');
    ToolsTexts(3) = uilabel('parent',ControlPointPanelSubGrid(3),...
                            'text','����������������',...
                            'horizontalalignment','left',...
                            'verticalalignment','center',...
                            'tooltip','Nd',...
                            'fontsize',zhfontsize);
    ToolsValues(3) = uieditfield(ControlPointPanelSubGrid(3),'numeric',...
                                 'fontsize',numfontsize,...
                                 'horizontalalignment','center',...
                                 'value',Nd_def,...
                                 'limit',[0 Inf],...
                                 'lowerlimitinclusive','off');
    ToolsTexts(4) = uilabel('parent',ControlPointPanelSubGrid(4),...
                            'text','������������',...
                            'horizontalalignment','left',...
                            'verticalalignment','center',...
                            'tooltip','ITER',...
                            'fontsize',zhfontsize);
    ToolsValues(4) = uieditfield(ControlPointPanelSubGrid(4),'numeric',...
                                 'fontsize',numfontsize,...
                                 'horizontalalignment','center',...
                                 'value',ITER_def,...
                                 'limit',[0 Inf],...
                                 'lowerlimitinclusive','off');
    
    % �ļ����
    FilePanel = uipanel(ToolsGrid);
    FilePanelGrid = uigridlayout(FilePanel,[2 4]);
    FilenameTexts = uilabel('parent',FilePanelGrid,...
                            'text','��Ŀ����',...
                            'horizontalalignment','left',...
                            'verticalalignment','center',...
                            'fontsize',zhfontsize);
    FilenameValues = uieditfield(FilePanelGrid,...
                                 'fontsize',numfontsize,...
                                 'horizontalalignment','center',...
                                 'value',filename_def);
    FilenameValues.Layout.Column = [2 4];
    LoadButton = uibutton('parent',FilePanelGrid,...
                          'text','����',...
                          'icon','load.png');
    SaveButton = uibutton('parent',FilePanelGrid,...
                          'text','����',...
                          'icon','save.png',...
                          'buttonpushedfcn',@Save);
    SaveCopyAsButton = uibutton('parent',FilePanelGrid,...
                                'text','���',...
                                'icon','savecopy.png');
    RunButton = uibutton('parent',FilePanelGrid,...
                         'text','����',...
                         'icon','run.png');
                     
    %% ���ûص�����
    % ����Specifications�༭��Ļص�����
    set(SpecificationsValues,'valuechangedfcn',{@ChangeSpecifications,...
                                                Dp_def,Dh_def,Dd_def});
    
    % ����Specifications��ѡ��Ļص�����
    set(FlagsValues(1),'valuechangedfcn',{@ChangeHub,Dh_def});
    set(FlagsValues(2),'valuechangedfcn',{@ChangeDuct,Dd_def,Cd_def,...
                                          TdoT_def,CDd_def});
    
    % ����BladeDesign��������ƵĻص�����
    set(Property1Values(1),'valuechangedfcn',{@ChangeNx,Meanline_x_items,...
                                              Thickness_x_items});
    
    % ����BladeDesign������ݿ��ƵĻص�����
    set(Property1Values(2),'valuechangedfcn',@ChangeMethod);
    set(FlagsValues(3),'valuechangedfcn',@ConstantMeanline);
    set(FlagsValues(4),'valuechangedfcn',@ConstantThickness);
    set(Meanline_x_in(1),'valuechangedfcn',@Change1stMeanline);
    set(Thickness_x_in(1),'valuechangedfcn',@Change1stThickness);
    
    % ����Inflow��������ƵĻص�����
    set(Property2Values(1),'valuechangedfcn',@ChangeNi);
    
    % ����Tools�༭��Ļص�����
    set(ToolsValues,'valuechangedfcn',@ChangeTools)
    
    % ����Tools��ť�Ļص�����
    set(LoadButton,'buttonpushedfcn',{@Load,OpenPropDirectory});
    set(SaveButton,'buttonpushedfcn',{@Save,OpenPropDirectory});
    set(SaveCopyAsButton,'buttonpushedfcn',@SaveCopyAs);
    set(RunButton,'buttonpushedfcn',@Execute);
    
end

%% ������������������漰���Ļص�����
% �༭��Ļص�����
function ChangeSpecifications(~,~,Dp_def,Dh_def,Dd_def)
    global Fig_Main SpecificationsValues;
    % ��������ֵ�������������������
    Z = round(get(SpecificationsValues(4),'value'));
    set(SpecificationsValues(4),'value',Z);
    Dp = get(SpecificationsValues(5),'value');
    Dh = get(SpecificationsValues(6),'value');
    Dd = get(SpecificationsValues(7),'value');
    % ��齰�ֱ����������ֱ��������ֵ����ǰ�ߴ��ں��ߣ��򵯳����ѿ�
    if Dh >= Dp
        message1 = sprintf('��ǰ���õĽ��ֱ������������ֱ�� \n ���������ã�');
        uialert(Fig_Main,message1,'�������ô���',...
                'closefcn',{@ResetDiameters1,SpecificationsValues,...
                            Dp_def,Dh_def});
    end
    % ��麭��ֱ����������ֱ��������ֵ����ǰ��С�ں��ߣ��򵯳����ѿ�
    if Dd < Dp
        message2 = sprintf('��ǰ���õĺ���ֱ��С��������ֱ�� \n ���������ã�');
        uialert(Fig_Main,message2,'�������ô���',...
                'closefcn',{@ResetDiameters2,SpecificationsValues,...
                            Dp_def,Dd_def});
    end
end
% ����ֱ�����ô������õĺ���
function ResetDiameters1(~,~,Dp_def,Dh_def)
    global SpecificationsValues;
    set(SpecificationsValues(5),'value',Dp_def);
    set(SpecificationsValues(6),'value',Dh_def);
end
function ResetDiameters2(~,~,Dp_def,Dd_def)
    global SpecificationsValues;
    set(SpecificationsValues(5),'value',Dp_def);
    set(SpecificationsValues(7),'value',Dd_def);
end

% ��ѡ��Ļص�����
function ChangeHub(hObject,~,Dh_def)
    global SpecificationsValues;
    if get(hObject,'value')
        % ��ѡ�����뽰챡�������༭���ֱ��
        set(SpecificationsValues(6),'enable','on');
        set(SpecificationsValues(6),'editable','on');
        set(SpecificationsValues(6),'value',Dh_def);
    else
        set(SpecificationsValues(6),'enable','off');
        set(SpecificationsValues(6),'editable','off');
        set(SpecificationsValues(6),'value',0);
    end
end

function ChangeDuct(hObject,~,Dd_def,Cd_def,TdoT_def,CDd_def)
    global SpecificationsValues DuctValues;
    if get(hObject,'value')
        % ��ѡ�����뺭����������༭������ؿؼ�
        set(SpecificationsValues(7),'enable','on');
        set(SpecificationsValues(7),'editable','on');
        set(SpecificationsValues(7),'value',Dd_def);
        set(SpecificationsValues(8),'enable','on');
        set(SpecificationsValues(8),'editable','on');
        set(SpecificationsValues(8),'value',Cd_def);
        set(DuctValues(1),'enable','on');
        set(DuctValues(1),'editable','on');
        set(DuctValues(1),'value',TdoT_def);
        set(DuctValues(2),'enable','on');
        set(DuctValues(2),'editable','on');
        set(DuctValues(2),'value',CDd_def);
    else
        set(SpecificationsValues(7),'enable','off');
        set(SpecificationsValues(7),'editable','off');
        set(SpecificationsValues(7),'value',0);
        set(SpecificationsValues(8),'enable','off');
        set(SpecificationsValues(8),'editable','off');
        set(SpecificationsValues(8),'value',0);
        set(DuctValues(1),'enable','off');
        set(DuctValues(1),'editable','off');
        set(DuctValues(1),'value',0);
        set(DuctValues(2),'enable','off');
        set(DuctValues(2),'editable','off');
        set(DuctValues(2),'value',0);
    end
end

%% ������ҶƬ������������漰���Ļص�����
% ���������
function ChangeNx(hObject,ED,Meanline_x_items,Thickness_x_items)
    global numfontsize;
    global Table1PanelGrid rx_in CDp_x_in Meanline_x_in Thickness_x_in...
           CL_x_in T0oDp_x_in Property1Values FlagsValues;
    % ��������ֵ�������������������  
    set(hObject,'value',round(get(hObject,'value')));
    Nx_new = get(hObject,'value');
    % ED�д����˺ͱ���ֵ������ص����ݣ������ǰ��ֵ
    Nx_previous = ED.PreviousValue;
    
    % ��GUI��ԭ��������ֵ����
    rx_previous = ones(Nx_previous,1);
    CDp_x_previous = ones(Nx_previous,1);
    CL_x_previous = ones(Nx_previous,1);
    T0oDp_x_previous = ones(Nx_previous,1);
    for index = 1 : Nx_previous
        rx_previous(index) = get(rx_in(index),'value');
        CDp_x_previous(index) = get(CDp_x_in(index),'value');
        CL_x_previous(index) = get(CL_x_in(index),'value');
        T0oDp_x_previous(index) = get(T0oDp_x_in(index),'value'); 
    end    
    % rx��ȱʡֵ��ȷ��˼·�����������˵�����һ�ε���ֵ�㣬����λ��ƽ������
    rx_new = zeros(Nx_new,1);
    for index = 1 : Nx_new-2
        rx_new(index) = 0.2+(index-1)*(1-0.2)/(Nx_new-1-1);
    end
    rx_new(Nx_new-1) = (1+rx_new(Nx_new-2))/2;
    rx_new(Nx_new) = 1;
    % CDp_x��CL_x��T0oDp_x��ȱʡֵ��ȷ��˼·��pchip��ֵ
    CDp_x_new = pchip(rx_previous,CDp_x_previous,rx_new);
    CL_x_new = pchip(rx_previous,CL_x_previous,rx_new);
    T0oDp_x_new = pchip(rx_previous,T0oDp_x_previous,rx_new);
    
    RowHeightCell = cell(1,Nx_new+1);
    for index = 1 : Nx_new+1
        RowHeightCell{index} = '1x';
    end    
    if Nx_new > Nx_previous
        % ����޸ĺ���������
        set(Table1PanelGrid,'rowheight',RowHeightCell);
        % �޸�ԭ�б������ֵ
        for index = 1 : Nx_previous
            set(rx_in(index),'value',rx_new(index));
            set(CDp_x_in(index),'value',CDp_x_new(index));
            set(CL_x_in(index),'value',CL_x_new(index));
            set(T0oDp_x_in(index),'value',T0oDp_x_new(index));
        end   
        % �����µı��
        for index = Nx_previous+1 : Nx_new
            rx_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                       'fontsize',numfontsize,...
                                       'horizontalalignment','center',...
                                       'value',rx_new(index),...
                                       'limit',[0 1]);
                               
            CDp_x_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                          'fontsize',numfontsize,...
                                          'horizontalalignment','center',...
                                          'value',CDp_x_new(index),...
                                          'limit',[0 Inf]);

            Meanline_x_in(index) = uidropdown('parent',Table1PanelGrid,...
                                              'items',Meanline_x_items,...
                                              'fontsize',numfontsize);

            Thickness_x_in(index) = uidropdown('parent',Table1PanelGrid,...
                                               'items',Thickness_x_items,...
                                               'fontsize',numfontsize);

            CL_x_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                         'fontsize',numfontsize,...
                                         'horizontalalignment','center',...
                                         'value',CL_x_new(index),...
                                         'limit',[0 Inf]);

            T0oDp_x_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                            'fontsize',numfontsize,...
                                            'horizontalalignment','center',...
                                            'value',T0oDp_x_new(index),...
                                            'limit',[0 Inf]);
            
            if strcmp(get(Property1Values(2),'value'),'CLmax')
                % ��Ʒ���Ϊ'CLmax'������༭����ϵ���ͺ񾶱�
                set(CL_x_in(index),'enable','on');
                set(CL_x_in(index),'editable','on');
                set(T0oDp_x_in(index),'enable','on');
                set(T0oDp_x_in(index),'editable','on');
            else
                set(CL_x_in(index),'enable','off');
                set(CL_x_in(index),'editable','off');
                set(T0oDp_x_in(index),'enable','off');
                set(T0oDp_x_in(index),'editable','off');
            end   
            
            if get(FlagsValues(3),'value')
                set(Meanline_x_in(index),'value',...
                    get(Meanline_x_in(1),'value'));
                set(Meanline_x_in(index),'enable','off');
            else    
                set(Meanline_x_in(index),'enable','on');
            end 
            
            if get(FlagsValues(4),'value')
                set(Thickness_x_in(index),'value',...
                    get(Thickness_x_in(1),'value'));
                set(Thickness_x_in(index),'enable','off');
            else
                set(Thickness_x_in(index),'enable','on');
            end    
        end
    elseif Nx_new < Nx_previous
        % ɾ������ı��
        for index = Nx_new+1 : Nx_previous
            delete(rx_in(index));
            delete(CDp_x_in(index));
            delete(Meanline_x_in(index));
            delete(Thickness_x_in(index));
            delete(CL_x_in(index));
            delete(T0oDp_x_in(index));
        end  
        set(Table1PanelGrid,'rowheight',RowHeightCell);
        % �޸�ʣ��������ֵ
        for index = 1 : Nx_new
            set(rx_in(index),'value',rx_new(index));
            set(CDp_x_in(index),'value',CDp_x_new(index));
            set(CL_x_in(index),'value',CL_x_new(index));
            set(T0oDp_x_in(index),'value',T0oDp_x_new(index));
        end   
    else
        return;
    end    
end

% ������ݿ���
function ChangeMethod(hObject,~)
    global Property1Values CL_x_in T0oDp_x_in;
    Nx = get(Property1Values(1),'value');
    method = get(hObject,'value');
    if strcmp(method,'CLmax')
        % ��Ʒ���Ϊ'CLmax'������༭����ϵ���ͺ񾶱�
        for index = 1 : Nx
            set(CL_x_in(index),'enable','on');
            set(CL_x_in(index),'editable','on');
            set(T0oDp_x_in(index),'enable','on');
            set(T0oDp_x_in(index),'editable','on');
        end    
    else
        for index = 1 : Nx
            set(CL_x_in(index),'enable','off');
            set(CL_x_in(index),'editable','off');
            set(T0oDp_x_in(index),'enable','off');
            set(T0oDp_x_in(index),'editable','off');
        end   
    end
end

function ConstantMeanline(hObject,~)
    global Property1Values Meanline_x_in;
    Nx = get(Property1Values(1),'value');
    if get(hObject,'value')
        % ��ѡ��͹Ե��ͬһ��������һ��������������ɱ༭
        for index = 2 : Nx
            set(Meanline_x_in(index),'enable','off');
        end    
    else
        for index = 2 : Nx
            set(Meanline_x_in(index),'enable','on');
        end
    end    
end

function ConstantThickness(hObject,~)
    global Property1Values Thickness_x_in;
    Nx = get(Property1Values(1),'value');
    if get(hObject,'value')
        % ��ѡ��Ҷ��ͬһ��������һ��������������ɱ༭
        for index = 2 : Nx
            set(Thickness_x_in(index),'enable','off');
        end    
    else
        for index = 2 : Nx
            set(Thickness_x_in(index),'enable','on');
        end
    end    
end

function Change1stMeanline(hObject,~)
    global Property1Values FlagsValues Meanline_x_in;
    Nx = get(Property1Values(1),'value');
    meanline = get(hObject,'value');
    if get(FlagsValues(3),'value')
        % ��ѡ��͹Ե��ͬһ������һ�и��ĺ���������һ������
        for index = 2 : Nx
            set(Meanline_x_in(index),'value',meanline);
        end
    else
        return
    end    
end

function Change1stThickness(hObject,~)
    global Property1Values FlagsValues Thickness_x_in;
    Nx = get(Property1Values(1),'value');
    thickness = get(hObject,'value');
    if get(FlagsValues(4),'value')
        % ��ѡ��Ҷ��ͬһ������һ�и��ĺ���������һ������
        for index = 2 : Nx
            set(Thickness_x_in(index),'value',thickness);
        end
    else
        return
    end    
end

%% �������ⲿ�������������漰���Ļص�����
function ChangeNi(hObject,ED)
    global numfontsize;
    global Table2PanelGrid ri_in VA_i_in VT_i_in;
    % ��������ֵ�������������������  
    set(hObject,'value',round(get(hObject,'value')));
    Ni_new = get(hObject,'value');
    % ED�д����˺ͱ���ֵ������ص����ݣ������ǰ��ֵ
    Ni_previous = ED.PreviousValue;
    
    % ��GUI��ԭ��������ֵ����
    ri_previous = ones(Ni_previous,1);
    VA_i_previous = ones(Ni_previous,1);
    VT_i_previous = ones(Ni_previous,1);
    for index = 1 : Ni_previous
        ri_previous(index) = get(ri_in(index),'value');
        VA_i_previous(index) = get(VA_i_in(index),'value');
        VT_i_previous(index) = get(VT_i_in(index),'value');
    end    
    % ri��ȱʡֵ��ȷ��˼·�����������˵�����һ�ε���ֵ�㣬����λ��ƽ������
    ri_new = zeros(Ni_new,1);
    for index = 1 : Ni_new-2
        ri_new(index) = 0.2+(index-1)*(1-0.2)/(Ni_new-1-1);
    end
    ri_new(Ni_new-1) = (1+ri_new(Ni_new-2))/2;
    ri_new(Ni_new) = 1;
    % VA_i��VT_i��ȱʡֵ��ȷ��˼·��pchip��ֵ
    VA_i_new = pchip(ri_previous,VA_i_previous,ri_new);
    VT_i_new = pchip(ri_previous,VT_i_previous,ri_new);
    
    RowHeightCell = cell(1,Ni_new+1);
    for index = 1 : Ni_new+1
        RowHeightCell{index} = '1x';
    end    
    if Ni_new > Ni_previous
        % ����޸ĺ���������
        set(Table2PanelGrid,'rowheight',RowHeightCell);
        % �޸�ԭ�б������ֵ
        for index = 1 : Ni_previous
            set(ri_in(index),'value',ri_new(index));
            set(VA_i_in(index),'value',VA_i_new(index));
            set(VT_i_in(index),'value',VT_i_new(index));
        end   
        % �����µı��
        for index = Ni_previous+1 : Ni_new
            ri_in(index) = uieditfield(Table2PanelGrid,'numeric',...
                                       'fontsize',numfontsize,...
                                       'horizontalalignment','center',...
                                       'value',ri_new(index),...
                                       'limit',[0 1]);
                               
            VA_i_in(index) = uieditfield(Table2PanelGrid,'numeric',...
                                         'fontsize',numfontsize,...
                                         'horizontalalignment','center',...
                                         'value',VA_i_new(index),...
                                         'limit',[0 Inf]);

            VT_i_in(index) = uieditfield(Table2PanelGrid,'numeric',...
                                         'fontsize',numfontsize,...
                                         'horizontalalignment','center',...
                                         'value',VT_i_new(index),...
                                         'limit',[0 Inf]);
        end
    elseif Ni_new < Ni_previous
        % ɾ������ı��
        for index = Ni_new+1 : Ni_previous
            delete(ri_in(index));
            delete(VA_i_in(index));
            delete(VT_i_in(index));
        end  
        set(Table2PanelGrid,'rowheight',RowHeightCell);
        % �޸�ʣ��������ֵ
        for index = 1 : Ni_new
            set(ri_in(index),'value',ri_new(index));
            set(VA_i_in(index),'value',VA_i_new(index));
            set(VT_i_in(index),'value',VT_i_new(index));
        end   
    else
        return;
    end    
end

%% ���������������빤�ߡ����漰���Ļص�����
% ���Ƶ�����ֵ����
function ChangeTools(hObject,~)
    % ��������ֵ�������������������
    input = round(get(hObject,'value'));
    set(hObject,'value',input);
end

% �ļ���ز����Ļص�����
function Load(~,~,OpenPropDirectory)
    global pt numfontsize;
    global SpecificationsValues FlagsValues Property1Values Table1PanelGrid...
           rx_in CDp_x_in Meanline_x_in Thickness_x_in CL_x_in T0oDp_x_in...
           Property2Values Table2PanelGrid ri_in VA_i_in VT_i_in DuctValues...
           ToolsValues FilenameValues;
    rest = pwd;
    LoadDirectory(rest,OpenPropDirectory);
    % uiload�������ڴ��ļ�ѡ��Ŀ¼
    uiload;
    
    % ��ȡ���µ��밴ťǰGUI�����е�Nx��Ni
    Nx_new = pt.input.Nx;
    Nx_previous = get(Property1Values(1),'value');
    Ni_new = pt.input.Ni;
    Ni_previous = get(Property2Values(1),'value');
    
    RowHeightCell = cell(1,Nx_new+1);
    for index = 1 : Nx_new+1
        RowHeightCell{index} = '1x';
    end    
    
    % �����ݵ��롰���������
    set(SpecificationsValues(1),'value',pt.input.N);
    set(SpecificationsValues(2),'value',pt.input.VS);
    set(SpecificationsValues(3),'value',pt.input.T);
    set(SpecificationsValues(4),'value',pt.input.Z);
    set(SpecificationsValues(5),'value',pt.input.Dp);
    set(FlagsValues(1),'value',pt.input.Hub_flag);
    set(SpecificationsValues(6),'value',pt.input.Dh);
    set(FlagsValues(2),'value',pt.input.Duct_flag);
    set(SpecificationsValues(7),'value',pt.input.Dd);
    set(SpecificationsValues(8),'value',pt.input.Cd);
    
    % �����ݵ��롰ҶƬ���������
    set(Property1Values(1),'value',pt.input.Nx);
    set(Property1Values(2),'value',pt.input.ChordMethod);
    set(FlagsValues(3),'value',pt.input.Meanline_flag);
    set(FlagsValues(4),'value',pt.input.Thickness_flag);
    if Nx_new > Nx_previous
        % �������ݵ�������������ԭGUI��������ʾ��
        set(Table1PanelGrid,'rowheight',RowHeightCell);
        % �޸�ԭ�б������ֵ
        for index = 1 : Nx_previous
            set(rx_in(index),'value',pt.input.rx(index));
            set(CDp_x_in(index),'value',pt.input.CDp_x(index));
            set(Meanline_x_in(index),'value',pt.input.Meanline_x{index});
            set(Thickness_x_in(index),'value',pt.input.Thickness_x{index});
            set(CL_x_in(index),'value',pt.input.CL_x(index));
            set(T0oDp_x_in(index),'value',pt.input.T0oDp_x(index));
            
            if strcmp(pt.input.ChordMethod,'CLmax')
                % ��Ʒ���Ϊ'CLmax'������༭����ϵ���ͺ񾶱�
                set(CL_x_in(index),'enable','on');
                set(CL_x_in(index),'editable','on');
                set(T0oDp_x_in(index),'enable','on');
                set(T0oDp_x_in(index),'editable','on');
            else
                set(CL_x_in(index),'enable','off');
                set(CL_x_in(index),'editable','off');
                set(T0oDp_x_in(index),'enable','off');
                set(T0oDp_x_in(index),'editable','off');
            end   
            
            if pt.input.Meanline_flag
                set(Meanline_x_in(index),'value',...
                    get(Meanline_x_in(1),'value'));
                set(Meanline_x_in(index),'enable','off');
            else    
                set(Meanline_x_in(index),'enable','on');
            end 
            
            if pt.input.Thickness_flag
                set(Thickness_x_in(index),'value',...
                    get(Thickness_x_in(1),'value'));
                set(Thickness_x_in(index),'enable','off');
            else
                set(Thickness_x_in(index),'enable','on');
            end    
        end   
        % �����µı��
        for index = Nx_previous+1 : Nx_new
            rx_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                       'fontsize',numfontsize,...
                                       'horizontalalignment','center',...
                                       'value',pt.input.rx(index),...
                                       'limit',[0 1]);
                               
            CDp_x_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                          'fontsize',numfontsize,...
                                          'horizontalalignment','center',...
                                          'value',pt.input.CDp_x(index),...
                                          'limit',[0 Inf]);

            Meanline_x_in(index) = uidropdown('parent',Table1PanelGrid,...
                                              'items',Meanline_x_items,...
                                              'value',pt.input.Meanline_x{index},...
                                              'fontsize',numfontsize);

            Thickness_x_in(index) = uidropdown('parent',Table1PanelGrid,...
                                               'items',Thickness_x_items,...
                                               'value',pt.input.Thickness_x{index},...
                                               'fontsize',numfontsize);

            CL_x_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                         'fontsize',numfontsize,...
                                         'horizontalalignment','center',...
                                         'value',pt.input.CL_x(index),...
                                         'limit',[0 Inf]);

            T0oDp_x_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                            'fontsize',numfontsize,...
                                            'horizontalalignment','center',...
                                            'value',pt.input.T0oDp_x(index),...
                                            'limit',[0 Inf]);
            
            if strcmp(pt.input.ChordMethod,'CLmax')
                % ��Ʒ���Ϊ'CLmax'������༭����ϵ���ͺ񾶱�
                set(CL_x_in(index),'enable','on');
                set(CL_x_in(index),'editable','on');
                set(T0oDp_x_in(index),'enable','on');
                set(T0oDp_x_in(index),'editable','on');
            else
                set(CL_x_in(index),'enable','off');
                set(CL_x_in(index),'editable','off');
                set(T0oDp_x_in(index),'enable','off');
                set(T0oDp_x_in(index),'editable','off');
            end   
            
            if pt.input.Meanline_flag
                set(Meanline_x_in(index),'value',...
                    get(Meanline_x_in(1),'value'));
                set(Meanline_x_in(index),'enable','off');
            else    
                set(Meanline_x_in(index),'enable','on');
            end 
            
            if pt.input.Thickness_flag
                set(Thickness_x_in(index),'value',...
                    get(Thickness_x_in(1),'value'));
                set(Thickness_x_in(index),'enable','off');
            else
                set(Thickness_x_in(index),'enable','on');
            end    
        end
    elseif Nx_new < Nx_previous
        % ɾ������ı��
        for index = Nx_new+1 : Nx_previous
            delete(rx_in(index));
            delete(CDp_x_in(index));
            delete(Meanline_x_in(index));
            delete(Thickness_x_in(index));
            delete(CL_x_in(index));
            delete(T0oDp_x_in(index));
        end  
        set(Table1PanelGrid,'rowheight',RowHeightCell);
        % �޸�ʣ��������ֵ
        for index = 1 : Nx_new
            set(rx_in(index),'value',pt.input.rx(index));
            set(CDp_x_in(index),'value',pt.input.CDp_x(index));
            set(Meanline_x_in(index),'value',pt.input.Meanline_x{index});
            set(Thickness_x_in(index),'value',pt.input.Thickness_x{index});
            set(CL_x_in(index),'value',pt.input.CL_x(index));
            set(T0oDp_x_in(index),'value',pt.input.T0oDp_x(index));
            
            if strcmp(pt.input.ChordMethod,'CLmax')
                % ��Ʒ���Ϊ'CLmax'������༭����ϵ���ͺ񾶱�
                set(CL_x_in(index),'enable','on');
                set(CL_x_in(index),'editable','on');
                set(T0oDp_x_in(index),'enable','on');
                set(T0oDp_x_in(index),'editable','on');
            else
                set(CL_x_in(index),'enable','off');
                set(CL_x_in(index),'editable','off');
                set(T0oDp_x_in(index),'enable','off');
                set(T0oDp_x_in(index),'editable','off');
            end   
            
            if pt.input.Meanline_flag
                set(Meanline_x_in(index),'value',...
                    get(Meanline_x_in(1),'value'));
                set(Meanline_x_in(index),'enable','off');
            else    
                set(Meanline_x_in(index),'enable','on');
            end 
            
            if pt.input.Thickness_flag
                set(Thickness_x_in(index),'value',...
                    get(Thickness_x_in(1),'value'));
                set(Thickness_x_in(index),'enable','off');
            else
                set(Thickness_x_in(index),'enable','on');
            end    
        end   
    else    
        return
    end    
    
    % �����ݵ��롰�ⲿ����������
    set(Property2Values(1),'value',pt.input.Ni);
    set(Property2Values(2),'value',pt.input.rho);
    if Ni_new > Ni_previous
        % ����޸ĺ���������
        set(Table2PanelGrid,'rowheight',RowHeightCell);
        % �޸�ԭ�б������ֵ
        for index = 1 : Ni_previous
            set(ri_in(index),'value',pt.input.ri(index));
            set(VA_i_in(index),'value',pt.input.VA_i(index));
            set(VT_i_in(index),'value',pt.input.VT_i(index));
        end   
        % �����µı��
        for index = Ni_previous+1 : Ni_new
            ri_in(index) = uieditfield(Table2PanelGrid,'numeric',...
                                       'fontsize',numfontsize,...
                                       'horizontalalignment','center',...
                                       'value',pt.input.ri(index),...
                                       'limit',[0 1]);
                               
            VA_i_in(index) = uieditfield(Table2PanelGrid,'numeric',...
                                         'fontsize',numfontsize,...
                                         'horizontalalignment','center',...
                                         'value',pt.input.VA_i(index),...
                                         'limit',[0 Inf]);

            VT_i_in(index) = uieditfield(Table2PanelGrid,'numeric',...
                                         'fontsize',numfontsize,...
                                         'horizontalalignment','center',...
                                         'value',pt.input.VT_i(index),...
                                         'limit',[0 Inf]);
        end
    elseif Ni_new < Ni_previous
        % ɾ������ı��
        for index = Ni_new+1 : Ni_previous
            delete(ri_in(index));
            delete(VA_i_in(index));
            delete(VT_i_in(index));
        end  
        set(Table2PanelGrid,'rowheight',RowHeightCell);
        % �޸�ʣ��������ֵ
        for index = 1 : Ni_new
            set(ri_in(index),'value',pt.input.ri(index));
            set(VA_i_in(index),'value',pt.input.VA_i(index));
            set(VT_i_in(index),'value',pt.input.VT_i(index));
        end   
    else
        return;
    end 
    
    % �����ݵ��롰������ز�����
    set(DuctValues(1),'value',pt.input.TdoT);
    set(DuctValues(2),'value',pt.input.CDd);
    
    % �����ݵ��롰���������빤�ߡ�
    set(FlagsValues(5),'value',pt.input.Analyze_flag);
    set(FlagsValues(6),'value',pt.input.Geometry_flag);
    set(FlagsValues(7),'value',pt.input.Coordinate_flag);
    set(FlagsValues(8),'value',pt.input.Printing_flag);
    set(ToolsValues(1),'value',pt.input.Mp);
    set(ToolsValues(2),'value',pt.input.Np);
    set(ToolsValues(3),'value',pt.input.Nd);
    set(ToolsValues(4),'value',pt.input.ITER);
    set(FilenameValues,'value',pt.filename);
end

function Save(~,~,OpenPropDirectory)
    global Fig_Main pt;
    global SpecificationsValues FlagsValues Property1Values rx_in...
           CDp_x_in Meanline_x_in Thickness_x_in CL_x_in T0oDp_x_in...
           Property2Values ri_in VA_i_in VT_i_in DuctValues ToolsValues...
           FilenameValues;
    filename = get(FilenameValues,'value');
    
    % existence��ʾ�Ƿ����ͬ����Ŀ
    existence = ChangeDirectory(OpenPropDirectory,filename);
    if existence
        message = sprintf('��ǰ��Ŀ���Ѿ����� \n ����������иò�����������ԭ���ļ� \n ȷ��������');
        SaveSelection = uiconfirm(Fig_Main,message,'���Ǿ���',...
                                  'options',{'����ԭ���ļ�','���Ϊ������','ȡ������'},...
                                  'defaultoption',2,...
                                  'canceloption',3,...
                                  'icon','warning');
        if strcmp(SaveSelection,'����ԭ���ļ�')
            % ѡ�񡰸���ԭʼ�ļ���,����Ҫ���в���
        elseif strcmp(SaveSelection,'���Ϊ������')
            % ѡ�����Ϊ�����ơ�������SaveCopyAs����
            SaveCopyAs;
        else
            % ѡ��ȡ�����������˳�
            return
        end    
    end
    
    % ����������������е�����ֵ
    N = get(SpecificationsValues(1),'value');       % ���ת��(RPM)
    VS = get(SpecificationsValues(2),'value');      % �н��ٶ�(m/s)
    T = get(SpecificationsValues(3),'value');       % ������(N)
    Z = get(SpecificationsValues(4),'value');       % ҶƬ����
    Dp = get(SpecificationsValues(5),'value');      % ҶƬֱ��(m)
    Hub_flag = get(FlagsValues(1),'value');         % ������Ƿ���뽰챣�
    Dh = get(SpecificationsValues(6),'value');      % ���ֱ��(m)
    Duct_flag = get(FlagsValues(2),'value');        % ������Ƿ���뺭����
    Dd = get(SpecificationsValues(7),'value');      % ����ֱ��(m)
    Cd = get(SpecificationsValues(8),'value');      % �����ҳ�(m)
    
    % ������ҶƬ����������е�����ֵ
    Nx = get(Property1Values(1),'value');           % ��������
    ChordMethod = get(Property1Values(2),'value');  % ��Ʒ���
    Meanline_flag = get(FlagsValues(3),'value');    % ͹Ե���Ƿ�ͬһ��
    Thickness_flag = get(FlagsValues(4),'value');   % Ҷ���Ƿ�ͬһ��
    rx = zeros(Nx,1);                               % ����λ��
    CDp_x = zeros(Nx,1);                            % ����ϵ��
    Meanline_x = cell(Nx,1);                        % ͹Ե��
    Thickness_x = cell(Nx,1);                       % Ҷ��
    CL_x = zeros(Nx,1);                             % ����ϵ��
    T0oDp_x = zeros(Nx,1);                          % �񾶱�
    for index = 1 : Nx
        rx(index) = get(rx_in(index),'value');
        CDp_x(index) = get(CDp_x_in(index),'value');
        Meanline_x{index} = get(Meanline_x_in(index),'value');
        Thickness_x{index} = get(Thickness_x_in(index),'value');
        CL_x(index) = get(CL_x_in(index),'value');
        T0oDp_x(index) = get(T0oDp_x_in(index),'value');
    end    
    
    % �������ⲿ�����������е�������
    Ni = get(Property2Values(1),'value');           % ��������
    rho = get(Property2Values(2),'value');          % �����ܶ�
    ri = zeros(Ni,1);                               % ����λ��
    VA_i = zeros(Ni,1);                             % ���������ٶ�
    VT_i = zeros(Ni,1);                             % ���������ٶ�
    for index = 1 : Ni
        ri(index) = get(ri_in(index),'value');
        VA_i(index) = get(VA_i_in(index),'value');
        VT_i(index) = get(VT_i_in(index),'value');
    end    
    
    % ������������ز������е�������
    TdoT = get(DuctValues(1),'value');              % ��������
    CDd = get(DuctValues(2),'value');               % ����ϵ��
    
    % ���������������빤�ߡ��е�������
    Analyze_flag = get(FlagsValues(5),'value');     % �Ƿ�����������ߣ�
    Geometry_flag = get(FlagsValues(6),'value');    % �Ƿ������ʾģ�ͣ�
    Coordinate_flag = get(FlagsValues(7),'value');  % �Ƿ���������ꣿ
    Printing_flag = get(FlagsValues(8),'value');    % �Ƿ����stl�ļ���
    Mp = get(ToolsValues(1),'value');               % ҶƬ���Ƶ�����
    Np = get(ToolsValues(2),'value');               % ҶƬ���������
    Nd = get(ToolsValues(3),'value');               % ��������������
    ITER = get(ToolsValues(4),'value');             % ����������
    
    % ��GUI���������ֵ�������ṹ������input
    input.part1 = '���������';
    input.N = N;
    input.VS = VS;
    input.T = T;
    input.Z = Z;
    input.Dp = Dp;
    input.Hub_flag = Hub_flag;
    input.Dh = Dh;
    input.Duct_flag = Duct_flag;
    input.Dd = Dd;
    input.Cd = Cd;
    
    input.part2 = 'ҶƬ�������';
    input.Nx = Nx;
    input.ChordMethod = ChordMethod;
    input.Meanline_flag = Meanline_flag;
    input.Thickness_flag = Thickness_flag;
    input.rx = rx;
    input.CDp_x = CDp_x;
    input.Meanline_x = Meanline_x;
    input.Thickness_x = Thickness_x;
    input.CL_x = CL_x;
    input.T0oDp_x = T0oDp_x;
    
    input.part3 = '�ⲿ��������';
    input.Ni = Ni;
    input.rho = rho;
    input.ri = ri;
    input.VA_i = VA_i;
    input.VT_i = VT_i;
    
    input.part4 = '������ز���';
    input.TdoT = TdoT;
    input.CDd = CDd;
    
    input.part5 = '���������빤��';
    input.Analyze_flag = Analyze_flag;
    input.Geometry_flag = Geometry_flag;
    input.Coordinate_flag = Coordinate_flag;
    input.Printing_flag = Printing_flag;
    input.Mp = Mp;
    input.Np = Np;
    input.Nd = Nd;
    input.ITER = ITER;
    
    % ��������Ŀ��Ϣ��װ���ṹ������pt��
    pt.filename = filename;
    pt.date = date;
    pt.input = input;
    pt.design = [];
    pt.geometry = [];
    pt.states = [];
    
    % ����
    save(filename,'pt');
end

function SaveCopyAs(~,~)
    global pt;
    global SpecificationsValues FlagsValues Property1Values rx_in...
           CDp_x_in Meanline_x_in Thickness_x_in CL_x_in T0oDp_x_in...
           Property2Values ri_in VA_i_in VT_i_in DuctValues ToolsValues...
           FilenameValues;
    filename = get(FilenameValues,'value');
    
    % ����������������е�����ֵ
    N = get(SpecificationsValues(1),'value');       % ���ת��(RPM)
    VS = get(SpecificationsValues(2),'value');      % �н��ٶ�(m/s)
    T = get(SpecificationsValues(3),'value');       % ������(N)
    Z = get(SpecificationsValues(4),'value');       % ҶƬ����
    Dp = get(SpecificationsValues(5),'value');      % ҶƬֱ��(m)
    Hub_flag = get(FlagsValues(1),'value');         % ������Ƿ���뽰챣�
    Dh = get(SpecificationsValues(6),'value');      % ���ֱ��(m)
    Duct_flag = get(FlagsValues(2),'value');        % ������Ƿ���뺭����
    Dd = get(SpecificationsValues(7),'value');      % ����ֱ��(m)
    Cd = get(SpecificationsValues(8),'value');      % �����ҳ�(m)
    
    % ������ҶƬ����������е�����ֵ
    Nx = get(Property1Values(1),'value');           % ��������
    ChordMethod = get(Property1Values(2),'value');  % ��Ʒ���
    Meanline_flag = get(FlagsValues(3),'value');    % ͹Ե���Ƿ�ͬһ��
    Thickness_flag = get(FlagsValues(4),'value');   % Ҷ���Ƿ�ͬһ��
    rx = zeros(Nx,1);                               % ����λ��
    CDp_x = zeros(Nx,1);                            % ����ϵ��
    Meanline_x = cell(Nx,1);                        % ͹Ե��
    Thickness_x = cell(Nx,1);                       % Ҷ��
    CL_x = zeros(Nx,1);                             % ����ϵ��
    T0oDp_x = zeros(Nx,1);                          % �񾶱�
    for index = 1 : Nx
        rx(index) = get(rx_in(index),'value');
        CDp_x(index) = get(CDp_x_in(index),'value');
        Meanline_x{index} = get(Meanline_x_in(index),'value');
        Thickness_x{index} = get(Thickness_x_in(index),'value');
        CL_x(index) = get(CL_x_in(index),'value');
        T0oDp_x(index) = get(T0oDp_x_in(index),'value');
    end    
    
    % �������ⲿ�����������е�������
    Ni = get(Property2Values(1),'value');           % ��������
    rho = get(Property2Values(2),'value');          % �����ܶ�
    ri = zeros(Ni,1);                               % ����λ��
    VA_i = zeros(Ni,1);                             % ���������ٶ�
    VT_i = zeros(Ni,1);                             % ���������ٶ�
    for index = 1 : Ni
        ri(index) = get(ri_in(index),'value');
        VA_i(index) = get(VA_i_in(index),'value');
        VT_i(index) = get(VT_i_in(index),'value');
    end    
    
    % ������������ز������е�������
    TdoT = get(DuctValues(1),'value');              % ��������
    CDd = get(DuctValues(2),'value');               % ����ϵ��
    
    % ���������������빤�ߡ��е�������
    Analyze_flag = get(FlagsValues(5),'value');     % �Ƿ�����������ߣ�
    Geometry_flag = get(FlagsValues(6),'value');    % �Ƿ������ʾģ�ͣ�
    Coordinate_flag = get(FlagsValues(7),'value');  % �Ƿ���������ꣿ
    Printing_flag = get(FlagsValues(8),'value');    % �Ƿ����stl�ļ���
    Mp = get(ToolsValues(1),'value');               % ҶƬ���Ƶ�����
    Np = get(ToolsValues(2),'value');               % ҶƬ���������
    Nd = get(ToolsValues(3),'value');               % ��������������
    ITER = get(ToolsValues(4),'value');             % ����������
    
    % ��GUI���������ֵ�������ṹ������input
    input.part1 = '���������';
    input.N = N;
    input.VS = VS;
    input.T = T;
    input.Z = Z;
    input.Dp = Dp;
    input.Hub_flag = Hub_flag;
    input.Dh = Dh;
    input.Duct_flag = Duct_flag;
    input.Dd = Dd;
    input.Cd = Cd;
    
    input.part2 = 'ҶƬ�������';
    input.Nx = Nx;
    input.ChordMethod = ChordMethod;
    input.Meanline_flag = Meanline_flag;
    input.Thickness_flag = Thickness_flag;
    input.rx = rx;
    input.CDp_x = CDp_x;
    input.Meanline_x = Meanline_x;
    input.Thickness_x = Thickness_x;
    input.CL_x = CL_x;
    input.T0oDp_x = T0oDp_x;
    
    input.part3 = '�ⲿ��������';
    input.Ni = Ni;
    input.rho = rho;
    input.ri = ri;
    input.VA_i = VA_i;
    input.VT_i = VT_i;
    
    input.part4 = '������ز���';
    input.TdoT = TdoT;
    input.CDd = CDd;
    
    input.part5 = '���������빤��';
    input.Analyze_flag = Analyze_flag;
    input.Geometry_flag = Geometry_flag;
    input.Coordinate_flag = Coordinate_flag;
    input.Printing_flag = Printing_flag;
    input.Mp = Mp;
    input.Np = Np;
    input.Nd = Nd;
    input.ITER = ITER;
    
    % ��������Ŀ��Ϣ��װ���ṹ������pt��
    pt.filename = filename;
    pt.date = date;
    pt.input = input;
    pt.design = [];
    pt.geometry = [];
    pt.states = [];
    
    % Ŀǰ��������Ϊ��ɺ���Ļ�����ص�matlab�����������Ҫ���
    uisave('pt',[filename,'-copy']);
end

