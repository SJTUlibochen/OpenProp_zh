% 该程序负责OpenProp中文版主界面GUI的生成
%% 主界面
function MainPage
    %% 前处理
    clear variables;
    clear global;
    global titlefontsize zhfontsize subcolumnfontsize numfontsize OpenPropDirectory;
    global Fig_Main Menu TabPageMenu Tabs SpecificationsValues FlagsValues...
           Property1Values Table1PanelGrid rx_in CDp_x_in Meanline_x_in...
           Thickness_x_in CL_x_in T0oDp_x_in Property2Values Table2PanelGrid...
           ri_in VA_i_in VT_i_in DuctValues ToolsValues FilenameValues;
    %% 设定各参数的缺省值
    % 螺旋桨规格    
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
    % 叶片切面参数
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
    % 外部环境参数
    Ni_def = 10;
    rho_def = 1.29;
    ri_def = [0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;0.95;1];
    VA_i_def = ones(Ni_def,1);
    VT_i_def = zeros(Ni_def,1);
    % 涵道相关参数
    TpoT_def = 1;
    TdoT_def = 1-TpoT_def;
    CDd_def = 0.008;
    % 其他参数与工具
    Analyze_flag_def = 1;
    Geometry_flag_def = 1;
    Coordinate_flag_def = 1;
    Printing_flag_def = 1;
    Mp_def = 20;
    Np_def = 20;
    Nd_def = 12;
    ITER_def = 50;
    filename_def = 'DefaultProject';
    %% GUI元素的参数
    % 字体大小
    titlefontsize = 60;
    zhfontsize = 15;
    subcolumnfontsize = 17;
    numfontsize = 13;
    % 主界面初始高度和宽度
    Windowht = 750;
    Window = 1530;
    % 设定文件目录
    OpenPropDirectory = 'OpenProp_zh';
    %% 主界面中的GUI元素
    % 生成主界面
    Fig_Main = uifigure('position',[5 55 Window Windowht],...
                        'numbertitle','off',...
                        'name','OpenProp');
    % 生成顶部菜单栏
    MenuStrings = {'文件' '设置' '标签页'};
    Menu = zeros(1,length(MenuStrings));
    for index = 1 : length(MenuStrings)
        Menu(index) = uimenu('parent',Fig_Main,...
                             'text',MenuStrings{index});
    end    
    % 生成子工具栏
    FileStrings = {'保存' '另存为' '导入' '运行'};
    FileMenu = zeros(1,length(FileStrings));
    for index = 1 : length(FileStrings)
        FileMenu(index) = uimenu('parent',Menu(1),...
                                 'text',FileStrings{index});
    end
    
    SettingStrings = {'切换语言' '夜晚模式'};
    SettingMenu = zeros(1,length(SettingStrings));
    for index = 1 : length(SettingStrings)
        SettingMenu(index) = uimenu('parent',Menu(2),...
                                    'text',SettingStrings{index});
    end    
    
    TabPageStrings = {'仅进行叶片设计' '研究参数的影响'};
    for index = 1 : length(TabPageStrings)
        TabPageMenu(index) = uimenu('parent',Menu(3),...
                                    'text',TabPageStrings{index},...
                                    'enable','off');
    end    
    % 生成主界面中的标签栏uitabgroup
    Tabs = uitabgroup('parent',Fig_Main,...
                      'position',[0 0 Window Windowht]);
    % 生成Tabs的第一个标签页“仅进行叶片设计”
    SingleDesign = uitab('parent',Tabs,...
                         'title','仅进行叶片设计');
    % 生成Tabs的第二个标签页“研究参数的影响”
    ParametricStudy = uitab('parent',Tabs,...
                            'title','研究参数的影响');
    %% 标签页“仅进行叶片设计”中的GUI元素
    % 生成SingleDesign中的网格空间
    SingleDesignGrid = uigridlayout(SingleDesign,[3 3],...
                                    'columnwidth',{'1x','3x','1.5x'},...
                                    'rowheight',{'4x','1x','1.1x'},...
                                    'padding',[10 10 10 10]);
                                
    % 螺旋桨规格 Propeller Specification
    Specifications = uipanel('parent',SingleDesignGrid,...
                             'title','螺旋桨规格',...
                             'titleposition','centertop',...
                             'fontsize',subcolumnfontsize,...
                             'fontweight','bold',...
                             'scrollable','off');
                         
    % 叶片切面参数 Blade Design Values
    BladeDesign = uipanel('parent',SingleDesignGrid,...
                          'title','叶片切面参数',...
                          'titleposition','centertop',...
                          'fontsize',subcolumnfontsize,...
                          'fontweight','bold',...
                          'scrollable','on');
    BladeDesign.Layout.Row = [1 2];
    
    % 外部环境参数 Inflow Profile Values
    Inflow = uipanel('parent',SingleDesignGrid,...
                     'title','外部环境参数',...
                     'titleposition','centertop',...
                     'fontsize',subcolumnfontsize,...
                     'fontweight','bold',...
                     'scrollable','on');
    Inflow.Layout.Row = [1 2];
    
    % 涵道相关参数 Duct Design Values
    Duct = uipanel('parent',SingleDesignGrid,...
                   'title','涵道相关参数',...
                   'titleposition','centertop',...
                   'fontsize',subcolumnfontsize,...
                   'fontweight','bold',...
                   'scrollable','off');
    Duct.Layout.Row = 2;
    Duct.Layout.Column = 1;
    
    % 标题 Title
    Title = uilabel(SingleDesignGrid,...
                    'text','OpenProp',...
                    'horizontalalignment','center',...
                    'fontname','Bauhaus 93',...
                    'fontsize',titlefontsize,...
                    'fontweight','normal');
     
    % 其他参数与工具 Other Desgin Values & Tools
    Tools = uipanel('parent',SingleDesignGrid,...
                    'title','其他参数与工具',...
                    'titleposition','centertop',...
                    'fontsize',subcolumnfontsize,...
                    'fontweight','bold',...
                    'scrollable','off');
    Tools.Layout.Column = [2 3];
    
    %% 标签页“研究参数的影响”中的GUI元素
    % 生成ParametricStudy中的网格
    ParametricStudyGrid = uigridlayout(ParametricStudy,[1 1]);
    Announcement = uilabel(ParametricStudyGrid,...
                           'text','该功能尚未完成，敬请期待',...
                           'horizontalalignment','center',...
                           'verticalalignment','center',...
                           'fontsize',zhfontsize);
                           
    %% 子栏“螺旋桨规格”中的GUI元素
    % 生成Specifications中的网格
    SpecificationsGrid = uigridlayout(Specifications,[4 1],...
                                      'rowheight',{'3x','2x','2x','3x'},...
                                      'columnspacing',5,...
                                      'rowspacing',5,...
                                      'padding',[5 5 5 5],...
                                      'scrollable','off');   
                                  
    % GUI元素
    SpecificationsStrings = {'电机转速：','行进速度：','总升力：','叶片数量：',...
                             '叶片直径：','桨毂直径：','涵道直径：','涵道弦长：'};
    SpecificationsTips = {'N','VS','T','Z','Dp','Dh','Dd','Cd'};
    SpecificationsValues_def = {N_def,VS_def,T_def,Z_def,Dp_def,Dh_def,Dd_def,Cd_def};
    
    % 装置面板
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
    
    % 叶片面板
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
    
    % 桨毂面板
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
                                'text','设计中是否加入桨毂？',...
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
                                      
    % 涵道面板
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
                                'text','设计中是否加入涵道？',...
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
    
    % 设置Specifications编辑框中的显示形式和单位
    set(SpecificationsValues(1),'valuedisplayformat','%.0f RPM');
    set(SpecificationsValues(2),'valuedisplayformat','%.2f m/s');
    set(SpecificationsValues(3),'valuedisplayformat','%.3f N');
    set(SpecificationsValues(4),'valuedisplayformat','%.0f');
    set(SpecificationsValues(5),'valuedisplayformat','%.3f m');
    set(SpecificationsValues(6),'valuedisplayformat','%.3f m');
    set(SpecificationsValues(7),'valuedisplayformat','%.3f m');
    set(SpecificationsValues(8),'valuedisplayformat','%.3f m');
    
    %% 子栏“叶片切面参数”中的GUI元素
    % 生成BladeDesign中的网格
    BladeDesignGrid = uigridlayout(BladeDesign,[2 1],...
                                   'rowheight',{'2x','21x'},...
                                   'columnspacing',5,...
                                   'rowspacing',5,...
                                   'padding',[5 5 5 5],...
                                   'scrollable','off');
                               
    % 属性1面板
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
                                'text','切面数量：',...
                                'horizontalalignment','left',...
                                'verticalalignment','center',...
                                'tooltip','Nx，控制输入的切面数量',...
                                'fontsize',zhfontsize);
    Property1Values(1) = uieditfield(Property1PanelSubGrid(1),'numeric',...
                                     'fontsize',numfontsize,...
                                     'horizontalalignment','center',...
                                     'value',Nx_def,...
                                     'limit',[0 Inf],...
                                     'lowerlimitinclusive','off');
    Property1Texts(2) = uilabel('parent',Property1PanelSubGrid(2),...
                                'text','设计方法：',...
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
                                'text','凸缘线是否同一？',...
                                'value',Meanline_flag_def,...
                                'tooltip','Meanline_flag',...
                                'fontsize',zhfontsize);
    FlagsValues(4) = uicheckbox('parent',Property1PanelGrid,...
                                'text','叶型是否同一？',...
                                'value',Thickness_flag_def,...
                                'tooltip','Thickness_flag',...
                                'fontsize',zhfontsize);
                                
    % 表格1面板
    Table1Panel = uipanel(BladeDesignGrid);
    Table1PanelGrid = uigridlayout(Table1Panel,[Nx_def+1 6],...
                                   'columnspacing',2,...
                                   'rowspacing',2,...
                                   'padding',[10 10 10 0]);
    Table1Strings = {'径向位置','阻力系数','凸缘线','叶型','升力系数','厚径比'};
    Table1Tips = {'rx','CDp_x','Meanline','Thickness',...
                  'CL_x','T0oDp_x，叶厚与叶片直径之比'};
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
    % 缺省状态下Meanline和Thickness下拉框只有第一行可用
    set(Meanline_x_in(1),'enable','on');
    set(Thickness_x_in(1),'enable','on');
    %% 子栏“外部环境参数”中的GUI元素
    % 生成Inflow中的网格
    InflowGrid = uigridlayout(Inflow,[2 1],...
                              'rowheight',{'2x','21x'},...
                              'columnspacing',5,...
                              'rowspacing',5,...
                              'padding',[5 5 5 5],...
                              'scrollable','off');
                          
    % 属性2面板
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
                                'text','坐标数量：',...
                                'horizontalalignment','left',...
                                'verticalalignment','center',...
                                'tooltip','Ni，控制输入的坐标数量',...
                                'fontsize',zhfontsize);
    Property2Values(1) = uieditfield(Property2PanelSubGrid(1),'numeric',...
                                     'fontsize',numfontsize,...
                                     'horizontalalignment','center',...
                                     'value',Ni_def,...
                                     'limit',[0 Inf],...
                                     'lowerlimitinclusive','off');
    Property2Texts(2) = uilabel('parent',Property2PanelSubGrid(2),...
                                'text','流体密度：',...
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
                                
    % 表格2面板
    Table2Panel = uipanel(InflowGrid);
    Table2PanelGrid = uigridlayout(Table2Panel,[Ni_def+1 3],...
                                   'columnspacing',2,...
                                   'rowspacing',2,...
                                   'padding',[10 10 10 0]);
    Table2Strings = {'径向位置','轴向来流速度','切向来流速度'};
    Table2Tips = {'ri','VA_i，轴向来流速度与行进速度的比值',...
                  'VT_i，切向来流速度与行进速度的比值'};
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
    %% 子栏“涵道相关参数”中的GUI元素
    % 生成Duct中的网格
    DuctGrid = uigridlayout(Duct,[2 2],...
                            'columnwidth',{'2x','3.5x'},...
                            'rowheight',{'1x','1x'},...
                            'columnspacing',10,...
                            'rowspacing',10,...
                            'padding',[15 10 15 10],...
                            'scrollable','off');
    DuctStrings = {'升力比重：','阻力系数：'};
    DuctTips = {'TdoT，涵道所提供升力占总升力的比值','CDd'};
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
                        
    %% 子栏“其他设计参数与工具”中的GUI元素
    % 生成Tools中的网格
    ToolsGrid = uigridlayout(Tools,[1,3],...
                             'columnwidth',{'1.2x','2x','1x'},...
                             'columnspacing',5,...
                             'rowspacing',5,...
                             'padding',[5 5 5 5],...
                             'scrollable','off');
    
    % 复选框面板
    FlagsPanel = uipanel(ToolsGrid);
    FlagsPanelGrid = uigridlayout(FlagsPanel,[2,2]);
    FlagsStrings = {'是否绘制性能曲线？','是否绘制演示模型？',...
                    '是否输出点坐标？','是否输出stl文件？'};
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
    % 控制点面板
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
                            'text','叶片控制点数量：',...
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
                            'text','叶片切向段数量：',...
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
                            'text','涵道涡流环数量：',...
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
                            'text','最大迭代次数：',...
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
    
    % 文件面板
    FilePanel = uipanel(ToolsGrid);
    FilePanelGrid = uigridlayout(FilePanel,[2 4]);
    FilenameTexts = uilabel('parent',FilePanelGrid,...
                            'text','项目名：',...
                            'horizontalalignment','left',...
                            'verticalalignment','center',...
                            'fontsize',zhfontsize);
    FilenameValues = uieditfield(FilePanelGrid,...
                                 'fontsize',numfontsize,...
                                 'horizontalalignment','center',...
                                 'value',filename_def);
    FilenameValues.Layout.Column = [2 4];
    LoadButton = uibutton('parent',FilePanelGrid,...
                          'text','导入',...
                          'icon','load.png');
    SaveButton = uibutton('parent',FilePanelGrid,...
                          'text','保存',...
                          'icon','save.png',...
                          'buttonpushedfcn',@Save);
    SaveCopyAsButton = uibutton('parent',FilePanelGrid,...
                                'text','另存',...
                                'icon','savecopy.png');
    RunButton = uibutton('parent',FilePanelGrid,...
                         'text','运行',...
                         'icon','run.png');
                     
    %% 设置回调函数
    % 设置菜单栏的回调函数
    set(FileMenu(1),'menuselectedfcn',@Save);
    set(FileMenu(2),'menuselectedfcn',@SaveCopyAs);
    set(FileMenu(3),'menuselectedfcn',@Load);
    set(FileMenu(4),'menuselectedfcn',@Execute);
    
    % 设置Specifications编辑框的回调函数
    set(SpecificationsValues,'valuechangedfcn',{@ChangeSpecifications,...
                                                Dp_def,Dh_def,Dd_def});
    
    % 设置Specifications复选框的回调函数
    set(FlagsValues(1),'valuechangedfcn',{@ChangeHub,Dh_def});
    set(FlagsValues(2),'valuechangedfcn',{@ChangeDuct,Dd_def,Cd_def,...
                                          TdoT_def,CDd_def});
    
    % 设置BladeDesign表格数控制的回调函数
    set(Property1Values(1),'valuechangedfcn',{@ChangeNx,Meanline_x_items,...
                                              Thickness_x_items});
    
    % 设置BladeDesign表格内容控制的回调函数
    set(Property1Values(2),'valuechangedfcn',@ChangeMethod);
    set(FlagsValues(3),'valuechangedfcn',@ConstantMeanline);
    set(FlagsValues(4),'valuechangedfcn',@ConstantThickness);
    set(Meanline_x_in(1),'valuechangedfcn',@Change1stMeanline);
    set(Thickness_x_in(1),'valuechangedfcn',@Change1stThickness);
    
    % 设置Inflow表格数控制的回调函数
    set(Property2Values(1),'valuechangedfcn',@ChangeNi);
    
    % 设置Tools编辑框的回调函数
    set(ToolsValues,'valuechangedfcn',@ChangeTools)
    
    % 设置Tools按钮的回调函数
    set(LoadButton,'buttonpushedfcn',{@Load,Meanline_x_items,Thickness_x_items});
    set(SaveButton,'buttonpushedfcn',@Save);
    set(SaveCopyAsButton,'buttonpushedfcn',@SaveCopyAs);
    set(RunButton,'buttonpushedfcn',@Execute);
    
end

%% 子栏“螺旋桨规格”中涉及到的回调函数
% 编辑框的回调函数
function ChangeSpecifications(~,~,Dp_def,Dh_def,Dd_def)
    global Fig_Main SpecificationsValues;
    % 调整输入值，避免出现无意义的情况
    Z = round(get(SpecificationsValues(4),'value'));
    set(SpecificationsValues(4),'value',Z);
    Dp = get(SpecificationsValues(5),'value');
    Dh = get(SpecificationsValues(6),'value');
    Dd = get(SpecificationsValues(7),'value');
    % 检查桨毂直径与螺旋桨直径的输入值，若前者大于后者，则弹出提醒框
    if Dh >= Dp
        message1 = sprintf('当前设置的桨毂直径大于螺旋桨直径 \n 请重新设置！');
        uialert(Fig_Main,message1,'参数设置错误',...
                'closefcn',{@ResetDiameters1,SpecificationsValues,...
                            Dp_def,Dh_def});
    end
    % 检查涵道直径与螺旋桨直径的输入值，若前者小于后者，则弹出提醒框
    if Dd < Dp
        message2 = sprintf('当前设置的涵道直径小于螺旋桨直径 \n 请重新设置！');
        uialert(Fig_Main,message2,'参数设置错误',...
                'closefcn',{@ResetDiameters2,SpecificationsValues,...
                            Dp_def,Dd_def});
    end
end
% 发生直径设置错误后调用的函数
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

% 复选框的回调函数
function ChangeHub(hObject,~,Dh_def)
    global SpecificationsValues;
    if get(hObject,'value')
        % 勾选“加入桨毂”，允许编辑桨毂直径
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
        % 勾选“加入涵道”，允许编辑涵道相关控件
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

%% 子栏“叶片切面参数”中涉及到的回调函数
% 表格数控制
function ChangeNx(hObject,ED,Meanline_x_items,Thickness_x_items)
    global numfontsize;
    global Table1PanelGrid rx_in CDp_x_in Meanline_x_in Thickness_x_in...
           CL_x_in T0oDp_x_in Property1Values FlagsValues;
    % 调整输入值，避免出现无意义的情况  
    set(hObject,'value',round(get(hObject,'value')));
    Nx_new = get(hObject,'value');
    % ED中储存了和本次值更改相关的数据，如更改前的值
    Nx_previous = ED.PreviousValue;
    
    % 将GUI中原本的输入值导出
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
    % rx新缺省值的确定思路：保留两个端点和最后一段的中值点，其余位置平均分配
    rx_new = zeros(Nx_new,1);
    for index = 1 : Nx_new-2
        rx_new(index) = 0.2+(index-1)*(1-0.2)/(Nx_new-1-1);
    end
    rx_new(Nx_new-1) = (1+rx_new(Nx_new-2))/2;
    rx_new(Nx_new) = 1;
    % CDp_x、CL_x、T0oDp_x新缺省值的确定思路：pchip插值
    CDp_x_new = pchip(rx_previous,CDp_x_previous,rx_new);
    CL_x_new = pchip(rx_previous,CL_x_previous,rx_new);
    T0oDp_x_new = pchip(rx_previous,T0oDp_x_previous,rx_new);
    
    RowHeightCell = cell(1,Nx_new+1);
    for index = 1 : Nx_new+1
        RowHeightCell{index} = '1x';
    end    
    if Nx_new > Nx_previous
        % 如果修改后表格数增加
        set(Table1PanelGrid,'rowheight',RowHeightCell);
        % 修改原有表格输入值
        for index = 1 : Nx_previous
            set(rx_in(index),'value',rx_new(index));
            set(CDp_x_in(index),'value',CDp_x_new(index));
            set(CL_x_in(index),'value',CL_x_new(index));
            set(T0oDp_x_in(index),'value',T0oDp_x_new(index));
        end   
        % 生成新的表格
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
                % 设计方法为'CLmax'，允许编辑升力系数和厚径比
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
        % 删除多余的表格
        for index = Nx_new+1 : Nx_previous
            delete(rx_in(index));
            delete(CDp_x_in(index));
            delete(Meanline_x_in(index));
            delete(Thickness_x_in(index));
            delete(CL_x_in(index));
            delete(T0oDp_x_in(index));
        end  
        set(Table1PanelGrid,'rowheight',RowHeightCell);
        % 修改剩余表格输入值
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

% 表格内容控制
function ChangeMethod(hObject,~)
    global Property1Values CL_x_in T0oDp_x_in;
    Nx = get(Property1Values(1),'value');
    method = get(hObject,'value');
    if strcmp(method,'CLmax')
        % 设计方法为'CLmax'，允许编辑升力系数和厚径比
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
        % 勾选“凸缘线同一”，除第一行外下拉框均不可编辑
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
        % 勾选“叶型同一”，除第一行外下拉框均不可编辑
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
        % 勾选“凸缘线同一”，第一行更改后其余诸行一并更改
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
        % 勾选“叶型同一”，第一行更改后其余诸行一并更改
        for index = 2 : Nx
            set(Thickness_x_in(index),'value',thickness);
        end
    else
        return
    end    
end

%% 子栏“外部环境参数”中涉及到的回调函数
function ChangeNi(hObject,ED)
    global numfontsize;
    global Table2PanelGrid ri_in VA_i_in VT_i_in;
    % 调整输入值，避免出现无意义的情况  
    set(hObject,'value',round(get(hObject,'value')));
    Ni_new = get(hObject,'value');
    % ED中储存了和本次值更改相关的数据，如更改前的值
    Ni_previous = ED.PreviousValue;
    
    % 将GUI中原本的输入值导出
    ri_previous = ones(Ni_previous,1);
    VA_i_previous = ones(Ni_previous,1);
    VT_i_previous = ones(Ni_previous,1);
    for index = 1 : Ni_previous
        ri_previous(index) = get(ri_in(index),'value');
        VA_i_previous(index) = get(VA_i_in(index),'value');
        VT_i_previous(index) = get(VT_i_in(index),'value');
    end    
    % ri新缺省值的确定思路：保留两个端点和最后一段的中值点，其余位置平均分配
    ri_new = zeros(Ni_new,1);
    for index = 1 : Ni_new-2
        ri_new(index) = 0.2+(index-1)*(1-0.2)/(Ni_new-1-1);
    end
    ri_new(Ni_new-1) = (1+ri_new(Ni_new-2))/2;
    ri_new(Ni_new) = 1;
    % VA_i、VT_i新缺省值的确定思路：pchip插值
    VA_i_new = pchip(ri_previous,VA_i_previous,ri_new);
    VT_i_new = pchip(ri_previous,VT_i_previous,ri_new);
    
    RowHeightCell = cell(1,Ni_new+1);
    for index = 1 : Ni_new+1
        RowHeightCell{index} = '1x';
    end    
    if Ni_new > Ni_previous
        % 如果修改后表格数增加
        set(Table2PanelGrid,'rowheight',RowHeightCell);
        % 修改原有表格输入值
        for index = 1 : Ni_previous
            set(ri_in(index),'value',ri_new(index));
            set(VA_i_in(index),'value',VA_i_new(index));
            set(VT_i_in(index),'value',VT_i_new(index));
        end   
        % 生成新的表格
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
        % 删除多余的表格
        for index = Ni_new+1 : Ni_previous
            delete(ri_in(index));
            delete(VA_i_in(index));
            delete(VT_i_in(index));
        end  
        set(Table2PanelGrid,'rowheight',RowHeightCell);
        % 修改剩余表格输入值
        for index = 1 : Ni_new
            set(ri_in(index),'value',ri_new(index));
            set(VA_i_in(index),'value',VA_i_new(index));
            set(VT_i_in(index),'value',VT_i_new(index));
        end   
    else
        return;
    end    
end

%% 子栏“其他参数与工具”中涉及到的回调函数
% 控制点输入值修正
function ChangeTools(hObject,~)
    % 调整输入值，避免出现无意义的情况
    input = round(get(hObject,'value'));
    set(hObject,'value',input);
end

% 文件相关操作的回调函数
function Load(~,~,Meanline_x_items,Thickness_x_items)
    global pt numfontsize OpenPropDirectory;
    global SpecificationsValues FlagsValues Property1Values Table1PanelGrid...
           rx_in CDp_x_in Meanline_x_in Thickness_x_in CL_x_in T0oDp_x_in...
           Property2Values Table2PanelGrid ri_in VA_i_in VT_i_in DuctValues...
           ToolsValues FilenameValues;
    rest = pwd;
    LoadDirectory(rest,OpenPropDirectory);
    % uiload函数用于打开文件选择目录
    uiload;
    
    % 获取按下导入按钮前GUI界面中的Nx和Ni
    Nx_new = pt.input.Nx;
    Nx_previous = get(Property1Values(1),'value');
    Ni_new = pt.input.Ni;
    Ni_previous = get(Property2Values(1),'value');
    
    RowHeightCell = cell(1,Nx_new+1);
    for index = 1 : Nx_new+1
        RowHeightCell{index} = '1x';
    end    
    
    % 将数据导入“螺旋桨规格”
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
    
    % 将数据导入“叶片切面参数”
    set(Property1Values(1),'value',pt.input.Nx);
    set(Property1Values(2),'value',pt.input.ChordMethod);
    set(FlagsValues(3),'value',pt.input.Meanline_flag);
    set(FlagsValues(4),'value',pt.input.Thickness_flag);
    if Nx_new > Nx_previous
        % 导入数据的切面数量多于原GUI界面上显示的
        set(Table1PanelGrid,'rowheight',RowHeightCell);
        % 修改原有表格输入值
        for index = 1 : Nx_previous
            set(rx_in(index),'value',pt.input.rx(index));
            set(CDp_x_in(index),'value',pt.input.CDp_x(index));
            set(Meanline_x_in(index),'value',pt.input.Meanline_x{index});
            set(Thickness_x_in(index),'value',pt.input.Thickness_x{index});
            set(CL_x_in(index),'value',pt.input.CL_x(index));
            set(T0oDp_x_in(index),'value',pt.input.T0oDp_x(index));
            
            if strcmp(pt.input.ChordMethod,'CLmax')
                % 设计方法为'CLmax'，允许编辑升力系数和厚径比
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
        % 始终保证下拉框第一行可以使用
        set(Meanline_x_in(1),'enable','on');
        set(Thickness_x_in(1),'enable','on');
        % 生成新的表格
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
                % 设计方法为'CLmax'，允许编辑升力系数和厚径比
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
        % 删除多余的表格
        for index = Nx_new+1 : Nx_previous
            delete(rx_in(index));
            delete(CDp_x_in(index));
            delete(Meanline_x_in(index));
            delete(Thickness_x_in(index));
            delete(CL_x_in(index));
            delete(T0oDp_x_in(index));
        end  
        set(Table1PanelGrid,'rowheight',RowHeightCell);
        % 修改剩余表格输入值
        for index = 1 : Nx_new
            set(rx_in(index),'value',pt.input.rx(index));
            set(CDp_x_in(index),'value',pt.input.CDp_x(index));
            set(Meanline_x_in(index),'value',pt.input.Meanline_x{index});
            set(Thickness_x_in(index),'value',pt.input.Thickness_x{index});
            set(CL_x_in(index),'value',pt.input.CL_x(index));
            set(T0oDp_x_in(index),'value',pt.input.T0oDp_x(index));
            
            if strcmp(pt.input.ChordMethod,'CLmax')
                % 设计方法为'CLmax'，允许编辑升力系数和厚径比
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
        % 始终保证下拉框第一行可以使用
        set(Meanline_x_in(1),'enable','on');
        set(Thickness_x_in(1),'enable','on');
    else    
        return
    end    
    
    % 将数据导入“外部环境参数”
    set(Property2Values(1),'value',pt.input.Ni);
    set(Property2Values(2),'value',pt.input.rho);
    if Ni_new > Ni_previous
        % 如果修改后表格数增加
        set(Table2PanelGrid,'rowheight',RowHeightCell);
        % 修改原有表格输入值
        for index = 1 : Ni_previous
            set(ri_in(index),'value',pt.input.ri(index));
            set(VA_i_in(index),'value',pt.input.VA_i(index));
            set(VT_i_in(index),'value',pt.input.VT_i(index));
        end   
        % 生成新的表格
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
        % 删除多余的表格
        for index = Ni_new+1 : Ni_previous
            delete(ri_in(index));
            delete(VA_i_in(index));
            delete(VT_i_in(index));
        end  
        set(Table2PanelGrid,'rowheight',RowHeightCell);
        % 修改剩余表格输入值
        for index = 1 : Ni_new
            set(ri_in(index),'value',pt.input.ri(index));
            set(VA_i_in(index),'value',pt.input.VA_i(index));
            set(VT_i_in(index),'value',pt.input.VT_i(index));
        end   
    else
        return;
    end 
    
    % 将数据导入“涵道相关参数”
    set(DuctValues(1),'value',pt.input.TdoT);
    set(DuctValues(2),'value',pt.input.CDd);
    
    % 将数据导入“其他参数与工具”
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

function Save(~,~)
    global pt OpenPropDirectory;
    global Fig_Main SpecificationsValues FlagsValues Property1Values rx_in...
           CDp_x_in Meanline_x_in Thickness_x_in CL_x_in T0oDp_x_in...
           Property2Values ri_in VA_i_in VT_i_in DuctValues ToolsValues...
           FilenameValues;
    filename = get(FilenameValues,'value');
    % existence表示是否存在同名项目
    existence = ChangeDirectory(OpenPropDirectory,filename);
    if existence
        message = sprintf('当前项目名已经存在 \n 如果继续进行该操作，将覆盖原本文件 \n 确定继续吗？');
        SaveSelection = uiconfirm(Fig_Main,message,'覆盖警告',...
                                  'options',{'覆盖原版文件','另存为新名称','取消操作'},...
                                  'defaultoption',2,...
                                  'canceloption',3,...
                                  'icon','warning');
        if strcmp(SaveSelection,'覆盖原版文件')
            % 选择“覆盖原始文件”,不需要进行操作
        elseif strcmp(SaveSelection,'另存为新名称')
            % 选择“另存为新名称”，调用SaveCopyAs函数
            SaveCopyAs;
        else
            % 选择“取消操作”，退出
            % 这里使用return返回到函数Execute后，仍会继续进行计算和绘图工作
            return
        end    
    end
    
    % 子栏“螺旋桨规格”中的输入值
    N = get(SpecificationsValues(1),'value');       % 电机转速(RPM)
    VS = get(SpecificationsValues(2),'value');      % 行进速度(m/s)
    T = get(SpecificationsValues(3),'value');       % 总升力(N)
    Z = get(SpecificationsValues(4),'value');       % 叶片数量
    Dp = get(SpecificationsValues(5),'value');      % 叶片直径(m)
    Hub_flag = get(FlagsValues(1),'value');         % 设计中是否加入桨毂？
    Dh = get(SpecificationsValues(6),'value');      % 桨毂直径(m)
    Duct_flag = get(FlagsValues(2),'value');        % 设计中是否加入涵道？
    Dd = get(SpecificationsValues(7),'value');      % 涵道直径(m)
    Cd = get(SpecificationsValues(8),'value');      % 涵道弦长(m)
    
    % 子栏“叶片切面参数”中的输入值
    Nx = get(Property1Values(1),'value');           % 切面数量
    ChordMethod = get(Property1Values(2),'value');  % 设计方法
    Meanline_flag = get(FlagsValues(3),'value');    % 凸缘线是否同一？
    Thickness_flag = get(FlagsValues(4),'value');   % 叶型是否同一？
    rx = zeros(Nx,1);                               % 径向位置
    CDp_x = zeros(Nx,1);                            % 阻力系数
    Meanline_x = cell(Nx,1);                        % 凸缘线
    Thickness_x = cell(Nx,1);                       % 叶型
    CL_x = zeros(Nx,1);                             % 升力系数
    T0oDp_x = zeros(Nx,1);                          % 厚径比
    for index = 1 : Nx
        rx(index) = get(rx_in(index),'value');
        CDp_x(index) = get(CDp_x_in(index),'value');
        Meanline_x{index} = get(Meanline_x_in(index),'value');
        Thickness_x{index} = get(Thickness_x_in(index),'value');
        CL_x(index) = get(CL_x_in(index),'value');
        T0oDp_x(index) = get(T0oDp_x_in(index),'value');
    end    
    
    % 子栏“外部环境参数”中的输入量
    Ni = get(Property2Values(1),'value');           % 坐标数量
    rho = get(Property2Values(2),'value');          % 流体密度
    ri = zeros(Ni,1);                               % 径向位置
    VA_i = zeros(Ni,1);                             % 轴向来流速度
    VT_i = zeros(Ni,1);                             % 切向来流速度
    for index = 1 : Ni
        ri(index) = get(ri_in(index),'value');
        VA_i(index) = get(VA_i_in(index),'value');
        VT_i(index) = get(VT_i_in(index),'value');
    end    
    
    % 子栏“涵道相关参数”中的输入量
    TdoT = get(DuctValues(1),'value');              % 升力比重
    CDd = get(DuctValues(2),'value');               % 阻力系数
    
    % 子栏“其他参数与工具”中的输入量
    Analyze_flag = get(FlagsValues(5),'value');     % 是否绘制性能曲线？
    Geometry_flag = get(FlagsValues(6),'value');    % 是否绘制演示模型？
    Coordinate_flag = get(FlagsValues(7),'value');  % 是否输出点坐标？
    Printing_flag = get(FlagsValues(8),'value');    % 是否输出stl文件？
    Mp = get(ToolsValues(1),'value');               % 叶片控制点数量
    Np = get(ToolsValues(2),'value');               % 叶片切向段数量
    Nd = get(ToolsValues(3),'value');               % 涵道涡流环数量
    ITER = get(ToolsValues(4),'value');             % 最大迭代次数
    
    % 将GUI界面的输入值保存至结构体数组input
    input.part1 = '螺旋桨规格';
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
    
    input.part2 = '叶片切面参数';
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
    
    input.part3 = '外部环境参数';
    input.Ni = Ni;
    input.rho = rho;
    input.ri = ri;
    input.VA_i = VA_i;
    input.VT_i = VT_i;
    
    input.part4 = '涵道相关参数';
    input.TdoT = TdoT;
    input.CDd = CDd;
    
    input.part5 = '其他参数与工具';
    input.Analyze_flag = Analyze_flag;
    input.Geometry_flag = Geometry_flag;
    input.Coordinate_flag = Coordinate_flag;
    input.Printing_flag = Printing_flag;
    input.Mp = Mp;
    input.Np = Np;
    input.Nd = Nd;
    input.ITER = ITER;
    
    % 将所有项目信息封装至结构体数组pt中
    pt.filename = filename;
    pt.date = date;
    pt.input = input;
    pt.design = [];
    pt.geometry = [];
    pt.states = [];
    
    % 保存
    save(filename,'pt');
end

function SaveCopyAs(~,~)
    global pt;
    global SpecificationsValues FlagsValues Property1Values rx_in...
           CDp_x_in Meanline_x_in Thickness_x_in CL_x_in T0oDp_x_in...
           Property2Values ri_in VA_i_in VT_i_in DuctValues ToolsValues...
           FilenameValues;
    filename = get(FilenameValues,'value');
    
    % 子栏“螺旋桨规格”中的输入值
    N = get(SpecificationsValues(1),'value');       % 电机转速(RPM)
    VS = get(SpecificationsValues(2),'value');      % 行进速度(m/s)
    T = get(SpecificationsValues(3),'value');       % 总升力(N)
    Z = get(SpecificationsValues(4),'value');       % 叶片数量
    Dp = get(SpecificationsValues(5),'value');      % 叶片直径(m)
    Hub_flag = get(FlagsValues(1),'value');         % 设计中是否加入桨毂？
    Dh = get(SpecificationsValues(6),'value');      % 桨毂直径(m)
    Duct_flag = get(FlagsValues(2),'value');        % 设计中是否加入涵道？
    Dd = get(SpecificationsValues(7),'value');      % 涵道直径(m)
    Cd = get(SpecificationsValues(8),'value');      % 涵道弦长(m)
    
    % 子栏“叶片切面参数”中的输入值
    Nx = get(Property1Values(1),'value');           % 切面数量
    ChordMethod = get(Property1Values(2),'value');  % 设计方法
    Meanline_flag = get(FlagsValues(3),'value');    % 凸缘线是否同一？
    Thickness_flag = get(FlagsValues(4),'value');   % 叶型是否同一？
    rx = zeros(Nx,1);                               % 径向位置
    CDp_x = zeros(Nx,1);                            % 阻力系数
    Meanline_x = cell(Nx,1);                        % 凸缘线
    Thickness_x = cell(Nx,1);                       % 叶型
    CL_x = zeros(Nx,1);                             % 升力系数
    T0oDp_x = zeros(Nx,1);                          % 厚径比
    for index = 1 : Nx
        rx(index) = get(rx_in(index),'value');
        CDp_x(index) = get(CDp_x_in(index),'value');
        Meanline_x{index} = get(Meanline_x_in(index),'value');
        Thickness_x{index} = get(Thickness_x_in(index),'value');
        CL_x(index) = get(CL_x_in(index),'value');
        T0oDp_x(index) = get(T0oDp_x_in(index),'value');
    end    
    
    % 子栏“外部环境参数”中的输入量
    Ni = get(Property2Values(1),'value');           % 坐标数量
    rho = get(Property2Values(2),'value');          % 流体密度
    ri = zeros(Ni,1);                               % 径向位置
    VA_i = zeros(Ni,1);                             % 轴向来流速度
    VT_i = zeros(Ni,1);                             % 切向来流速度
    for index = 1 : Ni
        ri(index) = get(ri_in(index),'value');
        VA_i(index) = get(VA_i_in(index),'value');
        VT_i(index) = get(VT_i_in(index),'value');
    end    
    
    % 子栏“涵道相关参数”中的输入量
    TdoT = get(DuctValues(1),'value');              % 升力比重
    CDd = get(DuctValues(2),'value');               % 阻力系数
    
    % 子栏“其他参数与工具”中的输入量
    Analyze_flag = get(FlagsValues(5),'value');     % 是否绘制性能曲线？
    Geometry_flag = get(FlagsValues(6),'value');    % 是否绘制演示模型？
    Coordinate_flag = get(FlagsValues(7),'value');  % 是否输出点坐标？
    Printing_flag = get(FlagsValues(8),'value');    % 是否输出stl文件？
    Mp = get(ToolsValues(1),'value');               % 叶片控制点数量
    Np = get(ToolsValues(2),'value');               % 叶片切向段数量
    Nd = get(ToolsValues(3),'value');               % 涵道涡流环数量
    ITER = get(ToolsValues(4),'value');             % 最大迭代次数
    
    % 将GUI界面的输入值保存至结构体数组input
    input.part1 = '螺旋桨规格';
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
    
    input.part2 = '叶片切面参数';
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
    
    input.part3 = '外部环境参数';
    input.Ni = Ni;
    input.rho = rho;
    input.ri = ri;
    input.VA_i = VA_i;
    input.VT_i = VT_i;
    
    input.part4 = '涵道相关参数';
    input.TdoT = TdoT;
    input.CDd = CDd;
    
    input.part5 = '其他参数与工具';
    input.Analyze_flag = Analyze_flag;
    input.Geometry_flag = Geometry_flag;
    input.Coordinate_flag = Coordinate_flag;
    input.Printing_flag = Printing_flag;
    input.Mp = Mp;
    input.Np = Np;
    input.Nd = Nd;
    input.ITER = ITER;
    
    % 将所有项目信息封装至结构体数组pt中
    pt.filename = filename;
    pt.date = date;
    pt.input = input;
    pt.design = [];
    pt.geometry = [];
    pt.states = [];
    
    % 目前导入和另存为完成后，屏幕焦点会回到matlab，这个问题需要解决
    uisave('pt',[filename,'-copy']);
end

function Execute(~,~)
    global pt;
    global NumbersValues;
    % 保存当前界面输入值
    Save;
    % 生成结果标签页
    ResultPage;
    
    % 进行计算
    pt.design = EppsOptimizer(pt.input);
    % 计算螺旋桨的扭矩，其中CQ为扭矩系数
    pt.design.Q = pt.design.CQ * 0.5 * pt.input.rho * pt.input.Vs^2 * pi*pt.input.Dp^2/4 * pt.input.Dp/2; % [Nm]  torque
    % 计算螺旋桨的转速
    omega = 2*pi*pt.input.N/60; % [rad/s]
    % 计算螺旋桨消耗的功率
    pt.design.P = pt.design.Q * omega;
    
    % 将计算结果导入数值结果的编辑框中
    NumbersValues_def = {pt.input.N,pt.input.VS,pt.input.Dp,pt.input.T,...
                         pt.design.Q,pt.design.P,pt.design.Js,pt.design.KT,...
                         pt.design.CT,pt.design.KQ,pt.design.CQ,pt.design.CP,...
                         pt.design.EFFY,pt.design.ADEFFY};
    for index = 1 : 14
        set(NumbersValues(index),'value',NumbersValues_def{index});
    end                 
    
    pt.geometry = Geometry(pt);
    
end

%% 生成结果界面
function ResultPage
    global pt;
    global zhfontsize subcolumnfontsize numfontsize;
    global Menu TabPageMenu Tabs FilenameValues;
    global NumbersValues PlotsValues PlotsPanel coordinate;
    
    global linecolor1_1 linecolor1_2 linecolor2_1 linecolor2_2...
           linecolor3_1 linecolor3_2 linecolor4_1 linecolor4_2;
    % 设置线条颜色
    linecolor1_1 = [0.047059 0.925490 0.866667];       % #0CECDD
    linecolor1_2 = [1.000000 0.952941 0.219608];       % #FFF338
    linecolor2_1 = [1.000000 0.403921 0.905882];       % #FF67E7
    linecolor2_2 = [0.768627 0.000000 1.000000];       % #C400FF
    linecolor3_1 = [0.627451 0.235294 0.470588];       % #A03C78
    linecolor3_2 = [0.929412 0.556863 0.486275];       % #ED8E7C
    linecolor4_1 = [0.803921 0.952941 0.635294];       % #CDF3A2
    linecolor4_2 = [0.576471 0.850980 0.639216];       % #93D9A3
    
    %% 生成新标签页
    filename = get(FilenameValues,'value');
    tabpagename = [filename,'的设计结果'];
    
    DesignResult = uitab('parent',Tabs,'title',tabpagename);
    % 获取当前标签页的数量
    tabpageamount = length(get(Tabs,'children'));
    % 切换至新生成的结果标签页
    set(Tabs,'selectedtab',DesignResult);
    
    % 将新生成的结果标签页添加至工具栏的“标签页”子栏下
    TabPageMenu(tabpageamount) = uimenu('parent',Menu(3),...
                                        'text',tabpagename,...
                                        'menuselectedfcn',@DeleteTab);
                                            
    %% 新标签页中的GUI元素
    % 生成DesignResult中的网格空间
    DesignResultGrid = uigridlayout(DesignResult,[1 2],...
                                    'columnwidth',{'1x','5x'});
    % 数值结果 Number Result
    Numbers = uipanel('parent',DesignResultGrid,...
                      'title','数值结果',...
                      'titleposition','centertop',...
                      'fontsize',subcolumnfontsize,...
                      'fontweight','bold');
                  
    % 图像结果 Plot Result
    Plots = uipanel('parent',DesignResultGrid,...
                    'title','图像结果',...
                    'titleposition','centertop',...
                    'fontsize',subcolumnfontsize,...
                    'fontweight','bold');
                
    %% 子栏“数值结果”中的GUI元素
    % 生成Numbers中的网格空间
    NumbersGrid = uigridlayout(Numbers,[14 2]);
    
    NumbersStrings = {'电机转速：','行进速度：','叶片直径：','总升力：',...
                      '总扭矩：','总功率：','进速系数：','推力系数(叶梢速度)：',...
                      '推力系数(行进速度)：','扭矩系数：','扭矩系数：',...
                      '功率系数：','EFFY：' 'ADEFFY：'};
    NumbersTips = {'N','VS','Dp','T','Q','P','Js','KT','CT','KQ','CQ',...
                   'CP','EFFY' 'ADEFFY'};
               
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
    
    % 设置Numbers编辑框中的显示形式和单位
    set(NumbersValues(1),'valuedisplayformat','%.1f RPM');
    set(NumbersValues(2),'valuedisplayformat','%.2f m/s');
    set(NumbersValues(3),'valuedisplayformat','%.2f m');
    set(NumbersValues(4),'valuedisplayformat','%.2f N');
    set(NumbersValues(5),'valuedisplayformat','%.2f N·m');
    set(NumbersValues(6),'valuedisplayformat','%.2f W');
    
    %% 子栏“图像结果”中的GUI元素
    % 生成Plots中的网格空间
    PlotsGrid = uigridlayout(Plots,[2 2],'rowheight',{'3.5x','9x'},...
                             'columnwidth',{'1x','3x'});
    
    PlotsStrings = {'升力系数分布曲线','厚径比分布曲线','流入速度分布曲线',...
                    '环量分布曲线','诱导速度分布曲线','相关角度的分布曲线',...
                    '伸张轮廓曲线','厚径比分布曲线','升力系数分布曲线',...
                    '性能曲线','叶切面二维轮廓','螺旋桨三维模型'};
                
    InputButtonBox = uibuttongroup('parent',PlotsGrid,...
                                   'title','来自输入',...
                                   'titleposition','centertop',...
                                   'fontsize',zhfontsize,...
                                   'fontweight','bold');
    InputButtonBoxGrid = uigridlayout(InputButtonBox,[3 1]);
                               
    PlotsPanel = uipanel('parent',PlotsGrid,...
                          'title','点击左侧按钮显示对应图像',...
                          'titleposition','centertop',...
                          'fontsize',zhfontsize,...
                          'fontweight','bold');
    PlotsPanel.Layout.Row = [1 2];
    PanelGrid = uigridlayout(PlotsPanel,[1 1]);
    % 绘制坐标系
    coordinate = uiaxes(PanelGrid);
    
    ResultButtonBox = uibuttongroup('parent',PlotsGrid,...
                                    'title','来自结果',...
                                    'titleposition','centertop',...
                                    'fontsize',zhfontsize,...
                                    'fontweight','bold');
    ResultButtonBox.Layout.Row = 2;
    ResultButtonBox.Layout.Column = 1;
    ResultButtonBoxGrid = uigridlayout(ResultButtonBox,[9 1]);
    
    for index = 1 : 3
        PlotsValues(index) = uibutton(InputButtonBoxGrid,'state',...
                                      'text',PlotsStrings{index},...
                                      'fontsize',zhfontsize);
    end
    for index = 4 : 12
        PlotsValues(index) = uibutton(ResultButtonBoxGrid,'state',...
                                      'text',PlotsStrings{index},...
                                      'fontsize',zhfontsize);
    end    
    
    % 设置输入按钮的回调函数
    set(PlotsValues(1),'valuechangedfcn',{@CLInputfcn,PlotsStrings});
    set(PlotsValues(2),'valuechangedfcn',{@ThicknessInputfcn,PlotsStrings});
    set(PlotsValues(3),'valuechangedfcn',{@InflowProfilefcn,PlotsStrings});
    
    if ~strcmp(pt.input.ChordMethod,'CLmax')
        % 设计方法不是CLmax时，升力系数和厚径比按钮不可用
        set(PlotsValues(1),'enable','off');
        set(PlotsValues(2),'enable','off');
    end    
    if pt.input.Analyze_flag == 0
        % 没有选择输出性能曲线时，性能曲线按钮不可用
        set(PlotsValues(10),'enable','off');
    end 
    if pt.input.Geometry_flag == 0
        % 没有选择输出几何图像时，二维轮廓和三维模型按钮不可用
        set(PlotsValues(11),'enable','off');
        set(PlotsValues(12),'enable','off');
    end       
    
    % 设置输出按钮的回调函数
    set(PlotsValues(4),'valuechangedfcn',{@Circulationfcn,PlotsStrings});
    set(PlotsValues(5),'valuechangedfcn',{@InducedVelocityfcn,PlotsStrings});
    set(PlotsValues(6),'valuechangedfcn',{@InflowAnglefcn,PlotsStrings});
    set(PlotsValues(7),'valuechangedfcn',{@ExpandedBladefcn,PlotsStrings});
    set(PlotsValues(8),'valuechangedfcn',{@ThicknessResultfcn,PlotsStrings});
    set(PlotsValues(9),'valuechangedfcn',{@CLResultfcn,PlotsStrings});
    set(PlotsValues(10),'valuechangedfcn',@Performancefcn);
end
    
%% 输入按钮的回调函数
% “升力系数分布曲线”的回调函数
function CLInputfcn(hObject,~,PlotsStrings)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    global CL_Input_line CL_Result_line CL_Input_point;
    global linecolor1_1 linecolor1_2;
    
    % 设置图像面板标题
    set(PlotsPanel,'title',PlotsStrings{1});
    
    % 查找被按下的按钮数量
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % 对按钮1进行操作后，若一个按钮也没有被按下，则所有按钮均可用
        set(PlotsValues,'enable','on');
        if pt.input.Analyze_flag == 0
            set(PlotsValues(10),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(11),'enable','off');
            set(PlotsValues(12),'enable','off');
        end
    else
        % 对按钮1进行操作后，若存在被按下的按钮，则除1、9外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(1),'enable','on');
        set(PlotsValues(9),'enable','on');
    end
    
    % 绘制图像
    if get(hObject,'value')
        % 避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        RhoRp = pt.input.Dh/pt.input.Dp;
        rx = pt.input.rx;
        CL_x = pt.input.CL_x;
        % rx_expanded是rx的插值序列，利用正弦函数拓展rx序列
        rx_expanded = RhoRp+(1-RhoRp)*(sin((0:60)*pi/(2*60)));
        % CL_x_expanded是CL_x基于pchip的插值序列
        CL_x_expanded = pchip(rx,CL_x,rx_expanded);
        
        % 如果此时被按下的按钮不为一个，则说明9也被按下
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        
        % 绘制升力系数曲线
        CL_Input_line = plot(coordinate,rx_expanded,CL_x_expanded,...
                             'linewidth',2,...
                             'color',linecolor1_1);
        hold(coordinate,'on');
        grid(coordinate,'on');
        box(coordinate,'on');
        
        % 绘制升力系数散点
        CL_Input_point = scatter(coordinate,rx,CL_x,...
                                 'marker','o',...
                                 'markeredgecolor',linecolor1_1,...
                                 'markerfacecolor',linecolor1_2,...
                                 'sizedata',25);
                             
        % 设置坐标系自动调整                         
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        
        if pushednum == 1
            legend(CL_Input_line,'CL');
        else
            legend([CL_Input_line,CL_Result_line],'Input:CL','Design:CL');
            set(PlotsPanel,'title','升力系数分布曲线(来自输入/来自结果)');
        end    
        
        % 设置坐标轴标签
        xlabel(coordinate,'rx','fontsize',16,'fontname','Times');
        ylabel(coordinate,'CL','fontsize',16,'fontname','Times');
        
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanel,'title','点击左侧按钮显示对应图像');
    else
        delete(CL_Input_line);
        delete(CL_Input_point);
        legend(CL_Result_line,'CL');
        set(PlotsPanel,'title',PlotsStrings{9});
    end
end

% “厚径比分布曲线”的回调函数
function ThicknessInputfcn(hObject,~,PlotsStrings)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    global T0oDp_Input_line T0oDp_Result_line T0oDp_Input_point;
    global linecolor2_1 linecolor2_2;
    
    % 设置图像面板标题
    set(PlotsPanel,'title',PlotsStrings{2});
    
    % 查找被按下的按钮数量
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % 对按钮2进行操作后，若一个按钮也没有被按下，则所有按钮均可用
        set(PlotsValues,'enable','on');
        if pt.input.Analyze_flag == 0
            set(PlotsValues(10),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(11),'enable','off');
            set(PlotsValues(12),'enable','off');
        end
    else
        % 对按钮2进行操作后，若存在被按下的按钮，则除2、8外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(2),'enable','on');
        set(PlotsValues(8),'enable','on');
    end
    
    % 绘制图像
    if get(hObject,'value')
        % 避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        RhoRp = pt.input.Dh/pt.input.Dp;
        rx = pt.input.rx;
        T0oDp_x = pt.input.T0oDp_x;
        % rx_expanded是rx的插值序列，利用正弦函数拓展rx序列
        rx_expanded = RhoRp+(1-RhoRp)*(sin((0:60)*pi/(2*60)));
        % T0oDp_x_expanded是T0oDp_x基于pchip的插值序列
        T0oDp_x_expanded = pchip(rx,T0oDp_x,rx_expanded);
        
        % 如果此时被按下的按钮不为一个，则说明8也被按下
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        
        % 绘制厚径比曲线
        T0oDp_Input_line = plot(coordinate,rx_expanded,T0oDp_x_expanded,...
                                'linewidth',2,...
                                'color',linecolor2_1);
        hold(coordinate,'on');
        grid(coordinate,'on');
        box(coordinate,'on');
        
        % 绘制厚径比散点
        T0oDp_Input_point = scatter(coordinate,rx,T0oDp_x,...
                                     'marker','o',...
                                     'markeredgecolor',linecolor2_1,...
                                     'markerfacecolor',linecolor2_2,...
                                     'sizedata',25);
                                 
        % 设置坐标系自动调整
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        
        if pushednum == 1
            legend(T0oDp_Input_line,'T0/Dp');
        else
            legend([T0oDp_Input_line,T0oDp_Result_line],'Input:T0/Dp','Design:T0/Dp');
            set(PlotsPanel,'title','厚径比分布曲线(来自输入/来自结果)');
        end    
        
        % 设置坐标轴标签
        xlabel(coordinate,'rx','fontsize',16,'fontname','Times');
        ylabel(coordinate,'T0/Dp','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanel,'title','点击左侧按钮显示对应图像');
    else
        delete(T0oDp_Input_line);
        delete(T0oDp_Input_point);
        legend(T0oDp_Result_line,'t0/Dp');
        set(PlotsPanel,'title',PlotsStrings{8});
    end
end

% “流入速度分布曲线”的回调函数
function InflowProfilefcn(hObject,~,PlotsStrings)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    global VA_line VT_line UASTAR_line UTSTAR_line;
    global linecolor1_1 linecolor2_1;
    
    % 设置图像面板标题
    set(PlotsPanel,'title',PlotsStrings{3});
    
    % 查找被按下的按钮数量
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % 对按钮3进行操作后，若一个按钮也没有被按下，则所有按钮均可用
        set(PlotsValues,'enable','on');
        if ~strcmp(pt.input.ChordMethod,'CLmax')
            set(PlotsValues(1),'enable','off');
            set(PlotsValues(2),'enable','off');
        end    
        if pt.input.Analyze_flag == 0
            set(PlotsValues(10),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(11),'enable','off');
            set(PlotsValues(12),'enable','off');
        end
    else
        % 对按钮3进行操作后，若存在被按下的按钮，则除3、5外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(3),'enable','on');
        set(PlotsValues(5),'enable','on');
    end
    
    % 绘制图像
    if get(hObject,'value')
        % 避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        ri = pt.input.ri;
        VA_i = pt.input.VA_i;
        VT_i = pt.input.VT_i;
        % 如果此时被按下的按钮不为一个，则说明5也被按下
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        
        % 绘制轴向来流速度曲线
        VA_line = plot(coordinate,ri,VA_i,...
                       'linewidth',2,...
                       'color',linecolor1_1);
        
        hold(coordinate,'on');
        grid(coordinate,'on');
        box(coordinate,'on');
        
        % 绘制切向来流速度曲线
        VT_line = plot(coordinate,ri,VT_i,...
                           'linewidth',2,...
                           'color',linecolor2_1);
                       
        % 设置坐标系自动调整
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        
        if pushednum == 1
            legend([VA_line,VT_line],'VA/Vs','VT/Vs');
        else
            legend([VA_line,VT_line,UASTAR_line,UTSTAR_line],...
                   'VA/Vs','VT/Vs','Ua*/Vs','Ut*/Vs');
            set(PlotsPanel,'title','速度分布曲线');
        end    
        
        % 设置坐标轴标签
        xlabel(coordinate,'ri','fontsize',16,'fontname','Times');
        ylabel(coordinate,'V/VS','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanel,'title','点击左侧按钮显示对应图像');
    else
        delete(VA_line);
        delete(VT_line);
        legend([UASTAR_line,UTSTAR_line],'Ua*/Vs','Ut*/Vs');
        set(PlotsPanel,'title',PlotsStrings{5});
    end    
end

%% 结果按钮的回调函数
% “环量分布曲线”的回调函数
function Circulationfcn(hObject,~,PlotsStrings)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    global Gamma_line Gamma_point;
    global linecolor3_1 linecolor3_2;
    
    % 设置图像面板标题
    set(PlotsPanel,'title',PlotsStrings{4});
    
    % 查找被按下的按钮数量
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % 对按钮4进行操作后，若一个按钮也没有被按下，则所有按钮均可用
        set(PlotsValues,'enable','on');
        if ~strcmp(pt.input.ChordMethod,'CLmax')
            set(PlotsValues(1),'enable','off');
            set(PlotsValues(2),'enable','off');
        end    
        if pt.input.Analyze_flag == 0
            set(PlotsValues(10),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(11),'enable','off');
            set(PlotsValues(12),'enable','off');
        end
    else
        % 对按钮4进行操作后，若存在被按下的按钮，则除4外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(4),'enable','on');
    end
    
    % 绘制曲线
    if get(hObject,'value')
        % 避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        RhoRp = pt.input.Dh/pt.input.Dp;
        rc = pt.design.RC;
        Gamma_c = pt.design.G';
        % rc_expanded是rc的插值序列，外拓了RhoRp和1两个边界值
        rc_expanded = [RhoRp,rc,1];
        % Gamma_c_expanded是Gamma_c基于pchip的插值序列
        Gamma_c_expanded = pchip(rc,Gamma_c,rc_expanded);
        % 清空图像
        cla(coordinate);
        
        % 绘制环量曲线
        Gamma_line = plot(coordinate,rc_expanded,Gamma_c_expanded,...
                          'linewidth',2,...
                          'color',linecolor3_1);
                      
        hold(coordinate,'on');
        grid(coordinate,'on');
        box(coordinate,'on');
        
        % 绘制环量散点
        Gamma_point = scatter(coordinate,rc,Gamma_c,...
                              'marker','o',...
                              'markeredgecolor',linecolor3_2,...
                              'markerfacecolor',linecolor3_1,...
                              'sizedata',25);
                          
        % 设置坐标系自动调整
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        
        % 设置坐标轴标签
        legend(Gamma_line,'τ/2π·Rp·VS');
        xlabel(coordinate,'rc','fontsize',16,'fontname','Times');
        ylabel(coordinate,'τ/2π·Rp·Vs','fontsize',16,'fontname','Times');
    else
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');
        legend(coordinate,'off');
        set(PlotsPanel,'title','点击左侧按钮显示对应图像');
    end
end

% “诱导速度分布曲线”的回调函数
function InducedVelocityfcn(hObject,~,PlotsStrings)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    global UASTAR_line UTSTAR_line VA_line VT_line UASTAR_point UTSTAR_point;
    global linecolor3_1 linecolor3_2 linecolor4_1 linecolor4_2;
     
    % 设置图像面板标题
    set(PlotsPanel,'title',PlotsStrings{5});
    
    %查找被按下的按钮数量
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % 对按钮5进行操作后，若一个按钮也没有被按下，则所有按钮均可用
        set(PlotsValues,'enable','on');
        if ~strcmp(pt.input.ChordMethod,'CLmax')
            set(PlotsValues(1),'enable','off');
            set(PlotsValues(2),'enable','off');
        end    
        if pt.input.Analyze_flag == 0
            set(PlotsValues(10),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(11),'enable','off');
            set(PlotsValues(12),'enable','off');
        end
    else
        % 对按钮5进行操作后，若存在被按下的按钮，则除3、5外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(3),'enable','on');
        set(PlotsValues(5),'enable','on');
    end
    
    % 绘制曲线
    if get(hObject,'value')
        % 避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        rc = pt.design.RC;
        UASTAR = pt.design.UASTAR;
        UTSTAR = pt.design.UTSTAR;
        % 如果此时被按下的按钮不为一个，则说明3也被按下
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        
        grid(coordinate,'on');
        box(coordinate,'on');
        % 绘制轴向诱导速度曲线
        UASTAR_line = plot(coordinate,rc,UASTAR,...
                           'linewidth',2,...
                           'color',linecolor3_1);
        hold(coordinate,'on');
        
        % 绘制轴向诱导速度散点
        UASTAR_point = scatter(coordinate,rc,UASTAR,...
                               'marker','o',...
                               'markeredgecolor',linecolor3_2,...
                               'markerfacecolor',linecolor3_1,...
                               'sizedata',25);
        hold(coordinate,'on');
        
        % 绘制切向诱导速度曲线
        UTSTAR_line = plot(coordinate,rc,UTSTAR,...
                           'linewidth',2,...
                           'color',linecolor4_1);
        hold(coordinate,'on');
        
        % 绘制切向诱导速度散点
        UTSTAR_point = scatter(coordinate,rc,UTSTAR,...
                               'marker','o',...
                               'markeredgecolor',linecolor4_2,...
                               'markerfacecolor',linecolor4_1,...
                               'sizedata',25);
                           
        % 设置坐标系自动调整
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        
        if pushednum == 1
            legend([UASTAR_line,UTSTAR_line],'Ua*/Vs','Ut*/Vs');
        else
            legend([VA_line,VT_line,UASTAR_line,UTSTAR_line],...
                   'VA/Vs','VT/Vs','Ua*/Vs','Ut*/Vs');
            set(PlotsPanel,'title','速度分布曲线');
        end    
        
        % 设置坐标轴标签
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'Velocity','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanel,'title','点击左侧按钮显示对应图像');
    else
        delete(UASTAR_line);
        delete(UTSTAR_line);
        delete(UASTAR_point);
        delete(UTSTAR_point);
        legend([VA_line,VT_line],'VA/Vs','VT/Vs');
        set(PlotsPanel,'title',PlotsStrings{3});
    end
end

% “相关角度的分布曲线”的回调函数
function InflowAnglefcn(hObject,~,PlotsStrings)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    global Beta_line BetaI_line Beta_point BetaI_point;
    global linecolor2_1 linecolor2_2 linecolor4_1 linecolor4_2;
    
    % 设置图像面板标题
    set(PlotsPanel,'title',PlotsStrings{6});
    
    % 查找被按下的按钮数量
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % 对按钮6进行操作后，若一个按钮也没有被按下，则所有按钮均可用
        set(PlotsValues,'enable','on');
        if ~strcmp(pt.input.ChordMethod,'CLmax')
            set(PlotsValues(1),'enable','off');
            set(PlotsValues(2),'enable','off');
        end    
        if pt.input.Analyze_flag == 0
            set(PlotsValues(10),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(11),'enable','off');
            set(PlotsValues(12),'enable','off');
        end
    else
        % 对按钮6进行操作后，若存在被按下的按钮，则除6外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(6),'enable','on');
    end
    
    % 绘制图像
    if get(hObject,'value')
        % 避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        RhoRp = pt.input.Dh/pt.input.Dp;
        rc = pt.design.RC;
        Beta_c = atand(pt.design.TANBC);
        BetaI_c = atand(pt.design.TANBIC);
        % rc_expanded是rc的插值序列，外拓了RhoRp和1两个边界值
        rc_expanded = [RhoRp,rc,1];
        % Beta_c_expanded和BetaI_c_expanded是Beta_c和BetaI_c基于spline的插值序列
        Beta_c_expanded = spline(rc,Beta_c,rc_expanded);
        BetaI_c_expanded = spline(rc,BetaI_c,rc_expanded);
        
        cla(coordinate);
        grid(coordinate,'on');
        box(coordinate,'on');
        
        % 绘制进角曲线
        Beta_line = plot(coordinate,rc_expanded,Beta_c_expanded,...
                            'linewidth',2,...
                            'color',linecolor2_1);
        hold(coordinate,'on');
        
        % 设置进角散点
        Beta_point = scatter(coordinate,rc,Beta_c,...
                             'marker','o',...
                             'markeredgecolor',linecolor2_1,...
                             'markerfacecolor',linecolor2_2,...
                             'sizedata',25);
        hold(coordinate,'on');
        
        % 设置水动力螺旋角曲线
        BetaI_line = plot(coordinate,rc_expanded,BetaI_c_expanded,...
                          'linewidth',2,...
                          'color',linecolor4_1);
        hold(coordinate,'on');
        
        % 设置水动力螺旋角散点
        BetaI_point = scatter(coordinate,rc,BetaI_c,...
                              'marker','o',...
                              'markeredgecolor',linecolor4_1,...
                              'markerfacecolor',linecolor4_2,...
                              'sizedata',25);
                          
        % 设置坐标系自动调整
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        
        % 设置坐标轴标签
        legend([Beta_line,BetaI_line],'β','βi');
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'Angle','fontsize',16,'fontname','Times');
    else
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');
        legend(coordinate,'off');
        set(PlotsPanel,'title','点击左侧按钮显示对应图像');
    end
end

% “伸张轮廓曲线”的回调函数
function ExpandedBladefcn(hObject,~,PlotsStrings)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    global CpoDp_line CpoDp_point;
    global linecolor3_1 linecolor3_2;
    
    % 设置图像面板标题
    set(PlotsPanel,'title',PlotsStrings{7});
    
    % 查找被按下的按钮数量
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % 对按钮7进行操作后，若一个按钮也没有被按下，则所有按钮均可用
        set(PlotsValues,'enable','on');
        if ~strcmp(pt.input.ChordMethod,'CLmax')
            set(PlotsValues(1),'enable','off');
            set(PlotsValues(2),'enable','off');
        end    
        if pt.input.Analyze_flag == 0
            set(PlotsValues(10),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(11),'enable','off');
            set(PlotsValues(12),'enable','off');
        end
    else
        % 对按钮7进行操作后，若存在被按下的按钮，则除7外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(7),'enable','on');
    end
    
    % 绘制图像
    if get(hObject,'value')
        % 避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        RhoRp = pt.input.Dh/pt.input.Dp;
        rc = pt.design.RC;
        CpoDp_c = pt.design.CpoDp;
        rc_expanded = RhoRp+(1-RhoRp)*(sin((0:60)*pi/(2*60)));
        CpoDp_c_expanded = InterpolateChord(rc,CpoDp_c,rc_expanded);
        % 清楚图像
        cla(coordinate);
        
        grid(coordinate,'on');
        box(coordinate,'on');
        % 绘制伸张轮廓曲线
        CpoDp_line = plot(coordinate,rc_expanded,CpoDp_c_expanded,...
                          'linewidth',2,...
                          'color',linecolor3_1);
        hold(coordinate,'on');
        
        % 绘制伸张轮廓散点
        CpoDp_point = scatter(coordinate,rc,CpoDp_c,...
                              'marker','o',...
                              'markeredgecolor',linecolor3_1,...
                              'markerfacecolor',linecolor3_2,...
                              'sizedata',25);
                          
        % 设置坐标系自动调整
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        
        % 设置坐标轴标签
        legend(CpoDp_line,'Cp/Dp');
        xlabel(coordinate,'rc','fontsize',16,'fontname','Times');
        ylabel(coordinate,'Cp/Dp','fontsize',16,'fontname','Times');
    else
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanel,'title','点击左侧按钮显示对应图像');
    end
end

% “厚径比分布曲线”的回调函数
function ThicknessResultfcn(hObject,~,PlotsStrings)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    global T0oDp_Input_line T0oDp_Result_line T0oDp_Result_point;
    global linecolor4_1 linecolor4_2;
    
    % 设置图像面板标题
    set(PlotsPanel,'title',PlotsStrings{8});
    
    % 查找被按下的按钮数量
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % 对按钮8进行操作后，若一个按钮也没有被按下，则所有按钮均可用
        set(PlotsValues,'enable','on');
        if ~strcmp(pt.input.ChordMethod,'CLmax')
            set(PlotsValues(1),'enable','off');
            set(PlotsValues(2),'enable','off');
        end  
        if pt.input.Analyze_flag == 0
            set(PlotsValues(10),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(11),'enable','off');
            set(PlotsValues(12),'enable','off');
        end
    else
        % 对按钮8进行操作后，若存在被按下的按钮，则除2、8外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(2),'enable','on');
        set(PlotsValues(8),'enable','on');
    end
    % 绘制图像
    if get(hObject,'value')
        % 避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        RhoRp = pt.input.Dh/pt.input.Dp;
        rc = pt.design.RC;
        T0oDp_c = pt.design.t0oDp;
        % rc_expanded是rc的插值序列，利用正弦函数拓展rx序列
        rc_expanded = RhoRp+(1-RhoRp)*(sin((0:60)*pi/(2*60)));
        % T0oDp_c_expanded是T0oDp_c基于pchip的插值序列
        T0oDp_c_expanded = pchip(rc,T0oDp_c,rc_expanded);
        % 如果此时被按下的按钮不为1，则说明2也被按下
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        
        grid(coordinate,'on');
        box(coordinate,'on');
        % 绘制厚径比曲线
        T0oDp_Result_line = plot(coordinate,rc_expanded,T0oDp_c_expanded,...
                                 'linewidth',2,...
                                 'color',linecolor4_1);
        hold(coordinate,'on');
        
        % 绘制厚径比散点
        T0oDp_Result_point = scatter(coordinate,rc,T0oDp_c,...
                                     'marker','o',...
                                     'markeredgecolor',linecolor4_1,...
                                     'markerfacecolor',linecolor4_2,...
                                     'sizedata',25);
                                 
        % 设置坐标系自动调整
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        
        if pushednum == 1
            legend(T0oDp_Result_line,'t0/Dp');
        else
            legend([T0oDp_Input_line,T0oDp_Result_line],'Input:T0/Dp','Design:T0/Dp');
            set(PlotsPanel,'title','厚径比分布曲线(来自输入/来自结果)');
        end    
        
        % 设置坐标轴标签
        xlabel(coordinate,'rc','fontsize',16,'fontname','Times');
        ylabel(coordinate,'T0/Dp','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanel,'title','点击左侧按钮显示对应图像');
    else
        delete(T0oDp_Result_line);
        delete(T0oDp_Result_point);
        legend(T0oDp_Input_line,'Cp/Dp');
        set(PlotsPanel,'title',PlotsStrings{2});
    end
end

% “升力系数分布曲线”的回调函数
function CLResultfcn(hObject,~,PlotsStrings)
    global pt;
    global Fig_Main PlotsValues PlotsPanel coordinate;
    global CL_Input_line CL_Result_line CL_Result_point;
    global linecolor1_1 linecolor1_2;
    
    % 设置图像面板标题
    set(PlotsPanel,'title',PlotsStrings{9});
    
    % 查找被按下的按钮数量
    [~,pushednum] = PushedFind;
    if pushednum == 0
        % 对按钮9进行操作后，若一个按钮也没有被按下，则所有按钮均可用
        set(PlotsValues,'enable','on');
        if ~strcmp(pt.input.ChordMethod,'CLmax')
            set(PlotsValues(1),'enable','off');
            set(PlotsValues(2),'enable','off');
        end  
        if pt.input.Analyze_flag == 0
            set(PlotsValues(10),'enable','off');
        end 
        if pt.input.Geometry_flag == 0
            set(PlotsValues(11),'enable','off');
            set(PlotsValues(12),'enable','off');
        end
    else
        % 对按钮9进行操作后，若存在被按下的按钮，则除1、9外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(1),'enable','on');
        set(PlotsValues(9),'enable','on');
    end
    
    % 绘制图像
    if get(hObject,'value')
        % 避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        RhoRp = pt.input.Dh/pt.input.Dp;
        rc = pt.design.RC;
        CL_c = pt.design.CL;
        rc_expanded = [RhoRp,rc,1];
        CL_c_expanded = spline(rc,CL_c,rc_expanded);
        
        % 如果此时被按下的按钮不为一个，则说明1也被按下
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        
        % 绘制升力系数曲线
        CL_Result_line = plot(coordinate,rc_expanded,CL_c_expanded,...
                              'linewidth',2,...
                              'color',linecolor1_1);
        hold(coordinate,'on');
        grid(coordinate,'on');
        box(coordinate,'on');
        
        % 绘制升力系数散点
        CL_Result_point = scatter(coordinate,rc,CL_c,...
                                  'marker','o',...
                                  'markeredgecolor',linecolor1_1,...
                                  'markerfacecolor',linecolor1_2,...
                                  'sizedata',25);
                              
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        
        xlabel(coordinate,'rc','fontsize',16,'fontname','Times');
        ylabel(coordinate,'CL','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');
        legend(coordinate,'off');
        set(PlotsPanel,'title','点击左侧按钮显示对应图像');
    else
        delete(CL_Result_line);
        delete(CL_Result_point);
        legend(CL_Input_line,'CL');
        set(PlotsPanel,'title',PlotsStrings{1});
    end
end

% “性能曲线”的回调函数
function Performancefcn(~,~)
    global pt;
    
    if pt.input.Analyze_flag
        pt.states = AnalyzeAuto(pt);
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotsPanel(15));
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
end

%% 菜单栏的回调函数
% 关闭标签页
function DeleteTab(hObject,~)
    global Fig_Main Tabs;
    tabpage = get(Tabs,'children');
    tabpageamount = length(tabpage);
    tabpagename = cell(tabpageamount,1);
    for index = 1 : tabpageamount
        tabpagename{index} = get(tabpage(index),'title');
    end 
    menutabpagename = get(hObject,'text');
    disp(tabpagename{3});
    disp(menutabpagename);
    if find(ismember(tabpagename,menutabpagename))
        % “标签页”菜单下有与标签页相对应的项
        confirmation = sprintf('确定要删除当前标签页吗？');
        DeleteSelection = uiconfirm(Fig_Main,confirmation,'删除确认',...
                                    'options',{'确定','取消'});
        if strcmp(DeleteSelection,'确定')
            index = ismember(tabpagename,menutabpagename);
            delete(tabpage(index));
            % 同时删除自身
            delete(hObject);
            % 切换至“仅进行叶片设计”
            set(Tabs,'selectedtab',tabpage(1));
        else
            return
        end    
    else
        error = sprintf('标签页不匹配，请检查！');
        uialert(Fig_Main,error,'删除失败');
    end    
end
