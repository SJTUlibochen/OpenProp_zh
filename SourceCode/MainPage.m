% 该程序负责OpenProp中文版主界面GUI的生成
%% 主界面
function MainPage
    %% 前处理
    clear variables;
    clear global;
    global titlefontsize zhfontsize subcolumnfontsize numfontsize;
    global Fig_Main;
    global SpecificationsValues FlagsValues Property1Values Table1PanelGrid...
           rx_in CDp_x_in Meanline_x_in Thickness_x_in CL_x_in T0oDp_x_in...
           Property2Values Table2PanelGrid ri_in VA_i_in VT_i_in DuctValues...
           ToolsValues FilenameValues
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
    FileStrings = {'保存' '另存为' '载入' '运行'};
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
    set(LoadButton,'buttonpushedfcn',{@Load,OpenPropDirectory});
    set(SaveButton,'buttonpushedfcn',{@Save,OpenPropDirectory});
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
function Load(~,~,OpenPropDirectory)
    global pt numfontsize;
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
        % 导入数据的切面数量多余原GUI界面上显示的
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

function Save(~,~,OpenPropDirectory)
    global Fig_Main pt;
    global SpecificationsValues FlagsValues Property1Values rx_in...
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

function Execute(hObject,ED)
    global PlotsPanels NumbersValues;
    global pt
    %% 保存界面数据
    SaveData;
    %% 生成新界面
    newPlots;
    %% 进行数据转换
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
    %% 通过EppsOptimizer函数得到pt.design的值
    pt.design = EppsOptimizer(pt.input);
    %计算螺旋桨的扭矩，其中CQ为扭矩系数
    pt.design.Q = pt.design.CQ * 0.5 * pt.input.rho * pt.input.Vs^2 * pi*pt.input.Dp^2/4 * pt.input.Dp/2; % [Nm]  torque
    %计算螺旋桨的转速
    omega = 2*pi*pt.input.N/60; % [rad/s]
    %计算螺旋桨消耗的功率
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
    %% 绘制结果报告
    Make_Reports;
    %% 绘制图像
    pt.geometry = Geometry(pt);
    %% 决定是否输出性能曲线
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

    %% 以下是关于参数研究的部分
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

%% 运行OpenProp后新界面几何元素的生成代码
function newPlots
    global pt;
    global Tabs TabPageStrings TabPageMenu Menu;
    global Filename NumbersValues PlotsValues PlotsStrings PlotsPanels coordinate;
    global zhfontsize subcolumnfontsize numfontsize;
    global lc1_1 lc1_2 lc2_1 lc2_2 lc3_1 lc3_2 lc4_1 lc4_2;
    %姜红[0.933333 0.721568 0.764705]；隐红灰[0.709803 0.596078 0.631372]
    %枇杷黄[0.988235 0.631372 0.023529]；蛋白石绿[0.341176 0.584313 0.447058]
    %星蓝[0.576470 0.709803 0.811764]；山梗紫[0.380392 0.392156 0.623529]
    lc1_1 = [0.047059 0.925490 0.866667];       %#0CECDD
    lc1_2 = [1.000000 0.952941 0.219608];       %#FFF338
    lc2_1 = [1.000000 0.403921 0.905882];       %#FF67E7
    lc2_2 = [0.768627 0.000000 1.000000];       %#C400FF
    lc3_1 = [0.627451 0.235294 0.470588];       %#A03C78
    lc3_2 = [0.929412 0.556863 0.486275];       %#ED8E7C
    lc4_1 = [0.803921 0.952941 0.635294];       %#CDF3A2
    lc4_2 = [0.576471 0.850980 0.639216];       %#93D9A3
    %% 生成主界面上的一个新标签页
    filename = get(Filename,'value');
    DesignResult = uitab('parent',Tabs,...
                         'title',[filename,'的设计结果']);
    %切换至新生成的结果标签页
    set(Tabs,'selectedtab',DesignResult);
    %将新生成的结果标签页添加至工具栏的“标签页”子栏下
    TabPageStrings{length(TabPageStrings)+1} = [filename,'的设计结果'];
    TabPageMenu(length(TabPageMenu)+1) = uimenu('parent',Menu(2),...
                                                'text',TabPageStrings{end});
    %% 生成DesignResult中的网格空间及空间中的各栏目
    %生成DesignResult中1*2的网格空间
    DesignResultGrid = uigridlayout(DesignResult,[1 2],...
                                    'columnwidth',{'2x','7x'});
    %生成“数值结果”栏
    Numbers = uipanel('parent',DesignResultGrid,...
                      'title','数值结果',...
                      'titleposition','centertop',...
                      'fontsize',subcolumnfontsize,...
                      'fontweight','bold');
    %生成“图像结果”栏
    Plots = uipanel('parent',DesignResultGrid,...
                    'title','图像结果',...
                    'titleposition','centertop',...
                    'fontsize',subcolumnfontsize,...
                    'fontweight','bold');
    %% “数值结果”子栏中的内容
    %生成Numbers中14*2的网格空间
    NumbersGrid = uigridlayout(Numbers,[14 2]);
    %Numbers中的参数
    NumbersStrings = {'装置行进速度：' '螺旋桨转速：' '螺旋桨直径：' '装置总推力：' '装置总扭矩：' '装置总功率：'...
                      '进速系数：' '推力系数(叶梢速度)：' '推力系数(行进速度)：' '扭矩系数：' '扭矩系数：' '功率系数：'...
                      'EFFY：' 'ADEFFY：'};
    NumbersTips = {'Vs' 'N' 'Dp' 'Thrust' 'Torque' 'Power'...
                   'Js' 'KT' 'CT' 'KQ' 'CQ' 'CP'...
                   'EFFY' 'ADEFFY'};
    %利用循环绘制14行
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
    %设置14行的编辑框中的显示形式和单位
    set(NumbersValues(1),'valuedisplayformat','%.2f m/s');
    set(NumbersValues(2),'valuedisplayformat','%.1f RPM');
    set(NumbersValues(3),'valuedisplayformat','%.2f m');
    set(NumbersValues(4),'valuedisplayformat','%.2f N');
    set(NumbersValues(5),'valuedisplayformat','%.2f N·m');
    set(NumbersValues(6),'valuedisplayformat','%.2f W');
    %% “图像结果”子栏中的内容
    %生成Plots中1*2的网格空间
    PlotsGrid = uigridlayout(Plots,[1 2],...
                             'columnwidth',{'1x','3x'});
    %选框和图像中需要使用的参数
    PlotsStrings = {'叶片伸张轮廓曲线(来自输入)'...     %原版'Expanded Blade (input blade)'
                    '叶片厚度分布曲线(来自输入)'...     %原版'Blade Thickness (input blade)'
                    '流入速度分布曲线'...       %原版'Inflow Profile'
                    '效率-螺旋桨直径'...       %原版'Efficiency vs Diameter'
                    '效率-螺旋桨转速'...       %原版'Efficiency vs Rotation Speed'
                    '环量分布曲线'...     %原版'Circulation Distribution'
                    '诱导速度分布曲线'...       %原版'Induced Velocity'
                    '相关角度的分布曲线'...      %原版'Inflow Angle'
                    '叶片伸张轮廓曲线(来自设计结果)'...       %原版'Expanded Blade (as designed)'
                    '叶片厚度分布曲线(来自设计结果)'...       %原版'Blade Thickness (as designed)'
                    '升力系数分布曲线'...       %原版'Lift Coefficient'
                    '性能曲线'...      %原版'Performance Curves'
                    '叶切面二维轮廓'...        %原版'2D Geometry'
                    '叶片三维模型'};      %原版'3D Geometry'
    %绘制14行单选按钮
    PlotsButtonBox = uibuttongroup('parent',PlotsGrid,...
                                   'title','选择显示的图像',...
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
    %设置各按钮回调函数
    set(PlotsValues(1),'valuechangedfcn',@ExpandedBladeInfcn);
    set(PlotsValues(2),'valuechangedfcn',@BladeThicknessInfcn);
    set(PlotsValues(3),'valuechangedfcn',@InflowProfilefcn);
    %“效率-螺旋桨直径”和“效率-螺旋桨转速”两曲线默认关闭
    set(PlotsValues(4),'enable','off');
    set(PlotsValues(5),'enable','off');
    if pt.input.Chord_flag ~= 0
        %选择弦长优化时，叶片伸张轮廓曲线(来自输入)不可用
        set(PlotsValues(1),'enable','off');
    end
    if pt.input.Analyze_flag == 0
        %没有选择输出性能曲线时，性能曲线不可用
        set(PlotsValues(12),'enable','off');
    end 
    if pt.input.Geometry_flag == 0
        %没有选择输出几何图像时，二维轮廓和三维模型不可用
        set(PlotsValues(13),'enable','off');
        set(PlotsValues(14),'enable','off');
    end       
    %绘制图像面板
    PlotsPanels = uipanel('parent',PlotsGrid,...
                          'title','点击左侧按钮显示对应图像',...
                          'titleposition','centertop',...
                          'fontsize',zhfontsize,...
                          'fontweight','bold');
    PanelGrid = uigridlayout(PlotsPanels,[1 1]);
    %绘制坐标系
    coordinate = uiaxes(PanelGrid);
end
%% 按下按钮后的回调函数
%按下“叶片伸张轮廓曲线(来自输入)”后的回调函数
function ExpandedBladeInfcn(hObject,ED)
    global PlotsPanels PlotsStrings
    global PlotsValues
    global pt;
    global Fig_Main coordinate Xr_XCpoDp_In_line Xr_XCpoDp_De_line Xr_XCpoDp_In_point;
    global lc1_1 lc1_2;
    %% 设置图像面板标题
    set(PlotsPanels,'title',PlotsStrings{1});
    %% 设置相应按钮的可用状态
    %查找被按下的按钮数量
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %对按钮1进行操作后，若一个按钮也没有被按下
        %则除4、5外的所有按钮均可用
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
        %对按钮1进行操作后，若存在被按下的按钮
        %则除1、9外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(1),'enable','on');
        set(PlotsValues(9),'enable','on');
    end
    %% 绘制曲线
    %按下“叶片伸张轮廓曲线(来自输入)”后的回调函数
    if get(PlotsValues(1),'value')
        %避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.input.Xr;
        XCpoDp = pt.input.XCpoDp;
        %XXr是Xr的一个插值序列，利用正弦函数拓展Xr序列
        XXr = Rhub_oR+(1-Rhub_oR)*(sin((0:60)*pi/(2*60)));
        %XXCpoDp为通过插值得到的XXr径向位置处的弦径比，该插值方法非matlab自带
        XXCpoDp = InterpolateChord(Xr,XCpoDp,XXr);
        %如果此时被按下的按钮不为1，则说明9也被按下
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
            set(PlotsPanels,'title','叶片伸张轮廓曲线(来自输入/来自设计结果)');
        end    
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'Cp/Dp','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanels,'title','点击左侧按钮显示对应图像');
    else
        delete(Xr_XCpoDp_In_line);
        delete(Xr_XCpoDp_In_point);
        legend(Xr_XCpoDp_De_line,'Cp/Dp');
        set(PlotsPanels,'title',PlotsStrings{9});
    end
end
%按下“叶片厚度分布曲线(来自输入)”后的回调函数
function BladeThicknessInfcn(hObject,ED)
    global PlotsPanels PlotsStrings
    global PlotsValues
    global pt;
    global Fig_Main coordinate Xr_Xt0oDp_In_line Xr_Xt0oDp_De_line Xr_Xt0oDp_In_point;
    global lc2_1 lc2_2;
    %% 设置图像面板标题
    set(PlotsPanels,'title',PlotsStrings{2});
    %% 设置相应按钮的可用状态
    %查找被按下的按钮数量
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %对按钮2进行操作后，若一个按钮也没有被按下
        %则除4、5外的所有按钮均可用
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
        %对按钮2进行操作后，若存在被按下的按钮
        %则除2、10外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(2),'enable','on');
        set(PlotsValues(10),'enable','on');
    end
    %% 绘制曲线
    if get(PlotsValues(2),'value')
        %避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.input.Xr;
        Xt0oDp = pt.input.Xt0oDp;
        %%XXr是Xr的一个插值序列，利用正弦函数拓展Xr序列
        XXr = Rhub_oR+(1-Rhub_oR)*(sin((0:60)*pi/(2*60)));
        %XXt0oDp为通过插值得到的XXr径向位置处的厚径比，该插值方法为pchip
        XXt0oDp = pchip(Xr,Xt0oDp,XXr);
        %如果此时被按下的按钮不为1，则说明10也被按下
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
            set(PlotsPanels,'title','叶片厚度分布曲线(来自输入/来自设计结果)');
        end    
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'t0/Dp','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanels,'title','点击左侧按钮显示对应图像');
    else
        delete(Xr_Xt0oDp_In_line);
        delete(Xr_Xt0oDp_In_point);
        legend(Xr_Xt0oDp_De_line,'t0/Dp');
        set(PlotsPanels,'title',PlotsStrings{10});
    end
end
%按下“流入速度分布曲线”后的回调函数
function InflowProfilefcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate Xr_VAI_line Xr_VTI_line Xr_UASTAR_line Xr_UTSTAR_line;
    global lc1_1 lc2_1;
    %% 设置图像面板标题
    set(PlotsPanels,'title',PlotsStrings{3});
    %% 设置相应按钮的可用状态
    %查找被按下的按钮数量
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %对按钮3进行操作后，若一个按钮也没有被按下
        %则除4、5外的所有按钮均可用
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
        %对按钮3进行操作后，若存在被按下的按钮
        %则除3、7外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(3),'enable','on');
        set(PlotsValues(7),'enable','on');
    end
    %% 绘制曲线
    if get(PlotsValues(3),'value')
        %避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        Xr = pt.input.ri;
        VAI = pt.input.VAI;
        VTI = pt.input.VTI;
        %如果此时被按下的按钮不为1，则说明7也被按下
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
            set(PlotsPanels,'title','速度分布曲线');
        end    
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'V/Vs','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanels,'title','点击左侧按钮显示对应图像');
    else
        delete(Xr_VAI_line);
        delete(Xr_VTI_line);
        legend([Xr_UASTAR_line,Xr_UTSTAR_line],'Ua*/Vs','Ut*/Vs');
        set(PlotsPanels,'title',PlotsStrings{7});
    end    
end

