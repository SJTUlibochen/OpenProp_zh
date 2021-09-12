% 该函数是新版GUI界面的测试程序
function NewGUITest
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
    Meanline_x_items = {'NACA a=0.8' 'NACA a=0.8 (modified)' 'Parabolic'};
    Thickness_x_items = {'NACA 65A010' 'NACA 65A010 (modified)' 'Elliptical' 'Parabolic' 'NACA 66 (DTRC modified)'};
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
    
    %% GUI元素的参数
    % 字体大小
    titlefontsize = 60;
    zhfontsize = 16;
    subcolumnfontsize = 18;
    numfontsize = 13;
    % 主界面初始高度和宽度
    Windowht = 750;
    Window = 1530;
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
                                    'rowheight',{'4x','1x','1x'},...
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
                    'verticalalignment','top',...
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
    FlagsValues(1) = uicheckbox('parent',HubPanelGrid,...
                                'text','设计中是否加入桨毂？',...
                                'value',Hub_flag_def,...
                                'tooltip','Hub_flag',...
                                'fontsize',zhfontsize);
    FlagsValues(1).Layout.Column = [1 2];
    SpecificationsTexts(6) = uilabel('parent',HubPanelGrid,...
                                     'text',SpecificationsStrings{index},...
                                     'horizontalalignment','left',...
                                     'verticalalignment','center',...
                                     'tooltip',SpecificationsTips{index},...
                                     'fontsize',zhfontsize);
    SpecificationsValues(6) = uieditfield(HubPanelGrid,'numeric',...
                                          'fontsize',numfontsize,...
                                          'horizontalalignment','center',...
                                          'value',SpecificationsValues_def{index},...
                                          'limit',[0 Inf],...
                                          'lowerlimitinclusive','off');
                                                          
    % 涵道面板
    DuctPanel = uipanel(SpecificationsGrid);
    DuctPanelGrid = uigridlayout(DuctPanel,[3 2],...
                                 'columnwidth',{'2x','3.5x'},...
                                 'columnspacing',0,...
                                 'rowspacing',10,...
                                 'padding',[10 10 10 10]);  
    FlagsValues(2) = uicheckbox('parent',DuctPanelGrid,...
                                'text','设计中是否加入涵道？',...
                                'value',Duct_flag_def,...
                                'tooltip','Duct_flag',...
                                'fontsize',zhfontsize);
    FlagsValues(2).Layout.Column = [1 2];
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
                                                  'lowerlimitinclusive','off');
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
                                    'items',{'CLmax','ConeyPLL','FAST2011dCTP','FAST2011dVAC','Brizzolara2007'},...
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
    Table1Tips = {'rx','CDp_x','Meanline','Thickness','CL_x','T0oDp_x，叶厚与叶片直径之比'};
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
                                          'fontsize',numfontsize);
                                      
        Thickness_x_in(index) = uidropdown('parent',Table1PanelGrid,...
                                           'items',Thickness_x_items,...
                                           'fontsize',numfontsize);
                                       
        CL_x_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                     'fontsize',numfontsize,...
                                     'horizontalalignment','center',...
                                     'value',CL_x_def(index),...
                                     'limit',[0 Inf]);
                                 
        T0oDp_x_in(index) = uieditfield(Table1PanelGrid,'numeric',...
                                       'fontsize',numfontsize,...
                                       'horizontalalignment','center',...
                                       'value',T0oDp_x_def(index),...
                                       'limit',[0 Inf]);
    end    
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
    Table2Tips = {'ri','VA_i，轴向来流速度与行进速度的比值','VT_i，切向来流速度与行进速度的比值'};
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
    ToolsGrid = uigridlayout(Tools,[2,2],...
                             'columnwidth',{'2x','1x'},...
                             'columnspacing',5,...
                             'rowspacing',5,...
                             'padding',[5 5 5 5],...
                             'scrollable','off');
    
    % 旗面板
    FlagsPanel = uipanel(ToolsGrid);
    FlagsPanelGrid = uigridlayout(FlagsPanel,[1,4]);
    FlagsStrings = {'是否绘制性能曲线？','是否输出点坐标？','是否输出stl文件',''}
    for index = 1 : 4
        FlagsValues(index+4) = uicheckbox('parent',DuctPanelGrid,...
                                          'text',,...
                                          'value',Duct_flag_def,...
                                          'tooltip','Duct_flag',...
                                          'fontsize',zhfontsize);
    end    
    % 控制点面板
    
    % 工具框
end