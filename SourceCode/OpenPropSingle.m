%该程序负责实现主界面和“仅进行叶片设计”功能
%关于数据结构的统一：统一使用char，需要索引时使用cell，不使用string
%关于函数命名的统一：关键词首字母大写，体现函数作用
%% 主界面部分
function OpenPropSingle
    clear variables;
    clear global;
    %% 设置全局变量
    global Fig_Main Tabs SingleDesign ParametricStudy TabPageStrings TabPageMenu Menu;
    global OpenPropDirectory Filename SpecificationsValues DuctValues FlagsValues FoilValues...
           N_R0 Xr_in Col_Label XCpoDp_in XCdp_in VAI_in VTI_in ri_in Xt0oDp_in skew0_in rake0_in...
           Meanline_cell Thickness_cell CalculatorsValues;
    global XCpoDp_values XCLmax_values XCdp_values;
    global zhfontsize numfontsize subcolumnfontsize buttonfontsize
    
    Meanline_cell = {'NACA a=0.8' 'NACA a=0.8 (modified)' 'Parabolic'};
    Thickness_cell = {'NACA 65A010' 'NACA 65A010 (modified)' 'Elliptical' 'Parabolic' 'NACA 66 (DTRC modified)'};
    
    %% “叶片设计参数”和“流入速度/装置行进速度”所需缺省值设定
    %根据叶切面半径与梢圆半径比值定义的径向位置
    Xr_def = [0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;0.95;1];
    N_R0 = length(Xr_def);      %径向位置的个数
    %定义径向位置处的弦径比，弦线长度与梢圆直径的比值
    XCpoDp_def = [0.1600;0.1812;0.2024;0.2196;0.2305;0.2311;0.2173;0.1807;0.1388;0.0010];
    %定义径向位置处的阻力系数(Drag Coefficient)，缺省值设为各处均为0.008
    XCdp_def = ones(N_R0,1).*0.008;
    %定义径向位置处的叶厚比δ，叶片的厚度与弦线长度比
    %t0oc0_def = [0.2055;0.1553;0.1180;0.09016;0.06960;0.05418;0.04206;0.03321;0.03228;0.03160];
    %定义径向位置处的厚径比，叶片的厚度与梢园直径的比值 Xt0oDp_def = t0oc0_def .* XCpoDp_def;
    Xt0oDp_def = [0.0329;0.0281;0.0239;0.0198;0.0160;0.0125;0.0091;0.0060;0.0045;0];
    %定义径向位置处的纵斜角
    skew0_def = zeros(N_R0,1);
    %定义径向位置处的侧斜角
    %实际上该变量最后是赋值给斜径比
    rake0_def = zeros(N_R0,1);
    %最大升力系数分布
    %下面的计算式是将叶梢到叶根的最大升力系数限制在0.2-0.5之间，沿径向方向线性分布
    XCLmax_def = 0.5 + (0.2-0.5)*(Xr_def-Xr_def(1))/(Xr_def(end)-Xr_def(1));
    %为变化回调设置默认值
    XCpoDp_values = XCpoDp_def;
    XCLmax_values = XCLmax_def;
    XCdp_values = XCdp_def;
    %定义径向位置处的轴向比值
    XVA_def = ones(N_R0,1);
    %定义径向位置处的切向比值
    XVT_def = zeros(N_R0,1);   
    %% “叶片规格”参数缺省值设定
    Z_def = 2;      %叶片数量
    N_def = 3500;        %螺旋桨转速(RPM)
    Dp_def = 0.4;      %螺旋桨直径(m)
    Thrust_def = 1.85;      %装置总推力要求(thrust,N)
    Vs_def = 5;     %装置行进速度(m/s)
    Dhub_def = 0.02;     %桨毂直径(m)
    rho_def = 1.29;     %流体密度(kg/m^3)
    Mp_def = 20;        %径向叶片分区数量(Number of vortex panels over the radius)
    Np_def = 20;        %切向叶片分区数量(Number of points over the chord)
    %% “无量纲参数”所需参数缺省值设定
    n_def = N_def/60;       %螺旋桨转速(r/s)
    lambda_def = n_def*2*pi*(Dp_def/2)/Vs_def;           %无量纲速度=叶梢圆周速度/装置行进速度(Non-dimensional)
    Js_def = Vs_def/(n_def*Dp_def);      %进速系数
    KT_def = Thrust_def/(rho_def*n_def^2*Dp_def^4);       %推力系数(Thrust Coeffcient)，以叶梢速度为特征速度
    CT_def = Thrust_def/(0.5*rho_def*Vs_def^2*pi*(Dp_def/2)^2);       %推力系数(Thrust Coefficient),以装置行进速度为特征速度    
    %% 导管螺旋桨相关参数Ducted Propeller variables
    TpoT_def = 1;        %螺旋桨推力占比，即螺旋桨提供的推力占总推力的比例(Thrust ratio)
    Cdd_def = 0.008;        %导管切面的阻力系数(Duct section drag coefficient)
    DdoDp_def = 1;        %导管直径/螺旋桨直径(duct D / prop D)
    filename = 'DefaultPropeller';      %设定文件名前缀的缺省值
    %% GUI中使用的元素尺寸的常量值
    %中文字体大小
    zhfontsize = 15;
    %编辑框内数字和单位大小
    numfontsize = 13;
    %叶片规格 & 叶片设计参数 & 流入速度/装置行进速度 & 设计选项 & 涵道螺旋桨 & 无量纲参数 & 工具
    %子栏标题字体大小
    subcolumnfontsize = 18;
    %按钮字体大小  Load & Save & Run OpenProp
    buttonfontsize = 16;
    %窗口的总高度和总宽度
    Windowht = 780;
    Window = 1530;
    %文件目录
    OpenPropDirectory = 'OpenProp_zh';
    %% 主界面绘制    
    close all;
    Fig_Main = uifigure('position',[5 40 Window Windowht],...
                        'resize','on',...
                        'autoresizechildren','on',...
                        'numbertitle','off',...
                        'name','OpenProp');
    %% 工具栏绘制
    %生成第一级工具栏
    MenuStrings = {'文件' '标签页'};
    MenuTips = {'进行文件相关操作' '点击对应标签页即可关闭'};
    Menu = zeros(1,length(MenuStrings));
    for index = 1 : length(MenuStrings)
        Menu(index) = uimenu('parent',Fig_Main,...
                             'text',MenuStrings{index},...
                             'tooltip',MenuTips{index});
    end    
    %生成子工具栏
    FileStrings = {'保存' '另存为' '载入' '运行'};
    FileMenu = zeros(1,length(FileStrings));
    for index = 1 : length(FileStrings)
        FileMenu(index) = uimenu('parent',Menu(1),'text',FileStrings{index});
    end    
    set(FileMenu(1),'menuselectedfcn',@SaveData);
    set(FileMenu(2),'menuselectedfcn',@SaveDataCopy);
    set(FileMenu(3),'menuselectedfcn',@LoadData);
    set(FileMenu(4),'menuselectedfcn',@Execute);
    TabPageStrings = {'仅进行叶片设计' '研究参数的影响'};
    for index = 1 : length(TabPageStrings)
        TabPageMenu(index) = uimenu('parent',Menu(2),...
                                    'text',TabPageStrings{index},...
                                    'enable','off');
    end    
    %% 生成主界面上的标签页
    %生成主界面中的标签栏uitabgroup
    Tabs = uitabgroup('parent',Fig_Main,...
                      'position',[0 0 Window Windowht],...
                      'selectionchangedfcn',@ToggleMode);
    %生成Tabs的第一个标签页“仅进行叶片设计”
    SingleDesign = uitab('parent',Tabs,...
                         'title','仅进行叶片设计');
    %生成Tabs的第二个标签页“研究参数的影响”
    ParametricStudy = uitab('parent',Tabs,...
                            'title','研究参数的影响');
    
    %% 生成SingleDesign中的网格空间及空间中的各栏目
    %生成SingleDesign中2*4的网格空间
    SingleDesignGrid = uigridlayout(SingleDesign,[2 4],...
                                    'columnwidth',{'2x','4x','2x','1.05x'},...
                                    'rowheight',{'3x','1x'});
    
    %生成“叶片规格”栏，原版'Specifications'
    Specifications = uipanel('parent',SingleDesignGrid,...
                             'title','叶片规格',...
                             'titleposition','centertop',...
                             'fontsize',subcolumnfontsize,...
                             'fontweight','bold',...
                             'scrollable','off');
    
    %生成“叶片设计参数”栏，原版'Blade Design Values'
    BladeDesign = uipanel('parent',SingleDesignGrid,...
                          'title','叶片设计参数',...
                          'titleposition','centertop',...
                          'fontsize',subcolumnfontsize,...
                          'fontweight','bold',...
                          'scrollable','on');
    
    %生成“流入速度/装置行进速度”栏，原版'Inflow Profile Values'
    Inflow = uipanel('parent',SingleDesignGrid,...
                     'title','流入速度/装置行进速度',...
                     'titleposition','centertop',...
                     'fontsize',subcolumnfontsize,...
                     'fontweight','bold',...
                     'scrollable','on');
    
    %生成“设计选项”栏，原版'Options'
    Flags = uipanel('parent',SingleDesignGrid,...
                    'title','设计选项',...
                    'titleposition','centertop',...
                    'fontsize',subcolumnfontsize,...
                    'fontweight','bold',...
                    'scrollable','on');
                         
    %生成“涵道螺旋桨”栏，原版'Ducted Propeller'
    Duct = uipanel('parent',SingleDesignGrid,...
                   'title','涵道螺旋桨',...
                   'titleposition','centertop',...
                   'fontsize',subcolumnfontsize,...
                   'fontweight','bold',...
                   'scrollable','off');

    %生成“无量纲参数”栏，原版'Non-dimensional Parameters'
    Calculators = uipanel('parent',SingleDesignGrid,...
                          'title','无量纲参数',...
                          'titleposition','centertop',...
                          'fontsize',subcolumnfontsize,...
                          'fontweight','bold',...
                          'scrollable','on');

    %生成“工具”栏，原版'Tools'
    Tools = uipanel('parent',SingleDesignGrid,...
                    'title','工具',...
                    'titleposition','centertop',...
                    'fontsize',subcolumnfontsize,...
                    'fontweight','bold',...
                    'scrollable','off');
    Tools.Layout.Column = [3 4];        %“工具”一栏横跨第二行的3、4列
       
    %% “叶片规格”子栏中的内容
    %生成Specifications中11*2的网格空间
    SpecificationsGrid = uigridlayout(Specifications,[11 2],...
                                      'scrollable','off');    
    %Specifications中的前9行参数   
    SpecificationsStrings = {'叶片数量：'...     %原版'Number of blades'
                             '螺旋桨转速：'...        %原版'Rotation speed'
                             '螺旋桨直径：'...        %原版'Rotor diameter'
                             '装置总推力要求：'...      %原版'Required thrust'
                             '装置行进速度：'...       %原版'Ship speed'
                             '桨毂直径：'...     %原版'Hub diameter'
                             '流体密度：'...     %原版'Fluid density'
                             '径向叶片分区数量：'...     %原版'radial panels'
                             '切向叶片分区数量：'};      %原版'chordwise panels'
    SpecificationsTips = {'Z' 'N' 'Dp' 'Thrust' 'Vs' 'Dhub' 'ρ' 'Mp' 'Np'};
    SpecificationsValues_def = {Z_def N_def Dp_def Thrust_def Vs_def Dhub_def rho_def Mp_def Np_def};
    %利用循环绘制前9行
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
    %设置前9行编辑框中的显示形式和单位
    set(SpecificationsValues(1),'valuedisplayformat','%.0f');
    set(SpecificationsValues(2),'valuedisplayformat','%.1f RPM');
    set(SpecificationsValues(3),'valuedisplayformat','%.2f m');
    set(SpecificationsValues(4),'valuedisplayformat','%.3f N');
    set(SpecificationsValues(5),'valuedisplayformat','%.2f m/s');
    set(SpecificationsValues(6),'valuedisplayformat','%.3f m');
    set(SpecificationsValues(7),'valuedisplayformat','%.3f kg/m^3');
    set(SpecificationsValues(8),'valuedisplayformat','%.0f');
    set(SpecificationsValues(9),'valuedisplayformat','%.0f');
    
    %Specifications中的后2行参数
    FoilStrings = {'凸缘线类型：'...      %原版'Meanline type'
                   '叶型类型：'};        %原版'Thickness type' 
    %利用循环绘制后2行
    for index = 1 : 2
        FoilTexts(index) = uilabel('parent',SpecificationsGrid,...
                                   'text',FoilStrings{index},...
                                   'horizontalalignment','left',...
                                   'verticalalignment','center',...
                                   'fontsize',zhfontsize);
        FoilValues(index) = uidropdown('parent',SpecificationsGrid,...
                                       'fontsize',numfontsize);
    end
    %设置后2行选框内容
    set(FoilValues(1),'items',Meanline_cell);
    set(FoilValues(2),'items',Thickness_cell);
    set(FoilValues(1),'itemsdata',[1 2 3]);
    set(FoilValues(2),'itemsdata',[1 2 3 4 5]);
    %% “涵道螺旋桨”子栏中的内容
    %生成Duct中3*2的网格空间
    DuctGrid = uigridlayout(Duct,[3 2],...
                            'columnwidth',{'1x','1x'},...
                            'rowheight',{'1x','1x','1x'},...
                            'scrollable','off');
    %Duct中的参数
    DuctStrings = {'螺旋桨推力占比：'...        %原版'Thrust Ratio'
                   '导管切面阻力系数：'...     %原版'Duct section drag (Cd)'
                   '导管直径/螺旋桨直径：'};      %原版'duct D / prop D'
    DuctTips = {'Tp/Thrust' 'Cdd' 'Dd/Dp'};
    DuctValues_def = {TpoT_def Cdd_def DdoDp_def};
    %利用循环绘制3行
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
    %设置参数范围
    set(DuctValues(1),'limit',[0 1]);
    set(DuctValues(1),'lowerlimitinclusive','off');
    set(DuctValues(2),'limit',[0 Inf]);
    set(DuctValues(3),'limit',[1 Inf]);
    %% “叶片设计参数”子栏中的内容
    %生成BladeDesign中11*6的网格空间
    BladeDesignGrid = uigridlayout(BladeDesign,[N_R0+1 6],...
                                   'columnspacing',0,...
                                   'rowspacing',0,...
                                   'scrollable','on');   
    %BladeDesign中的参数
    BladeDesignStrings = {'径向位置'...        %原版'r/R'
                          '弦径比'...     %原版'c/D'
                          '阻力系数'...        %原版'Cd'
                          '厚径比'...     %原版't0/D'
                          '侧倾角'...     %原版'Skew'
                          '斜径比'};      %原版'Xs/D'，应修正为'Zr/D'
    BladeDesignTips = {'r/Rp' 'Cp/Dp' 'Cdp' 't0/Dp' 'skew' 'Zr/Dp'};
    %利用循环绘制表格的题目
    for index = 1 : length(BladeDesignStrings)
        Col_Label(index) = uilabel('parent',BladeDesignGrid,...
                                   'text',BladeDesignStrings{index},...
                                   'horizontalalignment','center',...
                                   'verticalalignment','center',...
                                   'tooltip',BladeDesignTips{index},...
                                   'fontsize',zhfontsize);
    end
    %利用循环绘制表格的内容
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
    %% “流入速度/装置行进速度”子栏中的内容
    %生成Inflow中11*3的网格空间
    InflowGrid = uigridlayout(Inflow,[N_R0+1 3],...
                              'columnspacing',0,...
                              'rowspacing',0,...
                              'scrollable','on');  
    %Inflow中的参数
    InflowStrings = {'径向位置'...       %原版'r'，推测为'r/R'
                     '轴向比值'...       %原版'Va/Vs'
                     '切向比值'};        %原版'Vt/Vs'
    InflowTips = {'r/Rp' 'Va/Vs' 'Vt/Vs'};  
    %利用循环绘制表格
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
    %% “无量纲参数”子栏中的内容
    %生成Calculators中2*2的网格空间
    CalculatorsGrid = uigridlayout(Calculators,[3 4],...
                                  'scrollable','on');
    %Calculators中的参数
    CalculatorsStrings = {'进速系数：'...     %原版'Js = V/nD ='
                          '叶梢速度/行进速度：'...     %原版'L = omega*R/V ='
                          '推力系数(叶梢速度)：'...       %原版'KT = T/(rho*n^2*D^4) ='
                          '推力系数(行进速度)：'};      %原版'CT = T/(1/2*rho*V^2*pi*R^2) ='
    CalculatorsTips = {'Js = V/n*Dp' 'L = ω*Rp/V' 'KT = T/(ρ*n^2*Dp^4)' 'CT = T/(0.5πρ*V^2*Rp^2)'};
    CalculatorsValues_def = {Js_def lambda_def KT_def CT_def};  
    %利用循环绘制4个无量纲参数
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
    %% “设计选项”子栏中的内容
    %生成Flags中10*1的网格空间
    FlagsGrid = uigridlayout(Flags,[10 1],...
                             'rowspacing',15,...
                             'scrollable','on');
    %Flags中的参数
    SelectStrings = {'螺旋桨'...       %原版'Propeller'
                     '涡轮叶片'};       %原版'Turbine'
    FlagsStrings = {'是否输出桨毂'...     %原版'Hub'，该复选按钮用于判断输出图像中是否包括桨毂
                    '是否为涵道螺旋桨'...       %原版'Ducted'，该复选按钮用于判断是否生成导管螺旋桨
                    '是否进行弦长优化'...       %原版'Chord Optimization'，该复选按钮用于判断是否进行弦长优化
                    '是否考虑粘性阻力'...       %原版'Viscous forces'，该复选按钮用于判断是否考虑粘性力
                    '是否绘制性能曲线'...       %原版'Performance curve'，该复选按钮用于判断是否绘制性能曲线
                    '是否绘制叶片模型'...       %原版'Geometry plots'，该复选按钮用于判断是否输出2D和3D图形
                    '是否导出几何文档'...       %新增复选按钮，该复选按钮用于判断是否导出关于叶片几何信息的文档(txt)
                    '是否导出点坐标'...        %新增复选按钮，该复选按钮用于判断是否导出点坐标文档(csv)
                    '是否导出stl文件'};        %新增复选按钮，该复选按钮用于判断是否导出stl格式的可打印文件
    FlagsTips = {'Hub_flag，勾选后生成桨毂' 'Duct_flag，勾选后生成涵道' 'Chord_flag，勾选后进行弦长优化'...
                 'Viscous_flag，勾选后考虑粘性阻力' 'Analyze_flag，勾选后输出性能曲线' 'Geometry_flag，勾选后显示叶片模型'...
                 'txt_flag，勾选后输出txt格式叶片几何信息' 'csv_flag，勾选后输出csv格式坐标' 'stl_flag，勾选后输出stl格式模型'};
    FlagsValues_def = {1 1 0 1 0 1 1 1 0};    
    %绘制选框，使用'itemsdata'属性将选择文本分别与1和0对应，选择“螺旋桨”则返回值为1
    FlagsValues(1) = uidropdown('parent',FlagsGrid,...
                                'items',SelectStrings,...
                                'itemsdata',[1 0],...
                                'fontsize',zhfontsize,...
                                'valuechangedfcn',@PropTurb);
    %利用循环绘制9行复选框
    for index = 1 : length(FlagsStrings)
        FlagsValues(index+1) = uicheckbox('parent',FlagsGrid,...
                                          'text',FlagsStrings{index},...
                                          'value',FlagsValues_def{index},...
                                          'fontsize',zhfontsize,...
                                          'tooltip',FlagsTips{index});
    end
    %设置复选框的回调函数
    set(FlagsValues(3),'valuechangedfcn',@ChangeDuct);
    set(FlagsValues(4),'valuechangedfcn',@ChangeChord);
    set(FlagsValues(5),'valuechangedfcn',@ChangeViscous);
    %% “工具”子栏中的内容
    %生成Tools中2*4的网格空间
    ToolsGrid = uigridlayout(Tools,[2 4],...
                             'columnwidth',{'1x','1x','1x','1x'},...
                             'rowheight',{'2x','3x'},...
                             'columnspacing',25,...
                             'rowspacing',25,...
                             'scrollable','off');
    %绘制文件名文本框
    ToolsText = uilabel('parent',ToolsGrid,...
                        'text','当前项目名称：',...
                        'fontsize',zhfontsize,...
                        'horizontalalignment','center',...
                        'verticalalignment','center');
    %绘制编辑框
    Filename = uieditfield(ToolsGrid,...
                           'fontsize',numfontsize,...
                           'horizontalalignment','center',...
                           'value',filename);
    Filename.Layout.Column = [2 4];
    %绘制加载按钮
    LoadButton = uibutton('parent',ToolsGrid,...
                          'text','',...
                          'tooltip','将.mat文件载入GUI',...
                          'icon','load.png',...
                          'fontsize',buttonfontsize,...
                          'buttonpushedfcn',@LoadData);
    %绘制保存按钮
    SaveButton = uibutton('parent',ToolsGrid,...
                          'text','',...
                          'tooltip','将GUI保存至.mat文件',...
                          'icon','save.png',...
                          'fontsize',buttonfontsize,...
                          'buttonpushedfcn',@SaveData);
    %绘制另存为按钮
    SaveCopyButton = uibutton('parent',ToolsGrid,...
                              'text','',...
                              'tooltip','将GUI另存为.mat文件',...
                              'icon','savecopy.png',...
                              'fontsize',buttonfontsize,...
                              'buttonpushedfcn',@SaveDataCopy);
    %绘制运行按钮
    RunButton = uibutton('parent',ToolsGrid,...
                         'text','',...
                         'tooltip','运行OpenProp',...
                         'icon','run.png',...
                         'fontsize',buttonfontsize,...
                         'buttonpushedfcn',@Execute);
end
%% 切换“仅进行叶片设计”和“研究参数的影响”后回调的函数
function ToggleMode(hObject,ED)
    global Tabs SingleDesign ParametricStudy;
    global Filename;
    global SpecificationsValues DuctValues FlagsValues FoilValues CalculatorsValues...
           Xr_in XCpoDp_in XCdp_in VAI_in VTI_in ri_in Xt0oDp_in skew0_in rake0_in...
           Meanline_cell Thickness_cell XCpoDp_values XCLmax_values;
    if get(Tabs,'selectedtab') == SingleDesign
        %当前选择的模式为“仅进行叶片设计”，无需进行操作
    elseif get(Tabs,'selectedtab') == ParametricStudy
        %当前选择的模式为“研究参数的影响”，需要将“仅进行叶片设计”界面的数据暂存至临时文件中
        filename = get(Filename,'value');
        Z = get(SpecificationsValues(1),'value');       %叶片数量
        N = get(SpecificationsValues(2),'value');       %螺旋桨转速(RPM)
        Dp = get(SpecificationsValues(3),'value');      %螺旋桨直径(m)
        Thrust = get(SpecificationsValues(4),'value');      %装置总推力要求(N)
        Vs = get(SpecificationsValues(5),'value');      %装置行进速度(m/s)
        Dhub = get(SpecificationsValues(6),'value');        %桨毂直径(m)
        rho = get(SpecificationsValues(7),'value');     %流体密度(kg/m^3)
        Mp = get(SpecificationsValues(8),'value');      %径向叶片分区数量
        Np = get(SpecificationsValues(9),'value');      %切向叶片分区数量
        Meanline_index = get(FoilValues(1),'value');
        Meanline = char(Meanline_cell(Meanline_index));     %凸缘线类型
        Thickness_index	= get(FoilValues(2),'value');
        Thickness = char(Thickness_cell(Thickness_index));      %叶型类型
        TpoT = get(DuctValues(1),'value');      %螺旋桨推力占比
        Cdd = get(DuctValues(2),'value');       %导管切面阻力系数
        Rduct_oR = get(DuctValues(3),'value');     %导管直径/螺旋桨直径
        Js = get(CalculatorsValues(1),'value');      %进速系数
        L = get(CalculatorsValues(2),'value');       %叶梢速度/行进速度
        KT = get(CalculatorsValues(3),'value');      %推力系数(叶梢速度)
        CT = get(CalculatorsValues(4),'value');      %推力系数(行进速度)
        Propeller_flag = get(FlagsValues(1),'value');       %螺旋桨/涡轮叶片
        Hub_flag = get(FlagsValues(2),'value');     %是否输出桨毂
        Duct_flag = get(FlagsValues(3),'value');        %是否为涵道螺旋桨
        Chord_flag = get(FlagsValues(4),'value');       %是否进行弦长优化
        Viscous_flag = get(FlagsValues(5),'value');     %是否考虑粘性阻力
        Analyze_flag = get(FlagsValues(6),'value');     %是否绘制性能曲线
        Geometry_flag = get(FlagsValues(7),'value');        %是否绘制叶片模型
        txt_flag = get(FlagsValues(8),'value');     %是否导出几何文档
        csv_flag = get(FlagsValues(9),'value');     %是否导出点坐标
        stl_flag = get(FlagsValues(10),'value');     %是否导出stl文件   
        Xr = get(Xr_in,'value');        %径向位置
        XCpoDp = get(XCpoDp_in,'value');        %弦径比
        XCLmax = get(XCpoDp_in,'value');      %最大升力系数
        XCdp = get(XCdp_in,'value');      %阻力系数
        Xt0oDp = get(Xt0oDp_in,'value');      %厚径比
        skew0 = get(skew0_in,'value');      %侧倾角
        rake0 = get(rake0_in,'value');      %斜径比
        ri = get(ri_in,'value');        %径向位置
        VAI = get(VAI_in,'value');      %轴向比值
        VTI = get(VTI_in,'value');      %切向比值
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
%% “叶片规格”中参数改变后调用的函数
function UpdateSpecifications(hObject,ED)
    %更新“无量纲参数”中的变量值，与“叶片规格”中的参数保持一致
    global SpecificationsValues CalculatorsValues Fig_Main;
    Z = round(get(SpecificationsValues(1),'value'));        %叶片数量
    set(SpecificationsValues(1),'value',Z);
    N = get(SpecificationsValues(2),'value');      %螺旋桨转速(RPM)
    Dp = get(SpecificationsValues(3),'value');      %螺旋桨直径(m)
    Thrust = get(SpecificationsValues(4),'value');     %装置总推力要求(N)
    Vs = get(SpecificationsValues(5),'value');     %装置行进速度(m/s)
    Dhub = get(SpecificationsValues(6),'value');       %桨毂直径(m)
    rho = get(SpecificationsValues(7),'value');        %流体密度(kg/m^3)
    Mp = round(get(SpecificationsValues(8),'value'));      %径向叶片分区数量
    set(SpecificationsValues(8),'value',Mp);
    Np = round(get(SpecificationsValues(9),'value'));      %切向叶片分区数量
    set(SpecificationsValues(9),'value',Np);
    %检查桨毂直径与螺旋桨直径之间的关系，若前者大于后者，则弹出提醒框
    if Dhub >= Dp
        message = sprintf('当前设置的桨毂直径大于螺旋桨直径 \n 请重新设置！');
        uialert(Fig_Main,message,'参数设置错误','closefcn',@ResetDiameters);
    end
    n = N/60;       %螺旋桨转速(r/s)
    lambda = n*2*pi*(Dp/2)/Vs;       %无量纲速度=叶尖圆周速度/装置行进速度(Non-dimensional)
    Js = Vs/(n*Dp);      %进速系数
    KT = Thrust/(rho*n^2*Dp^4);      %推力系数(Thrust Coeffcient)，以叶梢速度为特征速度
    CT = Thrust/(0.5*rho*Vs^2*pi*(Dp/2)^2);      %推力系数(Thrust Coefficient)，以装置行进速度为特征速度
    set(CalculatorsValues(1),'value',Js);
    set(CalculatorsValues(2),'value',lambda);
    set(CalculatorsValues(3),'value',KT);
    set(CalculatorsValues(4),'value',CT);
end
%发生直径设置错误后调用的函数
function ResetDiameters(hObject,ED)
    global SpecificationsValues;
    set(SpecificationsValues(3),'value',0.4);
    set(SpecificationsValues(6),'value',0.02);
end
%% 修改“螺旋桨”/“涡轮叶片”按钮后回调的函数
function PropTurb(hObject,ED)
    global FlagsValues SpecificationsValues;
    if get(FlagsValues(1),'value')
        set(SpecificationsValues(4),'editable','on')
    else
        %选择“涡轮叶片”按钮后，“叶片规格”子栏下的“装置总推力要求”变为不可编辑
        set(SpecificationsValues(4),'editable','off')
    end
end
%% 修改“是否为涵道螺旋桨”按钮后回调的函数
function ChangeDuct(hObject,ED)
    global DuctValues FlagsValues;
    if get(FlagsValues(3),'value')
        set(DuctValues,'editable','on');
    else
        set(DuctValues,'editable','off');
    end
end
%% 修改“是否进行弦长优化”按钮后回调的函数
function ChangeChord(hObject,ED)
    global FlagsValues XCpoDp_in Col_Label;
    global XCpoDp_values XCLmax_values;
    if get(FlagsValues(4),'value')
        set(Col_Label(2),'text','最大升力系数');
        for index = 1:length(XCpoDp_in)
            XCpoDp_values(index) = get(XCpoDp_in(index),'value');
            set(XCpoDp_in(index),'value',XCLmax_values(index));
        end
    else
        set(Col_Label(2),'text','弦径比');
        for index = 1:length(XCpoDp_in)
            XCLmax_values(index) = get(XCpoDp_in(index),'value');
            set(XCpoDp_in(index),'value',XCpoDp_values(index));
        end        
    end
end
%% 修改“是否考虑粘性阻力”按钮后回调的函数
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
%% 点击“加载”后回调的函数
function LoadData(hObject,ED)
    global N_R0 OpenPropDirectory SpecificationsValues DuctValues FlagsValues FoilValues Filename...
           Xr_in XCpoDp_in XCdp_in VAI_in VTI_in ri_in Xt0oDp_in skew0_in rake0_in;
    global pt
    rest = pwd;
    LoadDir(rest,OpenPropDirectory);
    %uiload函数用于打开文件选择目录
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
%% 点击“保存”后回调的函数
function SaveData(hObject,ED)
    global OpenPropDirectory SpecificationsValues DuctValues CalculatorsValues FlagsValues FoilValues Filename...
           Xr_in XCpoDp_in XCdp_in VAI_in VTI_in ri_in Xt0oDp_in skew0_in rake0_in...
           Meanline_cell Thickness_cell;
    global Fig_Main pt;
    %% 从GUI界面中读取输入值
    filename = get(Filename,'value');       %当前项目名称
    %existence表示是否存在同名项目
    existence = ChangeDirectory(OpenPropDirectory,filename);
    if existence
        message = sprintf('当前项目名已经存在 \n 如果继续进行该操作，将覆盖原本文件 \n 确定继续吗？');
        SaveSelection = uiconfirm(Fig_Main,message,'覆盖警告',...
                                  'options',{'覆盖原版文件','另存为新名称','取消操作'},...
                                  'defaultoption',2,...
                                  'canceloption',3,...
                                  'icon','warning');
        if strcmp(SaveSelection,'覆盖原版文件')
            %选择“覆盖原始文件”,不需要进行操作
        elseif strcmp(SaveSelection,'另存为新名称')
            %选择“另存为新名称”，调用savecopy函数
            SaveDataCopy;
        else
            %选择“取消操作”，退出
            return
        end    
    end    
    Z = get(SpecificationsValues(1),'value');       %叶片数量
    N = get(SpecificationsValues(2),'value');       %螺旋桨转速(RPM)
    Dp = get(SpecificationsValues(3),'value');      %螺旋桨直径(m)
    Thrust = get(SpecificationsValues(4),'value');      %装置总推力要求(N)
    Vs = get(SpecificationsValues(5),'value');      %装置行进速度(m/s)
    Dhub = get(SpecificationsValues(6),'value');        %桨毂直径(m)
    rho = get(SpecificationsValues(7),'value');     %流体密度(kg/m^3)
    Mp = get(SpecificationsValues(8),'value');      %径向叶片分区数量
    Np = get(SpecificationsValues(9),'value');      %切向叶片分区数量
    Meanline_index = get(FoilValues(1),'value');
    Meanline = char(Meanline_cell(Meanline_index));     %凸缘线类型
    Thickness_index	= get(FoilValues(2),'value');
    Thickness = char(Thickness_cell(Thickness_index));      %叶型类型
    TpoT = get(DuctValues(1),'value');      %螺旋桨推力占比
    Cdd = get(DuctValues(2),'value');       %导管切面阻力系数
    Rduct_oR = get(DuctValues(3),'value');     %导管直径/螺旋桨直径
    Js = get(CalculatorsValues(1),'value');      %进速系数
    L = get(CalculatorsValues(2),'value');       %叶梢速度/行进速度
    KT = get(CalculatorsValues(3),'value');      %推力系数(叶梢速度)
    CT = get(CalculatorsValues(4),'value');      %推力系数(行进速度)
    Propeller_flag = get(FlagsValues(1),'value');       %螺旋桨/涡轮叶片
    Hub_flag = get(FlagsValues(2),'value');     %是否输出桨毂
    Duct_flag = get(FlagsValues(3),'value');        %是否为涵道螺旋桨
    Chord_flag = get(FlagsValues(4),'value');       %是否进行弦长优化
    Viscous_flag = get(FlagsValues(5),'value');     %是否考虑粘性阻力
    Analyze_flag = get(FlagsValues(6),'value');     %是否绘制性能曲线
    Geometry_flag = get(FlagsValues(7),'value');        %是否绘制叶片模型
    txt_flag = get(FlagsValues(8),'value');     %是否导出几何文档
    csv_flag = get(FlagsValues(9),'value');     %是否导出点坐标
    stl_flag = get(FlagsValues(10),'value');     %是否导出stl文件   
    Xr = cell2array(get(Xr_in,'value'));        %径向位置
    XCpoDp = cell2array(get(XCpoDp_in,'value'));        %弦径比
    XCLmax = cell2array(get(XCpoDp_in,'value'));      %最大升力系数
    XCdp = cell2array(get(XCdp_in,'value'));      %阻力系数
    Xt0oDp = cell2array(get(Xt0oDp_in,'value'));      %厚径比
    skew0 = cell2array(get(skew0_in,'value'));      %侧倾角
    rake0 = cell2array(get(rake0_in,'value'));      %斜径比
    ri = cell2array(get(ri_in,'value'));        %径向位置
    VAI = cell2array(get(VAI_in,'value'));      %轴向比值
    VTI = cell2array(get(VTI_in,'value'));      %切向比值
    ITER = 40;      %迭代次数
    Rhv = 0.5;      %hub vortex radius / hub radius
    %% 将GUI界面的输入值保存至结构体数组input
    input.part1 = '叶片规格';
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
    input.part2 = '叶片设计参数';
    input.Xr = Xr;
    input.XCdp = XCdp;
    input.XCLmax = XCLmax;
    input.XCpoDp = XCpoDp;
    input.Xt0oDp = Xt0oDp;
    input.skew0 = skew0;
    input.rake0 = rake0;
    input.part3 = '流入速度/装置行进速度';
    input.ri = ri;
    input.VAI = VAI;
    input.VTI = VTI;
    input.part4 = '设计选项';
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
    input.part5 = '涵道螺旋桨';
    input.TpoT = TpoT;
    input.Cdd = Cdd;
    input.Rduct_oR = Rduct_oR;
    input.part6 = '无量纲参数';
    input.Js = Js;
    input.L = L;
    input.KT = KT;
    input.CT = CT;
    input.part7 = '其他非输入参数';
    input.ITER = ITER;
    input.Rhv = Rhv;
    %% 将所有项目信息封装至结构体数组pt中
    pt.filename = filename;
    pt.date = date;
    pt.input = input;
    pt.design = [];     %结构体数组：design conditions
    pt.geometry = [];       %结构体数组：design geometry
    pt.states = [];     %结构体数组：off-design state analysis
    %% 保存
    save(filename,'pt');
end
%% 点击“另存为”后回调的函数
function SaveDataCopy(hObject,ED)
    global Fig_Main SpecificationsValues DuctValues CalculatorsValues FlagsValues FoilValues Filename...
           Xr_in XCpoDp_in XCdp_in VAI_in VTI_in ri_in Xt0oDp_in skew0_in rake0_in...
           Meanline_cell Thickness_cell;
    global pt;
    %% 从GUI界面中读取输入值
    filename = get(Filename,'value');       %当前项目名称
    Z = get(SpecificationsValues(1),'value');       %叶片数量
    N = get(SpecificationsValues(2),'value');       %螺旋桨转速(RPM)
    Dp = get(SpecificationsValues(3),'value');      %螺旋桨直径(m)
    Thrust = get(SpecificationsValues(4),'value');      %装置总推力要求(N)
    Vs = get(SpecificationsValues(5),'value');      %装置行进速度(m/s)
    Dhub = get(SpecificationsValues(6),'value');        %桨毂直径(m)
    rho = get(SpecificationsValues(7),'value');     %流体密度(kg/m^3)
    Mp = get(SpecificationsValues(8),'value');      %径向叶片分区数量
    Np = get(SpecificationsValues(9),'value');      %切向叶片分区数量
    Meanline_index = get(FoilValues(1),'value');
    Meanline = char(Meanline_cell(Meanline_index));     %凸缘线类型
    Thickness_index	= get(FoilValues(2),'value');
    Thickness = char(Thickness_cell(Thickness_index));      %叶型类型
    TpoT = get(DuctValues(1),'value');      %螺旋桨推力占比
    Cdd = get(DuctValues(2),'value');       %导管切面阻力系数
    Rduct_oR = get(DuctValues(3),'value');     %导管直径/螺旋桨直径
    Js = get(CalculatorsValues(1),'value');      %进速系数
    L = get(CalculatorsValues(2),'value');       %叶梢速度/行进速度
    KT = get(CalculatorsValues(3),'value');      %推力系数(叶梢速度)
    CT = get(CalculatorsValues(4),'value');      %推力系数(行进速度)
    Propeller_flag = get(FlagsValues(1),'value');       %螺旋桨/涡轮叶片
    Hub_flag = get(FlagsValues(2),'value');     %是否输出桨毂
    Duct_flag = get(FlagsValues(3),'value');        %是否为涵道螺旋桨
    Chord_flag = get(FlagsValues(4),'value');       %是否进行弦长优化
    Viscous_flag = get(FlagsValues(5),'value');     %是否考虑粘性阻力
    Analyze_flag = get(FlagsValues(6),'value');     %是否绘制性能曲线
    Geometry_flag = get(FlagsValues(7),'value');        %是否绘制叶片模型
    txt_flag = get(FlagsValues(8),'value');        %是否导出集合文档
    csv_flag = get(FlagsValues(9),'value');        %是否导出点坐标
    stl_flag = get(FlagsValues(10),'value');     %是否导出stl文件   
    Xr = cell2array(get(Xr_in,'value'));        %径向位置
    XCpoDp = cell2array(get(XCpoDp_in,'value'));        %弦径比
    XCLmax = cell2array(get(XCpoDp_in,'value'));      %最大升力系数
    XCdp = cell2array(get(XCdp_in,'value'));      %阻力系数
    Xt0oDp = cell2array(get(Xt0oDp_in,'value'));      %厚径比
    skew0 = cell2array(get(skew0_in,'value'));      %侧倾角
    rake0 = cell2array(get(rake0_in,'value'));      %斜径比
    ri = cell2array(get(ri_in,'value'));        %径向位置
    VAI = cell2array(get(VAI_in,'value'));      %轴向比值
    VTI = cell2array(get(VTI_in,'value'));      %切向比值
    ITER = 40;      %迭代次数
    Rhv = 0.5;      %hub vortex radius / hub radius
    %% 将GUI界面的输入值保存至结构体数组input
    input.part1 = '叶片规格';
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
    input.part2 = '叶片设计参数';
    input.Xr = Xr;
    input.XCdp = XCdp;
    input.XCLmax = XCLmax;
    input.XCpoDp = XCpoDp;
    input.Xt0oDp = Xt0oDp;
    input.skew0 = skew0;
    input.rake0 = rake0;
    input.part3 = '流入速度/装置行进速度';
    input.ri = ri;
    input.VAI = VAI;
    input.VTI = VTI;
    input.part4 = '设计选项';
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
    input.part5 = '涵道螺旋桨';
    input.TpoT = TpoT;
    input.Cdd = Cdd;
    input.Rduct_oR = Rduct_oR;
    input.part6 = '无量纲参数';
    input.Js = Js;
    input.L = L;
    input.KT = KT;
    input.CT = CT;
    input.part7 = '其他非输入参数';
    input.ITER = ITER;
    input.Rhv = Rhv;
    %% 将所有项目信息封装至结构体数组pt中
    pt.filename = [filename,'-copy'];
    pt.date = date;
    pt.input = input;
    pt.design = [];     %结构体数组：design conditions
    pt.geometry = [];       %结构体数组：design geometry
    pt.states = [];     %结构体数组：off-design state analysis
    %目前导入和另存为完成后，屏幕焦点会回到matlab，这个问题需要解决
    uisave('pt',[filename,'-copy']);
end
%将直接从GUI读取的cell转化为double数组的工具函数
function array = cell2array(cell)
    celllen = length(cell);
    array = zeros(celllen,1);
    for index = 1 : celllen
        array(index) = cell{index};
    end
end
%% 执行部分
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

