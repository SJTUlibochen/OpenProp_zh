%绘制在OpenPropSingle.m中没有设置的按钮回调函数，以及对应的曲线
function Make_Reports
    global PlotsValues;
    %设置按钮的回调函数
    set(PlotsValues(6),'valuechangedfcn',@Circulationfcn);
    set(PlotsValues(7),'valuechangedfcn',@InducedVelocityfcn);
    set(PlotsValues(8),'valuechangedfcn',@InflowAnglefcn);
    set(PlotsValues(9),'valuechangedfcn',@ExpandedBladeDefcn);
    set(PlotsValues(10),'valuechangedfcn',@BladeThicknessDefcn);
    set(PlotsValues(11),'valuechangedfcn',@Liftfcn);
    set(PlotsValues(12),'valuechangedfcn',@Performancefcn);
end
%按下“环量分布曲线”后的回调函数
function Circulationfcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate;
    global lc3_1 lc3_2;
    %% 设置图像面板标题
    set(PlotsPanels,'title',PlotsStrings{6});
    %% 设置相应按钮的可用状态
    %查找被按下的按钮数量
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %对按钮6进行操作后，若一个按钮也没有被按下
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
        %对按钮6进行操作后，若存在被按下的按钮
        %则除6外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(6),'enable','on');
    end
    %% 绘制曲线
    %按钮按下后才绘制图像
    if get(PlotsValues(6),'value')
        %避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.design.RC;
        XG = pt.design.G';
        %XXr是Xr的一个外拓序列，外拓了Rhub_oR和1两个边界值
        XXr = [Rhub_oR,Xr,1];
        %XXG为通过插值得到的XXr径向位置处的环量值，该插值方法为pchip
        XXG = pchip(Xr,XG,XXr);
        cla(coordinate);
        Xr_XG_line = plot(coordinate,XXr,XXG,...
                          'linewidth',2,...
                          'color',lc3_1);
        hold(coordinate,'on');
        grid(coordinate,'on');
        box(coordinate,'on');
        Xr_XG_point = scatter(coordinate,Xr,XG,...
                              'marker','o',...
                              'markeredgecolor',lc3_2,...
                              'markerfacecolor',lc3_1,...
                              'sizedata',25);
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        legend(Xr_XG_line,'τ/2πRpVs');
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'τ/2πRpVs','fontsize',16,'fontname','Times');
    else
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');
        legend(coordinate,'off');
        set(PlotsPanels,'title','点击左侧按钮显示对应图像');
    end
end
%按下“诱导速度分布曲线”后的回调函数
function InducedVelocityfcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate Xr_UASTAR_line Xr_UTSTAR_line Xr_VAI_line Xr_VTI_line Xr_UASTAR_point Xr_UTSTAR_point;
    global lc3_1 lc3_2 lc4_1 lc4_2;
    %% 设置图像面板标题
    set(PlotsPanels,'title',PlotsStrings{7});
    %% 设置相应按钮的可用状态
    %查找被按下的按钮数量
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %对按钮7进行操作后，若一个按钮也没有被按下
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
        %对按钮7进行操作后，若存在被按下的按钮
        %则除3、7外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(3),'enable','on');
        set(PlotsValues(7),'enable','on');
    end
    %% 绘制曲线
    %按钮按下后才绘制图像
    if get(PlotsValues(7),'value')
        %避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        Xr = pt.design.RC;
        UASTAR = pt.design.UASTAR;
        UTSTAR = pt.design.UTSTAR;
        %如果此时被按下的按钮不为1，则说明3也被按下
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        grid(coordinate,'on');
        box(coordinate,'on');
        Xr_UASTAR_line = plot(coordinate,Xr,UASTAR,...
                              'linewidth',2,...
                              'color',lc3_1);
        hold(coordinate,'on');
        Xr_UASTAR_point = scatter(coordinate,Xr,UASTAR,...
                                  'marker','o',...
                                  'markeredgecolor',lc3_2,...
                                  'markerfacecolor',lc3_1,...
                                  'sizedata',25);
        hold(coordinate,'on');
        Xr_UTSTAR_line = plot(coordinate,Xr,UTSTAR,...
                              'linewidth',2,...
                              'color',lc4_1);
        hold(coordinate,'on');
        Xr_UTSTAR_point = scatter(coordinate,Xr,UTSTAR,...
                                  'marker','o',...
                                  'markeredgecolor',lc4_2,...
                                  'markerfacecolor',lc4_1,...
                                  'sizedata',25);
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        if pushednum == 1
            legend([Xr_UASTAR_line,Xr_UTSTAR_line],'Ua*/Vs','Ut*/Vs');
        else
            legend([Xr_VAI_line,Xr_VTI_line,Xr_UASTAR_line,Xr_UTSTAR_line],...
                   'VA/Vs','VT/Vs','Ua*/Vs','Ut*/Vs');
            set(PlotsPanels,'title','速度分布曲线');
        end    
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'Velocity','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanels,'title','点击左侧按钮显示对应图像');
    else
        delete(Xr_UASTAR_line);
        delete(Xr_UTSTAR_line);
        delete(Xr_UASTAR_point);
        delete(Xr_UTSTAR_point);
        legend([Xr_VAI_line,Xr_VTI_line],'VA/Vs','VT/Vs');
        set(PlotsPanels,'title',PlotsStrings{3});
    end
end
%按下“相关角度的分布曲线”后的回调函数
function InflowAnglefcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate;
    global lc2_1 lc2_2 lc4_1 lc4_2;
    %% 设置图像面板标题
    set(PlotsPanels,'title',PlotsStrings{8});
    %% 设置相应按钮的可用状态
    %查找被按下的按钮数量
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %对按钮8进行操作后，若一个按钮也没有被按下
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
        %对按钮8进行操作后，若存在被按下的按钮
        %则除8外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(8),'enable','on');
    end
    %% 绘制曲线
    %按钮按下后才绘制图像
    if get(PlotsValues(8),'value')
        %避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.design.RC;
        XBeta = atand(pt.design.TANBC);
        XBetaI = atand(pt.design.TANBIC);
        %XXr是Xr的一个外拓序列，外拓了Rhub_oR和1两个边界值
        XXr = [Rhub_oR,Xr,1];
        %XXBeta和XXBetaI为通过插值得到的XXr径向位置处的角度值，该插值方法为spline
        XXBeta = spline(Xr,XBeta,XXr);
        XXBetaI = spline(Xr,XBetaI,XXr);
        cla(coordinate);
        grid(coordinate,'on');
        box(coordinate,'on');
        Xr_Beta_line = plot(coordinate,XXr,XXBeta,...
                            'linewidth',2,...
                            'color',lc2_1);
        hold(coordinate,'on');
        Xr_Beta_point = scatter(coordinate,Xr,XBeta,...
                                'marker','o',...
                                'markeredgecolor',lc2_1,...
                                'markerfacecolor',lc2_2,...
                                'sizedata',25);
        hold(coordinate,'on');
        Xr_XBetaI_line = plot(coordinate,XXr,XXBetaI,...
                              'linewidth',2,...
                              'color',lc4_1);
        hold(coordinate,'on');
        Xr_XBetaI_point = scatter(coordinate,Xr,XBetaI,...
                                  'marker','o',...
                                  'markeredgecolor',lc4_1,...
                                  'markerfacecolor',lc4_2,...
                                  'sizedata',25);
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        legend([Xr_Beta_line,Xr_XBetaI_line],'β','βi');
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'Angle','fontsize',16,'fontname','Times');
    else
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');
        legend(coordinate,'off');
        set(PlotsPanels,'title','点击左侧按钮显示对应图像');
    end
end
%按下“叶片伸张轮廓曲线(来自设计结果)”后的回调函数
function ExpandedBladeDefcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate Xr_XCpoDp_De_line Xr_XCpoDp_In_line Xr_XCpoDp_De_point;
    global lc3_1 lc3_2;
    %% 设置图像面板标题
    set(PlotsPanels,'title',PlotsStrings{9});
    %% 设置相应按钮的可用状态
    %查找被按下的按钮数量
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %对按钮9进行操作后，若一个按钮也没有被按下
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
        %对按钮9进行操作后，若存在被按下的按钮
        %则除1、9外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        if pt.input.Chord_flag == 0
            set(PlotsValues(1),'enable','on');
        end
        set(PlotsValues(9),'enable','on');
    end
    %% 绘制曲线
    %按钮按下后才绘制图像
    if get(PlotsValues(9),'value')
        %避免绘图时弹窗
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        set(0,'CurrentFigure',Fig_Main);
        Xr = pt.design.RC;
        XCpoDp = pt.design.CpoDp;
        XXr = Rhub_oR+(1-Rhub_oR)*(sin((0:60)*pi/(2*60)));
        XXCpoDp = InterpolateChord(Xr,XCpoDp,XXr);
        %如果此时被按下的按钮不为1，则说明1也被按下
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        grid(coordinate,'on');
        box(coordinate,'on');
        Xr_XCpoDp_De_line = plot(coordinate,XXr,XXCpoDp,...
                              'linewidth',2,...
                              'color',lc3_1);
        hold(coordinate,'on');
        Xr_XCpoDp_De_point = scatter(coordinate,Xr,XCpoDp,...
                                  'marker','o',...
                                  'markeredgecolor',lc3_1,...
                                  'markerfacecolor',lc3_2,...
                                  'sizedata',25);
        xlim(coordinate,'auto');
        ylim(coordinate,'auto')
        if pushednum == 1
            legend(Xr_XCpoDp_De_line,'Cp/Dp');
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
        delete(Xr_XCpoDp_De_line);
        delete(Xr_XCpoDp_De_point);
        legend(Xr_XCpoDp_In_line,'Cp/Dp');
        set(PlotsPanels,'title',PlotsStrings{1});
    end
end
%按下“叶片厚度分布曲线(来自设计结果)”后的回调函数
function BladeThicknessDefcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate Xr_Xt0oDp_De_line Xr_Xt0oDp_In_line Xr_Xt0oDp_De_point;
    global lc4_1 lc4_2;
    %% 设置图像面板标题
    set(PlotsPanels,'title',PlotsStrings{10});
    %% 设置相应按钮的可用状态
    %查找被按下的按钮数量
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %对按钮10进行操作后，若一个按钮也没有被按下
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
        %对按钮10进行操作后，若存在被按下的按钮
        %则除2、10外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(2),'enable','on');
        set(PlotsValues(10),'enable','on');
    end
    %% 绘制曲线
    %按钮按下后才绘制图像
    if get(PlotsValues(10),'value')
        %避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.design.RC;
        Xt0oDp = pt.design.t0oDp;
        %%XXr是Xr的一个插值序列，利用正弦函数拓展Xr序列
        XXr = Rhub_oR+(1-Rhub_oR)*(sin((0:60)*pi/(2*60)));
        %XXt0oDp为通过插值得到的XXr径向位置处的厚径比，该插值方法为pchip
        XXt0oDp = pchip(Xr,Xt0oDp,XXr);
        %如果此时被按下的按钮不为1，则说明2也被按下
        if pushednum == 1
            cla(coordinate);
        else
            hold(coordinate,'on');
        end
        grid(coordinate,'on');
        box(coordinate,'on');
        Xr_Xt0oDp_De_line = plot(coordinate,XXr,XXt0oDp,...
                              'linewidth',2,...
                              'color',lc4_1);
        hold(coordinate,'on');
        Xr_Xt0oDp_De_point = scatter(coordinate,Xr,Xt0oDp,...
                                  'marker','o',...
                                  'markeredgecolor',lc4_1,...
                                  'markerfacecolor',lc4_2,...
                                  'sizedata',25);
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        if pushednum == 1
            legend(Xr_Xt0oDp_De_line,'t0/Dp');
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
        delete(Xr_Xt0oDp_De_line);
        delete(Xr_Xt0oDp_De_point);
        legend(Xr_Xt0oDp_In_line,'Cp/Dp');
        set(PlotsPanels,'title',PlotsStrings{2});
    end
end
%按下“升力系数分布曲线”后的回调函数
function Liftfcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate;
    global lc1_1 lc1_2;
    %% 设置图像面板标题
    set(PlotsPanels,'title',PlotsStrings{11});
    %% 设置相应按钮的可用状态
    %查找被按下的按钮数量
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %对按钮11进行操作后，若一个按钮也没有被按下
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
        %对按钮11进行操作后，若存在被按下的按钮
        %则除11外的所有按钮均不可用
        set(PlotsValues,'enable','off');
        set(PlotsValues(11),'enable','on');
    end
    %% 绘制曲线
    %按钮按下后才绘制图像
    if get(PlotsValues(11),'value')
        %避免绘图时弹窗
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.design.RC;
        XCL = pt.design.CL;
        %XXr是Xr的一个外拓序列，外拓了Rhub_oR和1两个边界值
        XXr = [Rhub_oR,Xr,1];
        %XXG为通过插值得到的XXr径向位置处的环量值，该插值方法为spline
        XXCL = spline(Xr,XCL,XXr);
        cla(coordinate);
        Xr_CL_line = plot(coordinate,XXr,XXCL,...
                          'linewidth',2,...
                          'color',lc1_1);
        hold(coordinate,'on');
        grid(coordinate,'on');
        box(coordinate,'on');
        Xr_CL_point = scatter(coordinate,Xr,XCL,...
                              'marker','o',...
                              'markeredgecolor',lc1_1,...
                              'markerfacecolor',lc1_2,...
                              'sizedata',25);
        xlim(coordinate,'auto');
        ylim(coordinate,'auto');
        legend(Xr_CL_line,'CL');
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'CL','fontsize',16,'fontname','Times');
    else
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');
        legend(coordinate,'off');
        set(PlotsPanels,'title','点击左侧按钮显示对应图像');
    end
end
%按下“性能曲线”后的回调函数
function Performancefcn(hObject,ED)
    
end