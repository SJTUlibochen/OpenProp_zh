%������OpenPropSingle.m��û�����õİ�ť�ص��������Լ���Ӧ������
function Make_Reports
    global PlotsValues;
    %���ð�ť�Ļص�����
    set(PlotsValues(6),'valuechangedfcn',@Circulationfcn);
    set(PlotsValues(7),'valuechangedfcn',@InducedVelocityfcn);
    set(PlotsValues(8),'valuechangedfcn',@InflowAnglefcn);
    set(PlotsValues(9),'valuechangedfcn',@ExpandedBladeDefcn);
    set(PlotsValues(10),'valuechangedfcn',@BladeThicknessDefcn);
    set(PlotsValues(11),'valuechangedfcn',@Liftfcn);
    set(PlotsValues(12),'valuechangedfcn',@Performancefcn);
end
%���¡������ֲ����ߡ���Ļص�����
function Circulationfcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate;
    global lc3_1 lc3_2;
    %% ����ͼ��������
    set(PlotsPanels,'title',PlotsStrings{6});
    %% ������Ӧ��ť�Ŀ���״̬
    %���ұ����µİ�ť����
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %�԰�ť6���в�������һ����ťҲû�б�����
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
        %�԰�ť6���в����������ڱ����µİ�ť
        %���6������а�ť��������
        set(PlotsValues,'enable','off');
        set(PlotsValues(6),'enable','on');
    end
    %% ��������
    %��ť���º�Ż���ͼ��
    if get(PlotsValues(6),'value')
        %�����ͼʱ����
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.design.RC;
        XG = pt.design.G';
        %XXr��Xr��һ���������У�������Rhub_oR��1�����߽�ֵ
        XXr = [Rhub_oR,Xr,1];
        %XXGΪͨ����ֵ�õ���XXr����λ�ô��Ļ���ֵ���ò�ֵ����Ϊpchip
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
        legend(Xr_XG_line,'��/2��RpVs');
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'��/2��RpVs','fontsize',16,'fontname','Times');
    else
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');
        legend(coordinate,'off');
        set(PlotsPanels,'title','�����ఴť��ʾ��Ӧͼ��');
    end
end
%���¡��յ��ٶȷֲ����ߡ���Ļص�����
function InducedVelocityfcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate Xr_UASTAR_line Xr_UTSTAR_line Xr_VAI_line Xr_VTI_line Xr_UASTAR_point Xr_UTSTAR_point;
    global lc3_1 lc3_2 lc4_1 lc4_2;
    %% ����ͼ��������
    set(PlotsPanels,'title',PlotsStrings{7});
    %% ������Ӧ��ť�Ŀ���״̬
    %���ұ����µİ�ť����
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %�԰�ť7���в�������һ����ťҲû�б�����
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
        %�԰�ť7���в����������ڱ����µİ�ť
        %���3��7������а�ť��������
        set(PlotsValues,'enable','off');
        set(PlotsValues(3),'enable','on');
        set(PlotsValues(7),'enable','on');
    end
    %% ��������
    %��ť���º�Ż���ͼ��
    if get(PlotsValues(7),'value')
        %�����ͼʱ����
        set(0,'CurrentFigure',Fig_Main);
        Xr = pt.design.RC;
        UASTAR = pt.design.UASTAR;
        UTSTAR = pt.design.UTSTAR;
        %�����ʱ�����µİ�ť��Ϊ1����˵��3Ҳ������
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
            set(PlotsPanels,'title','�ٶȷֲ�����');
        end    
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'Velocity','fontsize',16,'fontname','Times');
    elseif pushednum == 0
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');     
        legend(coordinate,'off');
        set(PlotsPanels,'title','�����ఴť��ʾ��Ӧͼ��');
    else
        delete(Xr_UASTAR_line);
        delete(Xr_UTSTAR_line);
        delete(Xr_UASTAR_point);
        delete(Xr_UTSTAR_point);
        legend([Xr_VAI_line,Xr_VTI_line],'VA/Vs','VT/Vs');
        set(PlotsPanels,'title',PlotsStrings{3});
    end
end
%���¡���ؽǶȵķֲ����ߡ���Ļص�����
function InflowAnglefcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate;
    global lc2_1 lc2_2 lc4_1 lc4_2;
    %% ����ͼ��������
    set(PlotsPanels,'title',PlotsStrings{8});
    %% ������Ӧ��ť�Ŀ���״̬
    %���ұ����µİ�ť����
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %�԰�ť8���в�������һ����ťҲû�б�����
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
        %�԰�ť8���в����������ڱ����µİ�ť
        %���8������а�ť��������
        set(PlotsValues,'enable','off');
        set(PlotsValues(8),'enable','on');
    end
    %% ��������
    %��ť���º�Ż���ͼ��
    if get(PlotsValues(8),'value')
        %�����ͼʱ����
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.design.RC;
        XBeta = atand(pt.design.TANBC);
        XBetaI = atand(pt.design.TANBIC);
        %XXr��Xr��һ���������У�������Rhub_oR��1�����߽�ֵ
        XXr = [Rhub_oR,Xr,1];
        %XXBeta��XXBetaIΪͨ����ֵ�õ���XXr����λ�ô��ĽǶ�ֵ���ò�ֵ����Ϊspline
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
        legend([Xr_Beta_line,Xr_XBetaI_line],'��','��i');
        xlabel(coordinate,'r/Rp','fontsize',16,'fontname','Times');
        ylabel(coordinate,'Angle','fontsize',16,'fontname','Times');
    else
        cla(coordinate);
        xlabel(coordinate,'');
        ylabel(coordinate,'');
        legend(coordinate,'off');
        set(PlotsPanels,'title','�����ఴť��ʾ��Ӧͼ��');
    end
end
%���¡�ҶƬ������������(������ƽ��)����Ļص�����
function ExpandedBladeDefcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate Xr_XCpoDp_De_line Xr_XCpoDp_In_line Xr_XCpoDp_De_point;
    global lc3_1 lc3_2;
    %% ����ͼ��������
    set(PlotsPanels,'title',PlotsStrings{9});
    %% ������Ӧ��ť�Ŀ���״̬
    %���ұ����µİ�ť����
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %�԰�ť9���в�������һ����ťҲû�б�����
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
        %�԰�ť9���в����������ڱ����µİ�ť
        %���1��9������а�ť��������
        set(PlotsValues,'enable','off');
        if pt.input.Chord_flag == 0
            set(PlotsValues(1),'enable','on');
        end
        set(PlotsValues(9),'enable','on');
    end
    %% ��������
    %��ť���º�Ż���ͼ��
    if get(PlotsValues(9),'value')
        %�����ͼʱ����
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
        %�����ʱ�����µİ�ť��Ϊ1����˵��1Ҳ������
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
        delete(Xr_XCpoDp_De_line);
        delete(Xr_XCpoDp_De_point);
        legend(Xr_XCpoDp_In_line,'Cp/Dp');
        set(PlotsPanels,'title',PlotsStrings{1});
    end
end
%���¡�ҶƬ��ȷֲ�����(������ƽ��)����Ļص�����
function BladeThicknessDefcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate Xr_Xt0oDp_De_line Xr_Xt0oDp_In_line Xr_Xt0oDp_De_point;
    global lc4_1 lc4_2;
    %% ����ͼ��������
    set(PlotsPanels,'title',PlotsStrings{10});
    %% ������Ӧ��ť�Ŀ���״̬
    %���ұ����µİ�ť����
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %�԰�ť10���в�������һ����ťҲû�б�����
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
        %�԰�ť10���в����������ڱ����µİ�ť
        %���2��10������а�ť��������
        set(PlotsValues,'enable','off');
        set(PlotsValues(2),'enable','on');
        set(PlotsValues(10),'enable','on');
    end
    %% ��������
    %��ť���º�Ż���ͼ��
    if get(PlotsValues(10),'value')
        %�����ͼʱ����
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.design.RC;
        Xt0oDp = pt.design.t0oDp;
        %%XXr��Xr��һ����ֵ���У��������Һ�����չXr����
        XXr = Rhub_oR+(1-Rhub_oR)*(sin((0:60)*pi/(2*60)));
        %XXt0oDpΪͨ����ֵ�õ���XXr����λ�ô��ĺ񾶱ȣ��ò�ֵ����Ϊpchip
        XXt0oDp = pchip(Xr,Xt0oDp,XXr);
        %�����ʱ�����µİ�ť��Ϊ1����˵��2Ҳ������
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
        delete(Xr_Xt0oDp_De_line);
        delete(Xr_Xt0oDp_De_point);
        legend(Xr_Xt0oDp_In_line,'Cp/Dp');
        set(PlotsPanels,'title',PlotsStrings{2});
    end
end
%���¡�����ϵ���ֲ����ߡ���Ļص�����
function Liftfcn(hObject,ED)
    global PlotsPanels PlotsStrings;
    global PlotsValues;
    global pt;
    global Fig_Main coordinate;
    global lc1_1 lc1_2;
    %% ����ͼ��������
    set(PlotsPanels,'title',PlotsStrings{11});
    %% ������Ӧ��ť�Ŀ���״̬
    %���ұ����µİ�ť����
    [pushedlist,pushednum] = PushedFind;
    if pushednum == 0
        %�԰�ť11���в�������һ����ťҲû�б�����
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
        %�԰�ť11���в����������ڱ����µİ�ť
        %���11������а�ť��������
        set(PlotsValues,'enable','off');
        set(PlotsValues(11),'enable','on');
    end
    %% ��������
    %��ť���º�Ż���ͼ��
    if get(PlotsValues(11),'value')
        %�����ͼʱ����
        set(0,'CurrentFigure',Fig_Main);
        Dp = pt.input.Dp;
        Rp = Dp/2;
        Dhub = pt.input.Dhub;
        Rhub = Dhub/2;
        Rhub_oR = Rhub/Rp;
        Xr = pt.design.RC;
        XCL = pt.design.CL;
        %XXr��Xr��һ���������У�������Rhub_oR��1�����߽�ֵ
        XXr = [Rhub_oR,Xr,1];
        %XXGΪͨ����ֵ�õ���XXr����λ�ô��Ļ���ֵ���ò�ֵ����Ϊspline
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
        set(PlotsPanels,'title','�����ఴť��ʾ��Ӧͼ��');
    end
end
%���¡��������ߡ���Ļص�����
function Performancefcn(hObject,ED)
    
end