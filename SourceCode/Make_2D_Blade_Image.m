%该函数用于绘制2D叶切面轮廓
function [] = Make_2D_Blade_Image(RG,x2Dr,y2Dr)
    Mp = size(x2Dr,1)-1;
    Np = size(x2Dr,2)/2;
    if Geometry_flag 
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(13));
            axes(h);
    else    
        Fig2_S = figure('units','normalized','position',[0.31 .06 .4 .3],'name','Blade Image','numbertitle','off');
    end
                hold on;
                axis equal;     
                grid on;
                box on,
                title('2D Blade Image');  
                xlabel('X (2D) [m]');  
                ylabel('Y (2D) [m]');

                style      = ['r' 'g' 'b' 'm' 'k'];
                str_prefix = {'r/R = '};            

    flag = 1;

    plot(  [min(x2Dr),max(x2Dr)],0*[min(y2Dr),max(y2Dr)],'k')
    plot(0*[min(x2Dr),max(x2Dr)],  [min(y2Dr),max(y2Dr)],'k')

    for i = 1:ceil(Mp/5):Mp     % for five radial sections from root to tip
        handle_legend(flag) = plot(x2Dr(i,:),y2Dr(i,:),style(flag),'linewidth',2);


    %     plot(x2Dr(i,1:Np),y2Dr(i,1:Np),'g','linewidth',2)
    %     plot(x2Dr(i,Np+1:Np+Np),y2Dr(i,Np+1:Np+Np),'r','linewidth',2)


        plot([0.5*(x2Dr(i,1)+x2Dr(i,2*Np)),0.5*(x2Dr(i,Np)+x2Dr(i,Np+1))],[0.5*(y2Dr(i,1)+y2Dr(i,2*Np)),0.5*(y2Dr(i,Np)+y2Dr(i,Np+1))],style(flag),'linewidth',1);

        str_legend(flag) = strcat(str_prefix,num2str(RG(i)));

        flag = flag+1;
end

legend(handle_legend,str_legend,'location','northwest');
