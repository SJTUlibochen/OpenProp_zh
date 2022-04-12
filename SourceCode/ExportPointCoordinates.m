%叶片点坐标：
%X3D(i,j,k)     [m]，X坐标
%Y2D(i,j,k)     [m]，Y坐标
%Z3D(i,j,k)     [m]，Z坐标
%
%i = 1:Mp+1      径向分区线
%j = 1:2*Np      周向分区线
%k = 1:Z         叶片数量
function [] = ExportPointCoordinates(filename,Np,Mp,X3D,Y3D,Z3D)
    foldername = [filename,'_Coordinates'];
    mkdir(foldername);
    cd(foldername);
    %% 输出全部径向分区线
    %输出的径向分区线用于作为放样操作的截面
    for i = 1:Mp+1   
        section_curve_name = strcat('SectionCurve',num2str(i),'.csv');
        sc_fid = fopen(section_curve_name,'wt');
        for j = [1:Np,Np+2:2*Np,1]
            %Np和Np+1两点坐标相同，只保留Np的坐标
            %最后回到1点坐标，以形成闭环
            fprintf(sc_fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
        end
        fclose(sc_fid);
    end
    %% 输出部分周向分区线
    %输出的周向分区线用于作为放样操作的轨迹
    %输出周向分区线的选择标准：1 floor((Np+1)/2) Np Np+floor((Np+1)/2)
    n = 0;
    for j = [1 floor((Np+1)/2) Np Np+floor((Np+1)/2)]
        n = n + 1;
        guide_curve_name = strcat('GuideCurve',num2str(n),'.csv');
        gc_fid = fopen(guide_curve_name,'wt');
        for i = 1:Mp+1
            fprintf(gc_fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
        end
        fclose(gc_fid);
    end
end