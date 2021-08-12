%ҶƬ�����꣺
%X3D(i,j,k)     [m]��X����
%Y2D(i,j,k)     [m]��Y����
%Z3D(i,j,k)     [m]��Z����
%
%i = 1:Mp+1      ���������
%j = 1:2*Np      ���������
%k = 1:Z         ҶƬ����
function [] = ExportPointCoordinates(filename,Np,Mp,X3D,Y3D,Z3D)
    foldername = [filename,'_Coordinates'];
    mkdir(foldername);
    cd(foldername);
    %% ���ȫ�����������
    %����ľ��������������Ϊ���������Ľ���
    for i = 1:Mp+1   
        section_curve_name = strcat('SectionCurve',num2str(i),'.csv');
        sc_fid = fopen(section_curve_name,'wt');
        for j = [1:Np,Np+2:2*Np,1]
            %Np��Np+1����������ͬ��ֻ����Np������
            %���ص�1�����꣬���γɱջ�
            fprintf(sc_fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
        end
        fclose(sc_fid);
    end
    %% ����������������
    %��������������������Ϊ���������Ĺ켣
    %�����������ߵ�ѡ���׼��1 floor((Np+1)/2) Np Np+floor((Np+1)/2)
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