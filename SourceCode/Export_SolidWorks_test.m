function [] = Export_SolidWorks_test(filename_SolidWorks,Np,Mp,Z,X3D,Y3D,Z3D)
% 返回矩阵filename_SolidWorks的列数
filename_length = size(filename_SolidWorks, 2);
% 删去filename_SolidWorks后四位的'.txt'
foldername = filename_SolidWorks(1:filename_length-4);
% 新建以foldername为名的文件夹
mkdir(foldername);
% 目录切换至foldername中
cd(foldername);
%fid = fopen(filename_SolidWorks,'wt');  % 1-Aug-2013 BEPPS: 'w' changed to 'wt'  so newline character '\n' appears properly on Windows machines

% Prop Parameters at beginning of file
%fprintf(fid,'%g, ' ,Np); 
%fprintf(fid,'%g, ' ,Mp);
%fprintf(fid,'%g,\n',Z);

% Output curves defining each 2D section along the span
% for each section along the span
for i = 1:Mp+1   
    section_curve_name = strcat('SectionCurveTest',num2str(i),'.txt');
    sc_fid = fopen(section_curve_name,'wt');
    %fprintf(fid,'SectionCurve%g, \n',i);    

    % for each point along the suction and pressure surfaces
    % (trailing edge -> leading edge -> trailing edge, close the curve)
    for j = [1:Np,Np+2:2*Np-1,1] % (2*Np-1 points) does not double print the leading edge
        fprintf(sc_fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
    end
end