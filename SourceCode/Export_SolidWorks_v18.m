% -------------------------------------------------------------------------
% Created: Brenden Epps, 8/12/10
%
% Make SolidWorks.txt files, with coordinates for a single blade.
%
% Use this with SolidWorks macro v18.
%
% Blade geometry:
%   X3D(i,j,k) [m], X position in 3D space
%   Y2D(i,j,k) [m], Y position in 3D space
%   Z3D(i,j,k) [m], Z position in 3D space
%
%   i = 1:Mp+1      % for each section along the span
%   j = 1:2*Np      % for each point   along the upper and lower surfaces
%   k = 1:Z         % for each blade
%
% -------------------------------------------------------------------------

function [] = Export_SolidWorks_v18(filename_SolidWorks,Np,Mp,Z,X3D,Y3D,Z3D)

% 返回矩阵filename_SolidWorks的列数
filename_length = size(filename_SolidWorks, 2);
% 删去filename_SolidWorks后四位的'.txt'
foldername = filename_SolidWorks(1:filename_length-4);
% 新建以foldername为名的文件夹
mkdir(foldername);
% 目录切换至foldername中
cd(foldername);

% Output curves defining each 2D section along the span
% for each section along the span
for i = 1:Mp+1   
    section_curve_name = strcat('SectionCurve',num2str(i),'.txt');
    sc_fid = fopen(section_curve_name,'wt'); 
    % for each point along the suction and pressure surfaces
    % (trailing edge -> leading edge -> trailing edge, close the curve)
    for j = [1:Np,Np+2:2*Np-1,1] % (2*Np-1 points) does not double print the leading edge
        fprintf(sc_fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
    end
    fclose(sc_fid);
end

% Make guide curves
n = 0;   
% for 7 points along the chord
for j = [1 floor(Np/3) floor(2*Np/3) Np floor(4*Np/3) floor(5*Np/3) 2*Np-1]
    n = n + 1;
    guide_curve_name = strcat('GuideCurve',num2str(n),'.txt');
    gc_fid = fopen(guide_curve_name,'wt');
    %fprintf(fid,'GuideCurve%g, \n',n);
    for i = 1:Mp+1  % for each section along the span
        fprintf(gc_fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
    end
    fclose(gc_fid);
end

end % function