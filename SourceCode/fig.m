function H = fig(H)
backgroundcolor = [0.709803 0.596078 0.631372];
if nargin == 1   
    figure(H,'menubar','figure','toolbar','figure','resize','on','color',backgroundcolor);
    hold on, box on, grid on,    
else
    H = figure('menubar','figure','toolbar','figure','resize','on','color',backgroundcolor);
    hold on, box on, grid on,
end