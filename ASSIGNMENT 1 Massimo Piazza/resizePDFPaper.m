function [myFig] = resizePdfPaper(myFig)
% PROTOTYPE:
%   resizePdfPaper(myFig)
% 
% DESCRIPTION:
%   Returns a resized pdf image 
% 
% INPUT:
% Myfig: figure handle   
% 
% OUTPUT:
% plot
% 
% CALLED FUNCTIONS:
%   (none)

    %set(myFig, 'Units', 'Normalized', 'OuterPosition', [0.25 0.25 .5 1]);
    myFig.PaperPositionMode = 'auto';
    fig_pos = myFig.PaperPosition;
    myFig.PaperSize = [fig_pos(3) fig_pos(4)];
end