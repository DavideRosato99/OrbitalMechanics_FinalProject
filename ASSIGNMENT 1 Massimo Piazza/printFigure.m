function [] = printFigure(h, fileName, fileType, bitmapRes, ppi)
% PROTOTYPE:
%   printFigure(h, fileName, fileType, bitmapRes, ppi)
% 
% DESCRIPTION:
%   Automatic figure printing
% 
% INPUT:
% h: figure handle
% fileName
% fileType: {pdf, eps, svg, jpg, png}
% bitmapRes: 
% ppi:  pixel-per-inch
% 
% OUTPUT:
% plot
% 
% CALLED FUNCTIONS:
%   (none)

if nargin > 2
    switch fileType
        case 'pdf'
            fileTypeID = '-dpdf';
        case 'eps'
            fileTypeID = '-depsc';
        case 'svg'
            fileTypeID = '-dsvg';
        case 'jpg'
            fileTypeID = '-djpeg';
        case 'png'
            fileTypeID = '-dpng';
        otherwise
            fileType = 'pdf';
            fileTypeID = '-dpdf';
    end
else
    fileType = 'pdf';
    fileTypeID = '-dpdf';
end


% If bitmap resolution is specified (for non-vector file types)
if (nargin >= 5) && (fileType == "jpg" || fileType == "png")
    h.PaperPositionMode = 'auto';
    fig_pos = h.PaperPosition;
    set(h, 'PaperUnits', 'inches', 'PaperPosition', bitmapRes/ppi*[0 0 1 fig_pos(4)/fig_pos(3)]);
    print(h,'-dpng',sprintf('-r%d',ppi), sprintf('%s.%s', fileName, fileType));   
else % nominal behavior
    print(resizePDFPaper(h), sprintf('%s.%s', fileName, fileType), fileTypeID)

end

end

