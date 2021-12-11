function writeMP4Movie(myMovie, fileName)
% PROTOTYPE:
%   writeMP4Movie(myMovie, fileName)
% 
% DESCRIPTION:
%   Generates MP4 movie
% 
% INPUT:
% Mymovie
% fileName     
% 
% OUTPUT:
% movie
% 
% CALLED FUNCTIONS:
%   (none)

playAfterCapture = 0;
fps = 24;


if length(myMovie) > 1
        writerObj = VideoWriter( fileName, 'MPEG-4');
        writerObj.Quality = 100;
        writerObj.FrameRate = fps;
        
        open( writerObj );        
   
        writeVideo( writerObj, myMovie );
        disp(sprintf('%s was written', writerObj.Filename))
        close( writerObj );

        
        if playAfterCapture
            fig = figure; % create new figure for playing the movie
            movie(fig,myMovie,2)
        end
   
else
    warning('ERROR writing movie!')
end

end

