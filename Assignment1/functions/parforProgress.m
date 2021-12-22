function parforProgress(N)
% parfor_progress - Progress that works with parfor. It works by creating
% a file called parfor_progress.txt in your working directory, and then 
% keeping track of the parfor loop's progress within that file.
%
% PROTOTYPE
%   parfor_progress(N) initializes the progress monitor for a set of N
%   upcoming calculations.
%
%   parfor_progress updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   parfor_progress(0) deletes parfor_progress.txt and finalizes progress
%   bar.
%
%   dx=ode_2bp(t,x,mu,typeSim,date0) will perform a computation of the
%   orbit with J2 and Moon perturbations using the specified derivatives
%   equations in typeSim input.
%
% CALLED FUNCTIONS: -
%
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli Jason              10618465
%   Tagliati Alessia            10635119
%
% VERSIONS
%   2021-10-21: Release
%
% -------------------------------------------------------------------------

narginchk(0, 1);

if nargin < 1
    N = -1;
end

if N > 0
    f = fopen('parfor_progress.txt', 'w');
    if f<0
        error('Unable to open parfor_progress.txt');
    end
    fprintf(f, '%d\n', N); % Save N at the top of progress.txt
    fclose(f);
    
    if nargout == 0
        disp('0%');
    end
elseif N == 0
    delete('parfor_progress.txt');
else
    if ~exist('parfor_progress.txt', 'file')
        error('parfor_progress.txt not found. Run parfor_progress(N) before parfor_progress to initialize parfor_progress.txt.');
    end
    
    f = fopen('parfor_progress.txt', 'a');
    fprintf(f, '1\n');
    fclose(f);
    
    f = fopen('parfor_progress.txt', 'r');
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1)*100;
    
    if nargout == 0
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        disp(perc);
    end
end