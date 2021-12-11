function str = date2string(dateVec, format)
% PROTOTYPE:
%   str = date2string(dateVec, format)
% 
% DESCRIPTION:
%   Returns data in string form given the vector and the format chosen
% 
% INPUT:
% dateVec     vector of date and time: [year, month, day, hour, minutes, seconds]'
%                                  e.g. [2033, 05, 01, 00, 00, 00]'
% format     {date, time, datetime}
% 
% OUTPUT:
% str      date in sring form (YEAR MONTH DAY at HOUR:MINUTES:SECONDS) 
% 
% CALLED FUNCTIONS:
%   (none)

months = {'January', 'February', 'March', 'April', 'May', 'June', ...
	'July', 'August', 'September', 'October', 'November', 'December'};

if nargin < 2
    dateTime_str = sprintf('%i %s %i at %02i:%02i:%02i', dateVec(1), months{dateVec(2)}, dateVec(3), dateVec(4), dateVec(5), round(dateVec(6)));
    str = dateTime_str;
else
    switch format
        case 'date'
            str = sprintf('%i %s %i', dateVec(1), months{dateVec(2)}, dateVec(3));
        case 'time'
            str = sprintf('%02i:%02i:%02i', dateVec(4), dateVec(5), round(dateVec(6)));
        case 'datetime'
                str = dateTime_str;
        otherwise
                str = dateTime_str;
    end
end
end

