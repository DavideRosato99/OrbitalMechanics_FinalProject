function [mon, day, hr, minute, sec] = days2mdh(year, days)
% days2mdh - the function converts a fraction of year into month, day,
%            hours, minute and second
%
% INPUT:
%   year     double  [1x1]   year                                 [y]
%   days     double  [1x1]   day + fraction of the day            [d]
%
% OUTPUT:
%    mon     double  [1x1]   number of months integer greater or equal to 
%                            1 and lower or equal to 12.          [m]
%    day     double  [1x1]   number of days integer greater or equal to 1
%                            and lower or equal to the maximun days of the
%                            month considered.                    [d]
%    hr      double  [1x1]   Number of hours as integer greater or equal to
%                            0 and lower or equal to 23.          [h]
%    minute  double  [1x1]   Number of minutes as integer greater or equal 
%                            to 0 and lower or equal to 59.       [m]
%    sec     double  [1x1]   Number of seconds as a real greater or equal 
%                            to 0 and strictly lower than 60.     [s]    
% CONTRIBUTORS:
%   Rosato Davide               10618468
%   Saba Mohammadi Yengeje      10789462
%   Spinelli Jason              10618465
%   Tagliati Alessia            10635119
%
% VERSIONS
%   2021-12-17
% -------------------------------------------------------------------------

% Set up array od days in month
lMonth = [31 28 31 30 31 30 31 31 30 31 30 31];
dayOfYr = floor(days);

% Find month and day of the month
if rem(year - 1900, 4) == 0
    lMonth(2)= 29;
end

i = 1;
intTemp = 0;
while ( dayOfYr > intTemp + lMonth(i) ) && ( i < 12 )
    intTemp = intTemp + lMonth(i);
    i = i+1;
end

mon = i;
day = dayOfYr - intTemp;

% Find hours minutes and seconds 
temp = (days - dayOfYr )*24.0;
hr  = floor(temp);
temp = (temp - hr) * 60.0;
minute = fix(temp);
sec = floor((temp - minute) * 60.0);

end