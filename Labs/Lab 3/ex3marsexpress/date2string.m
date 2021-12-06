function str = date2string(dateVec)

months = {'January', 'February', 'March', 'April', 'May', 'June', ...
	'July', 'August', 'September', 'October', 'November', 'December'};
str = sprintf('%i %s %i at %02i:%02i:%02i', dateVec(1), months{dateVec(2)}, dateVec(3), dateVec(4), dateVec(5), round(dateVec(6)));

end

