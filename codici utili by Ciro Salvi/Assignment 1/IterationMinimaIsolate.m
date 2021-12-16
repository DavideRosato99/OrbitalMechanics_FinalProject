function DataMin = IterationMinimaIsolate(Data, PopSize)
% IterationMinimaIsolate isolates the minima of each ga iteration, given 
% the global population of data and the 'PopulationSize' setting used
%
% INPUT 
%	Data [nx4]          Matrix containing the whole ga population evaluated; 
%     Data(:, 1)        Mercury departure date  [mjd2000]
%     Data(:, 2)        Venus flyby date  [mjd2000]
%     Data(:, 3)        Jupiter arrival date  [mjd2000]
%     Data(:, 4)        Cost corresponding to the option [km/s]
%	PopSize[1]          'PopulationSize' setting used in ga simulator [-]
%
% OUTPUT
%   DataMin [n/PopSizex4] Matrix containing the best option of the population 
%                       evaluated at each iteration
% CONTRIBUTOR
%   Fabio Spada
% 
% VERSIONS
%   2021-02-11
%

nIt = size(Data, 1)/PopSize;
DataMin = zeros(nIt, 4);

for ii = 1:nIt
    DataIt = Data((ii-1)*PopSize+1 : ii*PopSize, :);
    idx = find(min(DataIt(:,4)));
    DataMin(ii, :) = DataIt(idx, :);
end
