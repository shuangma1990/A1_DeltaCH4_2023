% threshold = 1.1;
% [A10, B10, C10, emission_percentile10] = check_inv(filename4, sectorname1, threshold);


function [A, B, C, emission_percentile] = check_inv(filename, sectorname, threshold)
filedir='/Users/shuangma/RESEARCH/DATA/Basu_isotope_GOSAT_shared_with_Jworden/';
% read lat lon
lat_array = ncread([filedir 'calc_cov-inversion-full-apri.nc'],'latitude');
lon_array = ncread([filedir 'calc_cov-inversion-full-apri.nc'],'longitude');
% lat_matrix = repmat(lat_array,1,120); % 90*120
% lon_matrix = repmat(lon_array',90,1); % 90*120
GC=GEOSChem_xygrids_local('2x3'); % _local added 2x3 option
lat = reshape(GC.y,[],1); % vertical first
lon = reshape(GC.x,[],1); % vertical first
area= reshape(GC.area,[],1); % vertical first
% read data, unit is kg CH4/grid cell/second, no need to use area
data = ncread([filedir filename],sectorname);
size(data);
% convert data to a 90 by 120 by 240 by 100 matrix
data = permute(data, [2 1 3 4]);
size(data);
% prior cov should have 90 * 120 rows/columns before filtering out 2% of
% low emission pixels
% and that's how many vars I have
% my A B C variables should be 100 element long vectors
% % I'll name my variables as Vi, where i=1:[90 * 120]
% % I can use mat2array to get Vi, so what I need now is the 3D matrix and
% % know how the 3D elements are read


% take average for 240 months to get mean
data_tsmean = squeeze(nanmean(data,3));
% take ensemble mean to get 90*120 matrix
data_tsmean_enmean = nanmean(data_tsmean,3);
plotglobal(data_tsmean_enmean);
% convert to array
data = reshape(data_tsmean_enmean,[],1); % vertical first

% combine columns
A = [lat lon area data];
% total emission is
total_emission = sum(data);
% specify a minimum value of pixel level CH4 emission that want to ignore
TF = A(:,4) == 0;
TF = A(:,4) < threshold;
% remove rows TF
A(TF,:) = [];

criteria_amount = sum(A(:,4));
disp([sectorname filename]);
emission_percentile  = criteria_amount/total_emission
% now I know which (TF) gridcells (rows) to delete before calculating the final
% cov matrix
% get the arrays ready
size(data_tsmean);
data_tsmean_2D = reshape(data_tsmean, [], 100);
size(data_tsmean_2D);
% remove rows TF
data_tsmean_2D(TF,:) = [];
% data_tsmean_2D is the matrix I need for calculating cov
B = cov(data_tsmean_2D');
% B is 1349*1349
% see if B is invertible
C = inv(B);




%% sample for covariance matrix calculation, if final matrix has n columns as ensemble numbers, and m rows as number of pixels, need to traspose before using cov() function
% A = [3 1 4 5];
% B = [70 120 90 70];
% C = [5 2 1 4];
% D = [A;B;C]; % need to transpose the D matrix to get cov matrix;
% cov(D')
% % cov(D') is the way to calcualte cov matrix for multiple variables
% % should be made of these two var cov elements:
% cov(A,B)
% cov(B,C)
% cov(A,C)
% % D' is what I need, it's made of A B C, each as columns

%% % example of deleting rows if column value satisfy criterias
% A = [0 2 3 4 5 6 ; 11 0 0 0 15 16 ; 21 0 0 0 0 26 ;  31 0 33 34 35 36 ; 41 42 43 0 0 0]
% % Specify you conditions
% TF1 = A(:,1)==0 
% TF2 = all(A(:,2:5)==0,2) & A(:,6) ~= 0 
% TF6 = A(:,2) == 0 & A(:,3) ~= 0 
% % combine them
% TFall = TF1 & TF2 & TF6
% % remove
% A(TFall,:) = []
end
