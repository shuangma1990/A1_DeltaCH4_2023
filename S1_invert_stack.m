% this script generate cov matrix from Basu's ensemble of prior/posterior fluxes
% filedir='/Users/shuangma/RESEARCH/WORKFLOW/A1_DeltaCH4_2023/S1_output/';
% csvwrite([filedir,'calc_cov-inversion-full-apri_cov_fossil.csv'],B1);
% csvwrite([filedir,'calc_cov-inversion-full-apri_cov_microbial.csv'],B2);
% csvwrite([filedir,'calc_cov-inversion-full-apri_cov_pyrogenic.csv'],B3);
% 
% csvwrite([filedir,'calc_cov-inversion-full-apos_cov_fossil.csv'],B4);
% csvwrite([filedir,'calc_cov-inversion-full-apos_cov_microbial.csv'],B5);
% csvwrite([filedir,'calc_cov-inversion-full-apos_cov_pyrogenic.csv'],B6);

% The eight nc files are priors and posteriors of:
% CH?-only using Climatological prior
% CH? and delta13CH4 using Climatological prior
% CH?-only using year-specific prior
% CH? and delta13CH4 using year-specific prior

sectorname1 = 'fossil/emission';
sectorname2 = 'microbial/emission';
sectorname3 = 'pyrogenic/emission';

filename1 = 'calc_cov-inversion-full-apri.nc'; %'year' specific
filename2 = 'calc_cov-inversion-full-apos.nc';
filename3 = 'calc_cov-inversion-full-clim_prior-apri.nc';
filename4 = 'calc_cov-inversion-full-clim_prior-apos.nc';
threshold1 = 0.262;
threshold2 = 2;
[A1, B1, C1, C1_block] = check_inv_stack(filename1, sectorname1, threshold1, threshold2);
% [A1, B1, C1, emission_percentile1] = check_inv_stack(filename1, sectorname1, threshold);
% threshold = 0.75;
% [A2, B2, C2, emission_percentile2] = check_inv(filename1, sectorname2, threshold);
% threshold = 0.7;
% [A3, B3, C3, emission_percentile3] = check_inv(filename1, sectorname3, threshold);
% 
% threshold = 0.1;
% [A4, B4, C4, emission_percentile4] = check_inv(filename2, sectorname1, threshold);
% threshold = 0.85;
% [A5, B5, C5, emission_percentile5] = check_inv(filename2, sectorname2, threshold);
% threshold = 0.055;
% [A6, B6, C6, emission_percentile6] = check_inv(filename2, sectorname3, threshold);
% 
% 
% threshold = 0.26;% INF
% [A7, B7, C7, emission_percentile7] = check_inv(filename3, sectorname1, threshold);
% threshold = 0.75;
% [A8, B8, C8, emission_percentile8] = check_inv(filename3, sectorname2, threshold);
% threshold = 0.05;
% [A9, B9, C9, emission_percentile9] = check_inv(filename3, sectorname3, threshold);
% 
% 
% threshold = 1.1; % INF
% [A10, B10, C10, emission_percentile10] = check_inv(filename4, sectorname1, threshold);
% threshold = 0.75;
% [A11, B11, C11, emission_percentile11] = check_inv(filename4, sectorname2, threshold);
% threshold = 0.05;
% [A12, B12, C12, emission_percentile12] = check_inv(filename4, sectorname3, threshold);

%%
filename1 = 'calc_cov-inversion-full-ch4only-apri.nc'; %'year' specific
filename2 = 'calc_cov-inversion-full-ch4only-apos.nc';
filename3 = 'calc_cov-inversion-full-ch4only-clim_prior-apri.nc';
filename4 = 'calc_cov-inversion-full-ch4only-clim_prior-apos.nc';

sectorname1 = 'fossil/emission';
sectorname2 = 'microbial/emission';
sectorname3 = 'pyrogenic/emission';

threshold = 0.262;
[D1, E1, F1, emission_percentile1] = check_inv(filename1, sectorname1, threshold);
threshold = 0.75;
[D2, E2, F2, emission_percentile2] = check_inv(filename1, sectorname2, threshold);
threshold = 0.06;
[D3, E3, F3, emission_percentile3] = check_inv(filename1, sectorname3, threshold);


threshold = 3.3;%2;%3.3;
[D4, E4, F4, emission_percentile4] = check_inv(filename2, sectorname1, threshold);
threshold = 0.85;
[D5, E5, F5, emission_percentile5] = check_inv(filename2, sectorname2, threshold);
threshold = 0.055;
[D6, E6, F6, emission_percentile6] = check_inv(filename2, sectorname3, threshold);


threshold = 0.26;% INF
[D7, E7, F7, emission_percentile7] = check_inv(filename3, sectorname1, threshold);
threshold = 0.75;
[D8, E8, F8, emission_percentile8] = check_inv(filename3, sectorname2, threshold);
threshold = 0.05;
[D9, E9, F9, emission_percentile9] = check_inv(filename3, sectorname3, threshold);


threshold = 1.7; % 
[D10, E10, F10, emission_percentile10] = check_inv(filename4, sectorname1, threshold);
threshold = 0.75;
[D11, E11, F11, emission_percentile11] = check_inv(filename4, sectorname2, threshold);
threshold = 0.05;
[D12, E12, F12, emission_percentile12] = check_inv(filename4, sectorname3, threshold);

