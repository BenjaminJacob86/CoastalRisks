%function prepare_waterlevel
% This code extract integrated water levels from SCHISM-WWM output and to prepare input files for Xbeach
% Note that Xbeach x-axis is alwyas oriented towards the coast, y-axis is
% along the coast
% WWM data are produced on levante with extract4xbeach2.py
% author: Wei Chen (wei.chen@hereon.de)

clear all
close all

%  input data--------------------------------------------------------------
%  specify the directory, station name and extract period
input_dir = '/gpfs/work/chenw1/Xbeach_wei/jonswap_file/';
%input_dir = '/work/gg0028/g260126/Xbeach_wei/';
output_dir= '/gpfs/work/chenw1/Xbeach_wei/jonswap_file/'; 

% read parameters and write into xbeach input files
yr_id = '2090';
    var_name = {['Zeta','_',yr_id]};
    file_var_ID = fopen([input_dir,var_name{1},'.txt'],'r');
    indata{1} = textscan(file_var_ID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines',1);
    fclose(file_var_ID);

    Zeta = indata{1,1};
    time = 3600*(1:length(Zeta{1}))-3600;
   % ele_nr = [time', Zeta{4}, Zeta{2}];  % inside Ems
    ele_nr = [time', Zeta{8}, Zeta{13}];  % TI_31
   % ele_nr = [time', Zeta{14}, Zeta{18}];  % TI_41
   %     ele_nr = [time', Zeta{19}, Zeta{21}]; % Sci. Total Environ 5km-4km domain
   %     ele_nr = [time', Zeta{22}, Zeta{24}, Zeta{25}, Zeta{26}]; % Sci. Total Environ 12km-4km domain
            
    % use arcmap to see the location of the domain
    % 2 to 4 are for Xbeach domain inside Ems
    % 8 to 13 are for Xbeach domain 31
    % 14 to 18 are for Xbeach domain 41
    file_write_var_ID = fopen([output_dir, ['waterlevel_TI_31','_2p_',yr_id,'.txt']],'w');
    fprintf(file_write_var_ID, '%f %f %f\n', ele_nr.');
    fclose(file_write_var_ID);

