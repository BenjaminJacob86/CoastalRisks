%function prepare_jonswap
% This code extract integrated wave parameters (e.g. Hs Tp) from WWM output and to prepare input files for Xbeach
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

yr_id = '2090';
% read parameters and write into xbeach input files
for i = 1:4
    var_name = {['HS_',yr_id]; ['TP_',yr_id]; ['Dir_',yr_id]; ['DSP_',yr_id]};
    file_var_ID = fopen([input_dir,var_name{i},'.txt'],'r');
    indata{i} = textscan(file_var_ID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines',1);
    fclose(file_var_ID);
end
    % column 1 is time
    % column 14~18 are for Xbeach domain 41

    HS = indata{1,1};
    TP = indata{1,2};
    Dir= indata{1,3};
    DSP= indata{1,4};
    gammajsp = HS{1,1}*0+3.3;
    duration = HS{1,1}*0+3600;
    dtbc     = HS{1,1}*0+2;  %note it was 1. I change it to 2 after 1/19/2024
    
    

    sid = length(HS)-1;
    DSP_rad=[];S=[];
    skip_nr = 1;    % skip the column (1 is time, 2 ~ sid is station nr.)
   
for j = 1:sid
    DSP_rad = DSP{1,j+1}*pi./180;    % transform directional spreading in degree to radians;
    S        = 1./(DSP_rad.^2)-1;
    jonswap_nr{j} = [HS{j+skip_nr}, TP{j+skip_nr}, Dir{j+skip_nr}, gammajsp, S, duration, dtbc];
    file_write_var_ID = fopen([output_dir, ['jonswap',num2str(j),'_',yr_id,'.txt']],'w');
    fprintf(file_write_var_ID, '%f %f %f %f %f %f %f\n', jonswap_nr{1,j}.');
    fclose(file_write_var_ID);
end

1+1