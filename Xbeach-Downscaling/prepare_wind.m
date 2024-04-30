%function prepare_wind
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
yr_id = '2017';
for i = 1:2
    var_name = {['WinX','_',yr_id]; ['WinY','_',yr_id]};
    file_var_ID = fopen([input_dir,var_name{i},'.txt'],'r');
    indata{i} = textscan(file_var_ID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines',1);
    fclose(file_var_ID);
end
    winx = indata{1,1};
    winy = indata{1,2}; 
    sid = length(winx)-1;
    time = 3600*(1:length(winx{1,1}))-3600;
    
    win_vel=[];win_deg=[];
    skip_nr = 1;    % skip the column (1 is time, 2 ~ sid is station nr.)

for j = 1:sid
    win_vel = sqrt(winx{1,j+1}.^2 + winy{1,j+1}.^2);    % transform directional spreading in degree to radians;
    win_deg = 180*atan(winy{1,j+1}./winx{1,j+1})/pi;
    %win_deg(find(win_deg<=0))=win_deg(find(win_deg<=0))+360;
    wind_nr{j} = [time', win_vel, win_deg];
    file_write_var_ID = fopen([output_dir, ['wind',num2str(j),'_',yr_id,'.txt']],'w');
    fprintf(file_write_var_ID, '%f %f %f\n', wind_nr{1,j}.');
    fclose(file_write_var_ID);
end    
    
