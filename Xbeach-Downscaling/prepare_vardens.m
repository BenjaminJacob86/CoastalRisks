function prepare_vardens
% This code extract variance density spectra files from wavewatch3 output and to prepare input files for Xbeach
% Note that Xbeach x-axis is alwyas oriented towards the coast, y-axis is
% along the coast
% author: Wei Chen (wei.chen@hereon.de)

clear all
close all


%  input data--------------------------------------------------------------
%  specify the directory, staion name and extract period
input_dir = '/gpfs/work/chenw1/Xbeach_wei/xbeach-repo/Rest_Coast/wave_spectral';
input_dir = '/gpfs/work/chenw1/Xbeach_wei/xbeach-repo/Rest_Coast/2013';
output_dir= '/gpfs/work/chenw1/Xbeach_wei/xbeach-repo/Rest_Coast/2013/wavefiles';
%output_dir= '/gpfs/work/chenw1/Xbeach_wei/xbeach-repo/Rest_Coast/wave_spectral/wavefiles';
yr_id = '2013'
if strcmp(yr_id, '2017')
input_fn = 'ww3.201710_spec.nc'
else
input_fn = 'ww3_2013_11_12_spec.nc'
end
station_list = 'points.list'
station_name = 'Fino-1';
day = 38;                          %how many days you want to specify from ww3 data
alfa = +270;                         %angle of x-axis from north, when thetanaut = 1; 
%alfa = +90;                        %angle of x-axis from east, when thetanaut = 0;
dtbc = 2;                          %in seconds, timestep used to describe time series of wave energy and long wave flux at offshore boundary
rt   = 3600;                       %in seconds, duration of wave spectrum at offshore boundary, in morphological time
t_start      = 1;                  %start time of extracting



%--------------------------------------------------------------------------
t_steps      = day*24;                  %how many time steps you want to extract

% find the right station number
file_id = fopen([input_dir,'/',station_list]);
stat_list = textscan(file_id,'%s','delimiter','\n');
fclose(file_id);
temp = stat_list{1};
if strcmp(yr_id,'2017')
for stat_id = 1:length(temp)
    if strcmp(station_name, temp{stat_id}(27:end-1))
        station_nr = stat_id;
    end
end
else
    station_nr = 956;
end

% read data from input_fn
frequency = double(ncread([input_dir, '/', input_fn],'frequency'));
direction = double(ncread([input_dir, '/', input_fn],'direction'))+alfa; %note in Xbeach 0 degree is toward the land
direction(find(direction>=360))=direction(find(direction>=360))-360;
direction(find(direction<360&direction>180)) = direction(find(direction<360&direction>180)) -360; % make direction in range -180 ~ 180 degree

fileID = fopen([output_dir,'/','filelist.txt'],'wt');
fprintf(fileID,'FILELIST');
fclose(fileID);

for tid = 1:t_steps
    density   = squeeze(double(ncread([input_dir, '/', input_fn],'efth',[1,1,station_nr,tid],[Inf,Inf,1,1])));    
    
    %sort density in accordence of the direction in an accending order
    dir_den = [direction, density];
    dir_den_srt = sortrows(dir_den);
    dir_srt = dir_den_srt(:,1);
    den_srt = dir_den_srt(:,2:end)*pi/180;
    
    % make vardens.txt file
    file_var_ID = fopen([output_dir,'/',['vardens',num2str(tid),'.inp']],'w');
    fprintf(file_var_ID, '%.0f\n', length(frequency));
    fprintf(file_var_ID, '%7.3f\n', frequency);
    fprintf(file_var_ID, '%.0f\n', length(direction));
    fprintf(file_var_ID, '%7.3f\n', dir_srt); 
    fprintf(file_var_ID, [repmat('  %.7e',1,length(direction)),'\n'], den_srt);     
    fclose(file_var_ID);
    
    %dlmwrite([output_dir,'/',['vardens',num2str(tid),'.inp']],double(length(frequency)));
    %dlmwrite([output_dir,'/',['vardens',num2str(tid),'.inp']], double(frequency),'precision','%.3f','-append');
    %dlmwrite([output_dir,'/',['vardens',num2str(tid),'.inp']],double(length(direction)),'-append');
    %dlmwrite([output_dir,'/',['vardens',num2str(tid),'.inp']],double(dir_srt),'precision','%.3f','-append');
    %dlmwrite([output_dir,'/',['vardens',num2str(tid),'.inp']],double(den_srt),'delimiter',' ','-append');
    
    % make filelist.txt file
    fileID = fopen([output_dir,'/','filelist.txt'],'a');
    formatSpec = ['\n%d ' '%d ' 'vardens',num2str(tid),'.inp'];
    fprintf(fileID, formatSpec, [rt dtbc]);
    fclose(fileID);
end

exit
