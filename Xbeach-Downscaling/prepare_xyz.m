%% prepare bathymetry data (x.grd, y.grd z.grd)

input_dir = '/gpfs/work/chenw1/Xbeach_wei/Bathymetry/';
dikes       = 'dikes_clip.tif';
domain_name = 'Xbeach_new_coast_31.tif';

%output_dir:
output_dir = '/gpfs/work/chenw1/Xbeach_wei/';
output_folder = 'extract_domain';
dcut = 50;           % cut part of the domain side avoiding strange numerical errors 
alpha = 270*pi/180;  % rotate the domain to allow land boundary always on the right hand side


% read bathymetry
II = geotiffinfo([input_dir,domain_name]);
[xx,yy] = pixcenters(II);
[X, Y] = meshgrid(xx,yy);
Z = imread([input_dir,domain_name]);
% X, Y, Z are 10*10 meter grid;

% read non-erodable structures
II2 = geotiffinfo([input_dir,dikes]);
[xx2,yy2] = pixcenters(II2);
[X2, Y2] = meshgrid(xx2,yy2);
Z2 = imread([input_dir,dikes]);


%
%---------- outputs ----------------


% origin of the model domain relative to the world (lon, lat) coordinates:
xori = X(1,1);           %lon
yori = Y(1,1);           %lat

dlmwrite([output_dir,output_folder,'/','XBeach_tidal_inlet_xori_yori.txt'],[xori yori],'delimiter','\t');

nremove = size(Z,2);           

Xq_new = (X(1:end,1:end-dcut)-xori); 
Yq_new = (Y(1:end,1:end-dcut)-yori); 
Zq6_new= (Z(1:end,1:end-dcut)); 

X_xbeach = Xq_new*cos(alpha)+Yq_new*sin(alpha);
Y_xbeach = -sin(alpha)*Xq_new+Yq_new*cos(alpha);

Xq_dike= (X2(:,1:end)-xori); 
Yq_dike= (Y2(:,1:end)-yori); 
Z2_dike= Z2;

X_dike_xbeach = Xq_dike*cos(alpha)+Yq_dike*sin(alpha);
Y_dike_xbeach = -sin(alpha)*Xq_dike+Yq_dike*cos(alpha);

nelayer = Zq6_new*0+10;
nelayer(find(Zq6_new>=4))=0;


Xq_loc = X_dike_xbeach(find(Z2_dike>=-100));Yq_loc = Y_dike_xbeach(find(Z2_dike>=-100));

for i = 1:size(Xq_loc,1)
nelayer(find(abs(X_xbeach-Xq_loc(i))<=4&abs(Y_xbeach-Yq_loc(i))<=4))=0;
end

% make the lateral boundary nonerodable (to avoid unrealistic numerical errors)
nelayer(:,1:2)=0;nelayer(:,end-1:end)=0;

Zq6_new(find(Zq6_new<=-100))=10;
figure;contourf(Xq_new,Yq_new,Zq6_new,100);caxis([-25 20]);
figure;contourf(Xq_new,Yq_new,nelayer);

%%
dlmwrite([output_dir,output_folder,'/','x.grd'],X_xbeach','delimiter',' ');
dlmwrite([output_dir,output_folder,'/','y.grd'],Y_xbeach','delimiter',' ');
dlmwrite([output_dir,output_folder,'/','z_2090.grd'],-0.8+Zq6_new','delimiter',' ');
dlmwrite([output_dir,output_folder,'/','nelayer.grd'],nelayer','delimiter',' '); 