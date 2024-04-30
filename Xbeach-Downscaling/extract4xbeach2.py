"""
Extract SCHISM/WWM spatial fields for xbeach
"""

import os
import sys
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import* # import schism class to read grid structure
s=schism_setup()

s.init_node_tree(latlon=False) # build nearest neighbour serach tree


workdir=os.getcwd()
setupdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/' # where the schism grid is
ncdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs_all/' #schims ncdir

#setupdir='/work/gg0028/g260114/RUNS/GermanBight/xaver/GB_2013_xaver/'  # for the 2013 storm case min_stack = 
#ncdir='/work/gg0028/g260114/RUNS/GermanBight/xaver/GB_2013_xaver/output_xaver/outputs_all/'


# all wave outputs and sea level are contained within out2d files
#s.ds=schism_outputs_by_variable(ncdir=ncdir,varlist=['out2d'],max_stack=2).ds # load access to data test with only first stack
os.chdir(setupdir)
s.ds=schism_outputs_by_variable(ncdir=ncdir,varlist=['out2d'],min_stack=33,max_stack=37).ds # load access to data; for 2017: min_stack =299, max_stack=303, for 2013: min_stack =33, max_stack =37
os.chdir(workdir)

m=[[1., 366727, 5920100],
[2., 364628, 5920690],
[3., 364246, 5919332],
[4., 363864, 5917974],
[5., 363482, 5916616],
[6., 363100, 5915258],
[7., 389320, 5959870],
[8., 392150, 5959870],
[9., 394980, 5959870],
[10., 397810, 5959870],
[11., 400640, 5959870],
[12., 403470, 5959870],
[13., 386795, 5958743],
[14., 388905, 5958743],
[15., 391015, 5958743],
[16., 393125, 5958743],
[17., 395235, 5958743],
[18., 379270, 5961200],
[19., 381770, 5961200],
[20., 384000, 5961200],
[21., 379270, 5968200],
[22., 381770, 5968200],
[23., 384000, 5968200],
[24., 379270, 5955300],
[25., 384000, 5955300]]
m=np.asarray(m)

xq=m[:,1]
yq=m[:,2]

plt.ion()
s.plot_domain_boundaries(latlon=False)
#plt.plot(xq,yq,'ko')

nn=s.node_tree_xy.query(list(zip(xq,yq)))[1] #nearest neighbour grid nodes
x,y=np.asarray(s.x),np.asarray(s.y)
#plt.plot(x[nn],y[nn],'r+')


# redice xarray handle to nearest neighbourts
s.ds['out2d']=s.ds['out2d'].sel(nSCHISM_hgrid_node=nn)

Zeta=s.ds['out2d']['elevation'].values
HS=s.ds['out2d']['sigWaveHeight'].values
TP=s.ds['out2d']['peakPeriod'].values
DSP=s.ds['out2d']['meanDirSpreading'].values
Dir=s.ds['out2d']['dominantDirection'].values  # 1 !Peak (dominant) direction (degr) {dominantDirection}  2D
WinX=s.ds['out2d']['windSpeedX'].values
WinY=s.ds['out2d']['windSpeedY'].values

t=s.ds['out2d']['time'].values

# load reference tiem from param.nml
p=param()
reftime=dt.datetime(int(p.get_parameter('start_year')),
int(p.get_parameter('start_month')),
int(p.get_parameter('start_day')),
int(p.get_parameter('start_hour')),0,0)	 

header='time (seconds since {:s}) | Zeta'.format(str(reftime))
np.savetxt('Zeta_2013.txt',np.vstack((t,Zeta.T)).T,header=header)

header='time (seconds since {:s}) | HS'.format(str(reftime))
np.savetxt('HS_2013.txt',np.vstack((t,HS.T)).T,header=header)

header='time (seconds since {:s}) | TP'.format(str(reftime))
np.savetxt('TP_2013.txt',np.vstack((t,TP.T)).T,header=header)

header='time (seconds since {:s}) | DSP'.format(str(reftime))
np.savetxt('DSP_2013.txt',np.vstack((t,DSP.T)).T,header=header)

header='time (seconds since {:s}) | Dir'.format(str(reftime))
np.savetxt('Dir_2013.txt',np.vstack((t,Dir.T)).T,header=header)

header='time (seconds since {:s}) | WinX'.format(str(reftime))
np.savetxt('WinX_2013.txt',np.vstack((t,WinX.T)).T,header=header)

header='time (seconds since {:s}) | WinY'.format(str(reftime))
np.savetxt('WinY_2017.txt',np.vstack((t,WinY.T)).T,header=header)
