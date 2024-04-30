# Flood risk
import sys
import os
from glob import glob
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities')
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hzg-utilities/')
from schism import *
import dask
dask.delayed()
import xarray as xr
plt.ion()
import time

#HW
#Swash Regime Rhigh / Dhigh < 1 Runup confined to foreshore

#Overwash Regime
#Rhigh/Dhigh > 1 and Rlow/Dhigh < 1 Runup exceeds the elevation of the First line of defence

#Terrestrial Inunda-tion Regime Rhigh/Dhigh > 1 and Rlow/Dhigh > 1 Elevation of the base of the swash motion exceeds the el-
#evation of the First line of defence.


setupdir='/gpfs/work/villal/storage/schism/JadeB/'
ncdir='/gpfs/work/villal/storage/schism/JadeB/outputs/'
ncdir='/gpfs/work/jacobb/data/linkedfiles/'

setupdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/'
ncdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs_merged/'


os.chdir(setupdir)
s=schism_setup()

#ds=schism_output2(ncdir).nc


p=param(setupdir+'/param.nml')
reftime=dt.datetime(np.int(p.get_parameter('start_year')),
            np.int(p.get_parameter('start_month')),
            np.int(p.get_parameter('start_day')),
            np.int(p.get_parameter('start_hour')),0,0)
ndays_per_stack=1
def find_stack_for_date(date):
	""" return schism stacknr (1 based) contiaininge the searched date """
	return int((date-reftime)/dt.timedelta(days=ndays_per_stack))+1
date0=dt.datetime(2017,10,15)
date1=dt.datetime(2017,11,1)
stack0=find_stack_for_date(date0)
stack1=find_stack_for_date(date1)

acc=schism_outputs_by_variable(ncdir,varlist='out2d',min_stack=stack0,max_stack=stack1)

avg=acc.get('elevation')#.sel(nSCHISM_hgrid_node=nn)
dates=np.asarray(reftime,np.datetime64)+np.asarray(avg.time.values,np.timedelta64(1,'s'))

def read_time_selection(self):
	np.datetime64(self.exfrom.get())-self.dates[0]
	i0=np.int((np.datetime64(self.exfrom.get())-self.dates[0])/self.dt)
	i1=np.int((np.datetime64(self.exto.get())-self.dates[0])/self.dt)
	return i0,i1

# select
i00=24
jump=12
i0=i00
i1=i00+13

data=ds.elev[i0:,:]

ts1=ds.elev[i0:,0]

#DataArray.rolling(dim=None, min_periods=None, center=False, **window_kwargs)

#max=ds.elev.rolling(time=14,min_periods=13).max()
#max[:,0].values


wwndw=12
modus='trim'#,'exact'
data=avg
max=data.coarsen(time=wwndw,boundary=modus).max()
min=data.coarsen(time=wwndw,boundary=modus).min()

ts1=data[:,0]

torg=ts1.time.values
dtnew=1800  # closer to multiple of 0.24
tinterp=np.arange(torg[0],torg[-1]+dtnew,dtnew)

max=data.interp(time=tinterp).coarsen(time=wwndw*2,boundary=modus).max()
mmin=data.interp(time=tinterp).coarsen(time=wwndw*2,boundary=modus).min()




inode=21000
inode=31000
plt.clf()
ts=data.interp(time=tinterp)[:,inode]
dates=np.asarray(reftime,np.datetime64)+np.asarray(ts.time.values,np.timedelta64(1,'s'))
plt.plot(dates,ts,'.-')
window_max=max[:,inode].values
window_min=min[:,inode].values
inds=np.asarray([ np.where(ts==val)[0][0] for val in window_max])
inds2=np.asarray([ np.where(ts==val)[0][0] for val in window_min])
plt.plot(dates[inds],window_max,'ro')
plt.plot(dates[inds2],window_min,'ko')

keep1=np.hstack((True,np.diff(inds)>8))
keep2=np.hstack((True,np.diff(inds2)>8))

window_max=window_max[keep1]
window_min=window_min[keep2]
inds=inds[keep1]
inds2=inds2[keep2]

plt.clf()
ts=data[:,inode].values
plt.plot(dates,ts,'.-')
window_max=max[keep1,inode].values
window_min=min[keep2,inode].values
inds=np.asarray([ np.where(ts==val)[0][0] for val in window_max])
inds2=np.asarray([ np.where(ts==val)[0][0] for val in window_min])
plt.plot(dates[inds],window_max,'ro')
plt.plot(dates[inds2],window_min,'ko')
plt.gcf().autofmt_xdate()




Rhigh=data.max(axis=0).values
s.plotAtnodes(Rhigh)
plt.figure()
s.plotAtnodes(s.depths)


plt.figure()
Rhigh=window_max.max(axis=0)
s.plotAtnodes(Rhigh)

#fix overlaps - which should not be there?
idiff=np.where(np.diff(inds)<6)




lw=data.min(axis=0).values
hw=data.max(axis=0).values

#shift


s

Dhigh= # Dike height 8 ms over MHW ~ 8+ 1m  above NHN











