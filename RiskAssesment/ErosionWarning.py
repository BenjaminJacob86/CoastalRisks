# Early Warning for erosion
# compare seagrass and no seagrass scenarios
import os
import sys
sys.path.insert(0,'/work/gg0028/SCHISM/schism-hereon-utilities/')
from schism import*
plt.ion()
# Erosion Class

#Sed type  ;  Sd50 [m]  ;  Srho [kg/m3]  ;  Wsed [m/s]  ;  Erate ;  tau_ce   [Pa]  ;  ierosion
#    1  6.000000000000000E-005   2650.00000000000   2.159285844822444E-003  1.200000000000000E-003  0.118610510449951  0
#    1  7.000000000000001E-005   2650.00000000000   2.931197281148137E-003  1.200000000000000E-003  0.126236844675281  0
#    1  1.000000000000000E-004   2650.00000000000   5.901979093719125E-003  1.200000000000000E-003  0.143418736131942  0
#    1  1.250000000000000E-004   2650.00000000000   9.051553092008261E-003  1.200000000000000E-003  0.153956425595107  0
#    1  2.500000000000000E-004   2650.00000000000   2.978786443183705E-002  1.200000000000000E-003  0.190816725093894  0
#    1  5.000000000000000E-004   2650.00000000000   6.817372710987034E-002  1.200000000000000E-003  0.264568938546120  0
#    1  1.000000000000000E-003   2650.00000000000   0.116975861995752       1.200000000000000E-003  0.489441249602227  0
#    1  1.990000000000000E-003   2650.00000000000   0.176875255182802       1.200000000000000E-003   1.18393298911042  0
#



######### Load Setup #################
rundir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/Veg_REF/'
os.chdir(rundir)
s=schism_setup()

#region_limit=[8.0,9.2,53.38, 54.35]
#D=np.asarray(s.depths)
#s.plotAtnodesGeo(D,region_limit=region_limit)

ncdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/Veg_CNTRL/outputs02/'
ncdir='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/Veg_REF/outputs_all/'
ncdir2='/work/gg0028/g260114/RUNS/GermanBight/GB_2017_wave_sed/wwm_veg/Veg_max/outputs_all/'


# 3 day fore cast
s.ds=schism_outputs_by_variable(ncdir,varlist=['out2d']).ds
s.ds2=schism_outputs_by_variable(ncdir2,varlist=['out2d']).ds


#### Determine critical shear stress based on sedimnent compositions #############

# take average critical shear stress based on sediment composition
nsed=8
tau_ce=[0.118610510449951,0.126236844675281,0.143418736131942,0.153956425595107,0.190816725093894,0.264568938546120,0.489441249602227,1.18393298911042]
bed_frac=np.zeros((8,s.nnodes))
for i in range(1,9):
    name='bed_frac_{:d}.ic'.format(i)
    s.read_gr3(name)
    bed_frac[i-1,:]=s.gr3[name.split('.')[0]]
tau_ce_mean=np.sum([bed_frac[i,:]*tau_ce[i] for i in range(nsed)],axis=0)

d50=np.asarray([6.000000000000000E-005,7.000000000000001E-005,1.000000000000000E-004,1.250000000000000E-004,
2.500000000000000E-004,5.000000000000000E-004,1.000000000000000E-003,1.990000000000000E-003])*1000


fig, axis = plt.subplots(3,3)
flatAx=axis.flatten()
for i in range(8):
    f=bed_frac[i,:]
    ph,ch,ax=s.plotAtnodes(f*100,ax=flatAx[i])
    ch.set_label('fraction [%]')
    #ax.set_title('d50: {:.2f} mm | tau_ce {:.2f} Pa: '.format(d50[i],tau_ce[i]))
    ax.set_title('d50: {:.2f} mm '.format(d50[i]))

    
ph,ch,ax=s.plotAtnodes(tau_ce_mean,ax=flatAx[-1])
ch.set_label('w avg. tau_{crit} [Pa]')


#dcrit=[0.118610510449951,0.126236844675281,0.143418736131942,0.153956425595107,0.190816725093894,0.264568938546120,0.489441249602227
# 1.18393298911042].


######### Plot atmo conditions for refrence of storm / select period ########################
p=param()
reftime=dt.datetime(int(p.get_parameter('start_year')),int(p.get_parameter('start_month')),int(p.get_parameter('start_day')),int(p.get_parameter('start_hour')),0,0)
dates1=reftime+s.ds['out2d'].time.values*dt.timedelta(seconds=1)  #np.timedelta64(1,'s')

if False:
    coord=(7.352281806246741, 53.75857500576785)
    nn=s.find_nearest_node(coord[0],coord[1])
    elev=s.ds['out2d']['elevation'][:,nn].values
    time=s.ds['out2d']['time'].values
    wind=np.sqrt((s.ds['windSpeed']['windSpeed'][:,:,nn]**2).sum(dim='ivs')).values

    plt.figure()
    #plt.plot(time/86400,elev)
    plt.plot(time/86400,wind)

    #p1=263-265
    #p2=299-301
    plt.vlines(263,0,21,colors='blue')
    plt.vlines(265,0,21,colors='blue')

    plt.vlines(299,0,21,colors='r')
    plt.vlines(301,0,21,colors='r')


    day=time/86400
    inds1=[np.where(day==val)[0][0] for val in [263,265]]
    inds2=[np.where(day==val)[0][0] for val in [299,301]]

else:
    time=s.ds['out2d']['time'].values    
    day=time/86400
    inds1=[np.where(day==val)[0][0] for val in [263,265]]
    inds2=[np.where(day==val)[0][0] for val in [299,301]]
    
    dates2=reftime+s.ds2['out2d'].time.values*dt.timedelta(seconds=1)  #np.timedelta64(1,'s')
    
############ load bottom stress corresponding to time period #####################

tau_mag=np.sqrt((s.ds['bottomStress']['bottomStress']**2).sum(dim='ivs'))    #reference
tau_magB=np.sqrt((s.ds2['bottomStress']['bottomStress']**2).sum(dim='ivs'))  # NBS
#time=tau_mag.time.values

#tau_mag1=tau_mag.sel(time=slice(time[inds1[0]],time[inds1[1]]))
#tau_mag1=tau_mag.sel(time=slice(time[inds2[0]],time[inds2[1]]))
#tau_mag2=tau_mag.sel(time=slice(time[inds2[0]],time[inds2[1]]))

# reduce to time reference:
tau_mag=tau_mag.sel(time=slice(time[inds2[0]],time[inds2[1]]))
tau_magB=tau_magB.sel(time=slice(time[inds2[0]],time[inds2[1]]))  # NBS

#M1=tau_mag1.mean(axis=0).values
#M2=tau_mag2.mean(axis=0).values

#plt.figure()
#ph,ch,ax=s.plotAtnodes(M1,cmap=plt.cm.turbo)
#ch.set_label('<tau> [Pa]')
#ph.set_clim((0,1.2))

#plt.figure()
#ph,ch,ax=s.plotAtnodes(M2,cmap=plt.cm.turbo)
#ch.set_label('<tau> [Pa]')
#ph.set_clim((0,1.2))


# 72 h
# expand array
nt=len(tau_mag)
tau_bar=np.tile(tau_ce_mean,(nt,1)) # criticcal shear stress 

#a = xr.DataArray(
#    data=tau_bar,
#    dims=tau_mag.dims,
#    coords=tau_mag.coords
#    )
    
# load compare data sets    

# convert to xarray     
#a = xr.DataArray(
#    data=tau_bar,
#    dims=tau_mag1.dims,
#    coords=tau_mag1.coords
#    )
#
#b = xr.DataArray(
#    data=tau_bar,
#    dims=tau_mag2.dims,
#    coords=tau_mag2.coords
#    )
#

tau_bar = xr.DataArray(
    data=tau_bar,
    dims=tau_mag.dims,
    coords=tau_mag.coords
    )

##### Checke the critical share stress excidence ratio ###################
Fexceedence=tau_mag/tau_bar
#Fexceedence1.name='Factor_Tau_ce_day_263-265'
#Fexceedence1.to_netcdf(Fexceedence1.name+'.nc')

FexceedenceB=tau_magB/tau_bar

# save netcdf
#Fexceedence2=tau_bar = xr.DataArray(
#    data=tau_bar,
#    dims=tau_mag.dims,
#    coords=tau_mag.coords
#    )/b
#Fexceedence2.name='Factor_Tau_ce_day_299-301'
#Fexceedence2.to_netcdf(Fexceedence2.name+'.nc')
#    
#Fexceedence=tau_mag/a
#Fexceedence.name='Factor_Tau_ce'
#Fexceedence.to_netcdf(Fexceedence.name+'.nc')



##### Count Ration of critical shear stress acceedance #############
#A=Fexceedence.values
A1=Fexceedence.values
R1=(A1>1).sum(axis=0)/nt

A2=FexceedenceB.values
R2=(A2>1).sum(axis=0)/nt

########  Colors for it -> create new color map
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
green=np.array([0/256,128/256,0/256,1])
yellow=np.array([255/256,255/256,128/256,1])
orange=np.array([255/256,165/256,0/256,1])
red=np.array([256/256,0/256,0/256,1])
clrs=np.vstack((green,yellow,orange,red))
mymap = ListedColormap(clrs)  # reguster colormap for matplotlib
###########


#ph,ch,ax=s.plotAtnodes(R1,cmap=mymap)

#landcolor='darkslategray'
landcolor='teal'#'indigo'
landcolor='default'
ph,ch,ax=s.plotAtnodesGeo(R1,cmap=mymap,landcolor=landcolor) #,extend='None'
ph.set_clim((0,1)) #
#ch.set_ticks(np.linspace(0,1,5))
ch.set_ticks(np.linspace(0,.75,4))
#ch.set_ticklabels(['No','low','medium','high','severe'])
ch.set_ticklabels(['No','low','increased','high'])
#ch.set_label('Erosion Risk')
plt.sca(ax)
plt.title('Erosion Risk')
plt.tight_layout()


plt.figure()
ph,ch,ax=s.plotAtnodesGeo(R2,cmap=mymap,landcolor=landcolor) #,extend='None'
ph.set_clim((0,1)) #
#ch.set_ticks(np.linspace(0,1,5))
ch.set_ticks(np.linspace(0,.75,4))
#ch.set_ticklabels(['No','low','medium','high','severe'])
ch.set_ticklabels(['No','low','increased','high'])
#ch.set_label('Erosion Risk')
plt.sca(ax)
plt.title('Erosion Risk')
plt.tight_layout()


# bin into intervalls
#  below
bins=[0,0.25,.5,.75]
R1bins=R1.copy()
R2bins=R2.copy()

for lbin in bins:
    gt = R1 >  lbin
    R1bins[gt] =lbin

    gt = R2 >  lbin
    R2bins[gt] =lbin

# nomralize toi levels 1,2,3    
R1bins/=0.25  
R2bins/=0.25  


#
## depending on the colormap one of the two interpolation functions from chatgpt
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.colors as mcolors
#
#def create_segmented_colormap(base_cmap, num_segments):
#    # Get the existing colormap
#    cmap = plt.get_cmap(base_cmap)
#
#    # Determine the number of colors to sample from the existing colormap
#    num_colors = len(cmap.colors)
#    indices = np.linspace(0, num_colors - 1, num_segments, dtype=int)
#
#    # Sample colors from the existing colormap
#    colors = [cmap.colors[i] for i in indices]
#
#    # Create a segmented colormap using LinearSegmentedColormap
#    return mcolors.LinearSegmentedColormap.from_list(f'{base_cmap}_segmented', colors, N=num_segments)
#
## Example usage:
#base_cmap = 'viridis'  # Existing colormap
##base_cmap = 'seismic'  # Existing colormap
#base_cmap = 'jet'  # Existing colormap
#num_segments = 12  # Number of segments for the new colormap
#
#segmented_cmap = create_segmented_colormap(base_cmap, num_segments)
#
## Plot a colorbar to visualize the new segmented colormap
#plt.colorbar(plt.cm.ScalarMappable(cmap=segmented_cmap))
#plt.show()



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def create_segmented_colormap(base_cmap, num_segments):
    # Get the existing colormap
    cmap = plt.get_cmap(base_cmap)

    # Determine the number of colors to sample from the existing colormap
    num_colors = 256  # Number of colors in jet colormap
    indices = np.linspace(0, num_colors - 1, num_segments, dtype=int)

    # Sample colors from the existing colormap
    colors = [cmap(i) for i in np.linspace(0, 1, num_colors)]
    sampled_colors = [colors[i] for i in indices]

    # Create a segmented colormap using LinearSegmentedColormap
    return mcolors.LinearSegmentedColormap.from_list(f'{base_cmap}_segmented', sampled_colors, N=num_segments)

# Example usage:
#base_cmap = 'jet'  # Existing colormap
base_cmap = 'seismic'  # Existing colormap
num_segments = 7  # Number of segments for the new colormap

segmented_cmap = create_segmented_colormap(base_cmap, num_segments)

# Plot a colorbar to visualize the new segmented colormap
plt.colorbar(plt.cm.ScalarMappable(cmap=segmented_cmap))
plt.show()


os.chdir('/work/gg0028/g260114/RUNS/GermanBight/Pics/')

region_limit=[8.0,9.2,53.38, 54.35]


for region_limit in [None,[8.0,9.2,53.38, 54.35]]:

    if region_limit==None:
        regstr=''
    else:
        regstr='_'+'_'.join([str(i) for i in region_limit])

    plt.figure()
    ph,ch,ax=s.plotAtnodesGeo(R1,cmap=mymap,region_limit=region_limit,landcolor=landcolor) #,extend='None'
    ph.set_clim((0,1)) #
    #ch.set_ticks(np.linspace(0,1,5))
    ch.set_ticks(np.linspace(0,.75,4))
    #ch.set_ticklabels(['No','low','medium','high','severe'])
    ch.set_ticklabels(['No','low','increased','high'])
    #ch.set_label('Erosion Risk')
    plt.sca(ax)
    plt.title('Erosion Risk')
    plt.tight_layout()
    plt.savefig('ErosionRiskNoNBS'+regstr+'.png',dpi=300)

    plt.figure()
    ph,ch,ax=s.plotAtnodesGeo(R2,cmap=mymap,region_limit=region_limit,landcolor=landcolor) #,extend='None'
    ph.set_clim((0,1)) #
    #ch.set_ticks(np.linspace(0,1,5))
    ch.set_ticks(np.linspace(0,.75,4))
    #ch.set_ticklabels(['No','low','medium','high','severe'])
    ch.set_ticklabels(['No','low','increased','high'])
    #ch.set_label('Erosion Risk')
    plt.sca(ax)
    plt.title('Erosion Risk with NBS')
    plt.tight_layout()
    plt.savefig('ErosionRiskNBS'+regstr+'.png',dpi=300)

    plt.figure()
    new_cmap=segmented_cmap 
    ph,ch,ax=s.plotAtnodesGeo(R2bins-R1bins,cmap=new_cmap,region_limit=region_limit) #,extend='None'
    ph.set_clim((-3,3))
    plt.sca(ax)
    plt.title('Erosion Risk level change due to NBS')
    plt.tight_layout()
    plt.savefig('ErosionRiskReductionNBS'+regstr+'.png',dpi=300)
    
    plt.close('all')