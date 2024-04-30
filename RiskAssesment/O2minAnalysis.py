#Bio Geo Johannes
import os
import netCDF4
import sys
import csv
import matplotlib
from matplotlib import pyplot as plt
background=False
if background:
	matplotlib.use('Agg') # backend
else:
	plt.ion()	
import datetime as dt
import glob
from scipy.spatial import cKDTree
import datetime as dt
from netCDF4 import Dataset,MFDataset
import pickle
from scipy.interpolate import interp1d
# own and 3d party libraries

#levante
sys.path.insert(0,'/home/g/g260114/git/schism-hereon-utilities/')
sys.path.insert(0,'/home/g/g260114/git/schism-hereon-utilities/')
#strand
sys.path.insert(0,'/gpfs/work/ksddata/code/schism/scripts/schism-hzg-utilities/')
from schism import * # import schism functions
import cmocean	
import time

#setupdir='/work/bg1186/g260094/SNS/SNSE3D_01a_CMEMS/'
#ncdir='/work/bg1186/g260094/SNS/SNSE3D_01a_CMEMS/combined/'
#setupdir='/work/bg1186/g260094/Deutsche_Bucht/Elbe/ECO07d5_ref/'
#ncdir='/work/bg1186/g260094/Deutsche_Bucht/Elbe/ECO07d5_ref/combined/'
#

# RUN and output directory
setupdir='/gpfs/work/pein/data/GB_model/Elbe/ECO07d_5_CMEMS_2012_ref/'
ncdir=setupdir+'outputs_paper/'


# load schism grid and output
os.chdir(setupdir)
s=schism_setup()
s.ds=schism_output2(ncdir).nc

lon,lat=np.asarray(s.lon),np.asarray(s.lat)


# stack 120# SNS model crashes 201206 ~15

#hypoxia=60 #  mmol m-3  ~ 2mg/l       

# Define oxigen Weights and hypoxia concentrations
mass_o2=15.994*2 #g/mo
#mass_o2=15.994/1000 #g/mo
Lper_m3=1000
#mg_per_g=1000
#60* mass_o2  #mg_m-3
#60*mass_o2/Lper_m3/mg_per_g
mmolO2_perM3_to_mg_perL=mass_o2/Lper_m3

#6 mg/l 
Vlow_oxy=np.round(6/mmolO2_perM3_to_mg_perL) # mmol of 6g/l
Vhypoxia=np.round(2/mmolO2_perM3_to_mg_perL) # mmol of 2g/l


# Load River Discharge for dependence
Q=np.loadtxt(setupdir+'flux.th')
Qt=Q[:,0]
Q=Q[:,1]

t=s.ds['time'].values
Qt=t[0]-np.timedelta64(1,'h')+np.timedelta64(1,'s')*Qt

plt.plot(Qt,-Q)
plt.gcf().autofmt_xdate()
plt.xlim((t[0],t[-1]))

#mass_o2/Lper_m3

#Concentrations of oxygen below about 60 mmol m-3 (1 mmol m-3= 1 µM) are termed “hypoxic”; regions with oxygen persistently below this threshold are referred to as “dead zones”: normal respiration is severely
#oxy

 #Oxygen concentrations above 6mg/l are considered to support marine life with minimal problems, while concentrations less than 2mg/l (hypoxia, i.e., oxygen deficiency) are considered to cause severe problems (Best et al., 2007; Levin et al., 2009).
#180*mmolO2_perM3_to_mg_perL

 
ti=0
#idep=

ibtm=s.ds['node_bottom_index'][0,:].values-1
field=s.ds['ECO_oxy'][ti,:].values    # DO in mmol m-3

# calcoulate min
nt=len(s.ds['time'])
nodes=np.arange(len(s.x))

nnodes=s.nnodes
vmin=np.ones(nnodes)*999
imin=np.zeros(nnodes)
hypoxia=np.zeros(nnodes) # #2 g/l   
crit_ox=np.zeros(nnodes)  ##6 g/l   						


for ti in range(nt):
	print(ti)
	field=s.ds['ECO_oxy'][ti,:].values    # DO in mmol m-3
	bottom_filed=field[nodes,ibtm]

	iupdate=bottom_filed<vmin
	vmin[iupdate]=bottom_filed[iupdate]
	imin[iupdate]=ti

	# ctritical oxy
	ibelow=	bottom_filed <= Vlow_oxy #60
	crit_ox[ibelow]=crit_ox[ibelow]+1

	# hypxoia	
	ibelow=	bottom_filed <= Vhypoxia #60
	hypoxia[ibelow]=hypoxia[ibelow]+1


	
# Fokus Area
axlim=(9.832940049784687, 10.137623425894006, 53.4420865032227, 53.560149635446725)	 # ROI
inaxis=(axlim[0] <= lon) & (axlim[1] >= lon) & (axlim[2] <= lat) & (axlim[3] >= lat)
subind_max_hypo=np.argmax(hypoxia[inaxis])   # max hypo

# time Series of hypoxxy
ind_check=np.where(inaxis)[0][subind_max_hypo]

	

# loop 1  # fastes extraction? yes

# not used
t0=time.time()
TS=np.zeros(nt)
for ti in range(nt):
	field1=s.ds['ECO_no3'][ti,:,0].values
	TS[ti]=field1[ind_check]
dt1=time.time()-t0
print(dt1)

###
#
## loop2
#t0=time.time()
#TS=np.zeros(nt)
#for ti in range(nt):
#	ind0=ti*24
#	ind1=(ti+1)*24
#	field1=s.ds['ECO_no3'][ind0:ind1,:,0].values
#	TS[ind0:ind1]=field1[:,ind_check]
#dt2=time.time()-t0
#print(dt2)
#

###

# load nutirents
t0=time.time()
no3=np.zeros(nt)
nh4=np.zeros(nt)
pho=np.zeros(nt)
sil=np.zeros(nt)
for ti in range(nt):
	field1=s.ds['ECO_no3'][ti,:,0].values
	no3[ti]=field1[ind_check]

	field1=s.ds['ECO_nh4'][ti,:,0].values
	nh4[ti]=field1[ind_check]

	field1=s.ds['ECO_pho'][ti,:,0].values
	pho[ti]=field1[ind_check]

	field1=s.ds['ECO_sil'][ti,:,0].values
	sil[ti]=field1[ind_check]

ox=np.zeros(nt)	
for ti in range(nt):
	field1=s.ds['ECO_oxy'][ti,:,0].values
	ox[ti]=field1[ind_check]
	
 
# bioata
bio=dict.fromkeys(['ECO_fla','ECO_dia','ECO_bg', 'ECO_microzoo', 'ECO_mesozoo'])
for key in bio.keys():
	bio[key]=np.zeros(nt)

# loop over loop seems slow but I am to lazy to type	
# dont use dicts
for ti in range(nt):
	for key in bio.keys():
		field1=s.ds[key][ti,:,0].values
		bio[key][ti]=field1[ind_check]
bio_legend=[key.split('_')[1] for key in bio.keys()]
 
det=np.zeros(nt)
for ti in range(nt):
	field1=s.ds['ECO_det'][ti,:,0].values
	det[ti]=field1[ind_check]


	
#In Fig. 3 it would be good to use units in mmol or so for nutrients instead of mgC – which is rather unusual.
#For NO3, NH4 and SiO4 the factor to get from mgC to mmol is mg2mmol = 1/79.5. 
	
# unit conversion ->  mmol/m³
#mgC/m³
mg2mmol = 1/79.5
mg2mmol_po4 = 1/106/12

no3,nh4,sil=no3*mg2mmol,nh4*mg2mmol,sil*mg2mmol
pho=pho*mg2mmol_po4	




# custom clormarp
import matplotlib.colors as mcolors
colors1 = plt.cm.Reds_r(np.linspace(0., 1, 256))[32:]
colors2 = plt.cm.gray(np.linspace(0., 1, 256))[128:]
colors3 = plt.cm.summer_r(np.linspace(0., 1, 64))
colors = np.vstack((colors1, colors2, colors3))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
#
#fig,ax=plt.subplots(3,2,figsize=(10,10))
#ph,ch,ax0=s.plotAtnodes(vmin,ax=ax[0,0],cmap=mymap)	
#plt.axis(axlim)
#ph.set_clim((0,350)) # +60 mmol


	
# plotting	
#For PO4 it is mg2mmol = 1/106/12. 
plt.close('all')		
#fig,ax=plt.subplots(2,2,figsize=(10,6))
fig,ax=plt.subplots(3,2,figsize=(10,10))
plt.tight_layout()
# minimum
#ph,ch,ax0=s.plotAtnodes(vmin,ax=ax[0,0],cmap=cmocean.cm.oxy)	
ph,ch,ax0=s.plotAtnodes(vmin,ax=ax[0,0],cmap=mymap)	
s.plot_domain_boundaries(append=True)
plt.axis(axlim)
ch.set_ticks([0,60,100,150,200,250,300])
ch.set_label('$O_2$' + 'min bottom [mmol m' +'$^{-3}$'+']')
#ph.set_clim((0,300)) # +60 mmol
#ph.set_clim((0,400))
ph.set_clim((0,350)) # +60 mmol
ax[0,0].set_title('O$_2$'+'min')

# low oxy count
#ph,ch,ax0=s.plotAtnodes(hypoxia,cmap=plt.cm.Reds,ax=ax[0,1])	
ph,ch,ax0=s.plotAtnodes(crit_ox,cmap=plt.cm.Reds,ax=ax[0,1])	
s.plot_domain_boundaries(append=True)
plt.axis(axlim)
#ch.set_label('Hypoxia count [h]')
ch.set_label('Crit. O2 count [h]')
ax[0,1].set_ylabel('')
ax[0,1].set_yticks([])
#ax[0,1].set_title('hypoxia count')
ax[0,1].set_title('Crit. O2 count [h]')


axi=ax[0,1]
axi.plot(lon[ind_check],lat[ind_check],'k+')
#axi=ax[0,1]
axi.plot(lon[ind_check],lat[ind_check],'k+',markersize=12)

# Nutrients ##################
axi=ax[1,0]

for vari in [no3,nh4,pho,sil]:
	axi.plot(t,vari)
axi.legend(['no3','nh4','pho','sil'],ncol=4,frameon=False)	
#axi.set_ylabel('nutrients [mgC m' + '$^{-3}$'+ ']')
axi.set_ylabel('nutrients [mmol m' + '$^{-3}$'+ ']')
axi.grid()
axi.set_title('nutrients')

#ymin,ymax=axib.get_ylim()

# o2
def shade_hypox(axi):
	axib=axi.twinx()
	axib.plot(t,ox,'m',linewidth=0.3)
	axib.yaxis.label.set_color('m')
	axib.set_ylabel('O2 [mmol m'+'$^{-3}$'+']')
	#axib.hlines(60,t.min(),t.max(),color='m',linestyle='--')
	ymin,ymax=axib.get_ylim()
	axib.fill_between(t, ymin, ymax, where=ox <= 63,facecolor='darkred', alpha=0.4)
	axib.fill_between(t, ymin, ymax, where=ox <= Vlow_oxy,facecolor='red', alpha=0.3)
shade_hypox(axi)

# biology	
axi=ax[1,1]
for key in bio.keys():
	axi.plot(t,bio[key])
try:
	bio_legend[bio_legend.index('bg')]='cyano'
except:
	pass
axi.legend(bio_legend,ncol=3,frameon=False)
axi.set_ylabel('biomass [mgC m' + '$^{-3}$'+ ']')
axi.grid()
#bio_legend[bio_legend.index('Cyano')]='cyano'
axi.set_title('biota')
# o2
shade_hypox(axi)


# River Discharge
axi=ax[2,0]
#axi.axis('off')
axi.axis('on')
axi.plot(Qt,-Q)
axi.legend(['Q [m3/s]'],frameon=False)
axi.set_xlim(t[0],t[-1])
axi.set_ylim(0,3000)
ax[2,0].set_xlim(t[0],t[-1])
xt=axi.get_xticks()
axi.set_xticks(xt[::4])
shade_hypox(axi)
axi.set_title('River Discharge')

# detritus
axi=ax[2,1]
axi.set_ylabel('biomass [mgC m' + '$^{-3}$'+ ']')
axi.plot(t,det)
axi.legend(['detirtus'],frameon=False)
axi.grid()

# o2
shade_hypox(axi)
plt.tight_layout()
plt.savefig('HypoxiaWarningQ.png',dpi=300)



# adapt dateticks
axi=ax[2,1]

for axi in [ax[1,0],ax[1,1],ax[2,1]]:
	axi.set_xticklabels(axi.get_xticklabels(),rotation=45)
plt.tight_layout()

ax[2,0].set_xticklabels(axi.get_xticklabels(),rotation=45)
plt.savefig('HypoxiaPlots2.png',dpi=300)

ax[2,0].set_xticklabels(ax[2,0].get_xticklabels(),rotation=45)

ihypox=np.where(ox<=60)[0]   # less than 60
for axzom in ax.flatten()[2:]:
	axzom.set_xlim(t[ihypox[0]-24],t[ihypox[-1]+24])
plt.savefig('HypoxiaPlots_zoom.png',dpi=300)

	
ihypox=np.where(ox<=60)[0]
for axzom in ax.flatten()[2:]:
	axzom.set_xlim(t[0],t[-1])
	


plt.tight_layout()