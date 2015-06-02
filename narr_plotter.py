from netCDF4 import Dataset
import sys
from pylab import *
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from datetime import datetime
import scipy.ndimage

def regMap():
    figure(figsize=(10,8))
    m = Basemap(llcrnrlon=-110.14,llcrnrlat=26.28,urcrnrlon=-83.38,urcrnrlat=45.48,projection='lcc',lat_1=25.,lat_2=25.,lon_0=-95.,resolution ='l',area_thresh=1000.)
    m.drawcoastlines()
    m.drawcountries()
    m.drawstates()
    m.drawcounties()
    return m

def e2td(E):
    B = (np.log(E / 6.108)) / 17.27
    D = (237.3 * B) / (1 - B)
    return D

def cc(temp):
    e = 6.112 * np.exp((17.67*temp)/(temp + 243.5))
    return e

def q2w(w):
    return w/(1.+w) 

def w2e(w, p):
    eps = 0.622
    e = ((w/eps)*p)/(1 + (w/eps))
    return e


yyyymmdd = sys.argv[1]
hh = sys.argv[2]
type = sys.argv[3]

narr_path = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/narr-a/YYYYMM/YYYYMMDD/narr-a_221_YYYYMMDD_HH00_000.grb'
narr_path = narr_path.replace('YYYYMMDD', yyyymmdd)
narr_path = narr_path.replace('YYYYMM', yyyymmdd[:6])
narr_path = narr_path.replace('HH', hh)
print narr_path

d = Dataset(narr_path)
#for i in d.variables.keys():
print d.variables.keys()


ll = Dataset('narr_lat_lon.nc')
lat = ll.variables['lat'][:]
lon = ll.variables['lon'][:]
ll.close()

if type == 'sfc':
    dt = datetime.strptime(yyyymmdd+hh, '%Y%m%d%H')
    dt_str = datetime.strftime(dt, '%d %b %Y  %H UTC')
    mslp = d.variables['Mean_sea_level_pressure_ETA_model'][0,:,:]/100.
    m = regMap()
    
    title(dt_str + ' ' + 'Surface NARR-A', fontsize=15)
    sfc_temp = d.variables['Temperature_height_above_ground'][0,0,:,:]
    u_wind = d.variables['u_wind_height_above_ground'][0,0,:,:] * 1.94384
    v_wind = d.variables['v_wind_height_above_ground'][0,0,:,:] * 1.94384
    x,y = m(lon, lat)
    CS = m.contour(x,y, mslp, np.arange(940,1104,4), colors='k', linewidths=2)
    clabel(CS, CS.levels, fmt='%4.0f')
    stride = 5
    sfc_temp =((sfc_temp - 273.15)*1.8 + 32)
    cb = m.contourf(x,y,sfc_temp, np.arange(-40,132,2), cmap=get_cmap("RdYlBu_r"))
    fz = m.contour(x,y,sfc_temp, np.asarray([32]), colors='m', linestyles='--', linewidths=2)
    clabel(fz, fz.levels, fmt='%4.0f')
    barbs(x[::stride,::stride],y[::stride,::stride],u_wind[::stride,::stride], v_wind[::stride,::stride])
    divider = make_axes_locatable(gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = colorbar(cb, cax=cax)
    cb.set_label("Temperature [F]")
    tight_layout()
    #show()
    savefig(yyyymmdd + '.' + hh + '.sfc.png')

if type == 'sfccnt':
    dt = datetime.strptime(yyyymmdd+hh, '%Y%m%d%H')
    dt_str = datetime.strftime(dt, '%d %b %Y  %H UTC')
    mslp = d.variables['Mean_sea_level_pressure_ETA_model'][0,:,:]/100.
    m = regMap()
    
    title(dt_str + ' ' + 'Surface NARR-A', fontsize=15)
    sfc_temp = d.variables['Temperature_height_above_ground'][0,0,:,:]
    u_wind = d.variables['u_wind_height_above_ground'][0,0,:,:] * 1.94384
    v_wind = d.variables['v_wind_height_above_ground'][0,0,:,:] * 1.94384
    sh = d.variables['Specific_humidity_height_above_ground'][0,0,:,:]
    pres = d.variables['Pressure'][0,0,:,:]/100.
    e = w2e(q2w(sh), pres)
    dwpt = e2td(e)
    x,y = m(lon, lat)
    dwpt = dwpt*1.8 + 32.
    cb = m.contourf(x,y,dwpt, np.arange(50,82,2), cmap=get_cmap('Greens'))
    CS = m.contour(x,y, mslp, np.arange(940,1104,4), colors='k', linewidths=2)
    clabel(CS, CS.levels, fmt='%4.0f')
    
    sfc_temp =((sfc_temp - 273.15)*1.8 + 32)
    sfc_temp = scipy.ndimage.gaussian_filter(sfc_temp, 1.5)

    temp_levels = np.arange(-60,150,10)
    zero_level = np.where(temp_levels == 0)[0]
    temp_colors = np.repeat('r', len(temp_levels))
    temp_colors[zero_level] = 'm'
    temp_colors[:zero_level] = 'b'
    stride=5 
    tm = m.contour(x,y,sfc_temp, temp_levels, colors=temp_colors, linewidths=2, linestyles='--')
    clabel(tm, tm.levels, fmt='%4.0f')
    barbs(x[::stride,::stride],y[::stride,::stride],u_wind[::stride,::stride], v_wind[::stride,::stride])
    divider = make_axes_locatable(gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = colorbar(cb, cax=cax)
    cb.set_label("Dewpoint [F]")
    tight_layout()
    #show()
    savefig(yyyymmdd + '.' + hh + '.sfc.png')

if type == 'svr':
    m = regMap()
    mslp = d.variables['Mean_sea_level_pressure_ETA_model'][0,:,:]
    dt = datetime.strptime(yyyymmdd+hh, '%Y%m%d%H')
    dt_str = datetime.strftime(dt, '%d %b %Y  %H UTC')


    v_sfc = d.variables['v_wind_height_above_ground'][0,:,:] 
    u_sfc = d.variables['u_wind_height_above_ground'][0,:,:] 
    
    # 700 temperature:
    pres_idx = np.where(d.variables['isobaric'][:] == 700)[0]
    temp7 = d.variables['Temperature'][0,pres_idx, :,:]-273.15
    z7 = d.variables['Geopotential_height'][0, pres_idx,:,:]/1000.
    print z7, temp7

    # 500 temp:
    pres_idx = np.where(d.variables['isobaric'][:] == 500)[0]
    temp5 = d.variables['Temperature'][0,pres_idx, :,:]-273.15
    u5 = d.variables['u_wind'][0,pres_idx, :,:]
    v5 = d.variables['v_wind'][0,pres_idx, :,:]
    z5 = d.variables['Geopotential_height'][0, pres_idx,:,:]/1000.
    print z5, temp5

    cape = d.variables['Convective_available_potential_energy'][0,:,:][0]

    lr = ((temp7 - temp5) / (z7 - z5))[0]

    u_shear = 1.943 * (u5[0] - u_sfc)
    v_shear = 1.943 * (v5[0] - v_sfc)
    mag = np.sqrt(np.power(u_shear, 2) + np.power(v_shear, 2))
    u_shear = np.ma.masked_where(mag < 30, u_shear)[0]
    v_shear = np.ma.masked_where(mag < 30, v_shear)[0]

    sh = d.variables['Specific_humidity_height_above_ground'][0,0,:,:]
    pres = d.variables['Pressure'][0,0,:,:]/100.
    e = w2e(q2w(sh), pres)
    dwpt = e2td(e)
    x,y = m(lon, lat)
    dwpt = dwpt*1.8 + 32.

    sfc_dwpt = scipy.ndimage.gaussian_filter(dwpt, 1.5)
    cb = m.contour(x,y,sfc_dwpt, np.arange(50,82,2), colors='g')
    clabel(cb, cb.levels, fmt='%4.0f') 
    CS = m.contour(x,y, mslp, np.arange(940,1104,4), colors='k', linewidths=3)
    clabel(CS, CS.levels, fmt='%4.0f')
    
    title(dt_str + ' ' + 'Surface NARR-A', fontsize=15)
    stride=5 
    lr_levels = -1 * np.asarray([5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10])
    cape_levels = np.arange(500,6500,500)
    cb = m.contourf(x,y,cape, cape_levels, cmap='spring_r')
    barbs(x[::stride,::stride],y[::stride,::stride],u_shear[::stride,::stride], v_shear[::stride,::stride])
    divider = make_axes_locatable(gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = colorbar(cb, cax=cax)
    cb.set_label("CAPE [J/kg]")
    tight_layout()
    savefig(yyyymmdd + '.' + hh + '.svr.png')
    show()

def plotUA(level):
    dt = datetime.strptime(yyyymmdd+hh, '%Y%m%d%H')
    dt_str = datetime.strftime(dt, '%d %b %Y  %H UTC')
    pres_idx = np.where(d.variables['isobaric'][:] == level)[0]
    sfc_press = d.variables['Pressure'][0,0, :,:]
    mask = np.where(sfc_press/100. < level)
    
    temp = d.variables['Temperature'][0,pres_idx, :,:]-273.15
    rh = d.variables['Specific_humidity'][0, pres_idx, :, :]
    z = d.variables['Geopotential_height'][0, pres_idx,:,:]
    u = d.variables['u_wind'][0,pres_idx, :,:]*1.94384
    v = d.variables['v_wind'][0, pres_idx, :,:]*1.94384
    #print mask
    temp[0][mask] = np.nan
    rh[0][mask] = np.nan
    z[0][mask] = np.nan
    u[0][mask] = np.nan
    v[0][mask] = np.nan
    terrain = np.ones(v[0].shape)
    terrain[mask] = 0.
    #print temp.shape
    m = regMap()
    title(dt_str + ' ' + str(level) + ' mb NARR-A', fontsize=15)
    #print rh.shape, z.shape
    es = cc(temp)
    e = w2e(q2w(rh), level)
    rh = (e/es)*100.
    tight_layout()
    x,y = m(lon,lat)
    contourf(x,y,terrain, [0,.5,1], colors=['k','#FFFFFF'], alpha=.2)
    if level > 500:
        cb = m.contourf(x,y,rh[0], np.arange(70,105,5), cmap=get_cmap('Greens'))
    else:
        wind_spd = np.sqrt(np.power(u[0],2) + np.power(v[0], 2))
        cb = m.contourf(x,y,wind_spd, np.arange(60,240,10), cmap=get_cmap('Reds'), alpha=.8)
        m.barbs(x[::5,::5],y[::5,::5],u[0,::5,::5], v[0,::5,::5])
    #CS = m.contour(x,y,z[0], np.arange(5160-(60*10), 5880+(60*10),60), colors='k', linewidths=2)
    CS = m.contour(x,y,z[0], colors='k', linewidths=2)
    clabel(CS, CS.levels, fmt='%4.0f')
    temp_levels = np.arange(-80,82,2)
    zero_level = np.where(temp_levels == 0)[0]
    temp_colors = np.repeat('r', len(temp_levels))
    temp_colors[zero_level] = 'm'
    temp_colors[:zero_level] = 'b'
    CS = m.contour(x,y,temp[0], temp_levels, colors=temp_colors, linestyles='--', linewidths=1.5)
    clabel(CS, CS.levels, fmt='%4.0f')
    divider = make_axes_locatable(gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = colorbar(cb, cax=cax)
    cb.set_label("Relative Humidity [%]")
    savefig(yyyymmdd + '.' + hh + '.ua.' + str(level) + '.png') 
    #show()
print type
print int(type)
try:
    type = int(type)
    plotUA(type)
    d.close()
except Exception,e:
    print "Not an UA map."
    print e
    sys.exit()
    m = regMap()
    sfc_temp = d.variables['Temperature'][0,0,:,:]
    u_wind = d.variables['u_wind_height_above_ground'][0,0,:,:] * 1.94384
    v_wind = d.variables['v_wind_height_above_ground'][0,0,:,:] * 1.94384
    x,y = m(lon, lat)
    CS = m.contour(x,y, mslp, np.arange(940,1104,4), colors='k', linewidths=2)
    clabel(CS, CS.levels, fmt='%4.0f')
    stride = 5
    sfc_temp =((sfc_temp - 273.15)*1.8 + 32)
    cb = m.contourf(x,y,sfc_temp, np.arange(-40,132,2), cmap=get_cmap("jet"))
    barbs(x[::stride,::stride],y[::stride,::stride],u_wind[::stride,::stride], v_wind[::stride,::stride])
    colorbar(cb)
    tight_layout()
    savefig(yyyymmdd + '.' + hh + '.sfc.png')


    


if type == 'svr':
    mslp = d.variables['Mean_sea_level_pressure_ETA_model'][0,:,:]
    dt = datetime.strptime(yyyymmdd+hh, '%Y%m%d%H')
    dt_str = datetime.strftime(dt, '%d %b %Y  %H UTC')
    
    # 700 temperature:
    pres_idx = np.where(d.variables['isobaric'][:] == 700)[0]
    temp7 = d.variables['Temperature'][0,pres_idx, :,:]-273.15
    z7 = d.variables['Geopotential_height'][0, pres_idx,:,:]
    
    # 500 temp:
    pres_idx = np.where(d.variables['isobaric'][:] == 700)[0]
    temp5 = d.variables['Temperature'][0,pres_idx, :,:]-273.15
    z5 = d.variables['Geopotential_height'][0, pres_idx,:,:]
    

    lr = (temp7 - temp5) / (z7 - z5)
    print lr

#m = regMap()
#x,y = m(lon, lat)
#CS = m.contour(x,y, mslp/100., np.arange(940,1104,4), colors='k', linewidths=2)
#clabel(CS, CS.levels, fmt='%4.0f')
#m.contourf(x,y,cape, np.arange(0,6250,250), cmap=get_cmap('Reds'))
#colorbar()
#tight_layout()
#show()
