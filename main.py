import numpy as np 
from xarray_wrapper import wrapper as pp
import utilities as ut
import digest_input as input
from scipy.integrate import simps
import xarray as xr 

def main():
    gcm_output = 'diagfi80.nc'
    parfile = 'w43b.par'
    ## Initialisation of the planets and reading of parfile
    myp=input.Planet(parfile)
    myp.ini()
    rp = input.RunParameter(parfile)
    print("getting pressure")
    targetpress1D = np.logspace(np.log10(rp.p_lower),np.log10(rp.p_upper),rp.nlev)
    press = pp(file=gcm_output,var='p',x=rp.charx).getf()
    print("getting u and coordinates")
    u,xdim,ydim,zdim,tdim = pp(file=gcm_output,var='u',x=rp.charx).getfd() #zonal mean of zonal wind
    print("getting {}".format(rp.temperature_field))
    temp = pp(file=gcm_output,var=rp.temperature_field,x=rp.charx).getf() #zonal mean of temperature
    print("getting potential temperature")
    tpot = pp(file=gcm_output,var='teta',x=rp.charx).getf()#zonal mean of potential temperature
    print("getting v")
    v = pp(file=gcm_output,var='v',x=rp.charx).getf() #zonal mean of meridional wind
    if rp.is_omega:
        print("getting vertical wind")
        o = pp(file=gcm_output,var='w',x=rp.charx).getf() ##need to add a switch on 'w' because in dynamico the field is name "omega". #zonal mean of vertical wind
    print("Computing coupled terms")
    print("staru4D")
    staru4D = pp(gcm_output,var='u',compute='pert_x',x=rp.charx).getf() #u*
    print("starv4D")
    starv4D = pp(gcm_output,var='v',compute='pert_x',x=rp.charx).getf() #v*
    print("starT4D")
    starT4D = pp(gcm_output,var=rp.temperature_field,compute='pert_x',x=rp.charx).getf() #T*
    print("starTpot4D")
    starTpot4D = pp(gcm_output,var='teta',compute='pert_x',x=rp.charx).getf() #Tpot*
    if rp.is_omega:
        print("starw4D")
        starw4D = pp(gcm_output,var='w',compute='pert_x',x=rp.charx).getf() #w*
        wpup = ut.mean(starw4D*staru4D,axis=3) 
        wpTp = ut.mean(starw4D*starT4D,axis=3) 
        # del starw4D
    vpup = ut.mean(starv4D*staru4D,axis=3) 
    vpTp = ut.mean(starv4D*starT4D,axis=3) 
    upup = ut.mean(staru4D*staru4D,axis=3) 
    vpvp = ut.mean(starv4D*starv4D,axis=3) 
    
    # del staru4D,starv4D,starT4D
    print("Done with getting fields")
    print("Now interpolating stuff")
    u = ut.interpolate(targetpress1D,press,u,spline=rp.use_spline)
    tpot = ut.interpolate(targetpress1D,press,tpot,spline=rp.use_spline)
    v = ut.interpolate(targetpress1D,press,v,spline=rp.use_spline)
    temp =ut.interpolate(targetpress1D,press,temp,spline=rp.use_spline)
    vpup = ut.interpolate(targetpress1D,press,vpup,spline=rp.use_spline)
    vpTp = ut.interpolate(targetpress1D,press,vpTp,spline=rp.use_spline)
    if rp.is_omega:
        o = ut.interpolate(targetpress1D,press,o,spline=rp.use_spline)
        wpup = ut.interpolate(targetpress1D,press,wpup,spline=rp.use_spline)
        wpTp = ut.interpolate(targetpress1D,press,wpTp,spline=rp.use_spline)
    eke = ut.interpolate(targetpress1D,press,0.5*(vpvp+upup),spline=rp.use_spline)
    print("Done with interpolation")
    
    ### Dimensions
    nt,nz,nlat = u.shape 
    print(nt,nz,nlat)
    nlon = 1
    if rp.includels:
        day_per_year = np.ceil(myp.dayperyear())
        assert(day_per_year == rp.day_per_year)
        tdim = ut.fix_time_axis(tdim,day_per_year)
        lstab = np.zeros(nt)

    ## vertical coordinates 
    pseudoz = myp.pseudoz(targetpress1D,p0=targetpress1D[0]+1.)
    targetpress3D = np.tile(targetpress1D,(nt,1))
    targetpress3D = np.tile(np.transpose(targetpress3D),(nlat,1,1))
    targetpress3D = np.transpose(targetpress3D)
    print(targetpress3D.shape)
    
    ## Computing curvatures terms 
    lat2d = np.tile(ydim,(nz,1))
    acosphi2d = myp.acosphi(lat=lat2d)
    cosphi2d = acosphi2d/myp.a 
    latrad,lat2drad = ydim*np.pi/180.,lat2d*np.pi/180. #degrees to radians
    beta = myp.beta(lat=lat2d)
    f = myp.fcoriolis(lat=lat2d)
    tanphia = myp.tanphia(lat=lat2d)
    print("Computed coordinates and curvature terms")
    
    ## Angular momentum => this is zonal mean everywhere
    ## Careful => we use shallow atmosphere + hydrostatic balance approximation here
    dlat = np.abs(latrad[1]-latrad[0])
    dlon = 2*np.pi 
    dp = np.gradient(targetpress3D,axis=1,edge_order=2) #compute the gradient of pressure dp =>it's negative here
    dm = - myp.a*acosphi2d*dlon*dlat*dp/myp.g #this is the mass for each grid mesh. WARNING: here the gravity is constant with latitude
    wangmomperumass = myp.wangmom(u=u,lat=lat2d) #wind angular momentum per dm
    angmomperumass = myp.angmom(u=u,lat=lat2d) # total axial angular momentum per dm 
    ## Units below is kg.m2.s-1. E25
    ## normalisation to get rid of possible overflow issues
    angmom = dm*angmomperumass / 1.e25##
    wangmom = dm * wangmomperumass / 1.e25 #
    superindex = myp.superrot(u=u,lat=lat2d) # not sure I will use it but you never know
    
    ##Create output file and start to fill it with variables that will ALWAYS be outputted
    ds=ut.create_output_file(tdim,pseudoz,ydim)
    xar = xr.DataArray(targetpress3D,
                       coords={'Time':tdim,"pseudoalt":pseudoz,'lat':ydim},
                       dims=['Time','pseudoalt','lat'],
                       attrs={"units":"Pa","long_name":"Pressure"})
    ds['p'] = xar 
    xar = xr.DataArray(temp,
                       coords={'Time':tdim,"pseudoalt":pseudoz,'lat':ydim},
                       dims=['Time','pseudoalt','lat'],
                       attrs={"units":"K","long_name":"Temperature"})
    ds[rp.temperature_field]=xar 
    xar = xr.DataArray(tpot,
                       coords={'Time':tdim,"pseudoalt":pseudoz,'lat':ydim},
                       dims=['Time','pseudoalt','lat'],
                       attrs={"units":"K","long_name":"Potential Temperature"})
    ds['tpot']=xar 
    xar = xar = xr.DataArray(u,
                       coords={'Time':tdim,"pseudoalt":pseudoz,'lat':ydim},
                       dims=['Time','pseudoalt','lat'],
                       attrs={"units":"m/s","long_name":"Zonal wind"})
    ds['u']=xar
    xar = xr.DataArray(angmom,
                       coords={'Time':tdim,"pseudoalt":pseudoz,'lat':ydim},
                       dims=['Time','pseudoalt','lat'],
                       attrs={"units":"E25 kg.m2.s-1","long_name":"Axial Angular Momentum"})
    ds['angmom'] = xar 
    xar = xr.DataArray(wangmom,
                       coords={'Time':tdim,"pseudoalt":pseudoz,'lat':ydim},
                       dims=['Time','pseudoalt','lat'],
                       attrs={"units":"E25 kg.m2.s-1","long_name":"Axial Angular Momentum due to wind"})
    ds['wangmom'] = xar 
    
    
    if rp.compute_transport:
        ## Go into Mendonça 2020 diagnostics for transport
        ##1st, let's do Mean Meridional circulation : [Abar][Bbar], bar is temporal mean
        mmc_v_ang = ut.mean(v,axis=0) * ut.mean(angmom,axis=0) ## temporal mean of each fields (that where already zonal mean) and product
        mmc_w_ang = ut.mean(o,axis=0) * ut.mean(angmom,axis=0) ##same as line above for vertical wind
        mmc_v_tpot = ut.mean(v,axis=0) * ut.mean(tpot,axis=0)
        mmc_w_tpot = ut.mean(o,axis=0) * ut.mean(tpot,axis=0)
        ## I have a runtime warning, mean of empty slices apparently
        xar = xr.DataArray(mmc_v_ang,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":"E25 kg.m3.s-2","long_name":"Meridional transport of aam by mean circulation"})
        ds['mmc_v_ang']=xar
        xar = xr.DataArray(mmc_w_ang,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":"E25 kg.m3.s-2","long_name":"Vertical transport of aam by mean circulation"})
        ds['mmc_w_ang']=xar
        
        xar = xr.DataArray(mmc_v_tpot,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":" m.K.s-1","long_name":"Meridional transport of tpot by mean circulation"})
        ds['mmc_v_tpot']=xar
        xar = xr.DataArray(mmc_w_tpot,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":" m.K.s-1","long_name":"Vertical transport of tpot by mean circulation"})
        ds['mmc_w_tpot']=xar
        ## Now let's look at the stationnary waves: [Abar*Bbar*]
        ##need to recompute the angular momentum completely here and create starangmom4D => not sure how to do that
        ##for Tpot
        ## fields have not been interpolated on pressure grid here ! 
        sw_v_tpot = ut.mean(ut.mean(starv4D,axis=0) * ut.mean(starTpot4D,axis=0),axis=-1) #last axis is suppose to be longitude
        sw_w_tpot = ut.mean(ut.mean(starw4D,axis=0) * ut.mean(starTpot4D,axis=0),axis=-1) #last axis is suppose to be longitude
        starangmom4D = myp.angmom(u=staru4D,lat= np.tile(ydim[None,:,None],(nz,1,staru4D.shape[-1])))
        sw_v_ang = ut.mean(ut.mean(starv4D,axis=0) * ut.mean(starangmom4D,axis=0),axis=-1) #last axis is suppose to be longitude
        sw_w_ang = ut.mean(ut.mean(starw4D,axis=0) * ut.mean(starangmom4D,axis=0),axis=-1) #last axis is suppose to be longitude
        xar = xr.DataArray(sw_v_tpot,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":" m.K.s-1","long_name":"Meridional transport of tpot by stationnary waves"})
        ds['sw_v_tpot']=xar
        xar = xr.DataArray(sw_w_tpot,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":" m.K.s-1","long_name":"Vertical transport of tpot by stationnary waves"})
        ds['sw_w_tpot']=xar
        xar = xr.DataArray(sw_v_ang,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":" E25 kg.m2.s-2","long_name":"Meridional transport of aam by stationnary waves"})
        ds['sw_v_ang']=xar
        xar = xr.DataArray(sw_w_ang,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":" E25 kg.m2.s-2","long_name":"Vertical transport of aam by stationnary waves"})
        ds['sw_w_ang']=xar
        # ## Last but not least, the transient waves (eddies): [(A'B')bar] => the temporal mean is done on the product
        primeu4D = pp(gcm_output,var='u',compute='pert_t',x=rp.charx).getf() ##u'
        primev4D = pp(gcm_output,var='v',compute='pert_t',x=rp.charx).getf() ##v'
        primew4D = pp(gcm_output,var='w',compute='pert_t',x=rp.charx).getf() ##w'
        primeTpot4D = pp(gcm_output,var='teta',compute='pert_t',x=rp.charx).getf() ##Tpot'
        primeangmom4D = myp.angmom(u=primeu4D,lat=np.tile(ydim[None,:,None],(nz,1,primeu4D.shape[-1]))) #here we use a tiling to 3D dim since primeu4D is a 4D field and we want to use myp.angmom()
        tw_v_ang = ut.mean(primev4D*primeangmom4D,axis=(0,-1))
        tw_w_ang = ut.mean(primew4D*primeangmom4D,axis=(0,-1))
        tw_v_tpot = ut.mean(primev4D*primeTpot4D,axis=(0,-1))
        tw_w_tpot = ut.mean(primew4D*primeTpot4D,axis=(0,-1))
        xar = xr.DataArray(tw_v_ang,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":" E25 kg.m2.s-2","long_name":"Meridional transport of aam by transient waves"})
        ds['tw_v_ang']=xar 
        xar = xr.DataArray(tw_w_ang,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":" E25 kg.m2.s-2","long_name":"Vertical transport of aam by transient waves"})
        ds['tw_w_ang']=xar 
        xar = xr.DataArray(tw_v_tpot,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":" m.K.s-1","long_name":"Meridional transport of tpot by transient waves"})
        ds['tw_v_tpot']=xar 
        xar = xr.DataArray(tw_w_tpot,
                       coords={"pseudoalt":pseudoz,'lat':ydim},
                       dims=['pseudoalt','lat'],
                       attrs={"units":" m.K.s-1","long_name":"Vertical transport of tpot by transient waves"})
        ds['tw_w_tpot']=xar
    ds.to_netcdf("test.nc")
    ## We're done with transport diagnostics
    # print("Computing streamfunction now ")
    # psim = np.zeros((nt,nz,nlat))
    # alph = -2*np.pi*acosphi2d/myp.g
    # w = np.isnan(v)
    # v[w]=0. # so integrations will not fail later on
    # x = targetpress1D[:]
    # for tt in range(nt):
    #     for lla in range(nlat):
    #         y = v[tt,:,lla] #integrand is y
    #         # y = np.append(y,y[-1]) #we integrate from p=0 toward P
    #         for zz in range(0,nz):
    #             psim[tt,zz,lla] = -simps(y[zz:],x[zz:])*alph[0,lla] # Minus sign because x is ordered by decreasing values of P and integral is from 0 to P
    # print("We're done with scipy")
    # v[w]=np.nan #we put the NaN back in v
    # psim[w] = np.nan #we put the NaN in the streamfunction as well
    # print("Creating output file")
    
if __name__=="__main__":
    main()