import numpy as np 
from xarray_wrapper import wrapper as pp
import utilities as u 
import digest_input as input

def main():
    gcm_output='diagfi80.nc'
    parfile = 'w43b.par'
    myp=input.Planet(parfile)
    myp.ini()
    print("getting pressure")
    rp = input.RunParameter(parfile)
    press= pp(file=gcm_output,var='p',x=rp.charx).getf()
    print(press.shape)
    print("getting u and coordinates")
    u,xdim,ydim,zdim,tdim = pp(file=gcm_output,var='u',x=rp.charx).getfd()
    print(u.shape,xdim.size)
    print("getting {}".format(rp.temperature_field))
    temp = pp(file=gcm_output,var=rp.temperature_field,x=rp.charx).getf()
    print("getting potential temperature")
    tpot = pp(file=gcm_output,var='teta',x=rp.charx).getf()
    print("getting v")
    v = pp(file=gcm_output,var='teta',x=rp.charx).getf()
    if rp.is_omega:
        print("getting vertical wind")
        o=pp(file=gcm_output,var='w',x=rp.charx).getf() ##need to add a switch on 'w' because in dynamico the field is name "omega"
    print("Computing coupled terms")
    print("staru4D")
    staru4D=pp(gcm_output,var='u',compute='pert_x',x=rp.charx).getf()
    print("starv4D")
    starv4D=pp(gcm_output,var='v',compute='pert_x',x=rp.charx).getf()
    print("starT4D")
    starT4D=pp(gcm_output,var=rp.temperature_field,compute='pert_x',x=rp.charx).getf()
    if rp.is_omega:
        print("starw4D")
        starw4D=pp(gcm_output,var='w',compute='pert_x',x=rp.charx).getf()
        wpup=u.mean(starw4d*staru4D,axis=3) 
        wpTp=u.mean(starw4d*starT4D,axis=3) 
        del starw4D
    vpup=u.mean(starv4D*staru4D,axis=3) 
    vpTp=u.mean(starv4D*starT4D,axis=3) 
    upup=u.mean(staru4D*staru4D,axis=3) 
    vpvp=u.mean(starv4D*starv4D,axis=3) 
    
    del staru4D,starv4D,starT4D
    print("Done with getting fields")
    print("Now interpolating stuff")
    
if __name__=="__main__":
    main()