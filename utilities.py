import numpy as np 
from scipy.interpolate import interp1d


def interpolate(targetp1d,sourcep3d,fieldsource3d,spline=False):
    """targetp1d: 1D pressure field on which the interpolation is made
       sourcep3d: 3D pressure field from the GCM (dims are time,altitude,latitude)
       fieldsource3d: 3D field we want to interpolate on the 1D pressure grid (dims are time,altitude,latitude)
        spline (bool, optional): Wether to spline interpolate or not. Defaults to False.
    """
    nt,_,nlat = fieldsource3d.shape
    coordsource3d = -np.log(sourcep3d) #interpolatine in log(p)
    coordtarget1d = -np.log(targetp1d) #interpolate in log(p)
    nzt = coordtarget1d.size 
    fieldtarget3d = np.zeros((nt,nzt,nlat)) ##this is the results of the computation
    for lla in range(nlat):
        for tt in range(nt):
            xs = coordsource3d[tt,:,lla]
            ys = fieldsource3d[tt,:,lla]
            if not spline: #we use numpy interp function
                fieldtarget3d[tt,:,lla] = np.interp(coordtarget1d,xs,ys,left=np.nan,right=np.nan) #when interpolating, setting the border at NaN
            else: #we use scipy interpolation
                kind='linear' #could also be "cubic" or "quadratic"
                ff= interp1d(xs,ys,kind=kind,bounds_error=False)
                fieldtarget3d[tt,:,lla] = ff(coordtarget1d)
    return fieldtarget3d

def fix_time_axis(tdim,period):
    """This is only used to include Ls in the final diagfi. 
       Is this really useful ? Not sure
    """
    nt = tdim.size
    tdim = tdim % period 
    nperiod = 0
    corrected_tdim = np.empty(nt)
    corrected_tdim[0] = float(tdim[0]/period)
    for ii in range(1,nt):
        if tdim[ii]-tdim[ii-1] <0:
            nperiod +=1
        corrected_tdim[ii] = float(nperiod) + float(tdim[ii]/period)
    return corrected_ttdim

def correctnearzero(field):
    """Correct smal value to nan to get rid of numerical instabilities
    """
    negval = np.min(field[np.isfinite(field)]) #take min value where not Nan or Inf
    val = 2.
    removed = -val*negval
    w = np.where(np.abs(field)<=removed)
    field[w]=np.nan
    print("absolute values below this value are set to Nan", removed)
    return field

def mean(field,axis=None):
    """Just an elaborate way of taking the mean when there's Nan or masked_array

    Args:
        field (_type_): _description_
        axis (_type_, optional): _description_. Defaults to None.
    """
    if field is None:
        return None 
    if type(field).__name__=="MaskedArray":
        ##When using masked array, replace the mask with NaN
        field.set_fill_value(np.NaN)
        zout = np.ma.array(field).mean(axis=axis,dtype=np.float64)
        if axis is not None:
            zout.set_fill_value(np.NaN)
            return zout.filled()
        else:
            return zout
    else:
        return np.nanmean(field,axis=axis,dtype=np.float64)