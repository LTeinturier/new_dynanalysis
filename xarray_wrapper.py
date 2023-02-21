import xarray as xr
#################################################
#   This is modified from the fakepp.py script of
#   the original dynanalysis, from AS
#################################################
class wrapper():
    
    def __init__(self,
                 file=None,
                 var=None,
                 x=None,
                 y=None,
                 z=None,
                 t=None,
                 compute=None):
        self.file=file,
        self.var=var
        self.x=x
        self.y=y
        self.z=z
        self.t=t 
        self.compute = compute
        
    def open(self):
        ds = xr.open_dataset(self.file[0],decode_times=False)
        ##check the name of the coords
        if 'lon' in ds.coords:
            self.zonadim='lon'
        else:
            self.zonadim='longitude'
        if 'lat' in ds.coords:
            self.meridim='lat'
        else:
            self.meridim='latitude'
        if 'Time' in ds.coords:
            self.timedim='Time'
        else:
            self.timedim='time_counter'
        self.altidim='altitude'
        return ds
    
    def get(self):
        ##open dataset
        ds= self.open()
        #prepare stuff
        dalist = []
        dalistm = []
        if self.x is not None:
            if self.x in ["0,360","-180,180"]:
                dalistm.append( self.zonadim )
            else:
                dalist.append( (self.zonadim,int(self.x)) ) 
                
        if self.y is not None:
            dalist.append((self.meridim,int(self.y)))
        if self.z is not None:
            dalist.append((self.altidim,int(self.z)))
        if self.t is not None:
            dalist.append((self.timedim,int(self.t)))
        dsred = ds.sel(dict(dalist),method='nearest')[self.var]    ##fields value 
        ##averaging or computing deviation to zonal mean
        ## if list stayed [] before, nothing is done.
        if self.compute is not None:
            if "pert_x" in self.compute:
                print("computation X-[X]")
                mm=ds[self.var].mean(self.zonadim)
                dsred = dsred-mm
        else:
            dsred = dsred.mean(dalistm)
        return dsred
    
    def getf(self): 
        return self.get().values #return only the values of the field
    
    def getfd(self):
        ds=self.open()
        xd = ds.coords[self.zonadim]
        yd = ds.coords[self.meridim]
        zd = ds.coords[self.altidim]
        td = ds.coords[self.timedim]
        return self.getf(),xd.values,yd.values,zd.values,td.values ##returns value of the field + coordinates
    
    
if __name__=="__main__":
    gcm_output='diagfi80.nc'
    charx="-180,180"
    press=wrapper(file=gcm_output,var='p',x=charx).getf()
    print(press.shape)
#     import os 
#     path = '../wasp43b/'
#     file='diagfi80.nc'#os.path.join(path,'diagfi80.nc')
#     print(file)
#     charx="-180,180"
#     charx='12'
#     u=wrapper(file=file,var="u",x=charx,compute='pert_x').getf()
#     print(u.shape)
#     u=wrapper(file=file,var="u",x=charx).getf()
#     print(u.shape)
#     u=wrapper(file=file,var="u").getf()
#     print(u.shape)
#     u,x,y,z,t=wrapper(file=file,var="u").getfd()
#     print(u.shape)
    # print(x)
    # print(y)
    # print(z)
    # print(t)