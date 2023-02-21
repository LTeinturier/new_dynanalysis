import os 
import numpy as np 
import scipy.constants as cst 
from configobj import ConfigObj
import sys 
from distutils.util import strtobool

#-------------Basic physical constants-------------------
G = cst.G          # Gravitational constant
h = cst.h          # Planck's constant
sigma = cst.sigma  # Stefan-Boltzman constant
k = cst.k          # Boltzman thermodynamic constant
c = cst.c          # Speed of light
#-----------Thermodynamic constants----------------------
Avogadro = cst.Avogadro  # Avogadro's number
#Following will come out in J/(deg kmol), so
#that dividing Rstar by molecular weight gives
#gas constant appropriate for mks units
Rstarkilo = 1000.*k*Avogadro   #Universal gas constant
#-------------Useful planetary quantities----------------
astronomical_unit = cst.astronomical_unit       # astronomical unit in meters

class Planet():
    def __init__(self,file):
        desc={}
        self.file=file
        desc["file"]="Input file"
        self.name = None
        desc["name"] = "Name of the planet"
        self.a = None 
        desc["a"] = "Mean radius of planet (m)"
        self.mass = None 
        desc["mass"] = "Mass of planet (kg)"
        self.g = None 
        desc["g"] = "Surface gravitational acceleration (m/s**2)"
        self.L = None 
        desc["L"] = "Annual mean solar constant (current) (W/m**2)"
        self.albedo = None 
        desc["albedo"] = "Bond albedo (fraction)"
        self.rsm = None 
        desc["rsm"] = "Semi-major axis of orbit about Sun (m)"
        self.year = None 
        desc["year"] = "Sidereal length of year (s)"
        self.eccentricity = None 
        desc["eccentricity"] = "Eccentricity (unitless)"
        self.day = None 
        desc["day"] = "Mean tropical length of day (s)"
        self.obliquity = None 
        desc["obliquity"] = "Obliquity to orbit (degrees)"
        self.Lequinox = None 
        desc["Lequinox"] = "Longitude of equinox (degrees)"
        self.Tsbar = None 
        desc["Tsbar"] = "Mean surface temperature (K)"
        self.Tsmax = None 
        desc["Tsmax"] = "Maximum surface temperature (K)"
        self.M = None 
        desc["M"] = "Molecular weight (g/mol)"
        self.cp = None 
        desc["cp"] = "Specific heat capacity (J kg-1 K-1)"
        self.T0 = None 
        desc["T0"] = "Typical atmospheric temperature (K)"
        self.incl = None 
        desc["incl"] = "Orbit inclination (deg)"
        self.ascend = None 
        desc["ascend"] = "Longitude of ascending node (deg)"
        self.omeg = None 
        desc["omeg"] = "Argument of periapsis (deg)"
        self.date_peri = None 
        desc["date_peri"] = "Date of perihelion"
        self.date_equi = None 
        desc["date_equi"] = "Date of equinox"
        ## calculated
        self.R = None 
        desc["R"] = "planetary gas constant"
        self.dryadiab = None 
        desc["dryadiab"] = "dry adiabatic lapse rate"
        self.omega = None 
        desc["omega"] = "planetary rotation rate"
        self.density = None 
        desc["density"] = "density (kg m-3)"
        desc['desc']='dic'
        self.desc=desc
        
    def calculate(self):
        if self.M is not None:
            self.R = Rstarkilo/self.M
        else:
            self.R = None 
         # adiabatic lapse rate
        if self.cp is not None:
          self.dryadiab = self.g/self.cp
        else:
          self.dryadiab = None
        # planetary rotation rate
        self.omega = 2.*np.pi/self.day
        # density (assuming spherical shape)
        if self.mass is not None:
          self.density = self.mass / ((4./3.)*np.pi*(self.a**3))
        else:
          self.density = None
        print("Calculated remaining planetary parameters")
    def show(self):
        for k,v in list(vars(self).items()):
            print(k,v,self.desc[k])
            
    def fromfile(self):
        config=ConfigObj(self.file)
        for keys in config['Planet'].keys():
            if keys=='name':
                setattr(self,keys,config['Planet'][keys])
            else:
                setattr(self,keys,float(config['Planet'][keys]))
        print("Read the planetary parameter from {}".format(self.file))
        
    def ini(self):
        self.fromfile()
        # self.show()
        self.calculate()
    ###Computation used later
    def disk(self):
        """ Planetary disk area"""
        return np.pi*self.a*self.a
    
    def fcoriolis(self,lat=45.):
        """Coriolis parameter.

        Args:
            lat (float, degrees): latitude in deg. Defaults to 45..
        """
        return 2*self.omega*np.sin(np.deg2rad(lat))
    
    def Rossby(self,U=50.,L=None):
        """ Rossby number"""
        if L is None:
            L = self.a/2 #longueur caracteristique
        return U/(self.fcoriolis()*L)
        
    def H(self,T0=None,M=None):
        """scale height. M is in g/mol"""
        if Tp is None:
            T0=self.T0
        if M is None:
            RR=self.R
        else:
            RR=Rstarkilo/M
        return RR*T0/self.g
    
    def pseudoz(self,pressure,H=None,p0 = None):
        """pseudo-altitude (log(p) coordinates"""
        if H is None:
            H = self.H()
        if p0 is None:
            p0=1.e5 #1bar
        return H*np.log(p0/pressure)
    
    def acosphi(self,lat):
        return self.a*np.cos(np.deg2rad(lat))
    
    def tanphia(self,lat):
        return np.tan(np.deg2rad(lat))/self.a
    
    def beta(self,lat=None):
        """Variation of coriolis parameter with latitude (df/dy)"""
        if lat is None:
            lat=0.
        return 2*self.omega*np.cos(p.deg2rad(lat))/self.a
    
    def angmom(self,u=None,lat=None):
        """Axial Angular momentum
            if u and lat are None, computes omega*a**2
        """
        if lat is None:
            lat=0.
        if u is None:
            u=0.
        acosphi=self.acosphi(lat)
        return acosphi*((self.omega*acosphi)+u)
    
    def wangmon(self,u=None,lat=None):
        """Axial Angular momentum due to wind ONLY
        """
        if lat is None:
            lat=0.
        if u is None:
            u=0.
            return u*self.acosphi(lat)
    
    def superrot(self,u=None,lat=None):
        """Super-rotation index"""
        aam=self.angmom(u=u,lat=lat)/self.angmom
        return amm-1.
    
    def exner(self,p,p0=1.e5):
        return (p/p0)**(self.R/self.cp)
    
    def tpot(self,temp,p,p0=1.e5):
        """returns potential temperature from temperature"""
        return temp/self.exner(p,p0=p0)
    def invtpot(self,tpot,p,po = 1.e5):
        """returns temperature, using potential temperature"""
        return tpot*self.exner(p,p0=p0)
    
class RunParameter():
    def __init__(self,file):
        conf=ConfigObj(file)['Run_parameters']
        self.is_omega=bool(strtobool(conf['is_omega']))
        self.includels=bool(strtobool(conf['includels']))
        self.use_spline=bool(strtobool(conf['use_spline']))
        self.temperature_field=conf['temperature_field']
        self.p_upper=float(conf['p_upper'])
        self.p_lower = float(conf['p_lower'])
        self.nlev = float(conf['nlev'])
        self.day_per_year = float(conf['day_per_year']) 
        self.nopole= bool(strtobool(conf['nopole'])) 
        self.charx = conf['charx']
            
    def show(self):
        for k,v, in list(vars(self).items()):
            print(k,v)
if __name__=="__main__":
    pr = Planet('w43b.par')
    pr.ini()
    print(pr.name)
    print(pr.a)
    print(pr.omega)
    rp = RunParameter('w43b.par')
    print(rp.charx)
    rp.show()