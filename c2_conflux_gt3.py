# -*- coding: utf-8 -*-
"""
Created/Corrected:    25/09/09 09:09
Project:    GCP/GMB 2025 project
@author:    Dmitry Belikov
"""
import numpy as np
import xarray as xr
from datetime import datetime
from typing import Optional, Sequence, Tuple
import _set_case


class CnvFluxGt3(_set_case.SetCase):
    def __init__(self):
        super().__init__()
        # GT3 parameters
        self.in_imax = 360
        self.in_jmax = 180
        self.model_imax = 128
        self.model_jmax = 64
        self.spcnm = "FLUX"
        self.hunit = "UR4"
        self.kmax = 1

        # placeholders
        self.da = None
        self.times = None
        self.mlons = None
        self.mlats = None
        self.alon = None
        self.alat = None
        self.blat = None
        self.dlat = None
        self.gar_t42 = None
        self.gar_inp = None
        self.iloc = None
        self.jloc = None

        # Example usage
        self.convert("E1_gOng", "fch4")

    # -------------------------
    # Gaussian grid
    def setgrid_actm(self):
        def gauss(jmax: int, deltp: float = 1e-6, itrmax: int = 50) -> Tuple[np.ndarray, np.ndarray]:
            ctheta = np.zeros(jmax, float)
            gw = np.zeros(jmax, float)
            for jj in range(1, jmax//2+1):
                x0 = np.cos((jj-0.5)/jmax*np.pi)
                for _ in range(itrmax):
                    qpn = np.zeros(jmax+2, float)
                    qpn[0]=1.0
                    qpn[1]=np.sqrt(3.0)*x0
                    for n in range(2, jmax+2):
                        eps = np.sqrt(n**2/(4*n**2-1))
                        epsm = np.sqrt((n-1)**2/(4*(n-1)**2-1))
                        qpn[n] = (qpn[n-1]*x0 - qpn[n-2]*epsm)/eps
                    eps = np.sqrt(jmax**2/(4*jmax**2-1))
                    epsp = np.sqrt((jmax+1)**2/(4*(jmax+1)**2-1))
                    qdpn = qpn[jmax+1]*jmax*epsp - qpn[jmax-1]*(jmax+1)*eps
                    deltx = qpn[jmax]/qdpn*(1.0-x0**2)
                    x0 += deltx
                    if abs(deltx)<deltp: break
                ctheta[jj-1] = x0
                gw[jj-1] = (2*jmax-1)*(1-x0**2)/(jmax*qpn[jmax-1])**2
            for jj in range(jmax//2):
                ctheta[jmax-jj-1] = -ctheta[jj]
                gw[jmax-jj-1] = gw[jj]
            return ctheta, gw

        self.alon = 360.0*np.arange(self.model_imax)/self.model_imax
        ctheta, gw = gauss(self.model_jmax)
        self.alat = np.degrees(np.arcsin(ctheta))
        self.dlat = gw
        self.blat = np.zeros(self.model_jmax+1)
        self.blat[0] = 90.0
        self.blat[-1] = -90.0
        for j in range(1, self.model_jmax):
            self.blat[j] = 0.5*(self.alat[j]+self.alat[j-1])

    # -------------------------
    # Load NetCDF
    def load_netcdf(self, infile: str, varname: str):
        ds = xr.open_dataset(infile)
        da = ds[varname]
        self.mlons = np.mod(da.coords.get("lon", np.arange(0,360,360/self.in_imax)),360.0)
        self.mlats = da.coords.get("lat", np.linspace(-89.5,89.5,self.in_jmax))
        self.times = da.coords.get("time", np.array([np.datetime64("2000-01-01")]))
        self.da = da

    # -------------------------
    # Mapping indices
    def set_inout(self):
        if self.alon is None: self.setgrid_actm()
        mlon = np.where(self.mlons>0,self.mlons,self.mlons+360.0)
        pi = np.pi
        self.gar_t42 = 4.0*pi*self.R**2*(1.0/self.model_imax)*self.dlat
        self.gar_inp = np.zeros(self.in_jmax)
        dmdeg = 360/self.in_imax
        for j in range(1,self.in_jmax+1):
            angle = (-90.5 + j*dmdeg)*pi/180.0
            self.gar_inp[j-1] = np.cos(angle)*(dmdeg*self.R*np.pi/180.0)**2

        self.iloc = np.ones(self.in_imax,int)
        dlon = self.alon[1]-self.alon[0]
        for i in range(self.in_imax):
            amin = dlon*0.5
            chosen = 1
            for ii in range(self.model_imax):
                sub = abs(self.alon[ii]-mlon[i])
                if sub<amin:
                    amin=sub
                    chosen=ii+1
            self.iloc[i]=chosen

        self.jloc = np.zeros(self.in_jmax,int)
        for j in range(self.in_jmax):
            for jj in range(self.model_jmax):
                if self.mlats[j]-self.blat[jj]<=0 and self.mlats[j]-self.blat[jj+1]>=0:
                    self.jloc[j]=jj+1
                    break

    # -------------------------
    # GT3 header & record
    def gt3write_record(self, f, oval: np.ndarray, itime: int, jdate: Sequence[int]):
        head_fields = [' '*16]*64
        def put(idx,val,right=False):
            s=str(val)[:16]
            s = s.rjust(16) if right else s.ljust(16)
            head_fields[idx]=s
        put(0,9010,right=True)
        put(1,"GCP-Global")
        put(2,"QFLUX07")
        put(13,self.spcnm)
        put(15,self.hunit)
        put(24,itime,right=True)
        put(25,"MONTH")
        date_str = f"{jdate[0]:04d}{jdate[1]:02d}{jdate[2]:02d} {jdate[3]:02d}{jdate[4]:02d}{jdate[5]:02d}"
        put(26,date_str)
        put(28,f"GLON{self.model_imax}")
        put(31,f"GGLA{self.model_jmax}")
        put(34,"SFC1" if self.kmax==1 else f"CSIG{self.kmax:02d}")
        put(63,self.model_imax*self.model_jmax*self.kmax,right=True)
        for hf in head_fields: f.write(hf.encode('ascii'))

        arr_be = np.asfortranarray(oval).astype('>f4')
        f.write(arr_be.tobytes())

    # -------------------------
    # Top-level conversion
    def convert(self,infile:str,varname:str):
        self.load_netcdf(self.out_dir+infile+'.nc',varname)
        self.setgrid_actm()
        self.set_inout()
        outname=self.out_dir+infile+'.gt3'
        with open(outname,"wb") as f:
            arr = self.da.values
            if arr.ndim==2: arr=arr[np.newaxis,:,:]
            nt = arr.shape[0]
            iloc_py = self.iloc-1
            jloc_py = self.jloc-1
            for tindex in range(nt):
                budg = arr[tindex,:,:]
                gdz = np.zeros((self.model_imax,self.model_jmax))
                ndz = np.zeros((self.model_imax,self.model_jmax),int)
                for i_in in range(self.in_imax):
                    for j_in in range(self.in_jmax):
                        i_mod = iloc_py[i_in]
                        j_mod = jloc_py[j_in]
                        val = budg[j_in,i_in] if budg.shape==(self.in_jmax,self.in_imax) else budg[i_in,j_in]
                        gdz[i_mod,j_mod]+=float(val)
                        ndz[i_mod,j_mod]+=1
                mask=ndz>0
                gdz[mask]/=ndz[mask]
                rdz=np.where(gdz>=0,gdz,0.0).astype(np.float32)
                try:
                    tval=self.times[tindex]
                    py_dt=np.datetime64(tval).astype("datetime64[s]").astype(datetime)
                    jdate=[py_dt.year,py_dt.month,py_dt.day,py_dt.hour,py_dt.minute,py_dt.second]
                except:
                    jdate=[2000+tindex,1,1,0,0,0]
                oval=np.zeros((self.model_imax,self.model_jmax,self.kmax),np.float32)
                oval[:,:,0]=rdz
                self.gt3write_record(f,oval,tindex,jdate)
        print("GT3 conversion complete.")



