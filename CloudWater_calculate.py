# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 15:03:18 2019

@author: Administrator
fujian:23-29,115-121，0.25*0.25，25*25
era5
lat是从大到小的，和fnl的grib相反，lon正常
[27,104,20,112.5]

加上凝结-蒸发<0的情况
"""
from datetime import datetime
import numpy as np 
from netCDF4 import Dataset
import cftime
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
import sys
import multiprocessing
import os
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus'] = False

def calflux(u,v,q,clwc,isarea,isup,isdown,isleft,isright,zweight,area_grid,lat):
    [nt,nlev,ny,nx]=q.shape
    dt=3600
    qinflux=0 #kg
    qoutflux=0 #kg
    cinflux=0 #kg
    coutflux=0 #kg
    qsum_avg=0 #kg
    csum_avg=0 #kg
    qinflux_grid=np.zeros((ny,nx)) #kg
    qoutflux_grid=np.zeros((ny,nx)) #kg
    cinflux_grid=np.zeros((ny,nx)) #kg
    coutflux_grid=np.zeros((ny,nx)) #kg
    qsum_avg_grid=np.zeros((ny,nx)) #kg
    csum_avg_grid=np.zeros((ny,nx)) #kg
    uq=np.zeros((ny,nx)) #g/(s*cm)
    vq=np.zeros((ny,nx)) #g/(s*cm)
    uc=np.zeros((ny,nx)) #g/(s*cm)
    vc=np.zeros((ny,nx)) #g/(s*cm)
    
    #计算画图范围,周围一圈最后为0的行列,这里lat是从大到小
    islatsum=np.sum(isarea,axis=1)          
    islonsum=np.sum(isarea,axis=0)
    plat0=0;plon0=0;plat1=len(islatsum)-1;plon1=len(islonsum)-1;
    for i in range(len(islatsum)-1):
        if islatsum[i]==0 and islatsum[i+1]!=0:
            plat0=i
            break
    for i in range(len(islatsum)-1):
        if islatsum[len(islatsum)-1-i]==0 and islatsum[len(islatsum)-1-i-1]!=0:
            plat1=len(islatsum)-1-i
            break
    for i in range(len(islonsum)-1):
        if islonsum[i]==0 and islonsum[i+1]!=0:
            plon0=i
            break
    for i in range(len(islonsum)-1):
        if islonsum[len(islonsum)-1-i]==0 and islonsum[len(islonsum)-1-i-1]!=0:
            plon1=len(islonsum)-1-i
            break    
        
    for tt in range(nt):
        uarr=u[tt,:,:,:]
        varr=v[tt,:,:,:]
        uleft=uarr;uright=uarr;
        vup=varr;vdown=varr;
        carr=clwc[tt,:,:,:]
        qarr=q[tt,:,:,:]
        for ix in range(nx-1):
            uleft[:,:,ix+1]=(uarr[:,:,ix]+uarr[:,:,ix+1])/2
            uright[:,:,ix]=(uarr[:,:,ix]+uarr[:,:,ix+1])/2
        for iy in range(ny-1): #v的lat是从大到小，序号小的是上边界
            vup[:,iy+1,:]=(varr[:,iy,:]+varr[:,iy+1,:])/2
            vdown[:,iy,:]=(varr[:,iy,:]+varr[:,iy+1,:])/2
            
        qlevsum=np.zeros((ny,nx)) #g/cm2
        clevsum=np.zeros((ny,nx)) #g/cm2
        for ix in range(nx):
            for iy in range(ny): 
                if iy>=plat0 and iy<=plat1 and ix>=plon0 and ix<=plon1:
                    qlevsum[iy,ix]=np.sum(qarr[:,iy,ix]*zweight)/9.8*0.1
                    clevsum[iy,ix]=np.sum(carr[:,iy,ix]*zweight)/9.8*0.1
                    qsum_avg_grid[iy,ix]+=qlevsum[iy,ix]*area_grid[iy,ix]*10
                    csum_avg_grid[iy,ix]+=clevsum[iy,ix]*area_grid[iy,ix]*10
                    uq[iy,ix]+=np.sum(qarr[:,iy,ix]*uarr[:,iy,ix]*zweight)/9.8*10
                    vq[iy,ix]+=np.sum(qarr[:,iy,ix]*varr[:,iy,ix]*zweight)/9.8*10
                    uc[iy,ix]+=np.sum(carr[:,iy,ix]*uarr[:,iy,ix]*zweight)/9.8*10
                    vc[iy,ix]+=np.sum(carr[:,iy,ix]*varr[:,iy,ix]*zweight)/9.8*10
                    #格点边界输入输出
                    for iz in range(nlev):
                        #下边界,v>0输入，v<0输出
                        if (vdown[iz,iy,ix] > 0):
                            qinflux_grid[iy,ix]=qinflux_grid[iy,ix]+qarr[iz,iy,ix]*vdown[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                            cinflux_grid[iy,ix]=cinflux_grid[iy,ix]+carr[iz,iy,ix]*vdown[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                        else:
                            qoutflux_grid[iy,ix]=qoutflux_grid[iy,ix]-qarr[iz,iy,ix]*vdown[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                            coutflux_grid[iy,ix]=coutflux_grid[iy,ix]-carr[iz,iy,ix]*vdown[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                        #上边界,v>0输出，v<0输入
                        if (vup[iz,iy,ix] > 0):
                            qoutflux_grid[iy,ix]=qoutflux_grid[iy,ix]+qarr[iz,iy,ix]*vup[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                            coutflux_grid[iy,ix]=coutflux_grid[iy,ix]+carr[iz,iy,ix]*vup[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                        else:
                            qinflux_grid[iy,ix]=qinflux_grid[iy,ix]-qarr[iz,iy,ix]*vup[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                            cinflux_grid[iy,ix]=cinflux_grid[iy,ix]-carr[iz,iy,ix]*vup[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                        #左边界,u>0输入，u<0输出   
                        if (uleft[iz,iy,ix] > 0):
                            qinflux_grid[iy,ix]=qinflux_grid[iy,ix]+qarr[iz,iy,ix]*uleft[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25/9.8*dt
                            cinflux_grid[iy,ix]=cinflux_grid[iy,ix]+carr[iz,iy,ix]*uleft[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25/9.8*dt
                        else:
                            qoutflux_grid[iy,ix]=qoutflux_grid[iy,ix]-qarr[iz,iy,ix]*uleft[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25/9.8*dt
                            coutflux_grid[iy,ix]=coutflux_grid[iy,ix]-carr[iz,iy,ix]*uleft[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25/9.8*dt
                        #右边界,u>0输出，u<0输入    
                        if (uright[iz,iy,ix] > 0):
                            qoutflux_grid[iy,ix]=qoutflux_grid[iy,ix]+qarr[iz,iy,ix]*uright[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25/9.8*dt
                            coutflux_grid[iy,ix]=coutflux_grid[iy,ix]+carr[iz,iy,ix]*uright[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25/9.8*dt
                        else:
                            qinflux_grid[iy,ix]=qinflux_grid[iy,ix]-qarr[iz,iy,ix]*uright[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25/9.8*dt
                            cinflux_grid[iy,ix]=cinflux_grid[iy,ix]-carr[iz,iy,ix]*uright[iz,iy,ix]*zweight[iz]\
                            *6371000*2*np.pi/360*0.25/9.8*dt                 
                    
                        #边界输入输出分开叠加        
                        #下边界,v>0输入，v<0输出
                        if isdown[iy,ix]==1:
                            if (varr[iz,iy,ix] > 0):
                                qinflux=qinflux+qarr[iz,iy,ix]*varr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                                cinflux=cinflux+carr[iz,iy,ix]*varr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                            else:
                                qoutflux=qoutflux-qarr[iz,iy,ix]*varr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                                coutflux=coutflux-carr[iz,iy,ix]*varr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                        #上边界,v>0输出，v<0输入
                        if isup[iy,ix]==1:
                            if (varr[iz,iy,ix] > 0):
                                qoutflux=qoutflux+qarr[iz,iy,ix]*varr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                                coutflux=coutflux+carr[iz,iy,ix]*varr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                            else:
                                qinflux=qinflux-qarr[iz,iy,ix]*varr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                                cinflux=cinflux-carr[iz,iy,ix]*varr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25*np.cos(lat[iy,ix]/180*np.pi)/9.8*dt
                        #左边界,u>0输入，u<0输出
                        if isleft[iy,ix]==1:
                            if (uarr[iz,iy,ix] > 0):
                                qinflux=qinflux+qarr[iz,iy,ix]*uarr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25/9.8*dt
                                cinflux=cinflux+carr[iz,iy,ix]*uarr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25/9.8*dt
                            else:
                                qoutflux=qoutflux-qarr[iz,iy,ix]*uarr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25/9.8*dt
                                coutflux=coutflux-carr[iz,iy,ix]*uarr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25/9.8*dt
                        #右边界,u>0输出，u<0输入
                        if isright[iy,ix]==1:
                            if (uarr[iz,iy,ix] > 0):
                                qoutflux=qoutflux+qarr[iz,iy,ix]*uarr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25/9.8*dt
                                coutflux=coutflux+carr[iz,iy,ix]*uarr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25/9.8*dt
                            else:
                                qinflux=qinflux-qarr[iz,iy,ix]*uarr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25/9.8*dt
                                cinflux=cinflux-carr[iz,iy,ix]*uarr[iz,iy,ix]*zweight[iz]\
                                *6371000*2*np.pi/360*0.25/9.8*dt
                    
                    if isarea[iy,ix]==1:                           
                        qsum_avg+=qlevsum[iy,ix]*area_grid[iy,ix]*10
                        csum_avg+=clevsum[iy,ix]*area_grid[iy,ix]*10
                   
    return [qinflux,qoutflux,cinflux,coutflux,qsum_avg,csum_avg,\
            qinflux_grid,qoutflux_grid,cinflux_grid,coutflux_grid,qsum_avg_grid,csum_avg_grid,uq,vq,uc,vc]
    
    
if __name__ == '__main__':    
    multiprocessing.freeze_support()
    #参数文件
    argc = len(sys.argv)
    if argc!=2:
        print ('Need to input parafile')
        sys.exit(1)
    else:
        parafile=sys.argv[1]
    
    fpara=open(parafile,'r')
    para=fpara.readlines()
    fpara.close()
    tbstr=para[0].strip('\n').split('=')[1]
    testr=para[1].strip('\n').split('=')[1]
    ryears=para[2].strip('\n').split('=')[1]
    filedir=para[3].strip('\n').split('=')[1]
    outdir=para[4].strip('\n').split('=')[1]
    isboundaryfile=para[5].strip('\n').split('=')[1]
    pnum=int(para[6].strip('\n').split('=')[1])
    mapdir=para[7].strip('\n').split('=')[1]
    lonmin=float(para[8].strip('\n').split('=')[1])
    lonmax=float(para[9].strip('\n').split('=')[1])
    latmin=float(para[10].strip('\n').split('=')[1])
    latmax=float(para[11].strip('\n').split('=')[1])
    
    #一次时间过程里能跨年，但是不能超过12个月
    nyear=int(testr[0:4])-int(tbstr[0:4]) 
    if nyear==0:
        nmonth=int(testr[4:6])-int(tbstr[4:6])+1
    elif nyear==1:
        nmonth=13-int(tbstr[4:6])+int(testr[4:6])
        if nmonth>12:
            exit(1)
    else: 
        exit(1)
            
    mday=[31,28,31,30,31,30,31,31,30,31,30,31]
        
    zweight=np.zeros((24))
    zweight[0:13]=5000
    zweight[13]=3750
    zweight[14:23]=2500
    zweight[23]=1250
    
    totalqv=0
    totalqc=0
    totalrain=0
    totalyunshui=0
    timeqv=0
    timeqc=0
    percondenqv=0
    perrainqv=0
    perrainqc=0
    perrainall=0
    ny=int((latmax-latmin)/0.25+1)
    nx=int((lonmax-lonmin)/0.25+1)
    totalqv_grid=np.zeros((ny,nx))
    totalqc_grid=np.zeros((ny,nx))
    totalrain_grid=np.zeros((ny,nx))
    totalyunshui_grid=np.zeros((ny,nx))
    timeqv_grid=np.zeros((ny,nx))
    timeqc_grid=np.zeros((ny,nx))
    percondenqv_grid=np.zeros((ny,nx))
    perrainqv_grid=np.zeros((ny,nx))
    perrainqc_grid=np.zeros((ny,nx))
    perrainall_grid=np.zeros((ny,nx))
    uq1=np.zeros((ny,nx))
    vq1=np.zeros((ny,nx))
    uc1=np.zeros((ny,nx))
    vc1=np.zeros((ny,nx))
    
    ryears=int(ryears)
    for yy in range(ryears): #过去多年的时间过程平均
        tbstr=para[0].strip('\n').split('=')[1]
        testr=para[1].strip('\n').split('=')[1]
        tbegin=datetime(int(tbstr[0:4])-yy,int(tbstr[4:6]),1,0,0,0)
        tend=datetime(int(testr[0:4])-yy,int(testr[4:6]),mday[int(testr[4:6])-1],23,0,0)                
        delta=datetime(2017,1,1,1,0,0)-datetime(2017,1,1,0,0,0) #1hour
        ntotal=int((tend-tbegin)/delta)   
        tbstr=tbegin.strftime("%Y%m")
        testr=tend.strftime("%Y%m")
        
        mday_used=np.zeros((nmonth)) #每月天数，降水用
        mday_used[0]=mday[int(tbstr[4:6])-1]
        if int(tbstr[4:6])==2 and np.mod(int(tbstr[0:4]),4)==0:
            mday_used[0]=29
        #读数据文件
        nc_obj=Dataset(filedir+'\\era5_'+tbstr+'_quv.nc')
        nc_obj.set_auto_mask(False)
        lat=(nc_obj.variables['latitude'][:])  #这里是从大到小的
        lon=(nc_obj.variables['longitude'][:])
        lev=(nc_obj.variables['level'][:])
        q=(nc_obj.variables['q'][:])  #kg/kg
        nc_obj.close()
        nc_obj=Dataset(filedir+'\\era5_'+tbstr+'_qc.nc')
        nc_obj.set_auto_mask(False)
        clwc=(nc_obj.variables['clwc'][:])+(nc_obj.variables['ciwc'][:]) \
                +(nc_obj.variables['crwc'][:])+(nc_obj.variables['cswc'][:])  #kg/kg
        nc_obj.close()
        nc_obj=Dataset(filedir+'\\era5_'+tbstr+'_rain.nc')
        nc_obj.set_auto_mask(False)
        tpm=(nc_obj.variables['tp'][:])  #m
        evap=-(nc_obj.variables['e'][:]+nc_obj.variables['es'][:])  #地面蒸发，m
        lat0p1=(nc_obj.variables['latitude'][:])  #这里是从大到小的
        lon0p1=(nc_obj.variables['longitude'][:])
        nc_obj.close()
        
        flist=[]
        flist.append(filedir+'\\era5_'+tbstr)
        if nmonth>1:
            if nyear==0:        
                for i in range(nmonth-1):
                    tstr=tbstr[0:4]+'%02d'%(int(tbstr[4:6])+1+i)
                    mday_used[i+1]=mday[int(tstr[4:6])-1]
                    if int(tstr[4:6])==2 and np.mod(int(tstr[0:4]),4)==0:
                        mday_used[i+1]=29
                    flist.append(filedir+'\\era5_'+tstr)    
                    nc_obj=Dataset(filedir+'\\era5_'+tstr+'_rain.nc')
                    nc_obj.set_auto_mask(False)
                    temp=(nc_obj.variables['tp'][:])  #m
                    tpm=np.append(tpm, temp, axis=0)
                    temp=-(nc_obj.variables['e'][:]+nc_obj.variables['es'][:])  #m
                    evap=np.append(evap, temp, axis=0)
                    nc_obj.close()
            else: #跨年
                for i in range(12-int(tbstr[4:6])):
                    tstr=tbstr[0:4]+'%02d'%(int(tbstr[4:6])+1+i)
                    mday_used[i+1]=mday[int(tstr[4:6])-1]
                    if int(tstr[4:6])==2 and np.mod(int(tstr[0:4]),4)==0:
                        mday_used[i+1]=29
                    flist.append(filedir+'\\era5_'+tstr)
                    nc_obj=Dataset(filedir+'\\era5_'+tstr+'_rain.nc')
                    nc_obj.set_auto_mask(False)
                    temp=(nc_obj.variables['tp'][:])  #m
                    tpm=np.append(tpm, temp, axis=0)
                    temp=-(nc_obj.variables['e'][:]+nc_obj.variables['es'][:])  #m
                    evap=np.append(evap, temp, axis=0)
                    nc_obj.close()
                for i in range(int(testr[4:6])):
                    tstr=testr[0:4]+'%02d'%(1+i)
                    mday_used[i+13-int(tbstr[4:6])]=mday[int(tstr[4:6])-1]
                    if int(tstr[4:6])==2 and np.mod(int(tstr[0:4]),4)==0:
                        mday_used[i+1]=29
                    flist.append(filedir+'\\era5_'+tstr)
                    nc_obj=Dataset(filedir+'\\era5_'+tstr+'_rain.nc')
                    nc_obj.set_auto_mask(False)
                    temp=(nc_obj.variables['tp'][:])  #m
                    tpm=np.append(tpm, temp, axis=0)
                    temp=-(nc_obj.variables['e'][:]+nc_obj.variables['es'][:])  #m
                    evap=np.append(evap, temp, axis=0)
                    nc_obj.close()
            
        nx=len(lon)
        ny=len(lat)
        lon,lat=np.meshgrid(lon,lat)
        nlev=len(lev)
        nx0p1=len(lon0p1)
        ny0p1=len(lat0p1)
        lon0p1,lat0p1=np.meshgrid(lon0p1,lat0p1)
    
        
        #经度、纬度、是否区域内、上边界、下边界、左边界、右边界
        [dmp1,dmp2,isarea,isup,isdown,isleft,isright]=\
            np.loadtxt(isboundaryfile, delimiter='\t',unpack=True) 
        isarea=isarea.reshape(ny,nx)[::-1,:]
        isup=isup.reshape(ny,nx)[::-1,:]
        isdown=isdown.reshape(ny,nx)[::-1,:]
        isleft=isleft.reshape(ny,nx)[::-1,:]
        isright=isright.reshape(ny,nx)[::-1,:]
        
        #计算画图范围,周围一圈最后为0的行列,这里lat是从大到小        
        islatsum=np.sum(isarea,axis=1)          
        islonsum=np.sum(isarea,axis=0)
        plat0=0;plon0=0;plat1=len(islatsum)-1;plon1=len(islonsum)-1;
        for i in range(len(islatsum)-1):
            if islatsum[i]==0 and islatsum[i+1]!=0:
                plat0=i
                break
        for i in range(len(islatsum)-1):
            if islatsum[len(islatsum)-1-i]==0 and islatsum[len(islatsum)-1-i-1]!=0:
                plat1=len(islatsum)-1-i
                break
        for i in range(len(islonsum)-1):
            if islonsum[i]==0 and islonsum[i+1]!=0:
                plon0=i
                break
        for i in range(len(islonsum)-1):
            if islonsum[len(islonsum)-1-i]==0 and islonsum[len(islonsum)-1-i-1]!=0:
                plon1=len(islonsum)-1-i
                break    
        
        '''
        水汽和水凝物初值
        '''
        carr=clwc[0,:,:,:]
        qarr=q[0,:,:,:]
          
        qlevsum=np.zeros((ny,nx)) #g/cm2
        clevsum=np.zeros((ny,nx)) #g/cm2
        qsum_begin=0 #kg
        csum_begin=0 #kg  
        area=0 #m2
        qsum_begin_grid=np.zeros((ny,nx)) #kg
        csum_begin_grid=np.zeros((ny,nx)) #kg
        area_grid=np.zeros((ny,nx)) #m2
        
        for ix in range(nx):
            for iy in range(ny):
                qlevsum[iy,ix]=np.sum(qarr[:,iy,ix]*zweight)/9.8*0.1
                clevsum[iy,ix]=np.sum(carr[:,iy,ix]*zweight)/9.8*0.1 
                area_grid[iy,ix]=(6371000*2*np.pi/360*0.25)**2*np.cos(lat[iy,ix]/180*np.pi)
                if iy>=plat0 and iy<=plat1 and ix>=plon0 and ix<=plon1: 
                    qsum_begin_grid[iy,ix]=qlevsum[iy,ix]*area_grid[iy,ix]*10
                    csum_begin_grid[iy,ix]=clevsum[iy,ix]*area_grid[iy,ix]*10
                    if isarea[iy,ix]==1:                
                        area+=area_grid[iy,ix]                      
                        qsum_begin+=qlevsum[iy,ix]*area_grid[iy,ix]*10
                        csum_begin+=clevsum[iy,ix]*area_grid[iy,ix]*10
                                            
    
        '''
        地面降水和地面蒸发,月均值，0.1°
        '''    
        points=np.zeros((nx0p1*ny0p1,2))
        points[:,0]=lon0p1.reshape(nx0p1*ny0p1)
        points[:,1]=lat0p1.reshape(nx0p1*ny0p1)
        
        rain=0 #kg  
        rain_grid=np.zeros((ny,nx)) 
        evaporate=0 #kg  
        evaporate_grid=np.zeros((ny,nx)) 
        for tt in range(nmonth):
            mday1=mday_used[tt]
            tpnew=np.zeros((ny,nx))
            values=tpm[tt,:,:].reshape(nx0p1*ny0p1)
            tpnew = griddata(points, values, (lon,lat), method='nearest')
            enew=np.zeros((ny,nx))
            evalues=evap[tt,:,:].reshape(nx0p1*ny0p1)
            enew = griddata(points, evalues, (lon,lat), method='nearest')
            for ix in range(nx):
                for iy in range(ny):
                    rain_grid[iy,ix]+=tpnew[iy,ix]*area_grid[iy,ix]*1e3*mday1
                    evaporate_grid[iy,ix]+=enew[iy,ix]*area_grid[iy,ix]*1e3*mday1
                    if isarea[iy,ix]==1:
                        rain+=tpnew[iy,ix]*area_grid[iy,ix]*1e3*mday1
                        evaporate+=enew[iy,ix]*area_grid[iy,ix]*1e3*mday1
                        
                               
        '''
        水汽输入输出时间积分,1小时间隔
        '''             
        dt=3600
        qinflux=0 #kg
        qoutflux=0 #kg
        cinflux=0 #kg
        coutflux=0 #kg
        qsum_avg=0 #kg
        csum_avg=0 #kg
        qinflux_grid=np.zeros((ny,nx)) #kg
        qoutflux_grid=np.zeros((ny,nx)) #kg
        cinflux_grid=np.zeros((ny,nx)) #kg
        coutflux_grid=np.zeros((ny,nx)) #kg
        qsum_avg_grid=np.zeros((ny,nx)) #kg
        csum_avg_grid=np.zeros((ny,nx)) #kg
        uq=np.zeros((ny,nx)) #g/(s*cm)
        vq=np.zeros((ny,nx)) #g/(s*cm)
        uc=np.zeros((ny,nx)) #g/(s*cm)
        vc=np.zeros((ny,nx)) #g/(s*cm)
        
        #pnum = 4         
        pool = multiprocessing.Pool(processes=pnum)  # 开线程池
        result = []                    
        for ff in flist:
            nc_obj=Dataset(ff+'_quv.nc')
            nc_obj.set_auto_mask(False)
            q=(nc_obj.variables['q'][:])  #kg/kg
            u=(nc_obj.variables['u'][:])
            v=(nc_obj.variables['v'][:])
            nc_obj.close()
            nc_obj=Dataset(ff+'_qc.nc')
            nc_obj.set_auto_mask(False)
            clwc=(nc_obj.variables['clwc'][:])+(nc_obj.variables['ciwc'][:]) \
                    +(nc_obj.variables['crwc'][:])+(nc_obj.variables['cswc'][:])  #kg/kg
            nc_obj.close()        
            result.append(pool.apply_async(func=calflux, args=(u,v,q,clwc,\
                isarea,isup,isdown,isleft,isright,zweight,area_grid,lat,)))  
        pool.close()
        pool.join()
        # result里面0到n-1顺序保留n个进程的结果，每个里面有value对象保存，分别是[qinflux,qoutflux...]
        for res in result:
            ans = res.get()
            qinflux+=ans[0]
            qoutflux+=ans[1]
            cinflux+=ans[2]
            coutflux+=ans[3]
            qsum_avg+=ans[4]
            csum_avg+=ans[5]
            qinflux_grid+=ans[6]
            qoutflux_grid+=ans[7]
            cinflux_grid+=ans[8]
            coutflux_grid+=ans[9]
            qsum_avg_grid+=ans[10]
            csum_avg_grid+=ans[11]
            uq+=ans[12]
            vq+=ans[13]
            uc+=ans[14]
            vc+=ans[15]
        qsum_avg=qsum_avg/ntotal
        csum_avg=csum_avg/ntotal
        qsum_avg_grid=qsum_avg_grid/ntotal
        csum_avg_grid=csum_avg_grid/ntotal
        uq=uq/ntotal
        vq=vq/ntotal
        uc=uc/ntotal
        vc=vc/ntotal
        
        '''
        水汽和水凝物终值
        '''
        carr=clwc[-1,:,:,:]
        qarr=q[-1,:,:,:]
          
        qlevsum=np.zeros((ny,nx)) #g/cm2
        clevsum=np.zeros((ny,nx)) #g/cm2
        qsum_end=0 #kg
        csum_end=0 #kg 
        qsum_end_grid=np.zeros((ny,nx)) #kg
        csum_end_grid=np.zeros((ny,nx)) #kg  
        for ix in range(nx):
            for iy in range(ny):
                qlevsum[iy,ix]=np.sum(qarr[:,iy,ix]*zweight)/9.8*0.1
                clevsum[iy,ix]=np.sum(carr[:,iy,ix]*zweight)/9.8*0.1
                if iy>=plat0 and iy<=plat1 and ix>=plon0 and ix<=plon1: 
                    qsum_end_grid[iy,ix]=qlevsum[iy,ix]*area_grid[iy,ix]*10
                    csum_end_grid[iy,ix]=clevsum[iy,ix]*area_grid[iy,ix]*10
                    if isarea[iy,ix]==1:                       
                        qsum_end=qsum_end+qlevsum[iy,ix]*area_grid[iy,ix]*10
                        csum_end=csum_end+clevsum[iy,ix]*area_grid[iy,ix]*10
                                    
        '''
        计算最终变量(kg)
        '''
        days=np.sum(mday_used)
        conden=csum_end+coutflux+rain-csum_begin-cinflux
        totalrain+=rain
        if rain==0:
            norainyear+=1
        else:
            timeqv+=qsum_avg/rain*days
            timeqc+=csum_avg/rain*days*24
        if conden>=0: #凝结-蒸发>=0，用凝结
            totalqv+=(qsum_end+qoutflux+conden)
            totalqc+=(csum_begin+cinflux+conden)            
            totalyunshui+=(csum_begin+cinflux+conden-rain)            
            percondenqv+=100*(csum_begin+cinflux+conden)/(qsum_end+qoutflux+conden)
            perrainqv+=100*rain/(qsum_end+qoutflux+conden)
            perrainqc+=100*rain/(csum_begin+cinflux+conden)
            perrainall+=100*rain/(qsum_end+qoutflux+conden+csum_begin+cinflux+conden)
        else: #凝结-蒸发<0，用蒸发
            totalqv+=(qsum_begin+qinflux+conden+evaporate)
            totalqc+=(csum_end+coutflux+rain+conden)
            totalyunshui+=(csum_end+coutflux+conden)
            percondenqv+=100*(csum_end+coutflux+rain+conden)/(qsum_begin+qinflux+conden+evaporate)
            perrainqv+=100*rain/(qsum_begin+qinflux+conden+evaporate)
            perrainqc+=100*rain/(csum_end+coutflux+rain+conden)
            perrainall+=100*rain/(qsum_begin+qinflux+conden+evaporate+csum_end+coutflux+rain+conden)
        
        
        conden_grid=np.zeros((ny,nx))
        for ix in range(nx):
            for iy in range(ny):
                if iy>=plat0 and iy<=plat1 and ix>=plon0 and ix<=plon1: 
                    uq1[iy,ix]+=uq[iy,ix]
                    vq1[iy,ix]+=vq[iy,ix]
                    uc1[iy,ix]+=uc[iy,ix]
                    vc1[iy,ix]+=vc[iy,ix]
                    totalrain_grid[iy,ix]+=rain_grid[iy,ix]
                    if rain_grid[iy,ix]==0:
                        norainyear_grid[iy,ix]+=1
                    else:
                        timeqv_grid[iy,ix]+=qsum_avg_grid[iy,ix]/rain_grid[iy,ix]*days
                        timeqc_grid[iy,ix]+=csum_avg_grid[iy,ix]/rain_grid[iy,ix]*days*24
                    conden_grid[iy,ix]=(csum_end_grid[iy,ix]+coutflux_grid[iy,ix]+rain_grid[iy,ix]\
                               -csum_begin_grid[iy,ix]-cinflux_grid[iy,ix])
                    if conden_grid[iy,ix]>=0: #凝结-蒸发>=0，用凝结 
                        totalqv_grid[iy,ix]+=(qsum_end_grid[iy,ix]+qoutflux_grid[iy,ix]+conden_grid[iy,ix])
                        totalqc_grid[iy,ix]+=(csum_begin_grid[iy,ix]+cinflux_grid[iy,ix]+conden_grid[iy,ix])
                        totalyunshui_grid[iy,ix]+=(csum_begin_grid[iy,ix]+cinflux_grid[iy,ix]\
                                    +conden_grid[iy,ix]-rain_grid[iy,ix])               
                        percondenqv_grid[iy,ix]+=100*(csum_begin_grid[iy,ix]+cinflux_grid[iy,ix]+conden_grid[iy,ix])\
                                    /(qsum_end_grid[iy,ix]+qoutflux_grid[iy,ix]+conden_grid[iy,ix])
                        perrainqv_grid[iy,ix]+=100*rain_grid[iy,ix]\
                                    /(qsum_end_grid[iy,ix]+qoutflux_grid[iy,ix]+conden_grid[iy,ix])
                        perrainqc_grid[iy,ix]+=100*rain_grid[iy,ix]\
                                    /(csum_begin_grid[iy,ix]+cinflux_grid[iy,ix]+conden_grid[iy,ix])
                        perrainall_grid[iy,ix]+=100*rain_grid[iy,ix]/(qsum_end_grid[iy,ix]+qoutflux_grid[iy,ix]\
                                    +conden_grid[iy,ix]+csum_begin_grid[iy,ix]+cinflux_grid[iy,ix]+conden_grid[iy,ix])
                    else: #凝结-蒸发<0，用蒸发
                        conden_grid[iy,ix]=-conden_grid[iy,ix]
                        totalqv_grid[iy,ix]+=(qsum_begin_grid[iy,ix]+qinflux_grid[iy,ix]+conden_grid[iy,ix]+evaporate_grid[iy,ix])
                        totalqc_grid[iy,ix]+=(csum_end_grid[iy,ix]+coutflux_grid[iy,ix]+rain_grid[iy,ix]+conden_grid[iy,ix])
                        totalyunshui_grid[iy,ix]+=(csum_end_grid[iy,ix]+coutflux_grid[iy,ix]+conden_grid[iy,ix])                   
                        percondenqv_grid[iy,ix]+=100*(csum_end_grid[iy,ix]+coutflux_grid[iy,ix]+rain_grid[iy,ix]+conden_grid[iy,ix])\
                                    /(qsum_end_grid[iy,ix]+qoutflux_grid[iy,ix]+conden_grid[iy,ix])
                        perrainqv_grid[iy,ix]+=100*rain_grid[iy,ix]\
                                    /(qsum_begin_grid[iy,ix]+qinflux_grid[iy,ix]+conden_grid[iy,ix]+evaporate_grid[iy,ix])
                        perrainqc_grid[iy,ix]+=100*rain_grid[iy,ix]\
                                    /(csum_end_grid[iy,ix]+coutflux_grid[iy,ix]+rain_grid[iy,ix]+conden_grid[iy,ix])
                        perrainall_grid[iy,ix]+=100*rain_grid[iy,ix]\
                            /(qsum_begin_grid[iy,ix]+qinflux_grid[iy,ix]+conden_grid[iy,ix]+evaporate_grid[iy,ix]\
                                    +csum_end_grid[iy,ix]+coutflux_grid[iy,ix]+rain_grid[iy,ix]+conden_grid[iy,ix])
                               
    tbstr=para[0].strip('\n').split('=')[1]
    testr=para[1].strip('\n').split('=')[1] 
    
    outfilename=outdir+'/云水资源区域结果_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.txt'
    fout=open(outfilename,'w')
    fout.write('水汽总量 = %6.2f kg/m2\n'%(totalqv/ryears))
    fout.write('水凝物总量 = %6.2f kg/m2\n'%(totalqc/ryears))
    fout.write('降水总量 = %6.2f kg/m2\n'%(totalrain/ryears))   
    fout.write('云水总量 = %6.2f kg/m2\n'%(totalyunshui/ryears))
    fout.write('水汽更新周期 = %6.2f 天\n'%(timeqv/(ryears-norainyear)))
    fout.write('水凝物更新周期 = %6.2f 小时\n'%(timeqc/(ryears-norainyear)))
    fout.write('水汽凝结效率 = %6.2f %%\n'%(percondenqv/ryears))
    fout.write('水汽降水效率 = %6.2f %%\n'%(perrainqv/ryears))
    fout.write('水凝物降水效率 = %6.2f %%\n'%(perrainqc/ryears))
    fout.write('总水物质降水效率 = %6.2f %%\n'%(perrainall/ryears))
    fout.close()  
    
    totalqv_grid[totalqv_grid==0]=np.nan  
    totalqc_grid[totalqc_grid==0]=np.nan 
    totalrain_grid[totalrain_grid==0]=np.nan   
    totalyunshui_grid[totalyunshui_grid==0]=np.nan 
    timeqv_grid[timeqv_grid==0]=np.nan  
    timeqc_grid[timeqc_grid==0]=np.nan 
    percondenqv_grid[percondenqv_grid==0]=np.nan  
    perrainqv_grid[perrainqv_grid==0]=np.nan  
    perrainqc_grid[perrainqc_grid==0]=np.nan  
    perrainall_grid[perrainall_grid==0]=np.nan  
    uq1[uq1==0]=np.nan  
    vq1[vq1==0]=np.nan  
    uc1[uc1==0]=np.nan  
    vc1[vc1==0]=np.nan
    
    latmin=lat[plat1,plon1]
    latmax=lat[plat0,plon0]
    lonmin=lon[plat0,plon0]
    lonmax=lon[plat1,plon1]
    lon_num = np.arange(lonmin,lonmax+0.25,0.25)    
    lat_num =  np.arange(latmin,latmax+0.25,0.25)
    micaps_nx=len(lon_num)
    micaps_ny=len(lat_num)
    if ny>50 or nx>50:
        lon_num = np.arange(lonmin,lonmax+1,1)    
        lat_num =  np.arange(latmin,latmax+1,1)
    elif ny>20 or nx>20:
        lon_num = np.arange(lonmin,lonmax+0.5,0.5)    
        lat_num =  np.arange(latmin,latmax+0.5,0.5)
    if lon_num[-1]>lonmax:
        lon_num=lon_num[:-1]
    if lat_num[-1]>latmax:
        lat_num=lat_num[:-1]
    lon_label = lon_num
    lat_label = lat_num    
    
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    cs1=m.contourf(lon,lat,totalqv_grid/ryears,20,extend='both',latlon=True,cmap='jet')
    m.colorbar(cs1)
    plt.title('水汽总量 (kg/m2)')
    plt.savefig(outdir+'/水汽总量_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()
    
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    cs1=m.contourf(lon,lat,totalqc_grid/ryears,20,extend='both',latlon=True,cmap='jet')
    m.colorbar(cs1)  
    plt.title('水凝物总量 (kg/m2)')
    plt.savefig(outdir+'/水凝物总量_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()        
    
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    cs1=m.contourf(lon,lat,totalrain_grid/ryears,20,extend='both',latlon=True,cmap='jet')
    m.colorbar(cs1)     
    plt.title('降水总量 (kg/m2)')
    plt.savefig(outdir+'/降水总量_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()
    
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    cs1=m.contourf(lon,lat,totalyunshui_grid/ryears,20,extend='both',latlon=True,cmap='jet')
    m.colorbar(cs1)  
    plt.title('云水总量 (kg/m2)')
    plt.savefig(outdir+'/云水总量_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()        
                  
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    cs1=m.contourf(lon,lat,timeqv_grid/(ryears-norainyear_grid),20,extend='both',latlon=True,cmap='jet')
    m.colorbar(cs1)  
    plt.title('水汽更新周期 (天)')
    plt.savefig(outdir+'/水汽更新周期_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()  
    
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    cs1=m.contourf(lon,lat,timeqc_grid/(ryears-norainyear_grid),20,extend='both',latlon=True,cmap='jet')
    m.colorbar(cs1)  
    plt.title('水凝物更新周期 (小时)')
    plt.savefig(outdir+'/水凝物更新周期_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()  
    
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    cs1=m.contourf(lon,lat,percondenqv_grid/ryears,20,extend='both',latlon=True,cmap='jet')
    m.colorbar(cs1)  
    plt.title('水汽凝结效率 (%)')
    plt.savefig(outdir+'/水汽凝结效率_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()  
    
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    cs1=m.contourf(lon,lat,perrainqv_grid/ryears,20,extend='both',latlon=True,cmap='jet')
    m.colorbar(cs1) 
    plt.title('水汽降水效率 (%)')
    plt.savefig(outdir+'/水汽降水效率_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()  
    
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    cs1=m.contourf(lon,lat,perrainqc_grid/ryears,20,extend='both',latlon=True,cmap='jet')
    m.colorbar(cs1) 
    plt.title('水凝物降水效率 (%)')
    plt.savefig(outdir+'/水凝物降水效率_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()  
    
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    cs1=m.contourf(lon,lat,perrainall_grid/ryears,20,extend='both',latlon=True,cmap='jet')
    m.colorbar(cs1) 
    plt.title('总水物质降水效率 (%)')
    plt.savefig(outdir+'/总水物质降水效率_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()  
    
    if not os.path.exists(outdir+'/CWR'):
        os.mkdir(outdir+'/CWR')
    dirnames=['GQv','GQh','GR','GCWR','Tv','Th','P1','Ev','Eh','Em','UQv','UQh']
    for dnames in dirnames:
        if not os.path.exists(outdir+'/CWR/'+dnames):
            os.mkdir(outdir+'/CWR/'+dnames)

    fout=outdir+'/CWR/GQv/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 4 水汽总量(kg/m2)\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    cmin=np.nanmin(totalqv_grid/ryears)
    cmax=np.nanmax(totalqv_grid/ryears)
    clev=(cmax-cmin)/20
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i %10.6f %10.6f %10.6f 1.000000 0.000000\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny,clev,cmin,cmax))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(totalqv_grid[iy,ix]/ryears))
            f.write('\n\n')
    f.close()
    
    fout=outdir+'/CWR/GQh/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 4 水凝物总量(kg/m2)\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    cmin=np.nanmin(totalqc_grid/ryears)
    cmax=np.nanmax(totalqc_grid/ryears)
    clev=(cmax-cmin)/20
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i %10.6f %10.6f %10.6f 1.000000 0.000000\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny,clev,cmin,cmax))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(totalqc_grid[iy,ix]/ryears))
            f.write('\n\n')
    f.close()
    
    fout=outdir+'/CWR/GR/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 4 降水总量(kg/m2)\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    cmin=np.nanmin(totalrain_grid/ryears)
    cmax=np.nanmax(totalrain_grid/ryears)
    clev=(cmax-cmin)/20
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i %10.6f %10.6f %10.6f 1.000000 0.000000\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny,clev,cmin,cmax))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(totalrain_grid[iy,ix]/ryears))
            f.write('\n\n')
    f.close()
    
    fout=outdir+'/CWR/GCWR/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 4 云水总量(kg/m2)\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    cmin=np.nanmin(totalyunshui_grid/ryears)
    cmax=np.nanmax(totalyunshui_grid/ryears)
    clev=(cmax-cmin)/20
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i %10.6f %10.6f %10.6f 1.000000 0.000000\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny,clev,cmin,cmax))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(totalyunshui_grid[iy,ix]/ryears))
            f.write('\n\n')
    f.close()
    
    fout=outdir+'/CWR/Tv/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 4 水汽更新周期(天)\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    cmin=np.nanmin(timeqv_grid/(ryears-norainyear_grid))
    cmax=np.nanmax(timeqv_grid/(ryears-norainyear_grid))
    clev=(cmax-cmin)/20
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i %10.6f %10.6f %10.6f 1.000000 0.000000\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny,clev,cmin,cmax))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(timeqv_grid[iy,ix]/(ryears-norainyear_grid[iy,ix])))
            f.write('\n\n')
    f.close()
    
    
    fout=outdir+'/CWR/Th/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 4 水凝物更新周期(小时)\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    cmin=np.nanmin(timeqc_grid/(ryears-norainyear_grid))
    cmax=np.nanmax(timeqc_grid/(ryears-norainyear_grid))
    clev=(cmax-cmin)/20
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i %10.6f %10.6f %10.6f 1.000000 0.000000\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny,clev,cmin,cmax))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(timeqc_grid[iy,ix]/(ryears-norainyear_grid[iy,ix])))
            f.write('\n\n')
    f.close()
    
    fout=outdir+'/CWR/P1/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 4 水汽凝结效率(%)\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    cmin=np.nanmin(percondenqv_grid/ryears)
    cmax=np.nanmax(percondenqv_grid/ryears)
    clev=(cmax-cmin)/20
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i %10.6f %10.6f %10.6f 1.000000 0.000000\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny,clev,cmin,cmax))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(percondenqv_grid[iy,ix]/ryears))
            f.write('\n\n')
    f.close()
    
    fout=outdir+'/CWR/Ev/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 4 水汽降水效率(%)\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    cmin=np.nanmin(perrainqv_grid/ryears)
    cmax=np.nanmax(perrainqv_grid/ryears)
    clev=(cmax-cmin)/20
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i %10.6f %10.6f %10.6f 1.000000 0.000000\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny,clev,cmin,cmax))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(perrainqv_grid[iy,ix]/ryears))
            f.write('\n\n')
    f.close()
    
    fout=outdir+'/CWR/Eh/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 4 水凝物降水效率(%)\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    cmin=np.nanmin(perrainqc_grid/ryears)
    cmax=np.nanmax(perrainqc_grid/ryears)
    clev=(cmax-cmin)/20
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i %10.6f %10.6f %10.6f 1.000000 0.000000\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny,clev,cmin,cmax))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(perrainqc_grid[iy,ix]/ryears))
            f.write('\n\n')
    f.close()
    
    fout=outdir+'/CWR/Em/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 4 总水物质降水效率(%)\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    cmin=np.nanmin(perrainall_grid/ryears)
    cmax=np.nanmax(perrainall_grid/ryears)
    clev=(cmax-cmin)/20
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i %10.6f %10.6f %10.6f 1.000000 0.000000\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny,clev,cmin,cmax))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(perrainall_grid[iy,ix]/ryears))
            f.write('\n\n')
    f.close()
    
    
    ###水汽通量
    fout=outdir+'/CWR/UQv/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 11 水汽通量(g/(s*cm2))\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(uq[iy,ix]/ryears))
            f.write('\n\n')
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(vq[iy,ix]/ryears))
            f.write('\n\n')
    f.close()
    
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    wsp=np.sqrt((uq/ryears)*(uq/ryears)+(vq/ryears)*(vq/ryears))    
    wsp=int(np.mean(wsp))
    pu=uq/ryears
    pv=vq/ryears
    plon=lon
    plat=lat
    if ny>50 or nx>50:
        pu=pu[::2,::2]
        pv=pv[::2,::2]
        plon=plon[::2,::2]
        plat=plat[::2,::2]        
    m.quiver(plon,plat,pu,pv,\
                width=0.003,scale=wsp*3,scale_units='inches',latlon=True)
    m.quiver(lonmin+(lonmax-lonmin)*0.1,\
             latmin+(latmax-latmin)*0.1,\
             wsp,0,width=0.003,scale=wsp*3,scale_units='inches',color='r')
    plt.text(lonmin+(lonmax-lonmin)*0.1,\
             latmin+(latmax-latmin)*0.07,\
             str(wsp)+' g/(s*cm2)',color='r')
    plt.title('水汽通量 (g/(s*cm2))')
    plt.savefig(outdir+'/水汽通量_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()     
    
    fout=outdir+'/CWR/UQh/'+tbstr[2:6]+'0000.000'
    f=open(fout,'w')
    f.write('diamond 11 水凝物通量(g/(s*cm2))\n\n')
    f.write('%s %s 01 0 0 0 \n\n'%(tbstr[:4],tbstr[4:6]))
    f.write('0.250000 -0.250000 %10.6f %10.6f %10.6f %10.6f %i %i\n\n'\
            %(lonmin,lonmax,latmax,latmin,micaps_nx,micaps_ny))
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(uc[iy,ix]/ryears))
            f.write('\n\n')
    for iy in range(ny):
        if iy>=plat0 and iy<=plat1:
            for ix in range(nx):
                if ix>=plon0 and ix<=plon1:
                    f.write('%6.2f '%(vc[iy,ix]/ryears))
            f.write('\n\n')
    f.close()
    
    fig = plt.figure(figsize=(10,10),dpi=200)
    m = Basemap(projection = 'cyl', llcrnrlat=latmin,urcrnrlat=latmax,llcrnrlon=lonmin,urcrnrlon=lonmax)
    m.readshapefile(mapdir+'/gadm36_CHN_1','states',drawbounds=True,linewidth=2)
    m.readshapefile(mapdir+'/gadm36_CHN_2','states',drawbounds=True,linewidth=1)    
    plt.yticks(lat_num,lat_label)
    plt.xticks(lon_num,lon_label)
    wsp=np.sqrt((uc/ryears)*(uc/ryears)+(vc/ryears)*(vc/ryears)) 
    wsp=int(np.mean(wsp))
    pu=uc/ryears
    pv=vc/ryears
    plon=lon
    plat=lat
    if ny>50 or nx>50:
        pu=pu[::2,::2]
        pv=pv[::2,::2]
        plon=plon[::2,::2]
        plat=plat[::2,::2]        
    m.quiver(plon,plat,pu,pv,\
                width=0.003,scale=wsp*3,scale_units='inches',latlon=True)
    m.quiver(lonmin+(lonmax-lonmin)*0.1,\
             latmin+(latmax-latmin)*0.1,\
             wsp,0,width=0.003,scale=wsp*3,scale_units='inches',color='r')
    plt.text(lonmin+(lonmax-lonmin)*0.1,\
             latmin+(latmax-latmin)*0.07,\
             str(wsp)+' g/(s*cm2)',color='r')
    plt.title('水凝物通量 (g/(s*cm2))')
    plt.savefig(outdir+'/水凝物通量_month'+tbstr[4:6]+'-'+testr[4:6]+'_year'+tbegin.strftime("%Y")+'-'+tbstr[:4]+'.png')
    plt.close()
    
    