# -*- coding: utf-8 -*-
#"area": [30.5,113.5,24,118.5]
#year
#dd

import sys
import cdsapi
c = cdsapi.Client()

argc = len(sys.argv)  
if argc<9:
    print ('Need to input year mbegin mend latmin lonmin latmax lonmax outdir')
    sys.exit(1) 
    
year=sys.argv[1] 
mbegin=int(sys.argv[2])
mend=int(sys.argv[3])
latmin=float(sys.argv[4])
lonmin=float(sys.argv[5])
latmax=float(sys.argv[6])
lonmax=float(sys.argv[7])
outdir=sys.argv[8]
  
for dd in range(mbegin,mend+1):
    
    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type':'reanalysis',
            'format':'netcdf',
            'variable':[
                'specific_cloud_ice_water_content','specific_cloud_liquid_water_content',
                'specific_rain_water_content','specific_snow_water_content'
            ],
            'pressure_level':[
                '100','150',
                '200','250','300',
                '350','400','450',
                '500','550','600',
                '650','700','750',
                '775','800','825',
                '850','875','900',
                '925','950','975',
                '1000'
            ],
            'year':year,
            'month':'%02d'%(dd),
            'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31'
            ],
            'time':[
                '00:00','01:00','02:00',
                '03:00','04:00','05:00',
                '06:00','07:00','08:00',
                '09:00','10:00','11:00',
                '12:00','13:00','14:00',
                '15:00','16:00','17:00',
                '18:00','19:00','20:00',
                '21:00','22:00','23:00'
            ],
            "area": [latmax,lonmin,latmin,lonmax],
            
        },
        outdir+'\\era5_'+year+'%02d'%(dd)+'_qc.nc')
    
    #######
    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type':'reanalysis',
            'format':'netcdf',
            'variable':[
                'specific_humidity',
                'u_component_of_wind','v_component_of_wind'
            ],
            'pressure_level':[
                '100','150',
                '200','250','300',
                '350','400','450',
                '500','550','600',
                '650','700','750',
                '775','800','825',
                '850','875','900',
                '925','950','975',
                '1000'
            ],
            'year':year,
            'month':'%02d'%(dd),
            'day':[
            '01','02','03',
            '04','05','06',
            '07','08','09',
            '10','11','12',
            '13','14','15',
            '16','17','18',
            '19','20','21',
            '22','23','24',
            '25','26','27',
            '28','29','30',
            '31'
            ],
            'time':[
                '00:00','01:00','02:00',
                '03:00','04:00','05:00',
                '06:00','07:00','08:00',
                '09:00','10:00','11:00',
                '12:00','13:00','14:00',
                '15:00','16:00','17:00',
                '18:00','19:00','20:00',
                '21:00','22:00','23:00'
            ],
            "area": [latmax,lonmin,latmin,lonmax],
            
        },
        outdir+'\\era5_'+year+'%02d'%(dd)+'_quv.nc')
        
    #######
    c.retrieve(
            'reanalysis-era5-land-monthly-means',
            {
                'format':'netcdf',
                'variable':[
                    'total_evaporation','snow_evaporation','total_precipitation'
                ],
                'year':year,
                'month':'%02d'%(dd),
                "area": [latmax,lonmin,latmin,lonmax],
                'time':'00:00',
                'product_type':'monthly_averaged_reanalysis'
            },
            outdir+'\\era5_'+year+'%02d'%(dd)+'_rain.nc')
