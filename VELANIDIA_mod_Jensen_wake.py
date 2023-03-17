from datetime import datetime
from typing import List, Union, Any

import pandas as pd
import numpy as np
import xlsxwriter
import os
import math
import statistics
import kmeans1d
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
import csv

###Topographic info
##Altitude above sea level in m (altitude of the position +  hub
##height  %%for the moloch model it is only the altitude above sea level (add 50, 80 and 100m)

## issued according to point_coord names
##Moloch model points altitude and height agl

Alt_point=float(548)
Hmodel=np.array([50,80,100])

##WTG altitude and positions.
Alt_Velanidia=np.array([515,560,563,519])
WTG_coord_Velanidia=np.array([[39.471964,20.271951],[39.469746,20.276083],[39.470583,20.280561],[39.469556,20.286080]])

##Hub Height configurations in m.
H_hub=np.array([135,165,170])

##weather station altitude in m  -- given with the same sequence
##as the list station_name

H_wstation=float(77)
Agl_wstation=float(2.5)

##Read and write EAA and station data

date1=datetime(2020,8,4)
date2=datetime(2021,12,31)
date_generated=pd.date_range(start=date1 , end=date2)
# print(date_generated.strftime("%y%m%d"))

directory="C:/Users/plgeo/OneDrive/PC Desktop/FARIA/DATA EAA/Moloch_data"
directory_stations= "C:/Users/plgeo/OneDrive/PC Desktop/FARIA/DATA EAA/weather_station_data"
directory_to_save="C:/Users/plgeo/OneDrive/PC Desktop/FARIA/VELANIDIA results/"
# filetosave1=directory_to_save+"Chelona_old_met_mast_analysis_2022.xlsx"
# filetowrite1 = pd.ExcelWriter(filetosave1, engine='xlsxwriter')##excel file to save processed data

EAA_file1=directory+"/"+"points_EAA.xlsx"
print(EAA_file1)
point_directory=directory+"/l"

info_EAA= pd.read_excel(str(EAA_file1), sheet_name="Sheet1",header=0)
point_coords=list(info_EAA.iloc[:,0])##reads all rows of a specific column , e.g. [:,1] reads all rows of column 1
location_name=list(info_EAA.iloc[:,1])

def sign(x):
    if x>0:
       return 1
    elif x==0:
       return 0
    else:
       return -1

def distlatlon(lat1,lat2,lon1,lon2):
    lat1=math.radians(lat1)
    lat2=math.radians(lat2)
    lon1=math.radians(lon1)
    lon2 = math.radians(lon2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (math.sin(dlat / 2) ** 2 ) + math.cos(lat1) * math.cos(lat2) * (math.sin(dlon / 2) ** 2)
    c = 2 * math.asin(math.sqrt(a))

    ## earth radius in KM
    r = 6371
    D = c * r ## Shortest distance between points in km
    return D

##conversion of decimal degrees to egsa87 ( greek geodetic system of X,Y)
def latlon_to_xy(lat,lon):
    lat=math.radians(lat)
    lon=math.radians(lon)
    GE_WGS84_Alpha=6378137.000
    GE_WGS84_F_INV=298.257223563
    ##displace_geodetic_system
    e2=1-(1-1/GE_WGS84_F_INV)**2
    rad=GE_WGS84_Alpha/math.sqrt(1-e2*math.sin(lat)*math.sin(lat))
    ##spherical coordinates inverse transformation
    x=rad*math.cos(lat)*math.cos(lon)
    y=rad*math.cos(lat)*math.sin(lon)
    z=rad*(1-e2)*math.sin(lat)
    ##displace geocentric system
    x2=x+199.72
    y2=y-74.03
    z2=z-246.02
    ##ellipsoid xyz to philambda
    aradius=GE_WGS84_Alpha
    ##sphere xyz to philambda
    if abs(z2)<aradius:
       phi2=math.asin(z2/aradius)
    else:
        if z2>0:
           phi2=0.5*math.pi
        else:
           phi2=-0.5*math.pi
    if abs(x2)>0.001:
        lambda2=math.atan(y2/x2)
    else:
        if y2>0:
            lambda2=0.5*math.pi
        else:
            lambda2=-0.5*math.pi
    if x2<0:
       lambda2=math.pi-lambda2
    f=1/1/GE_WGS84_F_INV
    et2=e2/((1-f)*(1-f))
    phi2=math.atan(z2*(1+et2)/math.sqrt(x2**2+y2**2))
    acount=0
    aradius_old=10**30 ## a very large number
    while(abs(aradius-aradius_old)>0.00005  and acount<100):
        acount+=1
        aradius_old=aradius
        aradius=GE_WGS84_Alpha/math.sqrt(1-e2*math.sin(phi2)*math.sin(phi2)) ##ellipsoid_main_normal_section_radius(phi2)
        phi2=math.atan((z2+e2*aradius*math.sin(phi2))/math.sqrt(x2**2+y2**2))
    ##project philambda_to_xy
    kappa0=0.9996
    lambda0=24*math.pi
    lambda0=lambda0/180.00
    xoffset=500000
    yoffset=0.00
    dl=lambda2-lambda0
    t=math.tan(phi2)
    n2=(e2*math.cos(phi2)*math.cos(phi2))/(1-e2)
    L=dl*math.cos(phi2)
    Ni=GE_WGS84_Alpha/math.sqrt(1-e2*math.sin(phi2)*math.sin(phi2))
    ##Mi=ellipsoid_arc(phi2)
    e4=e2**2
    e6=e4*e2
    e8=e6*e2
    M0 = 1 + 0.75 * e2 + 0.703125 * e4 + 0.68359375 * e6 + 0.67291259765625 * e8
    M2 = 0.375 * e2 + 0.46875 * e4 + 0.5126953125 * e6 + 0.538330078125 * e8
    M4 = 0.05859375 * e4 + 0.1025390625 * e6 + 0.25 * e8
    M6 = 0.01139322916666667 * e6 + 0.025634765625 * e8
    M8 = 0.002408551504771226 * e8
    Mi = GE_WGS84_Alpha * (1 - e2) * (M0 * phi2 - M2 * math.sin(2 * phi2) + M4 * math.sin(4 * phi2) - M6 * math.sin(6 * phi2) + M8 * math.sin(8 * phi2))
    x = (((5 - 18 * t * t + t * t * t * t + 14 * n2 - 58 * t * t * n2) * L * L / 120.00 + (1 - t * t + n2) / 6.00) * L * L + 1) * L * kappa0 * Ni + xoffset
    y = Mi + (Ni * t / 2) * L * L + (Ni * t / 24) * (5 - t * t + 9 * n2 + 4 * n2 * n2) * L * L * L * L + (Ni * t / 720) * (61 - 58 * t * t) * L * L * L * L * L * L
    y = y * kappa0 + yoffset
    return x,y

##test coordinate conversion
# WTG_coord=WTG_coord_Iasmos
# WTG_coord_help=list(WTG_coord)
# for ii in range(len(WTG_coord_help)):
#     lat_wtg=WTG_coord[ii,0]
#     lon_wtg=WTG_coord[ii,1]
#     x,y=latlon_to_xy(lat_wtg,lon_wtg)
#     print(x)
#     print(y)

point=point_coords[1]
location=location_name[1]
Timestamp_50 = []
W50=np.zeros(shape=(len(date_generated),48),dtype=float); W80=np.zeros(shape=(len(date_generated),48),dtype=float);W100=np.zeros(shape=(len(date_generated),48),dtype=float)
D50=np.zeros(shape=(len(date_generated),48),dtype=float);D80=np.zeros(shape=(len(date_generated),48),dtype=float);D100=np.zeros(shape=(len(date_generated),48),dtype=float)
for j in range(len(date_generated)):
    datetoread = date_generated[j].strftime("%y%m%d")  ##converts dates to specific string format
    filetoread = point_directory + "_" + str(point) + "/" + "moloch_" + datetoread + "00.csv"
    if os.path.exists(filetoread):
        data = pd.read_csv(filetoread, sep=",")
        height_vector = list(data.iloc[:, 0])
        index_50 = [];
        index_80 = [];
        index_100 = []  ##equivalent code of find all indices of an array element in matlab
        for k in range(len(height_vector)):
            if height_vector[k] == 50:
                index_50.append(k)
            if height_vector[k] == 80:
                index_80.append(k)
            if height_vector[k] == 100:
                index_100.append(k)
        if len(index_50) <= 48:
            index_en50 = index_50[-1]
        else:
            index_en50 = 48
        timestamp_50 = list(data.iloc[0:index_en50, 1])
        u50 = np.array(data.iloc[0:index_en50, 2])
        v50 = np.array(data.iloc[0:index_en50, 3])
        if len(index_80) > 48:
            index_st80 = index_80[0]
            index_en80 = index_st80 + 47
        else:
            index_st80 = index_80[0]
            index_en80 = index_80[-1]

        timestamp_80 = list(data.iloc[index_st80:index_en80 + 1, 1])
        u80 = np.array(data.iloc[index_st80:index_en80 + 1, 2])
        v80 = np.array(data.iloc[index_st80:index_en80 + 1, 3])
        if len(index_100):
            if len(index_100) > 48:
                index_st100 = index_100[0]
                index_en100 = index_st100 + 47
            else:
                index_st100 = index_100[0]
                index_en100 = index_100[-1]

            timestamp_100 = list(data.iloc[index_st100:index_en100 + 1, 1])
            u100 = np.array(data.iloc[index_st100:index_en100 + 1, 2])
            v100 = np.array(data.iloc[index_st100:index_en100 + 1, 3])
        else:
            timestamp_100 = list(data.iloc[0:48, 1])
            u100 = np.array(data.iloc[0:48, 4])
            v100 = np.array(data.iloc[0:48, 5])

        for i in range(u50.size):
            uh50 = u50[i]
            vh50 = v50[i]
            W50[j, i] = math.sqrt((uh50 ** 2) + (vh50 ** 2))
        for jj in range(u80.size):
            uh80 = u80[jj]
            vh80 = v80[jj]
            W80[j, jj] = math.sqrt((uh80 ** 2) + (vh80 ** 2))
        for k in range(u100.size):
            uh100 = u100[k]
            vh100 = v100[k]
            W100[j, k] = math.sqrt((uh100 ** 2) + (vh100 ** 2))

        dir50 = np.zeros((u50.size, 1))
        u_50 = u50
        v_50 = v50
        del u50
        del v50
        for kk in range(u_50.size):
            u50 = u_50[kk]
            v50 = v_50[kk]
            if u50 == 0 and sign(v50) == 1:
                dr = 180
            elif u50 == 0 and sign(v50) == -1:
                dr = 0
            elif sign(u50) == 1 and sign(v50) == 1:
                dr = 270 - math.degrees(math.atan(abs(v50 / u50)))
            elif sign(u50) == 1 and sign(v50) == -1:
                dr = 270 + math.degrees(math.atan(abs(v50 / u50)))
            elif sign(u50) == -1 and sign(v50) == -1:
                dr = 90 - math.degrees(math.atan(abs(v50 / u50)))
            elif sign(u50) == -1 and sign(v50) == 1:
                dr = 90 + math.degrees(math.atan(abs(v50 / u50)))
            if dr < 0:
                dr = 0
            D50[j, kk] = dr

        u_80 = u80
        v_80 = v80
        del u80
        del v80
        del dr
        for ll in range(u_80.size):
            u80 = u_80[ll]
            v80 = v_80[ll]
            if u80 == 0 and sign(v80) == 1:
                dr = 180
            elif u80 == 0 and sign(v80) == -1:
                dr = 0
            elif sign(u80) == 1 and sign(v80) == 1:
                dr = 270 - math.degrees(math.atan(abs(v80 / u80)))
            elif sign(u80) == 1 and sign(v80) == -1:
                dr = 270 + math.degrees(math.atan(abs(v80 / u80)))
            elif sign(u80) == -1 and sign(v80) == -1:
                dr = 90 - math.degrees(math.atan(abs(v80 / u80)))
            elif sign(u80) == -1 and sign(v80) == 1:
                dr = 90 + math.degrees(math.atan(abs(v80 / u80)))
            if dr < 0:
                dr = 0
            D80[j, ll] = dr

        u_100 = u100
        v_100 = v100
        del u100
        del v100
        del dr
        for mm in range(u_100.size):
            u100 = u_100[mm]
            v100 = v_100[mm]
            if u100 == 0 and sign(v100) == 1:
                dr = 180
            elif u100 == 0 and sign(v100) == -1:
                dr = 0
            elif sign(u100) == 1 and sign(v100) == 1:
                dr = 270 - math.degrees(math.atan(abs(v100 / u100)))
            elif sign(u100) == 1 and sign(v100) == -1:
                dr = 270 + math.degrees(math.atan(abs(v100 / u100)))
            elif sign(u100) == -1 and sign(v100) == -1:
                dr = 90 - math.degrees(math.atan(abs(v100 / u100)))
            elif sign(u100) == -1 and sign(v100) == 1:
                dr = 90 + math.degrees(math.atan(abs(v100 / u100)))
            if dr < 0:
                dr = 0
            D100[j, mm] = dr

    Timestamp_50.append(timestamp_50)
    Timestamp_model = np.array((Timestamp_50))
Timestamp_model = Timestamp_model.flatten()
Timestamp_model = Timestamp_model.reshape(-1, 1)
W50 = W50.flatten()  # converts a 2-D array to a 1-D array - row vector
W_50 = W50.reshape(-1, 1)  # converts row vector to a column vector
W80 = W80.flatten();
W_80 = W80.reshape(-1, 1)
W100 = W100.flatten();
W_100 = W100.reshape(-1, 1)
D50 = D50.flatten();
D_50 = D50.reshape(-1, 1)
D80 = D80.flatten();
D_80 = D80.reshape(-1, 1)
D100 = D100.flatten();
D_100 = D100.reshape(-1, 1)
Model_data = np.column_stack([Timestamp_model, W_50, W_80, W_100, D_50, D_80,D_100])  # constructs a 2-D array where rows are the time frames and columns are the variable values

##weather station data

directory_stations="C:/Users/plgeo/OneDrive/PC Desktop/FARIA/DATA EAA/weather_station_data"
station_name=["igoumenitsa","stratoni","asprovalta","imeros"]
station_coord=[[39.541749,20.279882]];station_coord.append([40.515033,23.828907]);station_coord.append([40.724931,23.711809]);station_coord.append([40.955700,25.369020])

name = station_name[0]
stationfile = directory_stations + "/" + name + ".txt"
stationdata = open(stationfile, "r")
lines = stationdata.readlines()[1:]  ##returns all lines except first line
stationdata.close()
## read data line by line and split each line to create columns of variable values
date_10min = []
time_10min = []
Temp_10min = []
Pressure_10min = []
RH_10min = []
Wind_10min = []
Dir_10min = []
Wind_gust_10min = []

for line in lines:
    parts = line.split()
    if len(parts) > 1 and parts[2] != '---':
        date = str(parts[0])
        time = str(parts[1])
        temp_10min = float(parts[2])
        pressure_10min = float(parts[6])
        rh_10min = float(parts[5])
        wind_10min = float(parts[8])
        dir_10min = str(parts[9])
        wind_gust_10min = float(parts[10])
        Wind_10min.append(wind_10min)  ##stacks along columns
        Dir_10min.append(dir_10min)
        Wind_gust_10min.append(wind_gust_10min)
        Pressure_10min.append(pressure_10min)
        RH_10min.append(rh_10min)
        Temp_10min.append(temp_10min)
        date_10min.append(date)
        time_10min.append(time)
Timestamp = []
Temp = []
Wind = []
Pressure = []
Wind_gust = []
RH = []
Direction = []
for j in range(3, len(Wind_10min), 3):
    wind = statistics.mean(Wind_10min[j - 2:j + 1]) / 3.6
    wind_gust = max(Wind_gust_10min[j - 2:j + 1]) / 3.6
    direction = Dir_10min[j]
    pressure = statistics.mean(Pressure_10min[j - 2:j + 1])
    temp = statistics.mean(Temp_10min[j - 2:j + 1])
    rh = statistics.mean(RH_10min[j - 2:j + 1])
    date_var = str(date_10min[j])
    time_var = str(time_10min[j])
    time_var_help = time_var.split(':')
    hours_stamp = time_var_help[0]
    minutes = time_var_help[1]
    hours = int(hours_stamp)
    ##the code below ensures that during conversion from 10min to 30min, the timestamp will be either at :00 hour or :30 hour
    if hours_stamp != "23":
        if hours < 9:
            if minutes == "10":
                minutes = "00"
            elif minutes == "20" or minutes == "40":
                minutes = "30"
            elif minutes == "50":
                hours += 1
                minutes = "00"
                hours_stamp = "0" + str(hours)
        elif hours >= 9:
            if minutes == "10":
                minutes = "00"
            elif minutes == "20" or minutes == "40":
                minutes = "30"
            elif minutes == "50":
                hours += 1
                minutes = "00"
                hours_stamp = str(hours)
    elif hours_stamp == "23":
        if minutes == "10":
            minutes = "00"
        elif minutes == "20" or minutes == "40":
            minutes = "30"
        elif minutes == "50":
            hours = 0
            minutes = "00"
            hours_stamp = "0" + str(hours)
    time_var = hours_stamp + ":" + minutes
    timestamp = date_var + " " + time_var
    Pressure.append(pressure)
    Temp.append(temp)
    Wind.append(wind)
    Wind_gust.append(wind_gust)
    RH.append(rh)
    Direction.append(direction)
    Timestamp.append(timestamp)
Variable_names = ["timestamp", "Pressure", "Temperature", "Wind Speed", "Wind gust", "Relative Humidity","Wind Direction"]
Station_variables = np.dstack((np.array(Timestamp), np.array(Pressure), np.array(Temp), np.array(Wind),np.array(Wind_gust), np.array(RH),np.array(Direction)))  ##creates an array of variables columnwise
Station_variables = Station_variables[0, :, :]
Station_data = np.vstack((np.array(Variable_names), Station_variables))##concatenates variable names with variable values rowise
# print(Station_data)
del line

# ##read and process wind mast raw data - this is old masts data

del i
del j
del k
# del ii
# del kk
# del jj
del ll
del mm

#correlation between weather station and Moloch data
datestart = datetime.strptime("2020-08-04 00:00", "%Y-%m-%d %H:%M")
dateend = datetime.strptime("2021-12-31 23:30", "%Y-%m-%d %H:%M")
dates_generated = pd.date_range(datestart, dateend, freq='30T')
dates_generated=dates_generated.strftime("%Y-%m-%d %H:%M")
datestart=datestart.strftime("%Y-%m-%d %H:%M")##converts datetime to string with the specified format
dateend=dateend.strftime("%Y-%m-%d %H:%M")

date_s2020=datetime.strptime("2020-01-01 00:30","%Y-%m-%d %H:%M");date_s2020=date_s2020.strftime("%Y-%m-%d %H:%M")
date_e2020=datetime.strptime("2020-12-31 23:30","%Y-%m-%d %H:%M");date_e2020=date_e2020.strftime("%Y-%m-%d %H:%M")
date_s2021=datetime.strptime("2021-01-01 00:00","%Y-%m-%d %H:%M");date_s2021=date_s2021.strftime("%Y-%m-%d %H:%M")
date_e2021=datetime.strptime("2021-12-31 23:30","%Y-%m-%d %H:%M");date_e2021=date_e2021.strftime("%Y-%m-%d %H:%M")

hup=100##reference height agl
Md=28.9634*(10**(-3)) ## apparent molecular weight of dry air
Mw=18.01528*(10**(-3)) ## apparent molecular weight of water vapor
Rconst=8.31441 ## specific gas constant for dry air in J/Kg*K
G=9.80665 ## gravitational acceleration in m/s^2
L=0.0065 ## mean temperature lapse rate in lower troposphere in oC/m
Zc=0.99996 ## compressibility factor of the atmospheric air
Rad=6535*(10**3) ## earth's radius in m

moloch_coord=point.split('_')
moloch_coord=np.array(moloch_coord,dtype=float)
WTG_coord = WTG_coord_Velanidia
WTG_alt = Alt_Velanidia
station_info = Station_data
st_coord = station_coord[0]
st_alt = H_wstation
href = Agl_wstation
zo = float(0.15)#roughness length, provide 2 different versions if plausible, i.e. 0.15 and 0.30
lat_model = moloch_coord[0]
lon_model = moloch_coord[1]


lat_station = st_coord[0]
lon_station = st_coord[1]
moloch_data = Model_data
station_data = station_info[1:, :]  # this is going to skip first row, remember this is a numpy array
# print(station_data.shape)##check to see how many rows does each timeseries have ( how many timestamps) by checking no of rows and cols
station_timestamp = station_data[:, 0]  ##access of first column of numpy array
index_s2020=[];index_e2020=[];index_s2021=[];index_e2021=[];
index_start = [];
index_end = []
for j in range(station_timestamp.size):
    s_ts = station_timestamp[j]
    sts_minute = s_ts.split(':')
    sts_min = sts_minute[1]
    if s_ts == datestart:
       index_start.append(j)
    if s_ts == dateend:
       index_end.append(j)
    if s_ts==date_s2020:
       index_s2020.append(j)
    if s_ts==date_e2020:
       index_e2020.append(j)
    if s_ts==date_s2021:
       index_s2021.append(j)
    if s_ts==date_e2021:
       index_e2021.append(j)
index_s2020 = index_s2020[0]
index_start=index_start[0]
index_end=index_end[0]
    ##check if each list is empty and modify accordingly
if not index_e2020:
   index_e2020=17520
else:
   index_e2020=index_e2020[0]
if not index_s2021:
   index_s2021=17521
elif index_s2021:
   index_s2021=index_s2021[0]
if not index_e2021:
   index_e2021=len(station_timestamp)-48
else:
   index_e2021=index_e2021[0]

station_timestamp_end = station_timestamp[index_end - 1]
station_timestamp_help = station_timestamp[index_start:index_end]
station_w2020 = station_data[index_s2020:index_e2020 + 1, 3]
station_w2020 = station_w2020.astype('f')  ##convert array type from string to float
station_d2020 = station_data[index_s2020:index_e2020 + 1, 6]
station_w2021 = station_data[index_s2021:index_e2021 + 1, 3];
station_w2021 = station_w2021.astype('f')
station_d2021 = station_data[index_s2021:index_e2021 + 1, 6]
station_meanw2020 = np.mean(station_w2020, dtype=float)
station_sdvw2020 = np.std(station_w2020, dtype=float)
station_meanw2021 = np.mean(station_w2021, dtype=float)
station_sdvw2021 = np.std(station_w2021, dtype=float)

station_cl10m2020 = [station_meanw2020 * (math.log(10 / zo) / math.log(href / zo)), station_sdvw2020]  ##gives station climate  at 10m agl
station_cl10m2021 = [station_meanw2021 * (math.log(10 / zo) / math.log(href / zo)), station_sdvw2021]  ##gives station climate  at 10m agl

# print(station_cl10m2020)
# print(station_cl10m2021)
station_cl2020 = [station_meanw2020 * (math.log(hup / zo) / math.log(href / zo)), station_sdvw2020 * (math.log(hup / zo) / math.log(href / zo))]  ##gives station climate directly at 100m agl
station_cl2021=[station_meanw2021 * (math.log(hup / zo) / math.log(href / zo)), station_sdvw2021 * (math.log(hup / zo) / math.log(href / zo))]##gives station climate directly at 100m agl
Wmean_year1=station_cl2020[0]
Wmean_year2=station_cl2021[0]
Wmean_station_annual=statistics.mean([Wmean_year1,Wmean_year2])

##extract mean climate from the 3 old masts
directory_mast='C:/Users/plgeo/OneDrive/PC Desktop/FARIA/HPEIROS old masts results/'

mast1_file=directory_mast+'Varathron_old_met_mast_analysis_2022.xlsx'##VARATHRON is the reference mast.
mast2_file=directory_mast+'Chelona_old_met_mast_analysis_2022.xlsx'## CHELONA is the 2nd reference mast, highly correlated with VARATHRON in terms of wind speed and direction.
mast3_file=directory_mast+'Gymnoula_old_met_mast_analysis_2022.xlsx'## GYMNOULA is the 3rd reference mast, not very correlated with the other 2 masts.

dataframe1_year1=pd.read_excel(mast1_file,sheet_name='annual 2010',index_col=False)
dataframe1_year2=pd.read_excel(mast1_file,sheet_name='annual 2011',index_col=False)
dataframe2_year1=pd.read_excel(mast2_file,sheet_name='annual 2010',index_col=False)
dataframe2_year2=pd.read_excel(mast2_file,sheet_name='annual 2011',index_col=False)
dataframe3_year1=pd.read_excel(mast3_file,sheet_name='annual 2010',index_col=False)
dataframe3_year2=pd.read_excel(mast3_file,sheet_name='annual 2011',index_col=False)

height_list1=[];height_list2=[];height_list3=[]
Wold1=[];Wold2=[];Wold3=[]
## Obtain anemometer heights, annual mean wind speed at each height of interest
for h in range(3):
    if h!=1:## the intermediate anemometer is used as a control anemometer, it's output should not be taken into account
       hl_1=dataframe1_year1.iloc[0,h+1]
       hl_1=hl_1.split('m')
       hl1=float(hl_1[0])
       height_list1.append(hl1)
       w1_y1=float(dataframe1_year1.iloc[1,h+1])
       w1_y2=float(dataframe1_year2.iloc[1,h+1])
       w1_y=0.5*(w1_y1+w1_y2)## mean annual wind speed
       Wold1.append(w1_y)
       hl_2=dataframe2_year1.iloc[0,h+1]
       hl_2=hl_2.split('m')
       hl2=float(hl_2[0])
       height_list2.append(hl2)
       w2_y1=float(dataframe2_year1.iloc[1,h+1])
       w2_y2=float(dataframe2_year2.iloc[1,h+1])
       w2_y = 0.5 * (w2_y1 + w2_y2) ## mean annual wind speed
       Wold2.append(w2_y);
       hl_3=dataframe3_year1.iloc[0,h+1]
       hl_3=hl_3.split('m')
       hl3=float(hl_3[0])
       height_list3.append(hl3)
       w3_y1=float(dataframe3_year1.iloc[1,h+1])
       w3_y2=float(dataframe3_year2.iloc[1,h+1]) ## mean annual wind speed
       w3_y = 0.5 * (w3_y1 + w3_y2)  ## mean annual wind speed
       Wold3.append(w3_y)

##derivation of mean wind from 3 old mast mean data at the position of the weather station and at the position of moloch model , using spatial interpolation with inverse distance weighting
point_old1=[39.529389,20.303361];point_old2=[39.559258,20.338611];point_old3=[39.507648,20.421390]
lat_old1=point_old1[0];lon_old1=point_old1[1]
lat_old2=point_old2[0];lon_old2=point_old2[1]
lat_old3=point_old3[0];lon_old3=point_old3[1]
WeD1 = distlatlon(lat_old1, lat_station, lon_old1, lon_station)
WeD2 = distlatlon(lat_old2, lat_station, lon_old2, lon_station)
WeD3= distlatlon(lat_old3, lat_station, lon_old3, lon_station)
WeD1 = 1/(WeD1 ** 2)
WeD2 = 1/(WeD2 ** 2)
WeD3 = 1/(WeD3 ** 2)
WemD1 = distlatlon(lat_old1, lat_model, lon_old1, lon_model)
WemD2 = distlatlon(lat_old2, lat_model, lon_old2, lon_model)
WemD3= distlatlon(lat_old3, lat_model, lon_old3, lon_model)
WemD1 = 1/(WemD1 ** 2)
WemD2 = 1/(WemD2 ** 2)
WemD3 = 1/(WemD3 ** 2)
# print(WeD1/(WeD1+WeD2+WeD3));print(WeD2/(WeD1+WeD2+WeD3));print(WeD3/(WeD1+WeD2+WeD3)) ## prints the 3 factors of interpolation. no further action needed because the 3rd factor is way smaller.
## use the mean annual wind speed at reference height (30m agl)
wold1y=Wold1[1]
wold2y=Wold2[1]
wold3y=Wold3[1]
Wmean_annual=((WeD1 * wold1y)+(WeD2 * wold2y)+(WeD3 * wold3y)) / (WeD1 + WeD2 + WeD3)##mean wind climate of old mast @ 30m agl interpolated at weather station's position
Wmean_annual_mol=((WemD1 * wold1y)+(WemD2 * wold2y)+(WemD3 * wold3y)) / (WemD1 + WemD2 + WemD3)##mean wind climate of old mast @ 30m agl interpolated at moloch point's position

Wmean_50= Wmean_annual * (math.log(50 / zo) / math.log(30/zo))##mean wind climate of old mast @ 50m agl interpolated at weather station's position
Wmean_50_mol=Wmean_annual_mol*(math.log(50 / zo) / math.log(30/zo))##mean wind climate of old mast @ 50m agl interpolated at moloch point's position

# print(corr_factor)
# print(Wmean_annual)
# print(Wmean)

##obtain wind shear  and avg turbulence per dir sector from the closest to the station old mast. - the closest to the station old mast has the highest parameter WeD - inspection suggests VARATHRON - be careful this shear refers to the lowest portion of the PBL (20-30m agl)
## quality control implies that Chelona has more meaningful values of wind shear than VARATHRON
Dir_list=["N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW"]
Ashear_oldmast=[];TI_oldmast=[];S_sherror=[]; Std_corr_error=[]
for ash in range(len(Dir_list)):
    asheary1=float(dataframe2_year1.iloc[ash+13,2])
    asheary2=float(dataframe2_year2.iloc[ash+13,2])
    asheary=0.5*(asheary1+asheary2)##mean of the 2 years
    Ashear_oldmast.append(asheary)
    s_sherror1=float(dataframe1_year1.iloc[ash+13,4])
    s_sherror2=float(dataframe1_year2.iloc[ash+13,4])
    s_sherror=statistics.mean([s_sherror1,s_sherror2])
    S_sherror.append(s_sherror)
    tiy1=float(dataframe1_year1.iloc[ash+13,3])
    tiy2=float(dataframe1_year2.iloc[ash+13,3])
    tiy=0.5*(tiy1+tiy2)
    TI_oldmast.append(tiy)
S1=statistics.mean(Ashear_oldmast)## mean wind shear derived from old mast.
Wmean=Wmean_50*((hup/50)**S1)##mean wind climate of old mast @ 100m agl interpolated at weather station's position
Wmean_mol=Wmean_50_mol*((hup/50)**S1) ##mean wind climate of old mast @ 100m agl interpolated at moloch point's position
corr_factor_station=Wmean/Wmean_station_annual## correction factor of station wind speed at 100m agl based on old masts
Wmean2=Wmean_annual*((hup/30)**S1)## this proves that exponential law does not work except the lowest portion of PBL because we use a shear obtained for the 20-30m for the extr. to 100m . keep log law extrapolation.

Wmean_year1=Wmean_year1*corr_factor_station
Wmean_year2=Wmean_year2*corr_factor_station
moloch_timestamp_help=list(moloch_data[:,0])
# convert from one date format to another using datetime object.
moloch_timestamp = []
for x in moloch_timestamp_help:
    if len(x) > 16:
       moloch_date = datetime.strptime(x, "%Y-%m-%d %H:%M:%S")
       y = moloch_date.strftime("%Y-%m-%d %H:%M")
    else:
       y = x
    moloch_timestamp.append(y)
del moloch_timestamp_help
index_model_end=[]
for k in range(len(moloch_timestamp)):
    m_ts = moloch_timestamp[k]
    if m_ts==station_timestamp_end:
       index_model_end.append(k)
## correlation between station wind speed and model wind speed overall and per direction sector
index_model_end = index_model_end[0]
# print(index_model_end)
moloch_timestamp_help = moloch_timestamp[0:index_model_end + 1]
moloch_wspeed = moloch_data[0:index_model_end + 1, 1:4]
moloch_dir = moloch_data[0:index_model_end + 1, 4:7]
moloch_dir_help = moloch_dir[:, 0]
moloch_dir_help = list(moloch_dir_help)
# finds all elements of station_timestamps that exist in moloch_timestamps
##and creates a common timeseries for correlation and analysis. Also it upscales the station wind speed at 100m agl
station_timestamp_help = list(station_timestamp_help)
station_wspeed_for_common = station_data[index_start:index_end, 3]
station_wspeed_for_common = list(station_wspeed_for_common)
station_wcommon = [];
station_wcommon_up = [];
moloch_wcommon = [];
moloch_dcommon = []
index_common = []
for l in range(len(station_timestamp_help)):
    st_help = station_timestamp_help[l]
    if st_help in moloch_timestamp_help:
       ic = moloch_timestamp_help.index(st_help)
       if ic < len(station_wspeed_for_common):
         st_c = station_wspeed_for_common[ic]
         station_wcommon.append(st_c)
         w_up = np.array(st_c, dtype=float) *corr_factor_station* (math.log(hup / zo) / math.log(href / zo))  ##upscales station wind speed at 100m agl and corrects with old masts' spatially interpolated mean annual wind speed
         station_wcommon_up.append(w_up)
         mol_c = moloch_wspeed[ic, :]
         mold_c = moloch_dir_help[ic]
         moloch_wcommon.append(mol_c)
         moloch_dcommon.append(mold_c)
         index_common.append(ic)
station_wc = np.array(station_wcommon, dtype=float)
station_wc_up = np.array(station_wcommon_up, dtype=float)
station_common_up_mean = np.mean(station_wc_up)
moloch_wcommon = np.array(moloch_wcommon, dtype=float)

moloch_wcommon_50m = list(moloch_wcommon[:, 0])
moloch_wcommon_80m = list(moloch_wcommon[:, 1])
moloch_wcommon_100m = list(moloch_wcommon[:, 2])
Wind_mean_annual_moloch=statistics.mean(moloch_wcommon_100m)##mean moloch wind speed at 100m agl

corr_factor_moloch=Wmean_mol/Wind_mean_annual_moloch## geographical correction factor based on spatial interpolation from the 3 old masts.

input_dirm = []
labels=16
for zz in moloch_dcommon:
    z_help = float(zz)
    input_dirm.append(z_help)
clusters_m, centroids_m = kmeans1d.cluster(input_dirm,labels)  ##cluster wind direction of moloch model to 16 sectors with kmeans algorithm, centroids are the class centre and clusters are the sector number of each element, total 16
Station_dclass = [];R_sq = [];Linear_regression_model = [];Dclass_perc_station = [];Ashear_model=[]; S_shmodelerror=[]
for ii in range(len(centroids_m)):
    index_d = [];
    for jj in range(len(clusters_m)):
       class_d = clusters_m[jj]
       if class_d == ii:
          index_d.append(jj)
# print(index_d)
    moloch_wcommon_per_sector_50m = [];
    moloch_wcommon_per_sector_80m = [];
    moloch_wcommon_per_sector_100m = [];
    station_wcommon_per_sector = [];
    moloch_dcommon_per_sector = [];
    station_wcommon_up_per_sector = []
    for xx in index_d:
        moloch_wc_ps50 = moloch_wcommon_50m[xx]
        moloch_wc_ps80 = moloch_wcommon_80m[xx]
        moloch_wc_ps100 = moloch_wcommon_100m[xx]
        ## correct Moloch wind speed based on spatial interpolation from 3 old masts.
        moloch_wc_ps50=moloch_wc_ps50*corr_factor_moloch
        moloch_wc_ps80 = moloch_wc_ps80 * corr_factor_moloch
        moloch_wc_ps100 = moloch_wc_ps100 * corr_factor_moloch
        moloch_dc = input_dirm[xx]
        station_wc_ps = station_wcommon[xx]
        station_wc_ps_up = station_wcommon_up[xx]
        moloch_wcommon_per_sector_50m.append(moloch_wc_ps50)
        moloch_wcommon_per_sector_80m.append(moloch_wc_ps80)
        moloch_wcommon_per_sector_100m.append(moloch_wc_ps100)
        station_wcommon_per_sector.append(station_wc_ps)
        station_wcommon_up_per_sector.append(station_wc_ps_up)
        moloch_dcommon_per_sector.append(moloch_dc)
    dclass_perc = len(index_d) / len(station_wcommon)  ##percentage of wind direction per sector
    wcom50 = np.mean(moloch_wcommon_per_sector_50m, dtype=float)
    wcom80 = np.mean(moloch_wcommon_per_sector_80m, dtype=float)
    wcom100 = np.mean(moloch_wcommon_per_sector_100m, dtype=float)
    ## fit a linear regression model to the log of the moloch wind speed to derive wind shear per sector
    list_w50m_nn = [x1 + float(1.1) if x1 < 1.1 else x1 for x1 in moloch_wcommon_per_sector_50m]  ## exclude zeros or negative natural logarithms.
    list_w80m_nn = [x2 + float(1.1) if x2 < 1.1 else x2 for x2 in moloch_wcommon_per_sector_80m]
    list_w100m_nn = [x3 + float(1.1) if x3 < 1.1 else x3 for x3 in moloch_wcommon_per_sector_100m]
    log_w50 = np.log(list_w50m_nn, dtype=float)
    log_w80 = np.log(list_w80m_nn, dtype=float)
    log_w100 = np.log(list_w100m_nn, dtype=float)
    linear_log_model1 = LinearRegression(fit_intercept=True).fit(log_w50.reshape((-1, 1)), log_w80)
    linear_log_model2 = LinearRegression(fit_intercept=True).fit(log_w80.reshape((-1, 1)), log_w100)
    slope1 = linear_log_model1.coef_[0]
    slope2 = linear_log_model2.coef_[0]
    const1 = linear_log_model1.intercept_
    const2 = linear_log_model2.intercept_
    rsqm1 = linear_log_model1.score(log_w50.reshape((-1, 1)), log_w80)
    rsqm2 = linear_log_model2.score(log_w80.reshape((-1, 1)), log_w100)
    ashearmodel_1 = ((slope1 - 1) * np.mean(log_w50) + const1) / math.log(80 / 50)
    ashearmodel_2 = ((slope2 - 1) * np.mean(log_w80) + const2) / math.log(100 / 80)
    ashear_model = np.mean([ashearmodel_1, ashearmodel_2])
    Ashear_model.append(ashear_model)
    s_mlog = statistics.stdev(log_w100)
    s_shmodelerror1 = math.sqrt((1 - rsqm1) * (s_mlog ** 2))
    s_shmodelerror2 = math.sqrt((1 - rsqm2) * (s_mlog ** 2))
    s_shmodelerror = statistics.mean([s_shmodelerror1, s_shmodelerror2])  ## standard deviation of the shear factor fitting error.
    Dclass_perc_station.append(dclass_perc)  ##list of wind direction sector percentage of data
    S_shmodelerror.append(s_shmodelerror)
    Dclass_perc_station.append(dclass_perc)  ##list of wind direction sector percentage of data
    mean_dclass = min(moloch_dcommon_per_sector)  ##the minimum direction of each direction class, used in order to classify station direction
# print(mean_dclass)
    if mean_dclass < 348.5 and (mean_dclass > 326 or mean_dclass == 326):
       station_dclass = "NNW"
    elif mean_dclass < 326 and (mean_dclass > 303.5 or mean_dclass == 303.5):
       station_dclass = "NW"
    elif mean_dclass < 303.5 and (mean_dclass > 281 or mean_dclass == 281):
       station_dclass = "WNW"
    elif mean_dclass < 281 and (mean_dclass > 258.5 or mean_dclass == 258.5):
       station_dclass = "W"
    elif mean_dclass < 258.5 and (mean_dclass > 236 or mean_dclass == 236):
       station_dclass = "WSW"
    elif mean_dclass < 236 and (mean_dclass > 213.5 or mean_dclass == 213.5):
       station_dclass = "SW"
    elif mean_dclass < 213.5 and (mean_dclass > 191 or mean_dclass == 191):
       station_dclass = "SSW"
    elif mean_dclass < 191 and (mean_dclass > 168.5 or mean_dclass == 168.5):
       station_dclass = "S"
    elif mean_dclass < 168.5 and (mean_dclass > 146 or mean_dclass == 146):
       station_dclass = "SSE"
    elif mean_dclass < 146 and (mean_dclass > 123.5 or mean_dclass == 123.5):
       station_dclass = "SE"
    elif mean_dclass < 123.5 and (mean_dclass > 101 or mean_dclass == 101):
       station_dclass = "ESE";
    elif mean_dclass < 101 and (mean_dclass > 78.5 or mean_dclass == 78.5):
       station_dclass = "E";
    elif mean_dclass < 78.5 and (mean_dclass > 56 or mean_dclass == 56):
       station_dclass = "ENE";
    elif mean_dclass < 56 and (mean_dclass > 33.5 or mean_dclass == 33.5):
       station_dclass = "NE";
    elif mean_dclass<33.5 and (mean_dclass>11 or mean_dclass==11):
       station_dclass= "NNE"
    elif mean_dclass < 11 or (mean_dclass > 348.5 or mean_dclass == 348.5):
       station_dclass = "N";
    Station_dclass.append(station_dclass)  ##provides a list of the direction classes (16 sectors)
    ## relate upscaled wind speed from weather station with moloch wind speed at 100m agl per direction sector - linear regression model
    x_linear = np.array(station_wcommon_up_per_sector, dtype=float).reshape((-1, 1))
    y_linear = np.array(moloch_wcommon_per_sector_100m, dtype=float)
    linear_wind_model = LinearRegression(fit_intercept=True).fit(x_linear, y_linear)
    r_sq = linear_wind_model.score(x_linear, y_linear)  ##coefficient of determination or Rsquared
    R_sq.append(r_sq)
    s_y=statistics.stdev(y_linear)## standard deviation of the moloch model wind speed
    s_e=math.sqrt((1-r_sq)*(s_y**2))## standard deviation of the correlation error
    Std_corr_error.append(s_e)
    Linear_regression_model.append(linear_wind_model)

std_cerror=statistics.mean(Std_corr_error)## mean standard deviation of the correlation error

station_wspeed = list(station_w2020)
station_wspeed.extend(list(station_w2021))  ##creates a list with the full 2-year sequence of station wind speeds
station_dir = list(station_d2020)
station_dir.extend(list(station_d2021))  ##creates a list with the full 2-year sequence of station wind dirs
##upscale wind speed at weather station at 100m agl and corrects it with the spatially interpolated annual mean wind by the 3 old masts
station_wspeed_up = []
for ws in station_wspeed:
    ws = float(ws)
    ws_up = ws * corr_factor_station * (math.log(hup / zo) / math.log(href / zo))
    station_wspeed_up.append(ws_up)
del ws_up
del index_d
# print(statistics.mean(station_wspeed_up))
# create a complete 2 year time series of wind speed at 100m agl based on combined data from weather station and moloch model using the linear model per direction sector
Wind_model_up = [];Dir_model=[];Timestamp_m=[]
for kk in range(len(station_wspeed)):
    ws_up = station_wspeed_up[kk]
    station_d = station_dir[kk]
    timest=Timestamp[kk]
    if station_d!='---':
       for xxx in range(len(Station_dclass)):
           if station_d==Station_dclass[xxx]:
              d_class=xxx
    else:
       d_class=12## dominant direction class as derived from old masts.
    dir_station=centroids_m[d_class]
    doffset=11
    dir_station=dir_station-doffset
    Dir_model.append(dir_station)#converts verbal direction to deg. for the station
    exist_count = Station_dclass.count(station_d)  ##checks how many times the specific wind direction exists in the class of station directions created above
    if exist_count > 0:
      index_d = Station_dclass.index(station_d)
     # returns the index of the 1st time station_dir exists in the Station_dclass list
    else:
        index_d = R_sq.index(max(R_sq))
    # returns the index of the most well fitted direction class
    linear_wmodel = Linear_regression_model[index_d]
    wind_model = linear_wmodel.predict(np.array(ws_up).reshape(1, -1))
    # wind_model=ws_up
    wind_model = float(wind_model)
    Wind_model_up.append(wind_model)  ##complete 2 year time series of wind speed at 100m agl based on combined data from weather station and moloch model
    Timestamp_m.append(timest)
del exist_count
del index_d
del mean_dclass

moloch_wspeed_help = moloch_wspeed[:, 2];
moloch_wspeed_help = moloch_wspeed_help.astype('f')

##choose the timeseries with the highest mean wind speed in the common available time interval
mean_moloch_speed = statistics.mean(moloch_wspeed_help)
mean_recreated_speed = statistics.mean(Wind_model_up[-1 - len(moloch_wspeed_help) + 1:])

def myfunc1(moloch_wspeed_help):
    Wind_model_up[-1 - len(moloch_wspeed_help) + 1:] = list(moloch_wspeed_help)

def myfunc2():
    return

if mean_moloch_speed > mean_recreated_speed:
    myfunc1(moloch_wspeed_help)
else:
    myfunc2()

Dir_model[-1 - len(moloch_wspeed_help) + 1:] = list(moloch_dir_help)##replaces station wind dir with moloch wind dir for the existing Moloch model timestamps

Processed_data=np.dstack((np.vstack(Timestamp_m), np.vstack(Wind_model_up),np.vstack(Dir_model)))
csvfiletowrite1="C:/Users/plgeo/OneDrive/PC Desktop/FARIA/DATA EAA/Velanidia_synthetic_2years.csv"
with open (csvfiletowrite1,'w',encoding='UTF8', newline='') as file1:
     writer1=csv.writer(file1)
     for line in Processed_data:
         for l in line:
             writer1.writerow(l)
file1.close()

# std_measure=float(0.017)## anemometer measurement uncertainty  as an stdev in m/s
# d_wf=float(0.225) ## wind flow simulations uncertainty while considering wake losses in m/s
#
# station_pressure=station_data[:,1]
# station_temperature=station_data[:,2]
# station_RH=station_data[:,5]
# WTG_coord_help=list(WTG_coord)
# W1=station_wspeed_up
# W2=Wind_model_up
# Annual_mean_wind_speed_year1 = np.mean(Wind_model_up[0:17521]).astype('f')
# Annual_mean_wind_speed_year2 = np.mean(Wind_model_up[17521:]).astype('f')
# Annual_mean_model_speed=statistics.mean([Annual_mean_wind_speed_year1,Annual_mean_wind_speed_year2])
#
# Wind_WTG=np.zeros(shape=(len(W1),len(WTG_coord),len(H_hub)),dtype=float)
# Dens_WTG=np.zeros(shape=Wind_WTG.shape,dtype=float)
# Wind_100m=np.zeros(shape=(len(W1),len(WTG_coord)),dtype=float)
# Wind_uncertainty=np.zeros(shape=(len(W1),))
# for ii in range(len(WTG_coord_help)):
#     lat_wtg=WTG_coord[ii,0]
#     lon_wtg=WTG_coord[ii,1]
#     wtg_alt=WTG_alt[ii]
#     WeD1 = distlatlon(lat_station, lat_wtg, lon_station, lon_wtg)
#     WeD2 = distlatlon(lat_model, lat_wtg, lon_model, lon_wtg)
#     WeD1 = 1/(WeD1 ** 2)
#     WeD2 = 1/(WeD2 ** 2)
#     W_wtg=np.zeros(shape=(len(W1),1),dtype=float)
#     for ww in range(len(W1)):
#         w1=W1[ww]## upscaled station wind speed at corresp. timestamp
#         w2=W2[ww]## upscaled, mixture model wind speed at corresp. timestamp
#         # w_wtg = ((WeD1 * w1) + (WeD2 * w2)) / (WeD1 + WeD2)  ## the spatial interpolation using upscaled station data and model data is performed in all locations with WFs
#         w_wtg=w1
#         W_wtg[ww]=w_wtg
#     del w_wtg
#     for dd in range(W_wtg.size):
#         w_wtg = float(W_wtg[dd])
#         d1=float(Dir_model[dd])
#         mean_dclass=d1## just the mix model direction assigned to an old variable name( variable deleted) for convenience
#         if mean_dclass < 348.5 and (mean_dclass > 326 or mean_dclass == 326):
#             dsymb = "NNW"
#         elif mean_dclass < 326 and (mean_dclass > 303.5 or mean_dclass == 303.5):
#             dsymb = "NW"
#         elif mean_dclass < 303.5 and (mean_dclass > 281 or mean_dclass == 281):
#             dsymb = "WNW"
#         elif mean_dclass < 281 and (mean_dclass > 258.5 or mean_dclass == 258.5):
#             dsymb = "W"
#         elif mean_dclass < 258.5 and (mean_dclass > 236 or mean_dclass == 236):
#             dsymb = "WSW"
#         elif mean_dclass < 236 and (mean_dclass > 213.5 or mean_dclass == 213.5):
#             dsymb = "SW"
#         elif mean_dclass < 213.5 and (mean_dclass > 191 or mean_dclass == 191):
#             dsymb = "SSW"
#         elif mean_dclass < 191 and (mean_dclass > 168.5 or mean_dclass == 168.5):
#             dsymb = "S"
#         elif mean_dclass < 168.5 and (mean_dclass > 146 or mean_dclass == 146):
#             dsymb = "SSE"
#         elif mean_dclass < 146 and (mean_dclass > 123.5 or mean_dclass == 123.5):
#             dsymb = "SE"
#         elif mean_dclass < 123.5 and (mean_dclass > 101 or mean_dclass == 101):
#             dsymb = "ESE"
#         elif mean_dclass < 101 and (mean_dclass > 78.5 or mean_dclass == 78.5):
#             dsymb = "E"
#         elif mean_dclass < 78.5 and (mean_dclass > 56 or mean_dclass == 56):
#             dsymb = "ENE"
#         elif mean_dclass < 56 and (mean_dclass > 33.5 or mean_dclass == 33.5):
#             dsymb = "NE"
#         elif mean_dclass < 33.5 and (mean_dclass > 11 or mean_dclass == 11):
#             dsymb = "NNE"
#         elif mean_dclass < 11 or (mean_dclass > 348.5 or mean_dclass == 348.5):
#             dsymb = "N"
#         exist_count = Station_dclass.count(dsymb)  ##checks how many times the specific wind direction exists in the class of station directions created above , you may need to use "Dir_list" instead
#         if exist_count > 0:
#            index_d = Station_dclass.index(dsymb)  # returns the index of the 1st time station_dir exists in the Station_dclass list
#         else:
#            index_d = R_sq.index(max(R_sq))  # returns the index of the most well fitted direction class
#         # shear_fd = float(Ashear_model[index_d])
#         shear_fd = float(Ashear_oldmast[index_d])  ## wind shear derived from Chelona met mast.
#         ##modelling uncertainty per sector
#         s_ce = float(0.125)  ## vertical extrapolation uncertainty by using log law with annotated roughness length
#         s_sh = S_sherror[index_d]
#         s_unwind = s_ce + s_sh + std_measure + d_wf
#         Wind_uncertainty[dd] = s_unwind
#         P = float(station_pressure[dd])
#         P*=100 ##convert mb to Pa
#         T = float(station_temperature[dd])
#         T += 273.15 ##temperature in Kelvin
#         RH = float(station_RH[dd])
#         WDC1=[]
#         for ee in range(3):
#             h_hub=float(H_hub[ee])
#             wdc1= 0.5 / math.log(h_hub / zo)  # wake decay coefficient version 1
#             alt_hub = wtg_alt + h_hub ## altitude of the hub configuration
#             W_hubh = w_wtg * ((h_hub / 100) ** shear_fd) ## wind speed at hub height
#             # print(w_wtg)
#             # print(W_hubh)
#             Z = (Rad * st_alt) / (Rad + st_alt) ## geopotential altitude at weather station
#             Zhub = (Rad * alt_hub) / (Rad + alt_hub) ## geopotential altitude at WTG hubheight.
#             Thub = T - L * (Zhub - Z) ## temperature in Kelvin
#             Wind_WTG[dd, ii, ee] = W_hubh
#             Phub = P * (Thub / T) ** ((G * Md) / (Rconst * L)) ##%%atm pressure at WTG hub height.
#             f = 1.00070 + 3.113 * (10 ** (-8)) * Phub + 5.4 * (10 ** (-7)) * ((Thub - 273.15) ** 2)  #enhancement factor or ratio of effective saturation vapor pressure of water in moist air to saturation vapor pressure of water in moist air, empirical eq.
#             es = 1.7526 * (10 ** 11) * math.exp(-5315.56 / Thub) ##saturation vapor pressure of water , empirical eq.
#             phub = ((Phub * Md) / (Rconst * Thub * Zc)) * (1 - ((1 - (Mw / Md)) * (RH / 100) * f * es * (1 / Phub))) ## estimated air density at Phub
#             WDC1.append(wdc1)
#                    # print(phub)##checks if air density is calculated right
#             Dens_WTG[dd, ii, ee] = phub
#         Wind_100m[dd, ii] = w_wtg
#
# # print('The corresponding annual wind 2020 (mean, sdv) derived from each station is:', station_cl2020)##print them when comparison is going to take place
# # print('The corresponding mean annual wind 2021 (mean, sdv) derived from each station is:', station_cl2021)##print them when comparison is going to take place
# print('The annual wind climate for year 2020 , as derived from old masts is:', Wmean_year1)## this is mean climate at 100m agl
# print('The annual wind climate for year 2021 , as derived from old masts is:', Wmean_year2)## this is mean climate at 100m agl
# print('The annual wind climate for year 2020 , as derived from the mixture of models is:', Annual_mean_wind_speed_year1)## this is mean climate at 100m agl
# print('The annual wind climate for year 2021 , as derived from the mixture of models is:', Annual_mean_wind_speed_year2)## this is mean climate at 100m agl
# wind_climate=np.array([[Wmean_year1,Wmean_year2],[Annual_mean_wind_speed_year1,Annual_mean_wind_speed_year2]],ndmin=2)
# # # print(Dclass_perc_station)##checks dominant direction of the weather station
# filetosave=directory_to_save+"VELANIDIA_wake_model1_roughness_"+str(zo)+"_2022.xlsx"
# # filetosave=directory_to_save+"LOUTSA_wake_model1_without_masts_roughness_"+str(zo)+"_2022.xlsx"
# filetowrite = pd.ExcelWriter(filetosave, engine='xlsxwriter')##excel file to save processed data
# df_wind_climate=pd.DataFrame(wind_climate,index=["masts+station","EAA model+station"],columns=["wspeed @100m agl, 2020","wspeed @100m agl, 2021"])
# df_wind_climate.to_excel(filetowrite,sheet_name="annual wind speed",engine="xlsxwriter",startrow=0,startcol=0)
# # #
# #power output
#
# del WTG_coord
# del WTG_alt
# del WTG_coord_help
# del k
# del kk
# del dd
#
# D=170 ##rotor diameter in [m]
# RR=D/2 ## rotor radius in [m]
#
# ## filetosave='C:/Users/plgeo/OneDrive/Υπολογιστής/FARIA/Preliminary_wind_resource_assessment.xlsx'
#
# ## wind to power conversion - includes downstream wind speed derived from thrust curve and modified Jensen model
# directory_power_curve='C:/Users/plgeo/OneDrive/PC Desktop/FARIA/Wind turbine type/';
# file_curve=directory_power_curve+'SG_turbine_characteristics.xlsx'
# NUM_curve= pd.read_excel(file_curve,sheet_name='power curve')
# NUM_thrust=pd.read_excel(file_curve,sheet_name='Thrust coef')
#
# wind1=list(NUM_curve.iloc[1:,0])
# oper_dens=list(NUM_curve.iloc[0,1:])
# powerc1=[]
# Ct=[]
# pow_limit=float(5997)
# for y in range(len(wind1)):
#   powc1=list(NUM_curve.iloc[y+1,1:])
#   powc1=[pow_limit if powc>=6000 else powc for powc in powc1] #curtailment of power so that the nominal capacity is satisfied.
#   ct1=list(NUM_thrust.iloc[y+1,1:])
#   powerc1.append(powc1)
#   Ct.append(ct1)
# power_c1=np.array(powerc1,dtype=float).reshape(len(wind1),len(oper_dens))
# Ct=np.array(Ct,dtype=float).reshape(len(wind1),len(oper_dens))##thrust coeeficient curve.
# # print(power_c1)## prints the power curve matrix (rows= wind speed , cols = air density)
# electr_loss_factor=0.97
#
#
# WTG_coord = WTG_coord_Velanidia
# WTG_alt = Alt_Velanidia
# Density_wtg = Dens_WTG
# Wind_speed_wtg = Wind_WTG
# farm_name = location
# WTG_coord_help = list(WTG_coord)
#
# Distance = []
# for k in range(WTG_coord.shape[0] - 1):
#     lat_wtg0 = WTG_coord[k, 0]
#     lon_wtg0 = WTG_coord[k, 1]
#     lat_wtg1 = WTG_coord[k + 1, 0]
#     lon_wtg1 = WTG_coord[k + 1, 1]
#     dist_wtg = distlatlon(lat_wtg0, lat_wtg1, lon_wtg0, lon_wtg1)  ## finds distance between consecutive WTGs in km
#     Distance.append(dist_wtg)
# mean_wtg_distance = float(statistics.mean(Distance))  ## finds wind farm average distance between consecutive WTGs in km
#
# del lat_wtg0
# del lon_wtg0
# del lat_wtg1
# del lon_wtg1
#
# d_p=float(0.0004) # uncertainty in air density in kg/m^3
#
# ### power farm output without wake losses
# def Power_farm():
#     Power_wtg = np.zeros(shape=(Density_wtg.shape[0], Density_wtg.shape[1], Density_wtg.shape[2], 3), dtype=float)
#     for aa in range(Density_wtg.shape[0]):
#         d_ws = Wind_uncertainty[aa]
#         for bb in range(Density_wtg.shape[1]):
#             for cc in range(Density_wtg.shape[2]):
#                 densex = Density_wtg[aa, bb, cc]
#                 windex = Wind_speed_wtg[aa, bb, cc]
#                 densextr_help = np.array([densex - 1.68 * (d_p / math.sqrt(Density_wtg.shape[1] * Density_wtg.shape[2])), densex,densex + 1.68 * (d_p / math.sqrt(Density_wtg.shape[1] * Density_wtg.shape[2]))],dtype=float)  ##model air density uncertainty ,90% confidence interval
#                 windextr_help = np.array([windex - 1.68 * (d_ws / math.sqrt(Wind_speed_wtg.shape[1] * Wind_speed_wtg.shape[2])), windex,windex + 1.68 * (d_ws / math.sqrt(Wind_speed_wtg.shape[1] * Wind_speed_wtg.shape[2]))],dtype=float)  ##model wind speed uncertainty ,90% confidence interval
#                 ## extrapolation from wind speed and air density to power based on power curve of SG-6200 and conversion to half hour energy production
#                 for dd in range(densextr_help.size):
#                     windextr = windextr_help[dd]
#                     if windextr < 0:
#                        windextr = 0
#                     densextr = densextr_help[dd]
#                     indexw1 = [];
#                     indexw2 = [];
#                     indexw3 = []
#                     for kk in range(len(wind1)):
#                         w_c = wind1[kk]
#                         if w_c == windextr:
#                            indexw1.append(kk)
#                         elif w_c < windextr:
#                            indexw2.append(kk)
#                         elif w_c > windextr:
#                            indexw3.append(kk)
#                     indexd1 = [];
#                     indexd2 = [];
#                     indexd3 = []
#                     for ll in range(len(oper_dens)):
#                         d_c = oper_dens[ll]
#                         if d_c == densextr:
#                            indexd1.append(ll)
#                         elif d_c < densextr:
#                            indexd2.append(ll)
#                         elif d_c > densextr:
#                            indexd3.append(ll)
#                     if indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                        indexw1 = int(indexw1)
#                        indexd1 = int(indexd1)
#                        powerextr = powerc1[indexw1, indexd1]
#                     elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                        indexw1 = int(indexw1)
#                        indexds = int(indexd2[-1])
#                        indexdf = int(indexd3[0])
#                        ds = float(oper_dens[indexds])
#                        df = float(oper_dens[indexdf])
#                        pd_s = power_c1[indexw1, indexds]
#                        pd_f = power_c1[indexw1, indexdf]
#                        powerextr = pd_s + (pd_f - pd_s) * abs((densextr - ds) / (df - ds))
#                     elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
#                        indexw1 = int(indexw1)
#                        powerextr = power_c1[indexw1, 0]
#                     elif not indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                        indexd1 = int(indexd1)
#                        indexws = int(indexw2[-1])
#                        indexwf = int(indexw3[0])
#                        ws = float(wind1[indexws])
#                        wf = float(wind1[indexwf])
#                        ps = power_c1[indexws, indexd1]
#                        pf = power_c1[indexwf, indexd1]
#                        powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
#                     elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                        indexds = int(indexd2[-1])
#                        indexdf = int(indexd3[0])
#                        ds = float(oper_dens[indexds])
#                        df = float(oper_dens[indexdf])
#                        indexws = int(indexw2[-1])
#                        indexwf = int(indexw3[0])
#                        ws = float(wind1[indexws])
#                        wf = float(wind1[indexwf])
#                        pws1 = power_c1[indexws, indexds]
#                        pwf1 = power_c1[indexwf, indexds]
#                        pwextrs = pws1 + (pwf1 - pws1) * abs((windextr - ws) / (wf - ws))
#                        pws2 = power_c1[indexws, indexdf]
#                        pwf2 = power_c1[indexwf, indexdf]
#                        pwextrf = pws2 + (pwf2 - pws2) * abs((windextr - ws) / (wf - ws))
#                        powerextr = pwextrs + (pwextrf - pwextrs) * abs((densextr - ds) / (df - ds))
#                     elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
#                        indexws = int(indexw2[-1])
#                        indexwf = int(indexw3[0])
#                        ws = float(wind1[indexws])
#                        wf = float(wind1[indexwf])
#                        ps = power_c1[indexws, 0]
#                        pf = power_c1[indexwf, 0]
#                        powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
#                     elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[-1]) or densextr == float(oper_dens[-1])):
#                        indexws = int(indexw2[-1])
#                        indexwf = int(indexw3[0])
#                        ws = float(wind1[indexws])
#                        wf = float(wind1[indexwf])
#                        ps = power_c1[indexws, -1]
#                        pf = power_c1[indexwf, -1]
#                        powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
#                     else:
#                        powerextr = 0
#                     Power_wtg[aa, bb, cc, dd] = 0.5 * electr_loss_factor * (powerextr / 1000)  ## half hourly production in MWh
#     Annual_Power_wtg_year1 = Power_wtg[0:17521, :, :, :].sum(axis=0)
#     Annual_Power_wtg_year2 = Power_wtg[17521:, :, :, :].sum(axis=0)
#     AP1 = Annual_Power_wtg_year1.sum(axis=0)
#     AP2 = Annual_Power_wtg_year2.sum(axis=0)
#     return AP1, AP2
#
# Annual_Power_farm_year1, Annual_Power_farm_year2 = Power_farm()
# print(Annual_Power_farm_year1)
# print("\n")
#
# ## considering wake losses
# ##wake effect simulations
# mean_wtg_distance *= 1000
# distance_factor = mean_wtg_distance / D
# mix_coef=1-(1/distance_factor) #mixing coefficient -- a , used in EB or SSQ model for wake overlapping
# WDC1=statistics.mean(WDC1)
# WDC3=0.5*statistics.mean(TI_oldmast)#wake decay coefficient version 2 , this version seems to overestimate , while version 1 seems to underestimate compared to the constant value used in WaSP for onshore WFs, used under neutral stability conditions
# WDC2=0.075#wake decay coefficient version 3
# WDC=WDC1## switch between WDC3, WDC2 and WDC1, if it is WDC1 it is 3 values, each for each different hub height config.
# # print(WDC1);print(WDC2);print(WDC3)
# Ao=math.pi*(RR**2)## rotor swept area
#
# Power_wtg_we = np.zeros(shape=(Density_wtg.shape[0], Density_wtg.shape[1], Density_wtg.shape[2], 3), dtype=float)
# Angle_thres=[];Angle_dif=[];Lateral_dist=[];Check1=[];Check2=[];Theta1=[];Theta2=[];Overlap_factor=[]
# for aa in range(Density_wtg.shape[0]):
#     d_ws = Wind_uncertainty[aa]
#     wind_dir=float(Dir_model[aa])##computed wind direction from combo of station and model( mast cannot be used yet)
#     if wind_dir>=0.00 and wind_dir<270.00:
#        theta=wind_dir+90.00
#     elif wind_dir>=270.00 and wind_dir<=360:
#        theta=90.00-abs(360-wind_dir)
#     theta=math.radians(theta)
#     for bb in range(Density_wtg.shape[1]):
#            lat_wtgt = WTG_coord[bb, 0]
#            lon_wtgt = WTG_coord[bb, 1]
#            x_wtgt,y_wtgt=latlon_to_xy(lat_wtgt,lon_wtgt)
#            wake_red_array=np.empty(shape=(Density_wtg.shape[1],Density_wtg.shape[2]))##array of wake influence of target turbine from upwind turbines
#            for ee in range(Density_wtg.shape[1]):
#                if ee!=bb:
#                   lat_wtgup=WTG_coord[ee,0]
#                   lon_wtgup=WTG_coord[ee,1]
#                   x_wtgup,y_wtgup=latlon_to_xy(lat_wtgup,lon_wtgup)
#                   dist_up_down=math.sqrt(((x_wtgup-x_wtgt)**2)+((y_wtgup-y_wtgt)**2))## a better way to compute geographical distance because of calculation error in converting coordinates
#                   # dist_up_down=distlatlon(lat_wtgup,lat_wtgt,lon_wtgup,lon_wtgt)
#                   # dist_up_down=dist_up_down*1000 ## geographical distance between turbines in m
#                   ##determine cone of influence for upwind turbine wrt target turbine
#                   comp1=(x_wtgup-x_wtgt)*math.cos(theta)+(y_wtgup-y_wtgt)*math.sin(theta)+(RR/WDC)
#                   comp2=x_wtgup-x_wtgt+(RR/WDC)*math.cos(theta)
#                   comp3=y_wtgup-y_wtgt+(RR/WDC)*math.sin(theta)
#                   cone_of_infl=math.acos(comp1/math.sqrt((comp2**2)+(comp3**2)))## angle corresponding to the cone of wake influence of upwind to downwind turbine
#                   cone_of_infl=cone_of_infl*(180/math.pi)
#                   ##determination of the downwind (downstream)distance
#                   if theta>=0 and theta<(0.5*math.pi):
#                      compx = (x_wtgup - x_wtgt) * math.cos(theta)
#                      compy = (y_wtgup - y_wtgt) * math.sin(theta)
#                      d_down = math.sqrt((compx ** 2) + (compy ** 2))  ## downwind distance at wind_dir
#                   elif theta>=(0.5*math.pi) and theta<math.pi:
#                      phi=math.asin(abs(x_wtgup-x_wtgt)/dist_up_down)
#                      phi1=theta-phi-(0.5*math.pi)
#                      d_down=dist_up_down*math.cos(phi1)
#                   elif theta>=math.pi and theta<(1.5*math.pi):
#                      phi = math.asin(abs(x_wtgup - x_wtgt) / dist_up_down)
#                      phi1=(1.5*math.pi)-theta-phi
#                      d_down=dist_up_down*math.cos(phi1)
#                   else:
#                      phi = math.asin(abs(y_wtgup - y_wtgt) / dist_up_down)
#                      phi1=(2*math.pi)-theta-phi
#                      d_down=dist_up_down*math.cos(phi1)
#                   Rwake = RR + (WDC * d_down)  ## wake radius -- variable , function of wake decay constant and downwind distance
#                   angle_thres=math.atan(D/d_down)## calculates threshold of the cone of influence per wind direction in radians -- fradsen proposal on how to calculate threshold angle of wake cone.
#                   angle_thres=0.5*(angle_thres*(180/math.pi)+10)
#                   # angle_thres=math.atan((0.5*Rwake)/(RR/WDC))## alternative formulation of threshold angle based on personal geometrical interpretation
#                   # angle_thres=angle_thres*(180/math.pi)
#                   Angle_thres.append(angle_thres)
#                   ## if turbine bb is in the cone of wake influence of upwind turbine ee ( cone angle lower than or equal to angle threshold ,keep downwind distance at particular wind direction
#                   if cone_of_infl>0 and cone_of_infl<=angle_thres:
#                      D_down=d_down
#                      di_j = math.sqrt((dist_up_down ** 2) - (D_down ** 2))  ## this is the lateral component of the 2 turbines' geographical distance wrt wind direction
#                      Lateral_dist.append(di_j)
#                      # ## geometrical analytical computation of the wake overlapping area
#                      if di_j>=(Rwake-RR) and di_j<=(Rwake+RR):
#                         xdd1=((di_j ** 2) + (Rwake ** 2) - (RR ** 2)) / (2 * di_j)
#                         xdd2 = ((di_j ** 2) + (RR ** 2) - (Rwake ** 2)) / (2 * di_j)
#                         Check1.append(xdd1/Rwake)
#                         Check2.append(xdd2/RR)
#                         theta1=math.acos(xdd1/Rwake)
#                         theta2=math.acos(xdd2/RR)
#                         Theta1.append(math.degrees(theta1))
#                         Theta2.append(math.degrees(theta2))
#                         Atriangle1=0.5*(Rwake**2)*math.sin(theta1)
#                         Atriangle2=0.5*(RR**2)*math.sin(theta2)
#                         Acone1=0.5*(Rwake**2)*theta1
#                         Acone2=0.5*(RR**2)*theta2
#                         Ashadow =2*(Acone1+Acone2-Atriangle1-Atriangle2) ## overlapping area of wake effect of each upwind turbine - analytical solution assuming steady state
#                         # # Ashadow=0.5*(Rwake**2)*math.acos(0.5*di_j/Rwake)+0.5*(RR**2)*math.acos(di_j/D)-(0.25*di_j)*(math.sqrt((2*Rwake)**2-di_j**2)+math.sqrt(D**2-di_j**2))
#                         overlap_f = Ashadow / Ao  ## factor corresponding to overlapping area of wake effect of each upwind turbine to each downwind turbine
#                         print('WTG ' + str(bb + 1) + ' is in the wake influence of WTG ' + str(ee + 1))
#                      elif di_j<(Rwake-RR):
#                         xdd1=2*math.pi*Rwake
#                         xdd2=2*math.pi*RR
#                         Ashadow=Ao
#                         overlap_f=1
#                         print('WTG ' + str(bb + 1) + ' is in the wake influence of WTG ' + str(ee + 1))
#                      else:
#                         xdd1=0
#                         xdd2=0
#                         Ashadow=0
#                         overlap_f=0
#                         print('WTG ' + str(bb + 1) + ' is not in the wake influence of WTG ' + str(ee + 1))
#                   elif cone_of_infl==0:
#                      D_down=d_down
#                      Ashadow=Ao
#                      overlap_f=1
#                      print('WTG ' + str(bb + 1) + ' is in the wake influence of WTG ' + str(ee + 1))
#                   else:
#                      D_down=0
#                      Ashadow=0
#                      overlap_f=0
#                      print('WTG '+str(bb+1)+' is not in the wake influence of WTG '+str(ee+1))
#                else:
#                   D_down=0## this is the case of the upwind turbine and target turbine being the same.
#                   Ashadow=0
#                   overlap_f=0
#                   print('WTG ' + str(bb + 1) + ' is the same as OR they belong to different sites, WTG ' + str(ee + 1))
#                Overlap_factor.append(overlap_f)
#                for ff in range(Density_wtg.shape[2]):
#                    windextr = Wind_speed_wtg[aa, ee, ff]## upwind turbine wind velocity
#                    densextr=  Density_wtg[aa, ee, ff]
#                    indexw1 = [];
#                    indexw2 = [];
#                    indexw3 = []
#                    for kk in range(len(wind1)):
#                        w_c = wind1[kk]
#                        if w_c == windextr:
#                           indexw1.append(kk)
#                        elif w_c < windextr:
#                           indexw2.append(kk)
#                        elif w_c > windextr:
#                           indexw3.append(kk)
#                    indexd1 = [];
#                    indexd2 = [];
#                    indexd3 = []
#                    for ll in range(len(oper_dens)):
#                        d_c = oper_dens[ll]
#                        if d_c == densextr:
#                           indexd1.append(ll)
#                        elif d_c < densextr:
#                           indexd2.append(ll)
#                        elif d_c > densextr:
#                           indexd3.append(ll)
#                    ##determination of thrust coefficient from upwind velocity
#                    if indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                       indexw1 = int(indexw1)
#                       indexd1 = int(indexd1)
#                       thrustextr = Ct[indexw1, indexd1]
#                    elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                      indexw1 = int(indexw1)
#                      indexds = int(indexd2[-1])
#                      indexdf = int(indexd3[0])
#                      ds = float(oper_dens[indexds])
#                      df = float(oper_dens[indexdf])
#                      cd_s = Ct[indexw1, indexds]
#                      cd_f = Ct[indexw1, indexdf]
#                      thrustextr = cd_s + (cd_f - cd_s) * abs((densextr - ds) / (df - ds))
#                    elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
#                      indexw1 = int(indexw1)
#                      thrustextr = Ct[indexw1, 0]
#                    elif not indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                      indexd1 = int(indexd1)
#                      indexws = int(indexw2[-1])
#                      indexwf = int(indexw3[0])
#                      ws = float(wind1[indexws])
#                      wf = float(wind1[indexwf])
#                      cs = Ct[indexwf, indexd1]
#                      cf = Ct[indexwf, indexd1]
#                      thrustextr = cs + (cf - cs) * abs((windextr - ws) / (wf - ws))
#                    elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                      indexds = int(indexd2[-1])
#                      indexdf = int(indexd3[0])
#                      ds = float(oper_dens[indexds])
#                      df = float(oper_dens[indexdf])
#                      indexws = int(indexw2[-1])
#                      indexwf = int(indexw3[0])
#                      ws = float(wind1[indexws])
#                      wf = float(wind1[indexwf])
#                      cws1 = Ct[indexws, indexds]
#                      cwf1 = Ct[indexwf, indexds]
#                      cwextrs = cws1 + (cwf1 - cws1) * abs((windextr - ws) / (wf - ws))
#                      cws2 = Ct[indexws, indexdf]
#                      cwf2 = Ct[indexwf, indexdf]
#                      cwextrf = cws2 + (cwf2 - cws2) * abs((windextr - ws) / (wf - ws))
#                      thrustextr = cwextrs + (cwextrf - cwextrs) * abs((densextr - ds) / (df - ds))
#                    elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
#                      indexws = int(indexw2[-1])
#                      indexwf = int(indexw3[0])
#                      ws = float(wind1[indexws])
#                      wf = float(wind1[indexwf])
#                      cs = Ct[indexws, 0]
#                      cf = Ct[indexwf, 0]
#                      thrustextr = cs + (cf - cs) * abs((windextr - ws) / (wf - ws))
#                    elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[-1]) or densextr == float(oper_dens[-1])):
#                      indexws = int(indexw2[-1])
#                      indexwf = int(indexw3[0])
#                      ws = float(wind1[indexws])
#                      wf = float(wind1[indexwf])
#                      cs = Ct[indexws, -1]
#                      cf = Ct[indexwf, -1]
#                      thrustextr = cs + (cf - cs) * abs((windextr - ws) / (wf - ws))
#                    else:
#                      thrustextr = Ct[0,0]
#                    del indexw1
#                    del indexw2
#                    del indexw3
#                    del indexd1
#                    del indexd2
#                    del indexd3
#                    wind_up=windextr
#                    del windextr
#                    if D_down>0:
#                       w_red_f=(1-math.sqrt(1-thrustextr))*((RR/(RR+(WDC*D_down)))**2)*overlap_f ## wind reduction factor of target turbine from upwind turbine , due to wake effect, as derived from modified Jensen model.
#                    else:
#                        w_red_f=0
#                    wake_red_array[ee,ff]=w_red_f**2
#            wind_red_factor=wake_red_array.sum(axis=0)
#            # print(np.sqrt(wind_red_factor))
#            for cc in range(Density_wtg.shape[2]):
#                wi_red_f=wind_red_factor[cc]
#                wind_target=wind_up*(1-math.sqrt(wi_red_f))## wind speed of target turbine as derived from upwind turbine and wind reduction factor due to wake effect
#                # wind_target = wind_up * (1 - math.sqrt(mix_coef) * math.sqrt(wi_red_f))  ## wind speed of target turbine as derived from upwind turbine and wind reduction factor due to wake
#                windex=wind_target
#                # print(windex)
#                densex = Density_wtg[aa, bb, cc]
#                densextr_help = np.array([densex - 1.68 * (d_p / math.sqrt(Density_wtg.shape[1] * Density_wtg.shape[2])), densex,densex + 1.68 * (d_p / math.sqrt(Density_wtg.shape[1] * Density_wtg.shape[2]))],dtype=float)  ##model air density uncertainty ,90% confidence interval
#                windextr_help = np.array([windex - 1.68 * (d_ws / math.sqrt(Wind_speed_wtg.shape[1] * Wind_speed_wtg.shape[2])), windex,windex + 1.68 * (d_ws / math.sqrt(Wind_speed_wtg.shape[1] * Wind_speed_wtg.shape[2]))],dtype=float)  ##model wind speed uncertainty ,90% confidence interval
#                for dd in range(densextr_help.size):
#                    windextr = windextr_help[dd]
#                    if windextr < 0:
#                        windextr = 0
#                    densextr = densextr_help[dd]
#                    indexw1 = [];
#                    indexw2 = [];
#                    indexw3 = []
#                    for kk in range(len(wind1)):
#                        w_c = wind1[kk]
#                        if w_c == windextr:
#                           indexw1.append(kk)
#                        elif w_c < windextr:
#                           indexw2.append(kk)
#                        elif w_c > windextr:
#                           indexw3.append(kk)
#                    indexd1 = [];
#                    indexd2 = [];
#                    indexd3 = []
#                    for ll in range(len(oper_dens)):
#                        d_c = oper_dens[ll]
#                        if d_c == densextr:
#                           indexd1.append(ll)
#                        elif d_c < densextr:
#                           indexd2.append(ll)
#                        elif d_c > densextr:
#                           indexd3.append(ll)
#
#                 ## extrapolation from wind speed and air density to power based on power curve of SG-6200 and conversion to half hour energy production
#                    if indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                       indexw1 = int(indexw1)
#                       indexd1 = int(indexd1)
#                       powerextr = powerc1[indexw1, indexd1]
#                    elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                       indexw1 = int(indexw1)
#                       indexds = int(indexd2[-1])
#                       indexdf = int(indexd3[0])
#                       ds = float(oper_dens[indexds])
#                       df = float(oper_dens[indexdf])
#                       pd_s = power_c1[indexw1, indexds]
#                       pd_f = power_c1[indexw1, indexdf]
#                       powerextr = pd_s + (pd_f - pd_s) * abs((densextr - ds) / (df - ds))
#                    elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
#                       indexw1 = int(indexw1)
#                       powerextr = power_c1[indexw1, 0]
#                    elif not indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                       indexd1 = int(indexd1)
#                       indexws = int(indexw2[-1])
#                       indexwf = int(indexw3[0])
#                       ws = float(wind1[indexws])
#                       wf = float(wind1[indexwf])
#                       ps = power_c1[indexws, indexd1]
#                       pf = power_c1[indexwf, indexd1]
#                       powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
#                    elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
#                       indexds = int(indexd2[-1])
#                       indexdf = int(indexd3[0])
#                       ds = float(oper_dens[indexds])
#                       df = float(oper_dens[indexdf])
#                       indexws = int(indexw2[-1])
#                       indexwf = int(indexw3[0])
#                       ws = float(wind1[indexws])
#                       wf = float(wind1[indexwf])
#                       pws1 = power_c1[indexws, indexds]
#                       pwf1 = power_c1[indexwf, indexds]
#                       pwextrs = pws1 + (pwf1 - pws1) * abs((windextr - ws) / (wf - ws))
#                       pws2 = power_c1[indexws, indexdf]
#                       pwf2 = power_c1[indexwf, indexdf]
#                       pwextrf = pws2 + (pwf2 - pws2) * abs((windextr - ws) / (wf - ws))
#                       powerextr = pwextrs + (pwextrf - pwextrs) * abs((densextr - ds) / (df - ds))
#                    elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
#                       indexws = int(indexw2[-1])
#                       indexwf = int(indexw3[0])
#                       ws = float(wind1[indexws])
#                       wf = float(wind1[indexwf])
#                       ps = power_c1[indexws, 0]
#                       pf = power_c1[indexwf, 0]
#                       powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
#                    elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[-1]) or densextr == float(oper_dens[-1])):
#                       indexws = int(indexw2[-1])
#                       indexwf = int(indexw3[0])
#                       ws = float(wind1[indexws])
#                       wf = float(wind1[indexwf])
#                       ps = power_c1[indexws, -1]
#                       pf = power_c1[indexwf, -1]
#                       powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
#                    else:
#                       powerextr = 0
#                    Power_wtg_we[aa, bb, cc, dd] = 0.5 * electr_loss_factor * (powerextr / 1000)  ## half hourly production in MWh including wake effect
# Annual_Power_wtg_with_losses_year1=np.zeros(shape=(len(WTG_coord),Density_wtg.shape[2],densextr_help.size),dtype=float)
# Annual_Power_wtg_with_losses_year2 = np.zeros(shape=(len(WTG_coord), Density_wtg.shape[2], densextr_help.size), dtype=float)
# Annual_Power_wtg_with_losses_year1 = Power_wtg_we[0:17521, :, :, :].sum(axis=0)
# Annual_Power_wtg_with_losses_year2 = Power_wtg_we[17521:, :, :, :].sum(axis=0)
# Annual_PF_with_losses_year1=np.zeros(shape= (Density_wtg.shape[2], densextr_help.size), dtype=float)
# Annual_PF_with_losses_year2 = np.zeros(shape=(Density_wtg.shape[2], densextr_help.size), dtype=float)
# Annual_PF_with_losses_year1 = Annual_Power_wtg_with_losses_year1.sum(axis=0)
# Annual_PF_with_losses_year2 = Annual_Power_wtg_with_losses_year2.sum(axis=0)
# WDC=format(WDC,".3f")
# print('\n')
# print('WF ' + farm_name + ': ' + 'The array of 9 scenarios of APE of wake effect with WDC= ' + str(WDC) + ' for year 2020 in MWh is:\n',Annual_PF_with_losses_year1)  ## shows a matrix of 9 scenarios of annual APE including wake and electrical losses
# print('\n')
# print('WF ' + farm_name + ': ' + 'The array of 9 scenarios of APE of wake effect with WDC= ' + str(WDC) + ' for year 2021 in MWh is:\n',Annual_PF_with_losses_year2)
# print('\n')
#
# Wake_loss_year1 = np.zeros(shape=Annual_Power_farm_year1.shape, dtype=float)
# Wake_loss_year2 = np.zeros(shape=Annual_Power_farm_year2.shape, dtype=float)
# for sc in range(Annual_Power_farm_year1.shape[0]):
#     for ci in range(Annual_Power_farm_year1.shape[1]):
#         PF_no_loss1=Annual_Power_farm_year1[sc,ci]
#         PF_no_loss2=Annual_Power_farm_year2[sc,ci]
#         PF_loss1=Annual_PF_with_losses_year1[sc,ci]
#         PF_loss2=Annual_PF_with_losses_year2[sc,ci]
#         wake_loss_year1=(PF_no_loss1-PF_loss1)/PF_loss1
#         wake_loss_year2=(PF_no_loss2-PF_loss2)/PF_loss2
#         Wake_loss_year1[sc,ci]=wake_loss_year1
#         Wake_loss_year2[sc,ci]=wake_loss_year2
#
# df_2020=pd.DataFrame(Annual_PF_with_losses_year1,columns=["P95","P50","P5"])
# df_2021 = pd.DataFrame(Annual_PF_with_losses_year2, columns=["P95", "P50", "P5"])
# df_2020.to_excel(filetowrite,sheet_name="WDC="+str(WDC),startrow=0,startcol=0,index=False)
# df_2021.to_excel(filetowrite,sheet_name="WDC="+str(WDC),startrow=4,startcol=0,index=False)
# df_wl_2020=pd.DataFrame(Wake_loss_year1,columns=["P95","P50","P5"])
# df_wl_2021=pd.DataFrame(Wake_loss_year2,columns=["P95","P50","P5"])
# df_wl_2020.to_excel(filetowrite,sheet_name="wake loss @ zo="+str(zo)+", "+"WDC="+str(WDC),startrow=0,startcol=0,index=False)
# df_wl_2021.to_excel(filetowrite,sheet_name="wake loss @ zo="+str(zo)+", "+"WDC="+str(WDC),startrow=4,startcol=0,index=False)
# filetowrite.save()
# # # # #
# # # ##lines that check validity of geometrical calculations. (e.g. the def. set of sine and cosine)
# # # # print(min(Angle_thres))
# # # # print(max(Angle_thres))
# # # # print(min(Lateral_dist))
# # # # print(statistics.mean(Lateral_dist))
# # # # print(max(Lateral_dist))
# # # # print(min(Check2))
# # # # print(statistics.mean(Check2))
# # # # print(max(Check2))
# # # # print(min(Theta1))
# # # # print(max(Theta1))
# # # # print(min(Theta2))
# # # # print(max(Theta2))
# # # # print(min(Overlap_factor))
# # # # print(statistics.mean(Overlap_factor))
# # # # print(max(Overlap_factor))
# print(R_sq)
# print(S1)
# print(statistics.mean(Ashear_model))
# uncertainty_wind=statistics.mean(S_sherror)+s_ce+d_ws+std_measure
# print(uncertainty_wind)