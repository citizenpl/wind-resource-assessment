
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

###Topographic info
##Altitude above sea level in m (altitude of the position +  hub
##height  %%for the moloch model it is only the altitude above sea level (add 50, 80 and 100m)

## issued according to point_coord names
##mesoscale model points altitude and height agl

Alt_point=float(10)
Hmodel=np.array([50,80,100])

##WTG altitude and positions.
Alt_Iasmos=np.array([16,15,14,12,12,12])
WTG_coord_Iasmos=np.array([[41.118311,25.215096],[41.109337,25.214745],[41.104766,25.216400],[41.102801,25.221772],[41.097104,25.227500],[41.101425,25.235499]])

##Hub Height configurations in m.
H_hub=np.array([135,165,170])

##weather station altitude in m  -- given with the same sequence
##as the list station_name

H_wstation=float(10)
Agl_wstation=float(5)

##Read and write mesoscale and station data

date1=datetime(2020,8,4)
date2=datetime(2021,12,31)
date_generated=pd.date_range(start=date1 , end=date2)
# print(date_generated.strftime("%y%m%d"))

directory="C:/Users/
directory_stations= "C:/Users/weather_station_data"
directory_to_save="C:/Users/plgeo/OneDrive/PC Desktop/"

EAA_file1=directory+"/"+"points_EAA.xlsx"
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

point=point_coords[5]
location=location_name[5]
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

directory_stations="C:/Users/plgeo/OneDrive/PC Desktop/weather_station_data"
station_name=["igoumenitsa","stratoni","asprovalta","imeros"]
station_coord=[[39.541749,20.279882]];station_coord.append([40.515033,23.828907]);station_coord.append([40.724931,23.711809]);station_coord.append([40.955700,25.369020])

name = station_name[3]
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

##wind mast altitude info and location
####Vegetation: sparce trees / terrain: cray soil
Alt_mast=float(10)# altitude in m.
H_mast=float(80)## height agl in m. , main measurement , 'C4'
H_mast_compl=np.array([75.8,60.4,59.1],dtype='f')## height agl in m. , complementary measurements ('C3','C2','C1'), Α2 is wind direction and A7 is temperature.
WTG_mast=np.array([41.098639,25.208333],dtype='f')# location of met mast in decimal degs.

##read and process wind mast raw data

directory_mast='C:/Users/PC Desktop/ new mast data/'

##find and use all folders containing mast raw data
folderlabel='_RAWDATA_'
folderstamps=os.listdir(directory_mast)
folderstamps=[fold for fold in folderstamps if fold[0:len(folderlabel)]==folderlabel]

Timestamp_mast=[];Datestamp_mast=[];W1_av=[];W1_sdv=[];W1_max=[];W2_av=[];W2_sdv=[];W2_max=[];W3_av=[];W3_sdv=[];W3_max=[];W4_av=[];W4_sdv=[];W4_max=[];Dir_mast=[];Temp_mast=[]
for  f in range(len(folderstamps)):
     folderstamp=folderstamps[f]
     folder=directory_mast+folderstamp
     dir_list=os.listdir(folder)
     for file in dir_list:
         directory_file=os.path.join(folder,file)
         label=file.split('.')
         recogn=label[1]##αναγνωριστικό αρχείου τύπου '000' , τα '001' δε θα ανοιγονται
         if recogn=='000':
            filecontent=open(directory_file,"r")
            ##check the raw files to see how many lines to skip till the actual data, by finding in which line the sentence 'All vane values are referenced to North' is
            ##start writing from 2 lines after this line for the first file( in order to have the line with the variable names) and from 3 lines after this line for the rest of the files
            datalines_pre=filecontent.readlines()
            filecontent.close()
            print('File is opened')
            del filecontent
            for line_no in range(len(datalines_pre)):
                line_pre=datalines_pre[line_no]
                parse=line_pre.strip()
                exp_to_check='All vane values are referenced to North.'
                if parse==exp_to_check:
                   index_of_lines_to_skip=line_no
            filecontent = open(directory_file, "r")
            datalines=filecontent.readlines()[index_of_lines_to_skip+2:]
            filecontent.close()
            varnames=datalines[0].split()
            if f==0 and file==dir_list[0]:
              pos_W1av=[varnames.index(pos) for pos in varnames if pos=='C1-av']
              pos_W1sdv=[varnames.index(pos) for pos in varnames if pos=='C1-sdv']
              pos_W1max=[varnames.index(pos) for pos in varnames if pos=='C1-max']
              pos_W2av=[varnames.index(pos) for pos in varnames if pos=='C2-av']
              pos_W2sdv=[varnames.index(pos) for pos in varnames if pos=='C2-sdv']
              pos_W2max = [varnames.index(pos) for pos in varnames if pos == 'C2-max']
              pos_W3av=[varnames.index(pos) for pos in varnames if pos=='C3-av']
              pos_W3sdv=[varnames.index(pos) for pos in varnames if pos=='C3-sdv']
              pos_W3max = [varnames.index(pos) for pos in varnames if pos == 'C3-max']
              pos_W4av=[varnames.index(pos) for pos in varnames if pos=='C4-av']
              pos_W4sdv=[varnames.index(pos) for pos in varnames if pos=='C4-sdv']
              pos_W4max = [varnames.index(pos) for pos in varnames if pos == 'C4-max']
              pos_Dir=[varnames.index(pos) for pos in varnames if pos=='A2-av']
              pos_Temp=[varnames.index(pos) for pos in varnames if pos=='A7-av']
              print(pos_W1av)
            List_of_var_pos=[pos_W1av[0]+2,pos_W1sdv[0]+2,pos_W1max[0]+2,pos_W2av[0]+2,pos_W2sdv[0]+2,pos_W2max[0]+2,pos_W3av[0]+2,pos_W3sdv[0]+2,pos_W3max[0]+2,pos_W4av[0]+2,pos_W4sdv[0]+2,pos_W4max[0]+2,pos_Dir[0]+2,pos_Temp[0]+2]
            datalines=datalines[1:]
            for line in datalines:
                raw=line.split()
                datestamp=raw[0]
                timestamp=raw[0]+" "+raw[1]
                w1_av=raw[List_of_var_pos[0]]
                w1_sdv=raw[List_of_var_pos[1]]
                w1_max=raw[List_of_var_pos[2]]
                w2_av=raw[List_of_var_pos[3]]
                w2_sdv=raw[List_of_var_pos[4]]
                w2_max = raw[List_of_var_pos[5]]
                w3_av=raw[List_of_var_pos[6]]
                w3_sdv=raw[List_of_var_pos[7]]
                w3_max = raw[List_of_var_pos[8]]
                w4_av=raw[List_of_var_pos[9]]
                w4_sdv=raw[List_of_var_pos[10]]
                w4_max = raw[List_of_var_pos[11]]
                dir=raw[List_of_var_pos[12]]
                temp=raw[List_of_var_pos[13]]
                W1_av.append(w1_av);W1_sdv.append(w1_sdv);W1_max.append(w1_max);W2_av.append(w2_av);W2_sdv.append(w2_sdv);W2_max.append(w2_max)
                W3_av.append(w3_av);W3_sdv.append(w3_sdv);W3_max.append(w3_max);W4_av.append(w4_av);W4_sdv.append(w4_sdv);W4_max.append(w3_max)
                Dir_mast.append(dir);Temp_mast.append(temp);Timestamp_mast.append(timestamp);Datestamp_mast.append(datestamp)
         else:
             print('File should not be opened')
Mast_variables=np.dstack((np.array(Timestamp_mast),np.array(W1_av),np.array(W2_av),np.array(W3_av),np.array(W4_av),np.array(W1_sdv),np.array(W2_sdv),np.array(W3_sdv),np.array(W4_sdv),np.array(W1_max),np.array(W2_max),np.array(W3_max),np.array(W4_max),np.array(Dir_mast),np.array(Temp_mast)))

##wind mast processing

##derivation of wind climate per month and overall
Period_label=[]
for dateframe in Datestamp_mast:
    month_no=dateframe[2:4]
    year_no=dateframe[4:]
    monthtime_object = datetime.strptime(month_no, "%m")
    month_word = monthtime_object.strftime("%B")#finds the name of the month
    period_label=month_word+" "+year_no
    Period_label.append(period_label)

def unique(list1):
    # insert the list to the set
    list_set = set(list1)
    # convert the set to the list
    unique_list = list(list_set)
    return unique_list

List_of_periods=unique(Period_label)## list that returns the period corresponding to each month.
del month_no
del year_no
del month_word
del monthtime_object

h1=H_mast_compl[2]
h2=H_mast_compl[1]
h3=H_mast_compl[0]
h4=H_mast
Height_list=[str(h1)+"m",str(h2)+"m",str(h3)+"m",str(h4)+"m"]
Dir_list=["N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW"]
for jj in range(len(List_of_periods)):
    period_of_report=List_of_periods[jj]
    w1mean=0;w2mean=0;w3mean=0;w4mean=0;Wg1=[];Wg2=[];Wg3=[];Wg4=[];Ws1=[];Ws2=[];Ws3=[];Ws4=[];Temper_mast=[];TEMPER_mast=[];Ti1=[];TI1=[];Ti2=[];TI2=[];Ti3=[];TI3=[];Ti4=[];TI4=[]
    count1=0
    W1mean=0;W2mean=0;W3mean=0;W4mean=0;WG1=[];WG2=[];WG3=[];WG4=[];WS1=[];WS2=[];WS3=[];WS4=[]
    count=0
    count_d1 = 0;count_d2 = 0;count_d3 = 0;count_d4 = 0;count_d5 = 0;count_d6 = 0;count_d7 = 0;count_d8 = 0;
    count_d9 = 0;count_d10 = 0;count_d11 = 0;count_d12 = 0;count_d13 = 0;count_d14 = 0;count_d15 = 0;count_d16 = 0
    w_d1 = 0;w_d2 = 0;w_d3 = 0;w_d4 = 0;w_d5 = 0;w_d6 = 0;w_d7 = 0;w_d8 = 0;
    w_d9 = 0;w_d10 = 0;w_d11 = 0;w_d12 = 0;w_d13 = 0;w_d14 = 0;w_d15 = 0;w_d16 = 0
    countd1=0;countd2=0;countd3=0;countd4=0;countd5=0;countd6=0;countd7=0;countd8=0;
    countd9 = 0;countd10 = 0;countd11 = 0;countd12 = 0;countd13 = 0;countd14 = 0;countd15 = 0;countd16 = 0
    Dir_cluster=[];DIR_cluster=[]
    for kk in range(len(Datestamp_mast)):
        datestamp_mast=Datestamp_mast[kk]
        month_no = datestamp_mast[2:4]
        year_no = datestamp_mast[4:]
        monthtime_object = datetime.strptime(month_no, "%m")
        month_word = monthtime_object.strftime("%B")  # finds the name of the month
        period_label = month_word + " " + year_no
        w1=float(W1_av[kk]);w2=float(W2_av[kk]);w3=float(W3_av[kk]);w4=float(W4_av[kk])
        wg1=float(W1_max[kk]);wg2=float(W2_max[kk]);wg3=float(W3_max[kk]);wg4=float(W4_max[kk])
        ws1=float(W1_sdv[kk]);ws2=float(W2_sdv[kk]);ws3=float(W3_sdv[kk]);ws4=float(W4_sdv[kk])
        dir_mast=float(Dir_mast[kk]);temp_mast=float(Temp_mast[kk])
        if dir_mast < 348.5 and (dir_mast > 326 or dir_mast == 326):
           dmast_symbol = "NNW"
           w_d1+=w4
           count_d1+=1
        elif dir_mast < 326 and (dir_mast > 303.5 or dir_mast == 303.5):
           dmast_symbol = "NW"
           w_d2 += w4
           count_d2 += 1
        elif dir_mast < 303.5 and (dir_mast > 281 or dir_mast == 281):
           dmast_symbol = "WNW"
           w_d3 += w4
           count_d3 += 1
        elif dir_mast < 281 and (dir_mast > 258.5 or dir_mast == 258.5):
           dmast_symbol = "W"
           w_d4 += w4
           count_d4 += 1
        elif dir_mast < 258.5 and (dir_mast> 236 or dir_mast == 236):
           dmast_symbol = "WSW"
           w_d5+= w4
           count_d5 += 1
        elif dir_mast < 236 and (dir_mast > 213.5 or dir_mast == 213.5):
           dmast_symbol = "SW"
           w_d6 += w4
           count_d6 += 1
        elif dir_mast < 213.5 and (dir_mast > 191 or dir_mast == 191):
           dmast_symbol = "SSW"
           w_d7 += w4
           count_d7 += 1
        elif dir_mast < 191 and (dir_mast > 168.5 or dir_mast == 168.5):
           dmast_symbol = "S"
           w_d8 += w4
           count_d8 += 1
        elif dir_mast < 168.5 and (dir_mast > 146 or dir_mast == 146):
           dmast_symbol = "SSE"
           w_d9 += w4
           count_d9 += 1
        elif dir_mast< 146 and (dir_mast > 123.5 or dir_mast == 123.5):
           dmast_symbol= "SE"
           w_d10 += w4
           count_d10 += 1
        elif dir_mast < 123.5 and (dir_mast > 101 or dir_mast == 101):
           dmast_symbol = "ESE"
           w_d11 += w4
           count_d11+= 1
        elif dir_mast < 101 and (dir_mast > 78.5 or dir_mast == 78.5):
           dmast_symbol = "E"
           w_d12 += w4
           count_d12 += 1
        elif dir_mast < 78.5 and (dir_mast > 56 or dir_mast == 56):
           dmast_symbol = "ENE"
           w_d13 += w4
           count_d13 += 1
        elif dir_mast < 56 and (dir_mast > 33.5 or dir_mast == 33.5):
           dmast_symbol = "NE"
           w_d14 += w4
           count_d14 += 1
        elif dir_mast<33.5 and(dir_mast>11 or dir_mast==11):
           dmast_symbol = "NNE"
           w_d15 += w4
           count_d15 += 1
        elif dir_mast < 11 or (dir_mast > 348.5 or dir_mast == 348.5):
           dmast_symbol = "N"
           w_d16 += w4
           count_d16 += 1
        DIR_cluster.append(dmast_symbol)
        DCLASS_perc=[count_d16/len(DIR_cluster),count_d15/len(DIR_cluster),count_d14/len(DIR_cluster),count_d13/len(DIR_cluster),count_d12/len(DIR_cluster),count_d11/len(DIR_cluster),count_d10/len(DIR_cluster),count_d9/len(DIR_cluster),count_d8/len(DIR_cluster),count_d7/len(DIR_cluster),count_d6/len(DIR_cluster),count_d5/len(DIR_cluster),count_d4/len(DIR_cluster),count_d3/len(DIR_cluster),count_d2/len(DIR_cluster),count_d1/len(DIR_cluster)]
        count+=1
        W1mean+=w1
        W2mean+=w2
        W3mean+=w3
        W4mean+=w4

        if w1!=0:
           ti1=ws1/w1
        else:
           ti1=0
        if w2!=0:
           ti2=ws2/w2
        else:
           ti2=0
        if w3!=0:
           ti3=ws3/w3
        else:
           ti3=0
        if w4!=0:
           ti4=ws4/w4
        else:
           ti4=0
        WG1.append(wg1);WG2.append(wg2);WG3.append(wg3);WG4.append(wg4)
        WS1.append(ws1);WS2.append(ws2);WS3.append(ws3);WS4.append(ws4)
        TI1.append(ti1);TI2.append(ti2);TI3.append(ti3);TI4.append(ti4)
        TEMPER_mast.append(temp_mast)
        if period_label==period_of_report:
           count1+=1
           w1mean+=w1
           w2mean+=w2
           w3mean+=w3
           w4mean+=w4
           Temper_mast.append(temp_mast)
           Wg1.append(wg1);Wg2.append(wg2);Wg3.append(wg3);Wg4.append(wg4)
           Ws1.append(ws1);Ws2.append(ws2);Ws3.append(ws3);Ws4.append(ws4)
           Ti1.append(ti1);Ti2.append(ti2);Ti3.append(ti3);Ti4.append(ti4)
           if dir_mast < 348.5 and (dir_mast > 326 or dir_mast == 326):
               dmastsymbol = "NNW"
               countd1 += 1
           elif dir_mast < 326 and (dir_mast > 303.5 or dir_mast == 303.5):
               dmastsymbol = "NW"
               countd2 += 1
           elif dir_mast < 303.5 and (dir_mast > 281 or dir_mast == 281):
               dmastsymbol = "WNW"
               countd3 += 1
           elif dir_mast < 281 and (dir_mast > 258.5 or dir_mast == 258.5):
               dmastsymbol = "W"
               countd4 += 1
           elif dir_mast < 258.5 and (dir_mast > 236 or dir_mast == 236):
               dmastsymbol = "WSW"
               countd5 += 1
           elif dir_mast < 236 and (dir_mast > 213.5 or dir_mast == 213.5):
               dmastsymbol = "SW"
               countd6 += 1
           elif dir_mast < 213.5 and (dir_mast > 191 or dir_mast == 191):
               dmastsymbol = "SSW"
               countd7 += 1
           elif dir_mast < 191 and (dir_mast > 168.5 or dir_mast == 168.5):
               dmastsymbol = "S"
               countd8 += 1
           elif dir_mast < 168.5 and (dir_mast > 146 or dir_mast == 146):
               dmastsymbol = "SSE"
               countd9 += 1
           elif dir_mast < 146 and (dir_mast > 123.5 or dir_mast == 123.5):
               dmastsymbol = "SE"
               countd10 += 1
           elif dir_mast < 123.5 and (dir_mast > 101 or dir_mast == 101):
               dmastsymbol = "ESE"
               countd11 += 1
           elif dir_mast < 101 and (dir_mast > 78.5 or dir_mast == 78.5):
               dmastsymbol = "E"
               countd12 += 1
           elif dir_mast < 78.5 and (dir_mast > 56 or dir_mast == 56):
               dmastsymbol = "ENE"
               countd13 += 1
           elif dir_mast < 56 and (dir_mast > 33.5 or dir_mast == 33.5):
               dmastsymbol = "NE"
               countd14 += 1
           elif dir_mast < 33.5 and (dir_mast > 11 or dir_mast == 11):
               dmastsymbol = "NNE"
               countd15 += 1
           elif dir_mast < 11 or (dir_mast > 348.5 or dir_mast == 348.5):
               dmastsymbol = "N"
               countd16 += 1
           Dir_cluster.append(dmast_symbol)
           Dclass_perc = [countd16/len(Dir_cluster), countd15/len(Dir_cluster), countd14/len(Dir_cluster), countd13/len(Dir_cluster), countd12/len(Dir_cluster), countd11/len(Dir_cluster), countd10/len(Dir_cluster), countd9/len(Dir_cluster),countd8/len(Dir_cluster), countd7/len(Dir_cluster), countd6/len(Dir_cluster), countd5/len(Dir_cluster), countd4/len(Dir_cluster), countd3/len(Dir_cluster), countd2/len(Dir_cluster), countd1/len(Dir_cluster)]
    ###investigation of wind shear and turbulence intensity per direction sector
    Ashear = [];TI_MEAN = [];DCLASS_per=[];S_sherror=[]
    for ii in range(len(Dir_list)):
        class_dir=Dir_list[ii]
        index_d=[]
        dclass=DCLASS_perc[ii]
        dclass*=100
        dclass_str=str(float("{:.2f}".format(dclass)))+"%"
        DCLASS_per.append(dclass_str)
        for jj in range(len(DIR_cluster)):
            class_d = DIR_cluster[jj]
            if class_d==class_dir:
               index_d.append(jj)
            wmast_per_sector1=[];wmast_per_sector2=[];wmast_per_sector3=[];wmast_per_sector4=[]
            wsdv_per_sector1=[];wsdv_per_sector2=[];wsdv_per_sector3=[];wsdv_per_sector4=[];TI_ps=[]
        for xx in index_d:
            wmast_ps1=float(W1_av[xx])
            wmast_ps2=float(W2_av[xx])
            wmast_ps3=float(W3_av[xx])
            wmast_ps4=float(W4_av[xx])
            wsdv_ps1=float(W1_sdv[xx])
            wsdv_ps2=float(W2_sdv[xx])
            wsdv_ps3 = float(W3_sdv[xx])
            wsdv_ps4 = float(W4_sdv[xx])
            ti_ps4 = wsdv_ps4 / wmast_ps4
            ti_ps=ti_ps4
            wmast_per_sector1.append(wmast_ps1);wmast_per_sector2.append(wmast_ps2);wmast_per_sector3.append(wmast_ps3);wmast_per_sector4.append(wmast_ps4)
            wsdv_per_sector1.append(wsdv_ps1);wsdv_per_sector2.append(wsdv_ps2); wsdv_per_sector3.append(wsdv_ps3);wsdv_per_sector4.append(wsdv_ps4)
            TI_ps.append(ti_ps)#mean turbulence intensity per sector
        if TI_ps:
           TI_ps_mean=float("{:.2f}".format(statistics.mean(TI_ps)))
        else:
           TI_ps_mean=0
        ## log-linear fit to find wind shear per sector.
        list_w1_nn = [x1 + float(1.1) if x1 < 1.1 else x1 for x1 in wmast_per_sector1]  ## exclude zeros or negative natural logarithms.
        list_w2_nn = [x2 + float(1.1) if x2 < 1.1 else x2 for x2 in wmast_per_sector2]
        list_w4_nn = [x3 + float(1.1) if x3 < 1.1 else x3 for x3 in wmast_per_sector4]
        log_w1 = np.log(list_w1_nn, dtype=float)
        log_w2 = np.log(list_w2_nn, dtype=float)
        log_w4 = np.log(list_w4_nn, dtype=float)
        linear_log_model1 = LinearRegression(fit_intercept=True).fit(log_w1.reshape((-1, 1)), log_w4)
        linear_log_model2 = LinearRegression(fit_intercept=True).fit(log_w2.reshape((-1, 1)), log_w4)
        slope1 = linear_log_model1.coef_[0]
        slope2 = linear_log_model2.coef_[0]
        const1 = linear_log_model1.intercept_
        const2 = linear_log_model2.intercept_
        rsq1=linear_log_model1.score(log_w1.reshape((-1, 1)), log_w4)
        rsq2=linear_log_model2.score(log_w2.reshape((-1, 1)), log_w4)
        s_log=statistics.stdev(log_w4)
        s_sherror1=math.sqrt((1-rsq1)*(s_log**2))
        s_sherror2=math.sqrt((1-rsq2)*(s_log**2))
        s_sherror=statistics.mean([s_sherror1,s_sherror2])## standard deviation of the shear factor fitting error.
        ashear1 = ((slope1 - 1) * np.mean(log_w1) + const1) / math.log(h4 / h1)
        ashear2 = ((slope2 - 1) * np.mean(log_w2) + const2) / math.log(h4 / h2)
        ashear= float("{:.2f}".format(statistics.mean([ashear1, ashear2])))
        Ashear.append(ashear)
        wcom1=np.mean(wmast_per_sector1,dtype=float)
        wcom2=np.mean(wmast_per_sector2,dtype=float)
        wcom3=np.mean(wmast_per_sector3,dtype=float)
        wcom4=np.mean(wmast_per_sector4, dtype=float)
        S_sherror.append(s_sherror)
        TI_MEAN.append(TI_ps_mean)
    w1mean = w1mean / count1
    w2mean = w2mean / count1
    w3mean = w3mean / count1
    w4mean = w4mean / count1
    w1mean=float("{:.2f}".format(w1mean))##present result with 2 decimal places
    w2mean = float("{:.2f}".format(w2mean))
    w3mean = float("{:.2f}".format(w3mean))
    w4mean = float("{:.2f}".format(w4mean))
    wmax1=max(Wg1);wmax2=max(Wg2);wmax3=max(Wg3);wmax4=max(Wg4)
    wsdv1=float("{:.2f}".format(statistics.mean(Ws1)));wsdv2=float("{:.2f}".format(statistics.mean(Ws2)));wsdv3=float("{:.2f}".format(statistics.mean(Ws3)));wsdv4=float("{:.2f}".format(statistics.mean(Ws4)));
    Ti1=statistics.mean(Ti1);Ti2=statistics.mean(Ti2);Ti3=statistics.mean(Ti3);Ti4=statistics.mean(Ti4)
    Mean_list = np.array([w1mean, w2mean, w3mean, w4mean])
    Max_list=np.array([wmax1,wmax2,wmax3,wmax4])
    Sdv_list=np.array([wsdv1,wsdv2,wsdv3,wsdv4])
    Ti_list=np.array([float("{:.2f}".format(Ti1)),float("{:.2f}".format(Ti2)),float("{:.2f}".format(Ti3)),float("{:.2f}".format(Ti4))])
    Wind_list=np.vstack((Mean_list,Max_list,Sdv_list,Ti_list))
    Temperature_list=np.array([float("{:.1f}".format(min(Temper_mast))),float("{:.1f}".format(statistics.mean(Temper_mast))),float("{:.1f}".format(max(Temper_mast)))])
W1mean=W1mean/count
W2mean=W2mean/count
W3mean=W3mean/count
W4mean=W4mean/count
W1mean=float("{:.2f}".format(W1mean))##present result with 2 decimal places
W2mean = float("{:.2f}".format(W2mean))
W3mean = float("{:.2f}".format(W3mean))
W4mean = float("{:.2f}".format(W4mean))
WMEAN_cluster_help=[w_d16/count_d16,w_d15/count_d15,w_d14/count_d14,w_d13/count_d13,w_d12/count_d12,w_d11/count_d11,w_d10/count_d10,w_d9/count_d9,w_d8/count_d8,w_d7/count_d7,w_d6/count_d6,w_d5/count_d5,w_d4/count_d4,w_d3/count_d3,w_d2/count_d2,w_d1/count_d1]
WMEAN_cluster=[float("{:.2f}".format(wmean)) for wmean in WMEAN_cluster_help]
WMAX1=max(WG1);WMAX2=max(WG2);WMAX3=max(WG3);WMAX4=max(WG4)
Wsdv1=float("{:.2f}".format(statistics.mean(WS1)));Wsdv2=float("{:.2f}".format(statistics.mean(WS2)));Wsdv3=float("{:.2f}".format(statistics.mean(WS3)));Wsdv4=float("{:.2f}".format(statistics.mean(WS4)))
TI1=statistics.mean(TI1);TI2=statistics.mean(TI2);TI3=statistics.mean(TI3);TI4=statistics.mean(TI4)
MEAN_list = np.array([W1mean, W2mean, W3mean, W4mean])
MAX_list = np.array([WMAX1, WMAX2, WMAX3, WMAX4])
SDV_list=np.array([Wsdv1,Wsdv2,Wsdv3,Wsdv4])
TI_list=np.array([float("{:.2f}".format(TI1)),float("{:.2f}".format(TI2)),float("{:.2f}".format(TI3)),float("{:.2f}".format(TI4))])
WIND_list = np.vstack((MEAN_list, MAX_list,SDV_list,TI_list))
TEMP_list=np.array([float("{:.1f}".format(min(TEMPER_mast))),float("{:.1f}".format(statistics.mean(TEMPER_mast))),float("{:.1f}".format(max(TEMPER_mast)))])


del i
del j
del k
del ii
del kk
del jj
del ll
del mm

##correlation between weather station and Moloch data
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

mol_coord=point.split('_')
mol_coord=np.array(mol_coord,dtype=float)
WTG_coord = WTG_coord_Iasmos
WTG_alt = Alt_Iasmos
station_info = Station_data
st_coord = station_coord[3]
st_alt = H_wstation
href = Agl_wstation
zo = float(0.05)#roughness length, provide 2 different versions , i.e. 0.05 and 0.10
lat_model = mol_coord[0]
lon_model = mol_coord[1]

##derivation of mean wind from 3 old mast mean data at the position of the moloch model using spatial interpolation with inverse distance weighting
point_old1=[41.11245,25.04727];point_old2=[41.11563,25.11530];point_old3=[41.11503,25.01637]
lat_old1=point_old1[0];lon_old1=point_old1[1]
lat_old2=point_old2[0];lon_old2=point_old2[1]
lat_old3=point_old3[0];lon_old3=point_old3[1]
WeD1 = distlatlon(lat_old1, lat_model, lon_old1, lon_model)
WeD2 = distlatlon(lat_old2, lat_model, lon_old2, lon_model)
WeD3= distlatlon(lat_old3, lat_model, lon_old3, lon_model)
WeD1 = 1/(WeD1 ** 2)
WeD2 = 1/(WeD2 ** 2)
WeD3 = 1/(WeD3 ** 2)
wold1y1=4.075;wold1y2=4.150
wold2y1=4.375;wold2y2=4.600
wold3y1=3.225;wold3y2=3.025
Wmean_year1=((WeD1 * wold1y1) + (WeD2 * wold2y1)+(WeD3 * wold3y1)) / (WeD1 + WeD2 + WeD3)##mean wind climate @ 10m agl
Wmean_year2=((WeD1 * wold1y2) + (WeD2 * wold2y2)+(WeD3 * wold3y2)) / (WeD1 + WeD2 + WeD3)##mean wind climate @ 10m agl
Wmean_year1= Wmean_year1 * (math.log(hup / zo) / math.log(10 / zo))##mean wind climate @ 100m agl
Wmean_year2= Wmean_year2 * (math.log(hup / zo) / math.log(10 / zo))##mean wind climate @ 100m agl
Wmean=statistics.mean([Wmean_year1,Wmean_year2])

lat_station = st_coord[0]
lon_station = st_coord[1]
mol_data = Model_data
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
station_cl2020 = [station_meanw2020, station_sdvw2020]  ##gives station climate  at 10m agl
station_cl2021 = [station_meanw2021, station_sdvw2021]  ##gives station climate  at 10m agl
mol_timestamp_help=list(mol_data[:,0])
# convert from one date format to another using datetime object.
mol_timestamp = []
for x in mol_timestamp_help:
    if len(x) > 16:
       mol_date = datetime.strptime(x, "%Y-%m-%d %H:%M:%S")
       y = mol_date.strftime("%Y-%m-%d %H:%M")
    else:
       y = x
    mol_timestamp.append(y)
del mol_timestamp_help
index_model_end=[]
for k in range(len(mol_timestamp)):
    m_ts = mol_timestamp[k]
    if m_ts==station_timestamp_end:
       index_model_end.append(k)
## correlation between station wind speed and model wind speed overall and per direction sector
index_model_end = index_model_end[0]
# print(index_model_end)
mol_timestamp_help = mol_timestamp[0:index_model_end + 1]
mol_wspeed = mol_data[0:index_model_end + 1, 1:4]
mol_dir = mol_data[0:index_model_end + 1, 4:7]
mol_dir_help = mol_dir[:, 0]
mol_dir_help = list(mol_dir_help)
# finds all elements of station_timestamps that exist in moloch_timestamps
##and creates a common timeseries for correlation and analysis. Also it upscales the station wind speed at 100m agl
station_timestamp_help = list(station_timestamp_help)
station_wspeed_for_common = station_data[index_start:index_end, 3]
station_wspeed_for_common = list(station_wspeed_for_common)
station_wcommon = [];
station_wcommon_up = [];
wcommon = [];
dcommon = []
index_common = []
for l in range(len(station_timestamp_help)):
    st_help = station_timestamp_help[l]
    if st_help in mol_timestamp_help:
       ic = mol_timestamp_help.index(st_help)
       if ic < len(station_wspeed_for_common):
         st_c = station_wspeed_for_common[ic]
         station_wcommon.append(st_c)
         w_up = np.array(st_c, dtype=float) * (math.log(hup / zo) / math.log(href / zo))  ##upscales station wind speed at 100m agl
         station_wcommon_up.append(w_up)
         mol_c = wspeed[ic, :]
         mold_c = dir_help[ic]
         wcommon.append(mol_c)
         dcommon.append(mold_c)
         index_common.append(ic)
station_wc = np.array(station_wcommon, dtype=float)
station_wc_up = np.array(station_wcommon_up, dtype=float)
station_common_up_mean = np.mean(station_wc_up)
wcommon = np.array(wcommon, dtype=float)

wcommon_50m = list(mol_wcommon[:, 0])
wcommon_80m = list(mol_wcommon[:, 1])
wcommon_100m = list(mol_wcommon[:, 2])

input_dirm = []
labels=16
for zz in mol_dcommon:
    z_help = float(zz)
    input_dirm.append(z_help)
clusters_m, centroids_m = kmeans1d.cluster(input_dirm,labels)  ##cluster wind direction of moloch model to 16 sectors with kmeans algorithm, centroids are the class centre and clusters are the sector number of each element, total 16
Station_dclass = [];R_sq = []; Ashear_model=[]; Linear_regression_model = [];Dclass_perc_station = [];Std_corr_error=[]; S_shmodelerror=[]
for ii in range(len(centroids_m)):
    index_d = [];
    for jj in range(len(clusters_m)):
       class_d = clusters_m[jj]
       if class_d == ii:
          index_d.append(jj)
# print(index_d)
    mol_wcommon_per_sector_50m = [];
    mol_wcommon_per_sector_80m = [];
    mol_wcommon_per_sector_100m = [];
    station_wcommon_per_sector = [];
    mol_dcommon_per_sector = [];
    station_wcommon_up_per_sector = []
    for xx in index_d:
        mol_wc_ps50 = mol_wcommon_50m[xx]
        mol_wc_ps80 = mol_wcommon_80m[xx]
        mol_wc_ps100 = mol_wcommon_100m[xx]
        mol_dc = input_dirm[xx]
        station_wc_ps = station_wcommon[xx]
        station_wc_ps_up = station_wcommon_up[xx]
        mol_wcommon_per_sector_50m.append(mol_wc_ps50)
        mol_wcommon_per_sector_80m.append(mol_wc_ps80)
        mol_wcommon_per_sector_100m.append(mol_wc_ps100)
        station_wcommon_per_sector.append(station_wc_ps)
        station_wcommon_up_per_sector.append(station_wc_ps_up)
        mol_dcommon_per_sector.append(mol_dc)
    dclass_perc = len(index_d) / len(station_wcommon)  ##percentage of wind direction per sector
    ## fit a linear regression model to the log of the mol wind speed to derive wind shear per sector
    list_w50m_nn = [x1 + float(1.1) if x1 < 1.1 else x1 for x1 in mol_wcommon_per_sector_50m]  ## exclude zeros or negative natural logarithms.
    list_w80m_nn = [x2 + float(1.1) if x2 < 1.1 else x2 for x2 in mol_wcommon_per_sector_80m]
    list_w100m_nn = [x3 + float(1.1) if x3 < 1.1 else x3 for x3 in mol_wcommon_per_sector_100m]
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
    wcom50 = np.mean(mol_wcommon_per_sector_50m, dtype=float)
    wcom80 = np.mean(mol_wcommon_per_sector_80m, dtype=float)
    wcom100 = np.mean(mol_wcommon_per_sector_100m, dtype=float)
    Dclass_perc_station.append(dclass_perc)  ##list of wind direction sector percentage of data
    mean_dclass = min(moloch_dcommon_per_sector)  ##the minimum direction of each direction class, used in order to classify station direction
    S_shmodelerror.append(s_shmodelerror)
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
##upscale wind speed at weather station at 100m agl
station_wspeed_up = []
for ws in station_wspeed:
    ws = float(ws)
    ws_up = ws * (math.log(hup / zo) / math.log(href / zo))
    station_wspeed_up.append(ws_up)
del ws_up
del index_d
# print(statistics.mean(station_wspeed_up))
# create a complete 2 year time series of wind speed at 100m agl based on combined data from weather station and moloch model using the linear model per direction sector
Wind_model_up = [];Dir_model=[]
for kk in range(len(station_wspeed)):
    ws_up = station_wspeed_up[kk]
    station_d = station_dir[kk]
    if station_d!='---':
       for xxx in range(len(Station_dclass)):
           if station_d==Station_dclass[xxx]:
              d_class=xxx
    else:
       d_class=0## dominant direction class as derived from new mast.
    dir_station=centroids_m[d_class]
    doffset = 11
    dir_station = dir_station - doffset
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
    wind_model = float(wind_model)
    Wind_model_up.append(wind_model)  ##complete 2 year time series of wind speed at 100m agl based on combined data from weather station and moloch model
del exist_count
del index_d
del mean_dclass

mol_wspeed_help = mol_wspeed[:, 2];
mol_wspeed_help = mol_wspeed_help.astype('f')

##choose the timeseries with the highest mean wind speed in the common available time interval
mean_mol_speed = statistics.mean(mol_wspeed_help)
mean_recreated_speed = statistics.mean(Wind_model_up[-1 - len(mol_wspeed_help) + 1:])

def myfunc1(mol_wspeed_help):
    Wind_model_up[-1 - len(mol_wspeed_help) + 1:] = list(mol_wspeed_help)

def myfunc2():
    return

if mean_mol_speed > mean_recreated_speed:
    myfunc1(mol_wspeed_help)
else:
    myfunc2()

Dir_model[-1 - len(mol_wspeed_help) + 1:] = list(mol_dir_help)##replaces station wind dir with mol wind dir for the existing Mol model timestamps

std_measure=float(0.017)## anemometer measurement uncertainty  as an stdev in m/s

station_pressure=station_data[:,1]
station_temperature=station_data[:,2]
station_RH=station_data[:,5]
WTG_coord_help=list(WTG_coord)
W1=station_wspeed_up
W2=Wind_model_up
Annual_mean_wind_speed_year1 = np.mean(Wind_model_up[0:17521])
Annual_mean_wind_speed_year2 = np.mean(Wind_model_up[17521:])
Annual_mean_model_speed=statistics.mean([Annual_mean_wind_speed_year1,Annual_mean_wind_speed_year2])
Wind_WTG=np.zeros(shape=(len(W1),len(WTG_coord),len(H_hub)),dtype=float)
Dens_WTG=np.zeros(shape=Wind_WTG.shape,dtype=float)
Wind_100m=np.zeros(shape=(len(W1),len(WTG_coord)),dtype=float)
Wind_uncertainty=np.zeros(shape=(len(W1),))
for ii in range(len(WTG_coord_help)):
    lat_wtg=WTG_coord[ii,0]
    lon_wtg=WTG_coord[ii,1]
    wtg_alt=WTG_alt[ii]
    # WeD1 = distlatlon(lat_station, lat_wtg, lon_station, lon_wtg)
    # WeD2 = distlatlon(lat_model, lat_wtg, lon_model, lon_wtg)
    # WeD1 = 1/ (WeD1 ** 2)
    # WeD2 = 1/ ( WeD2 ** 2)
    W_wtg=np.zeros(shape=(len(W1),1),dtype=float)
    for ww in range(len(W1)):
        w1=W1[ww]## upscaled station wind speed at corresp. timestamp
        w2=W2[ww]## upscaled, complete model wind speed at corresp. timestamp
        w_wtg=w2*(Wmean/Annual_mean_model_speed)##upscale wind speed according to old mast's mean wind speed interpolation @ moloch model point
        W_wtg[ww]=w_wtg
    del w_wtg
    for dd in range(W_wtg.size):
        w_wtg=float(W_wtg[dd])
        d1=float(Dir_model[dd])
        mean_dclass = d1  ## just the mix model direction assigned to an old variable name( variable deleted) for convenience
        if mean_dclass < 348.5 and (mean_dclass > 326 or mean_dclass == 326):
            dsymb = "NNW"
        elif mean_dclass < 326 and (mean_dclass > 303.5 or mean_dclass == 303.5):
            dsymb = "NW"
        elif mean_dclass < 303.5 and (mean_dclass > 281 or mean_dclass == 281):
            dsymb = "WNW"
        elif mean_dclass < 281 and (mean_dclass > 258.5 or mean_dclass == 258.5):
            dsymb = "W"
        elif mean_dclass < 258.5 and (mean_dclass > 236 or mean_dclass == 236):
            dsymb = "WSW"
        elif mean_dclass < 236 and (mean_dclass > 213.5 or mean_dclass == 213.5):
            dsymb = "SW"
        elif mean_dclass < 213.5 and (mean_dclass > 191 or mean_dclass == 191):
            dsymb = "SSW"
        elif mean_dclass < 191 and (mean_dclass > 168.5 or mean_dclass == 168.5):
            dsymb = "S"
        elif mean_dclass < 168.5 and (mean_dclass > 146 or mean_dclass == 146):
            dsymb = "SSE"
        elif mean_dclass < 146 and (mean_dclass > 123.5 or mean_dclass == 123.5):
            dsymb = "SE"
        elif mean_dclass < 123.5 and (mean_dclass > 101 or mean_dclass == 101):
            dsymb = "ESE"
        elif mean_dclass < 101 and (mean_dclass > 78.5 or mean_dclass == 78.5):
            dsymb = "E"
        elif mean_dclass < 78.5 and (mean_dclass > 56 or mean_dclass == 56):
            dsymb = "ENE"
        elif mean_dclass < 56 and (mean_dclass > 33.5 or mean_dclass == 33.5):
            dsymb = "NE"
        elif mean_dclass < 33.5 and (mean_dclass > 11 or mean_dclass == 11):
            dsymb = "NNE"
        elif mean_dclass < 11 or (mean_dclass > 348.5 or mean_dclass == 348.5):
            dsymb = "N"
        exist_count = Dir_list.count(dsymb)  ##checks how many times the specific wind direction exists in the class of mast directions created above
        if exist_count > 0:
           index_dd = Dir_list.index(dsymb)  # returns the index of the 1st time station_dir exists in the Station_dclass list
        else:
           index_dd = R_sq.index(max(R_sq))  # returns the index of the most well fitted direction class
        shear_fd = Ashear[index_dd]  ## use the wind shear derived from the new mast.
        shear_fd = float(shear_fd)
        ##modelling uncertainty per sector
        s_ce = Std_corr_error[index_dd]
        s_sh = S_sherror[index_dd]
        s_unwind = s_ce + s_sh + std_measure
        Wind_uncertainty[dd] = s_unwind
        P = float(station_pressure[dd])
        P*=100 ##convert mb to Pa
        T = float(station_temperature[dd])
        T += 273.15 ##temperature in Kelvin
        RH = float(station_RH[dd])
        WDC1=[]
        for ee in range(3):
            h_hub=float(H_hub[ee])
            wdc1= 0.5 / math.log(h_hub / zo)  # wake decay coefficient version 1
            alt_hub = wtg_alt + h_hub ## altitude of the hub configuration
            W_hubh = w_wtg * ((h_hub / 100) ** shear_fd) ## wind speed at hub height
            # print(w_wtg)
            # print(W_hubh)
            Z = (Rad * st_alt) / (Rad + st_alt) ## geopotential altitude at weather station
            Zhub = (Rad * alt_hub) / (Rad + alt_hub) ## geopotential altitude at WTG hubheight.
            Thub = T - L * (Zhub - Z) ## temperature in Kelvin
            Wind_WTG[dd, ii, ee] = W_hubh
            Phub = P * (Thub / T) ** ((G * Md) / (Rconst * L)) ##%%atm pressure at WTG hub height.
            f = 1.00070 + 3.113 * (10 ** (-8)) * Phub + 5.4 * (10 ** (-7)) * ((Thub - 273.15) ** 2)  #enhancement factor or ratio of effective saturation vapor pressure of water in moist air to saturation vapor pressure of water in moist air, empirical eq.
            es = 1.7526 * (10 ** 11) * math.exp(-5315.56 / Thub) ##saturation vapor pressure of water , empirical eq.
            phub = ((Phub * Md) / (Rconst * Thub * Zc)) * (1 - ((1 - (Mw / Md)) * (RH / 100) * f * es * (1 / Phub))) ## estimated air density at Phub
            WDC1.append(wdc1)
                   # print(phub)##checks if air density is calculated right
            Dens_WTG[dd, ii, ee] = phub
        Wind_100m[dd, ii] = w_wtg

# print('The corresponding annual wind 2020 (mean, sdv) derived from each station is:', station_cl2020)##print them when comparison is going to take place
# print('The corresponding mean annual wind 2021 (mean, sdv) derived from each station is:', station_cl2021)##print them when comparison is going to take place
Wmod_year1=Annual_mean_wind_speed_year1*(Wmean/Annual_mean_model_speed)## correction of mixture of models with old masts' mean annual wind speed.
Wmod_year2=Annual_mean_wind_speed_year2*(Wmean/Annual_mean_model_speed)
print('The annual wind climate for year 2020 , as derived from met masts is:', Wmod_year1)## this is mean climate at 100m agl
print('The annual wind climate for year 2021 , as derived from met masts is:', Wmod_year2)## this is mean climate at 100m agl
print('The annual wind climate for year 2020 , as derived from the mixture of models is:', Annual_mean_wind_speed_year1)## this is mean climate at 100m agl
print('The annual wind climate for year 2021 , as derived from the mixture of models is:', Annual_mean_wind_speed_year2)## this is mean climate at 100m agl
wind_climate=np.array([[Wmod_year1,Wmod_year2],[Annual_mean_wind_speed_year1,Annual_mean_wind_speed_year2]],ndmin=2)
# print(Dclass_perc_station)##checks dominant direction of the weather station
filetosave=directory_to_save+"_wake_model1_roughness_"+str(zo)+"_2022.xlsx"
filetowrite = pd.ExcelWriter(filetosave, engine='xlsxwriter')##excel file to save processed data
# # #
# #power output
#
 del WTG_coord
 del WTG_alt
 del WTG_coord_help
 del k
 del kk
 del dd
#
 D=170 ##rotor diameter in [m]
 RR=D/2 ## rotor radius in [m]
#
# ## filetosave='C:/Users/plgeo/OneDrive/Υπολογιστής/FARIA/Preliminary_wind_resource_assessment.xlsx'
#
# ## wind to power conversion - includes downstream wind speed derived from thrust curve and modified Jensen model
directory_power_curve='C:/Users/plgeo/OneDrive/PC Desktop/Wind turbine type/';
file_curve=directory_power_curve+'SG_turbine_characteristics.xlsx'
NUM_curve= pd.read_excel(file_curve,sheet_name='power curve')
NUM_thrust=pd.read_excel(file_curve,sheet_name='Thrust coef')
#
wind1=list(NUM_curve.iloc[1:,0])
oper_dens=list(NUM_curve.iloc[0,1:])
powerc1=[]
Ct=[]
pow_limit=float(5997)
for y in range(len(wind1)):
   powc1=list(NUM_curve.iloc[y+1,1:])
   powc1=[pow_limit if powc>=6000 else powc for powc in powc1] #curtailment of power so that the nominal capacity is satisfied.
   ct1=list(NUM_thrust.iloc[y+1,1:])
   powerc1.append(powc1)
   Ct.append(ct1)
power_c1=np.array(powerc1,dtype=float).reshape(len(wind1),len(oper_dens))
Ct=np.array(Ct,dtype=float).reshape(len(wind1),len(oper_dens))##thrust coeeficient curve.
# # print(power_c1)## prints the power curve matrix (rows= wind speed , cols = air density)
electr_loss_factor=0.97
#
#
WTG_coord = WTG_coord_Velanidia
WTG_alt = Alt_Velanidia
Density_wtg = Dens_WTG
Wind_speed_wtg = Wind_WTG
farm_name = location
WTG_coord_help = list(WTG_coord)

#
# ### power farm output without wake losses
Distance = []
for k in range(WTG_coord.shape[0] - 1):
    lat_wtg0 = WTG_coord[k, 0]
    lon_wtg0 = WTG_coord[k, 1]
    lat_wtg1 = WTG_coord[k + 1, 0]
    lon_wtg1 = WTG_coord[k + 1, 1]
    dist_wtg = distlatlon(lat_wtg0, lat_wtg1, lon_wtg0, lon_wtg1)  ## finds distance between consecutive WTGs in km
    Distance.append(dist_wtg)
mean_wtg_distance = float(statistics.mean(Distance))  ## finds wind farm average distance between consecutive WTGs in km

del lat_wtg0
del lon_wtg0
del lat_wtg1
del lon_wtg1

d_p=float(0.0004) # uncertainty in air density in kg/m^3

### power farm output without wake losses
def Power_farm():
    Power_wtg = np.zeros(shape=(Density_wtg.shape[0], Density_wtg.shape[1], Density_wtg.shape[2], 3), dtype=float)
    for aa in range(Density_wtg.shape[0]):
        d_ws = Wind_uncertainty[aa]
        for bb in range(Density_wtg.shape[1]):
            for cc in range(Density_wtg.shape[2]):
                densex = Density_wtg[aa, bb, cc]
                windex = Wind_speed_wtg[aa, bb, cc]
                densextr_help = np.array([densex - 1.68 * (d_p / math.sqrt(Density_wtg.shape[1] * Density_wtg.shape[2])), densex,densex + 1.68 * (d_p / math.sqrt(Density_wtg.shape[1] * Density_wtg.shape[2]))],dtype=float)  ##model air density uncertainty ,90% confidence interval
                windextr_help = np.array([windex - 1.68 * (d_ws / math.sqrt(Wind_speed_wtg.shape[1] * Wind_speed_wtg.shape[2])), windex,windex + 1.68 * (d_ws / math.sqrt(Wind_speed_wtg.shape[1] * Wind_speed_wtg.shape[2]))],dtype=float)  ##model wind speed uncertainty ,90% confidence interval
                ## extrapolation from wind speed and air density to power based on power curve of SG-6200 and conversion to half hour energy production
                for dd in range(densextr_help.size):
                    windextr = windextr_help[dd]
                    if windextr < 0:
                       windextr = 0
                    densextr = densextr_help[dd]
                    indexw1 = [];
                    indexw2 = [];
                    indexw3 = []
                    for kk in range(len(wind1)):
                        w_c = wind1[kk]
                        if w_c == windextr:
                           indexw1.append(kk)
                        elif w_c < windextr:
                           indexw2.append(kk)
                        elif w_c > windextr:
                           indexw3.append(kk)
                    indexd1 = [];
                    indexd2 = [];
                    indexd3 = []
                    for ll in range(len(oper_dens)):
                        d_c = oper_dens[ll]
                        if d_c == densextr:
                           indexd1.append(ll)
                        elif d_c < densextr:
                           indexd2.append(ll)
                        elif d_c > densextr:
                           indexd3.append(ll)
                    if indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                       indexw1 = int(indexw1)
                       indexd1 = int(indexd1)
                       powerextr = powerc1[indexw1, indexd1]
                    elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                       indexw1 = int(indexw1)
                       indexds = int(indexd2[-1])
                       indexdf = int(indexd3[0])
                       ds = float(oper_dens[indexds])
                       df = float(oper_dens[indexdf])
                       pd_s = power_c1[indexw1, indexds]
                       pd_f = power_c1[indexw1, indexdf]
                       powerextr = pd_s + (pd_f - pd_s) * abs((densextr - ds) / (df - ds))
                    elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
                       indexw1 = int(indexw1)
                       powerextr = power_c1[indexw1, 0]
                    elif not indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                       indexd1 = int(indexd1)
                       indexws = int(indexw2[-1])
                       indexwf = int(indexw3[0])
                       ws = float(wind1[indexws])
                       wf = float(wind1[indexwf])
                       ps = power_c1[indexws, indexd1]
                       pf = power_c1[indexwf, indexd1]
                       powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
                    elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                       indexds = int(indexd2[-1])
                       indexdf = int(indexd3[0])
                       ds = float(oper_dens[indexds])
                       df = float(oper_dens[indexdf])
                       indexws = int(indexw2[-1])
                       indexwf = int(indexw3[0])
                       ws = float(wind1[indexws])
                       wf = float(wind1[indexwf])
                       pws1 = power_c1[indexws, indexds]
                       pwf1 = power_c1[indexwf, indexds]
                       pwextrs = pws1 + (pwf1 - pws1) * abs((windextr - ws) / (wf - ws))
                       pws2 = power_c1[indexws, indexdf]
                       pwf2 = power_c1[indexwf, indexdf]
                       pwextrf = pws2 + (pwf2 - pws2) * abs((windextr - ws) / (wf - ws))
                       powerextr = pwextrs + (pwextrf - pwextrs) * abs((densextr - ds) / (df - ds))
                    elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
                       indexws = int(indexw2[-1])
                       indexwf = int(indexw3[0])
                       ws = float(wind1[indexws])
                       wf = float(wind1[indexwf])
                       ps = power_c1[indexws, 0]
                       pf = power_c1[indexwf, 0]
                       powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
                    elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[-1]) or densextr == float(oper_dens[-1])):
                       indexws = int(indexw2[-1])
                       indexwf = int(indexw3[0])
                       ws = float(wind1[indexws])
                       wf = float(wind1[indexwf])
                       ps = power_c1[indexws, -1]
                       pf = power_c1[indexwf, -1]
                       powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
                    else:
                       powerextr = 0
                    Power_wtg[aa, bb, cc, dd] = 0.5 * electr_loss_factor * (powerextr / 1000)  ## half hourly production in MWh
    Annual_Power_wtg_year1 = Power_wtg[0:17521, :, :, :].sum(axis=0)
    Annual_Power_wtg_year2 = Power_wtg[17521:, :, :, :].sum(axis=0)
    AP1 = Annual_Power_wtg_year1.sum(axis=0)
    AP2 = Annual_Power_wtg_year2.sum(axis=0)
    return AP1, AP2

Annual_Power_farm_year1, Annual_Power_farm_year2 = Power_farm()
print(Annual_Power_farm_year1)
print("\n")

## considering wake losses
##wake effect simulations
mean_wtg_distance *= 1000
distance_factor = mean_wtg_distance / D
mix_coef=1-(1/distance_factor) #mixing coefficient -- a , used in EB or SSQ model for wake overlapping
WDC1=statistics.mean(WDC1)
WDC3=0.5*statistics.mean(TI_oldmast)#wake decay coefficient version 2 , this version seems to overestimate , while version 1 seems to underestimate compared to the constant value used in WaSP for onshore WFs, used under neutral stability conditions
WDC2=0.075#wake decay coefficient version 3
WDC=WDC1## switch between WDC3, WDC2 and WDC1, if it is WDC1 it is 3 values, each for each different hub height config.
# print(WDC1);print(WDC2);print(WDC3)
Ao=math.pi*(RR**2)## rotor swept area

Power_wtg_we = np.zeros(shape=(Density_wtg.shape[0], Density_wtg.shape[1], Density_wtg.shape[2], 3), dtype=float)
Angle_thres=[];Angle_dif=[];Lateral_dist=[];Check1=[];Check2=[];Theta1=[];Theta2=[];Overlap_factor=[]
for aa in range(Density_wtg.shape[0]):
    d_ws = Wind_uncertainty[aa]
    wind_dir=float(Dir_model[aa])##computed wind direction from combo of station and model( mast cannot be used yet)
    if wind_dir>=0.00 and wind_dir<270.00:
       theta=wind_dir+90.00
    elif wind_dir>=270.00 and wind_dir<=360:
       theta=90.00-abs(360-wind_dir)
    theta=math.radians(theta)
    for bb in range(Density_wtg.shape[1]):
           lat_wtgt = WTG_coord[bb, 0]
           lon_wtgt = WTG_coord[bb, 1]
           x_wtgt,y_wtgt=latlon_to_xy(lat_wtgt,lon_wtgt)
           wake_red_array=np.empty(shape=(Density_wtg.shape[1],Density_wtg.shape[2]))##array of wake influence of target turbine from upwind turbines
           for ee in range(Density_wtg.shape[1]):
               if ee!=bb:
                  lat_wtgup=WTG_coord[ee,0]
                  lon_wtgup=WTG_coord[ee,1]
                  x_wtgup,y_wtgup=latlon_to_xy(lat_wtgup,lon_wtgup)
                  dist_up_down=math.sqrt(((x_wtgup-x_wtgt)**2)+((y_wtgup-y_wtgt)**2))## a better way to compute geographical distance because of calculation error in converting coordinates
                  # dist_up_down=distlatlon(lat_wtgup,lat_wtgt,lon_wtgup,lon_wtgt)
                  # dist_up_down=dist_up_down*1000 ## geographical distance between turbines in m
                  ##determine cone of influence for upwind turbine wrt target turbine
                  comp1=(x_wtgup-x_wtgt)*math.cos(theta)+(y_wtgup-y_wtgt)*math.sin(theta)+(RR/WDC)
                  comp2=x_wtgup-x_wtgt+(RR/WDC)*math.cos(theta)
                  comp3=y_wtgup-y_wtgt+(RR/WDC)*math.sin(theta)
                  cone_of_infl=math.acos(comp1/math.sqrt((comp2**2)+(comp3**2)))## angle corresponding to the cone of wake influence of upwind to downwind turbine
                  cone_of_infl=cone_of_infl*(180/math.pi)
                  ##determination of the downwind (downstream)distance
                  if theta>=0 and theta<(0.5*math.pi):
                     compx = (x_wtgup - x_wtgt) * math.cos(theta)
                     compy = (y_wtgup - y_wtgt) * math.sin(theta)
                     d_down = math.sqrt((compx ** 2) + (compy ** 2))  ## downwind distance at wind_dir
                  elif theta>=(0.5*math.pi) and theta<math.pi:
                     phi=math.asin(abs(x_wtgup-x_wtgt)/dist_up_down)
                     phi1=theta-phi-(0.5*math.pi)
                     d_down=dist_up_down*math.cos(phi1)
                  elif theta>=math.pi and theta<(1.5*math.pi):
                     phi = math.asin(abs(x_wtgup - x_wtgt) / dist_up_down)
                     phi1=(1.5*math.pi)-theta-phi
                     d_down=dist_up_down*math.cos(phi1)
                  else:
                     phi = math.asin(abs(y_wtgup - y_wtgt) / dist_up_down)
                     phi1=(2*math.pi)-theta-phi
                     d_down=dist_up_down*math.cos(phi1)
                  Rwake = RR + (WDC * d_down)  ## wake radius -- variable , function of wake decay constant and downwind distance
                  angle_thres=math.atan(D/d_down)## calculates threshold of the cone of influence per wind direction in radians -- fradsen proposal on how to calculate threshold angle of wake cone.
                  angle_thres=0.5*(angle_thres*(180/math.pi)+10)
                  # angle_thres=math.atan((0.5*Rwake)/(RR/WDC))## alternative formulation of threshold angle based on personal geometrical interpretation
                  # angle_thres=angle_thres*(180/math.pi)
                  Angle_thres.append(angle_thres)
                  ## if turbine bb is in the cone of wake influence of upwind turbine ee ( cone angle lower than or equal to angle threshold ,keep downwind distance at particular wind direction
                  if cone_of_infl>0 and cone_of_infl<=angle_thres:
                     D_down=d_down
                     di_j = math.sqrt((dist_up_down ** 2) - (D_down ** 2))  ## this is the lateral component of the 2 turbines' geographical distance wrt wind direction
                     Lateral_dist.append(di_j)
                     # ## geometrical analytical computation of the wake overlapping area
                     if di_j>=(Rwake-RR) and di_j<=(Rwake+RR):
                        xdd1=((di_j ** 2) + (Rwake ** 2) - (RR ** 2)) / (2 * di_j)
                        xdd2 = ((di_j ** 2) + (RR ** 2) - (Rwake ** 2)) / (2 * di_j)
                        Check1.append(xdd1/Rwake)
                        Check2.append(xdd2/RR)
                        theta1=math.acos(xdd1/Rwake)
                        theta2=math.acos(xdd2/RR)
                        Theta1.append(math.degrees(theta1))
                        Theta2.append(math.degrees(theta2))
                        Atriangle1=0.5*(Rwake**2)*math.sin(theta1)
                        Atriangle2=0.5*(RR**2)*math.sin(theta2)
                        Acone1=0.5*(Rwake**2)*theta1
                        Acone2=0.5*(RR**2)*theta2
                        Ashadow =2*(Acone1+Acone2-Atriangle1-Atriangle2) ## overlapping area of wake effect of each upwind turbine - analytical solution assuming steady state
                        # # Ashadow=0.5*(Rwake**2)*math.acos(0.5*di_j/Rwake)+0.5*(RR**2)*math.acos(di_j/D)-(0.25*di_j)*(math.sqrt((2*Rwake)**2-di_j**2)+math.sqrt(D**2-di_j**2))
                        overlap_f = Ashadow / Ao  ## factor corresponding to overlapping area of wake effect of each upwind turbine to each downwind turbine
                        print('WTG ' + str(bb + 1) + ' is in the wake influence of WTG ' + str(ee + 1))
                     elif di_j<(Rwake-RR):
                        xdd1=2*math.pi*Rwake
                        xdd2=2*math.pi*RR
                        Ashadow=Ao
                        overlap_f=1
                        print('WTG ' + str(bb + 1) + ' is in the wake influence of WTG ' + str(ee + 1))
                     else:
                        xdd1=0
                        xdd2=0
                        Ashadow=0
                        overlap_f=0
                        print('WTG ' + str(bb + 1) + ' is not in the wake influence of WTG ' + str(ee + 1))
                  elif cone_of_infl==0:
                     D_down=d_down
                     Ashadow=Ao
                     overlap_f=1
                     print('WTG ' + str(bb + 1) + ' is in the wake influence of WTG ' + str(ee + 1))
                  else:
                     D_down=0
                     Ashadow=0
                     overlap_f=0
                     print('WTG '+str(bb+1)+' is not in the wake influence of WTG '+str(ee+1))
               else:
                  D_down=0## this is the case of the upwind turbine and target turbine being the same.
                  Ashadow=0
                  overlap_f=0
                  print('WTG ' + str(bb + 1) + ' is the same as OR they belong to different sites, WTG ' + str(ee + 1))
               Overlap_factor.append(overlap_f)
               for ff in range(Density_wtg.shape[2]):
                   windextr = Wind_speed_wtg[aa, ee, ff]## upwind turbine wind velocity
                   densextr=  Density_wtg[aa, ee, ff]
                   indexw1 = [];
                   indexw2 = [];
                   indexw3 = []
                   for kk in range(len(wind1)):
                       w_c = wind1[kk]
                       if w_c == windextr:
                          indexw1.append(kk)
                       elif w_c < windextr:
                          indexw2.append(kk)
                       elif w_c > windextr:
                          indexw3.append(kk)
                   indexd1 = [];
                   indexd2 = [];
                   indexd3 = []
                   for ll in range(len(oper_dens)):
                       d_c = oper_dens[ll]
                       if d_c == densextr:
                          indexd1.append(ll)
                       elif d_c < densextr:
                          indexd2.append(ll)
                       elif d_c > densextr:
                          indexd3.append(ll)
                   ##determination of thrust coefficient from upwind velocity
                   if indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                      indexw1 = int(indexw1)
                      indexd1 = int(indexd1)
                      thrustextr = Ct[indexw1, indexd1]
                   elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                     indexw1 = int(indexw1)
                     indexds = int(indexd2[-1])
                     indexdf = int(indexd3[0])
                     ds = float(oper_dens[indexds])
                     df = float(oper_dens[indexdf])
                     cd_s = Ct[indexw1, indexds]
                     cd_f = Ct[indexw1, indexdf]
                     thrustextr = cd_s + (cd_f - cd_s) * abs((densextr - ds) / (df - ds))
                   elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
                     indexw1 = int(indexw1)
                     thrustextr = Ct[indexw1, 0]
                   elif not indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                     indexd1 = int(indexd1)
                     indexws = int(indexw2[-1])
                     indexwf = int(indexw3[0])
                     ws = float(wind1[indexws])
                     wf = float(wind1[indexwf])
                     cs = Ct[indexwf, indexd1]
                     cf = Ct[indexwf, indexd1]
                     thrustextr = cs + (cf - cs) * abs((windextr - ws) / (wf - ws))
                   elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                     indexds = int(indexd2[-1])
                     indexdf = int(indexd3[0])
                     ds = float(oper_dens[indexds])
                     df = float(oper_dens[indexdf])
                     indexws = int(indexw2[-1])
                     indexwf = int(indexw3[0])
                     ws = float(wind1[indexws])
                     wf = float(wind1[indexwf])
                     cws1 = Ct[indexws, indexds]
                     cwf1 = Ct[indexwf, indexds]
                     cwextrs = cws1 + (cwf1 - cws1) * abs((windextr - ws) / (wf - ws))
                     cws2 = Ct[indexws, indexdf]
                     cwf2 = Ct[indexwf, indexdf]
                     cwextrf = cws2 + (cwf2 - cws2) * abs((windextr - ws) / (wf - ws))
                     thrustextr = cwextrs + (cwextrf - cwextrs) * abs((densextr - ds) / (df - ds))
                   elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
                     indexws = int(indexw2[-1])
                     indexwf = int(indexw3[0])
                     ws = float(wind1[indexws])
                     wf = float(wind1[indexwf])
                     cs = Ct[indexws, 0]
                     cf = Ct[indexwf, 0]
                     thrustextr = cs + (cf - cs) * abs((windextr - ws) / (wf - ws))
                   elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[-1]) or densextr == float(oper_dens[-1])):
                     indexws = int(indexw2[-1])
                     indexwf = int(indexw3[0])
                     ws = float(wind1[indexws])
                     wf = float(wind1[indexwf])
                     cs = Ct[indexws, -1]
                     cf = Ct[indexwf, -1]
                     thrustextr = cs + (cf - cs) * abs((windextr - ws) / (wf - ws))
                   else:
                     thrustextr = Ct[0,0]
                   del indexw1
                   del indexw2
                   del indexw3
                   del indexd1
                   del indexd2
                   del indexd3
                   wind_up=windextr
                   del windextr
                   if D_down>0:
                      w_red_f=(1-math.sqrt(1-thrustextr))*((RR/(RR+(WDC*D_down)))**2)*overlap_f ## wind reduction factor of target turbine from upwind turbine , due to wake effect, as derived from modified Jensen model.
                   else:
                       w_red_f=0
                   wake_red_array[ee,ff]=w_red_f**2
           wind_red_factor=wake_red_array.sum(axis=0)
           # print(np.sqrt(wind_red_factor))
           for cc in range(Density_wtg.shape[2]):
               wi_red_f=wind_red_factor[cc]
               wind_target=wind_up*(1-math.sqrt(wi_red_f))## wind speed of target turbine as derived from upwind turbine and wind reduction factor due to wake effect
               # wind_target = wind_up * (1 - math.sqrt(mix_coef) * math.sqrt(wi_red_f))  ## wind speed of target turbine as derived from upwind turbine and wind reduction factor due to wake
               windex=wind_target
               # print(windex)
               densex = Density_wtg[aa, bb, cc]
               densextr_help = np.array([densex - 1.68 * (d_p / math.sqrt(Density_wtg.shape[1] * Density_wtg.shape[2])), densex,densex + 1.68 * (d_p / math.sqrt(Density_wtg.shape[1] * Density_wtg.shape[2]))],dtype=float)  ##model air density uncertainty ,90% confidence interval
               windextr_help = np.array([windex - 1.68 * (d_ws / math.sqrt(Wind_speed_wtg.shape[1] * Wind_speed_wtg.shape[2])), windex,windex + 1.68 * (d_ws / math.sqrt(Wind_speed_wtg.shape[1] * Wind_speed_wtg.shape[2]))],dtype=float)  ##model wind speed uncertainty ,90% confidence interval
               for dd in range(densextr_help.size):
                   windextr = windextr_help[dd]
                   if windextr < 0:
                       windextr = 0
                   densextr = densextr_help[dd]
                   indexw1 = [];
                   indexw2 = [];
                   indexw3 = []
                   for kk in range(len(wind1)):
                       w_c = wind1[kk]
                       if w_c == windextr:
                          indexw1.append(kk)
                       elif w_c < windextr:
                          indexw2.append(kk)
                       elif w_c > windextr:
                          indexw3.append(kk)
                   indexd1 = [];
                   indexd2 = [];
                   indexd3 = []
                   for ll in range(len(oper_dens)):
                       d_c = oper_dens[ll]
                       if d_c == densextr:
                          indexd1.append(ll)
                       elif d_c < densextr:
                          indexd2.append(ll)
                       elif d_c > densextr:
                          indexd3.append(ll)

                ## extrapolation from wind speed and air density to power based on power curve of SG-6200 and conversion to half hour energy production
                   if indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                      indexw1 = int(indexw1)
                      indexd1 = int(indexd1)
                      powerextr = powerc1[indexw1, indexd1]
                   elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                      indexw1 = int(indexw1)
                      indexds = int(indexd2[-1])
                      indexdf = int(indexd3[0])
                      ds = float(oper_dens[indexds])
                      df = float(oper_dens[indexdf])
                      pd_s = power_c1[indexw1, indexds]
                      pd_f = power_c1[indexw1, indexdf]
                      powerextr = pd_s + (pd_f - pd_s) * abs((densextr - ds) / (df - ds))
                   elif indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
                      indexw1 = int(indexw1)
                      powerextr = power_c1[indexw1, 0]
                   elif not indexw1 and indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                      indexd1 = int(indexd1)
                      indexws = int(indexw2[-1])
                      indexwf = int(indexw3[0])
                      ws = float(wind1[indexws])
                      wf = float(wind1[indexwf])
                      ps = power_c1[indexws, indexd1]
                      pf = power_c1[indexwf, indexd1]
                      powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
                   elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[0]) and densextr < float(oper_dens[-1])):
                      indexds = int(indexd2[-1])
                      indexdf = int(indexd3[0])
                      ds = float(oper_dens[indexds])
                      df = float(oper_dens[indexdf])
                      indexws = int(indexw2[-1])
                      indexwf = int(indexw3[0])
                      ws = float(wind1[indexws])
                      wf = float(wind1[indexwf])
                      pws1 = power_c1[indexws, indexds]
                      pwf1 = power_c1[indexwf, indexds]
                      pwextrs = pws1 + (pwf1 - pws1) * abs((windextr - ws) / (wf - ws))
                      pws2 = power_c1[indexws, indexdf]
                      pwf2 = power_c1[indexwf, indexdf]
                      pwextrf = pws2 + (pwf2 - pws2) * abs((windextr - ws) / (wf - ws))
                      powerextr = pwextrs + (pwextrf - pwextrs) * abs((densextr - ds) / (df - ds))
                   elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr < float(oper_dens[0]) or densextr == float(oper_dens[0])):
                      indexws = int(indexw2[-1])
                      indexwf = int(indexw3[0])
                      ws = float(wind1[indexws])
                      wf = float(wind1[indexwf])
                      ps = power_c1[indexws, 0]
                      pf = power_c1[indexwf, 0]
                      powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
                   elif not indexw1 and not indexd1 and (windextr > float(wind1[1]) and windextr < float(wind1[-1])) and (densextr > float(oper_dens[-1]) or densextr == float(oper_dens[-1])):
                      indexws = int(indexw2[-1])
                      indexwf = int(indexw3[0])
                      ws = float(wind1[indexws])
                      wf = float(wind1[indexwf])
                      ps = power_c1[indexws, -1]
                      pf = power_c1[indexwf, -1]
                      powerextr = ps + (pf - ps) * abs((windextr - ws) / (wf - ws))
                   else:
                      powerextr = 0
                   Power_wtg_we[aa, bb, cc, dd] = 0.5 * electr_loss_factor * (powerextr / 1000)  ## half hourly production in MWh including wake effect
Annual_Power_wtg_with_losses_year1=np.zeros(shape=(len(WTG_coord),Density_wtg.shape[2],densextr_help.size),dtype=float)
Annual_Power_wtg_with_losses_year2 = np.zeros(shape=(len(WTG_coord), Density_wtg.shape[2], densextr_help.size), dtype=float)
Annual_Power_wtg_with_losses_year1 = Power_wtg_we[0:17521, :, :, :].sum(axis=0)
Annual_Power_wtg_with_losses_year2 = Power_wtg_we[17521:, :, :, :].sum(axis=0)
Annual_PF_with_losses_year1=np.zeros(shape= (Density_wtg.shape[2], densextr_help.size), dtype=float)
Annual_PF_with_losses_year2 = np.zeros(shape=(Density_wtg.shape[2], densextr_help.size), dtype=float)
Annual_PF_with_losses_year1 = Annual_Power_wtg_with_losses_year1.sum(axis=0)
Annual_PF_with_losses_year2 = Annual_Power_wtg_with_losses_year2.sum(axis=0)
WDC=format(WDC,".3f")
print('\n')
print('WF ' + farm_name + ': ' + 'The array of 9 scenarios of APE of wake effect with WDC= ' + str(WDC) + ' for year 2020 in MWh is:\n',Annual_PF_with_losses_year1)  ## shows a matrix of 9 scenarios of annual APE including wake and electrical losses
print('\n')
print('WF ' + farm_name + ': ' + 'The array of 9 scenarios of APE of wake effect with WDC= ' + str(WDC) + ' for year 2021 in MWh is:\n',Annual_PF_with_losses_year2)
print('\n')

Wake_loss_year1 = np.zeros(shape=Annual_Power_farm_year1.shape, dtype=float)
Wake_loss_year2 = np.zeros(shape=Annual_Power_farm_year2.shape, dtype=float)
for sc in range(Annual_Power_farm_year1.shape[0]):
    for ci in range(Annual_Power_farm_year1.shape[1]):
        PF_no_loss1=Annual_Power_farm_year1[sc,ci]
        PF_no_loss2=Annual_Power_farm_year2[sc,ci]
        PF_loss1=Annual_PF_with_losses_year1[sc,ci]
        PF_loss2=Annual_PF_with_losses_year2[sc,ci]
        wake_loss_year1=(PF_no_loss1-PF_loss1)/PF_loss1
        wake_loss_year2=(PF_no_loss2-PF_loss2)/PF_loss2
        Wake_loss_year1[sc,ci]=wake_loss_year1
        Wake_loss_year2[sc,ci]=wake_loss_year2

df_2020=pd.DataFrame(Annual_PF_with_losses_year1,columns=["P95","P50","P5"])
df_2021 = pd.DataFrame(Annual_PF_with_losses_year2, columns=["P95", "P50", "P5"])
df_2020.to_excel(filetowrite,sheet_name="WDC="+str(WDC),startrow=0,startcol=0,index=False)
df_2021.to_excel(filetowrite,sheet_name="WDC="+str(WDC),startrow=4,startcol=0,index=False)
df_wl_2020=pd.DataFrame(Wake_loss_year1,columns=["P95","P50","P5"])
df_wl_2021=pd.DataFrame(Wake_loss_year2,columns=["P95","P50","P5"])
df_wl_2020.to_excel(filetowrite,sheet_name="wake loss @ zo="+str(zo)+", "+"WDC="+str(WDC),startrow=0,startcol=0,index=False)
df_wl_2021.to_excel(filetowrite,sheet_name="wake loss @ zo="+str(zo)+", "+"WDC="+str(WDC),startrow=4,startcol=0,index=False)
filetowrite.save()
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
