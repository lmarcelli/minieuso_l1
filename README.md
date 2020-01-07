# minieuso_l1
Mini-EUSO software to create Level1 data starting from raw data. 
You need "TLE" file containing TLE, Moon and Sun info. It is necessary to run a python script to obtain it. 

Processing Level0 files to Level1:
- the behaviour of the two "strange" PMTs is corrected;
- TLE, Moon and Sun info are added to Level0;
- CPU time offset info is added: to each timestamp_unix"/"gtu_time" is associated a new "timestamp_unix_off"/"gtu_time_off" abtained from Lech "city correspondance lock" or from fit estrapolation fro run 1,3, 8).

How to run Level1 software
There are two possibility:
1) On a single file, running the executable directly (useful if you want to process a single file or few files):
run the exe by command line passing as arguments: a)level0 file path; b)"TLE" file (with complete path) to be used; c)cpu_time_offset (daily offset) to be added to timestamp_unix (optional, if not passed the software calculate it automatically). 
2) Automatically on all level0 file in a folder, running bash script "automatic_l1_v2.sh":
run the bash script by command line passing as arguments: a)the path of the folder containing Level0 files ; b)"TLE" file (with complete path) to be used; c)time_offset to be added to timestamp_unix/gtu_time (optional, if not passed the software calculate it automatically). 

Here below the structures for TLE, Moon and Sun info:

    
    typedef struct {
      Double_t abstime;
      Double_t timeunix;
      Double_t lat;
      Double_t lon;
      Double_t alt;
      Double_t x;
      Double_t y;
      Double_t z;
      Double_t vx;
      Double_t vy;
      Double_t vz;
      Double_t CdLon;
      Double_t CdLat;
      Double_t X_geo;				
      Double_t Y_geo;				
      Double_t Z_geo;							
    } positionDiss;
    
    typedef struct {			
      Double_t Sun_from_ISS_isDay;		
      Double_t Sun_from_ISS_alt;		
      Double_t Sun_from_ISS_az;		
      Double_t Sun_from_ISS_dist;		
      Double_t Sun_from_Earth_isDay;		
      Double_t Sun_from_Earth_alt;		
      Double_t Sun_from_Earth_az;		
      Double_t Sun_from_Earth_dist;
    } positionDSun;
    
    typedef struct {		
      Double_t Moon_from_ISS_isDay;		
      Double_t Moon_from_ISS_alt;		
      Double_t Moon_from_ISS_az;		
      Double_t Moon_from_ISS_dist;       	
      Double_t Moon_from_Earth_isDay;		
      Double_t Moon_from_Earth_alt;		
      Double_t Moon_from_Earth_az;		
      Double_t Moon_from_Earth_dist;
      Double_t FractionMoon;
    } positionDMoon;
    
    
    
   
