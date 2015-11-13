
#Two phases transfer (main.mod)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This AMPL model transcribe the NEP low-thrust OCP transfer from planet
#to planet via an intermediate fly-by into an NLP using the Hermite-Simpson
#Formula. The model is made of a series of files, main.mod, definitions.inc, data.dat, eph.inc,
#iniguess.inc etc.
#
#

#--------------------------------------------------------------------------
#Grid options
param n1:=25;				#Number of grid points in the first phase
param n2:=35;				#Number of grid points in the second phase
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Departure and Arrival Choice
param depart:='ast85' symbolic;
param flyby :='ast141' symbolic;
param target:='ast141' symbolic;
#--------------------------------------------------------------------------

include include/definitions.inc;

#--------------------------------------------------------------------------
#Fly-by
param rmin   :=   (6871)   / R;		#Minimum fly-by altitude (km)
param mupla  :=   (398600) / MU;	#Gravitational parameter of the planet (km^3/sec^2)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Spacecraft Design
param M	     :=   (1700);	
param Isp    :=   (3000)   /T;		#Engine specific impulse (sec)
param Tmax   :=   (0.15)   /A/M/1000;	#Maximum Thrust level 	 (N)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Mission Design
param C3dep     :=   (0.001**2)  /V**2;	#Initial Launch C3 (km^2/sec^2)
param C3arr     :=   (0.001**2)  /V**2;	#Arrival C3 allowed (km^2/sec^2)

#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Times (used as initial guesses)
param tI1	     :=   12098;	#Departure time (MJD2000)
param tT1	     :=   495;		#Transfer time (MJD2000)

param tI2	     :=   tI1 + tT1;	#Departure time (MJD2000)
param tT2	     :=   475;		#Transfer time (MJD2000)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Solver options
options snopt_options "outlev=2 Major_iterations=100";
options conopt_options "outlev=3";
#option log_file 'out/log.txt';
options solver snopt;
#--------------------------------------------------------------------------

include include/ephcalc1.inc;    #ephemerides calculations for phase1
include include/ephcalc2.inc;    #ephemerides calculations for phase2


#--------------------------------------------------------------------------
#Print the non dimensional units
printf "%17.16e, %17.16e, %17.16e \n ",R,V,M > out/units.out;

#------------------------------------------------------------------------ 

options snopt_options "Superbasics_limit=3000 outlev=2 Major_iterations=2000";
include include/guesslambertph1.inc;
include include/guesslambertph2.inc;		
include ocp.mod;

