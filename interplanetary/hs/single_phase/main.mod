
#Single phase transfer (main.mod)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This AMPL model transcribe the NEP low-thrust OCP transfer from planet
#to planet into an NLP using the Hermite-Simpson Formula. The model is made 
#of a series of files, main.mod, definitions.inc, data.dat, ephcalc.inc,
#iniguess.inc
#
#

#--------------------------------------------------------------------------
#Grid options
param n:=20;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Departure and Arrival Choice
param depart:='earth' symbolic;
param target:='mars' symbolic;
#--------------------------------------------------------------------------

include include/definitions.inc;

#--------------------------------------------------------------------------
#Spacecraft Design
param M	     :=   (1000);	
param Isp    :=   (2500)   /T;		#Engine specific impulse (sec)
param Tmax   :=   (0.33)   /A/M/1000;	#Maximum Thrust level 	 (N)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Mission Design
param C3dep     :=   (0.001**2)  /V**2;	#Initial Launch C3 (km^2/sec^2)
param C3arr     :=   (0.001**2)  /V**2;	#Arrival C3 allowed (km^2/sec^2)

#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Initial and final time (used as initial guesses)
 
param tI	     :=   2000;		#Departure time (MJD2000)
param tT	     :=   400;		#Arrival time (MJD2000)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Solver options
options snopt_options "outlev=1 Major_iterations=1500";
options conopt_options "outlev=3";
#option log_file 'out/log.txt';
#--------------------------------------------------------------------------

options solver snopt;
include include/ephcalc.inc;
solve;

options solver snopt;
include include/guesstangential2.inc;

options snopt_options "Superbasics_limit=3000 outlev=2 Major_iterations=30000";


options solver snopt;		
include ocp.mod;
