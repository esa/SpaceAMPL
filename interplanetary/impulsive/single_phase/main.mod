#--------------------------------------------------------------------------
#Single phase transfer (main.mod)
#Written by Dario Izzo (February 2009)
#
#This AMPL model transcribe the NEP low-thrust OCP problem from planet
#to planet into an NLP using an impulsive thrust model. Care is taken to
#make the correct substitutions as to obtain a small dimensional problem.
#Each impulse is collocated on each grid node and final position is enforced
#as a matching constraint.
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
#Grid options: as a fixed grid is used, n is the only parameter used
param n:=30;
shell "cd include; python writeequations.py 30; cd ..";	
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Departure and Arrival Planets. List of available objects is defined in
#definitions.inc, in the Set Planet
param depart:='earth' symbolic;
param target:='venus' symbolic;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Initial and final time (used as initial guesses)
param tI	     :=   3000.5;		#Departure time (MJD2000)
param tT	     :=   600;			#Transfer time (MJD2000)
#--------------------------------------------------------------------------

include include/definitions.inc;

#--------------------------------------------------------------------------
#Spacecraft Design
param M	     :=   (1500);	
param Isp    :=   (3800)   /T;			#Engine specific impulse (sec)
param Tmax   :=   (0.33)   /A/M/1000;	#Maximum Thrust level 	 (N)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Boundary Conditions
param C3dep     :=   (0.001**2)  /V**2;	#Initial Launch C3 (km^2/sec^2)
param C3arr     :=   (0.001**2)  /V**2;	#Arrival C3 allowed (km^2/sec^2)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Default Solver options
options solver snopt;
options snopt_options "outlev=2 Major_iterations=500";
options conopt_options "outlev=3";
#option log_file 'out/log.txt';
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Variables are loaded, equation of motion are set and ephemerides evaluated
include include/equations.inc;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#An Initial guess is built and its solution written in out/InitialGuess.Out
options solver snopt;
fix timod;
fix tfmod;
include include/guess.inc;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#The real problem is built and its solution written in out/solution.Out
options snopt_options "outlev=2 Major_iterations=30000";
options solver snopt;

include ocp.mod;
#--------------------------------------------------------------------------
