#####################################################################################
#
# Problem:  Optimal Control Problem (OCP)
# Dynamics: Variable mass point in 3D
# Transcription: Impulsive
#
# Author: Dario Izzo (Nov 2012)
#
#####################################################################################


#Parameters---------------------------------------------------
#Generic 
	param n:=20;				#Numbers of nodes
	param g0:=9.81;				#[m/s^2] Earth gravity constant at sea level
	param gmoon:=1.6229;			#[m/s^2] Moon gravity constant

#Spacecraft
	param Isp:=311;				#[s] Specific Impulse of the engine
	param maxthrust:=45760		;	#[N] max Thrust

#Initial Conditions
	param x0:=0;				#[m] Initial x
	param z0:=2300;				#[m] Initial z (height)
	param vx0:=150;                		#[m/s] Initial velocity in x 
	param vz0:=-44;			 	#[m/s] Initial Velocity in z
	param m0:=9472.06;             		#[kg] initial total mass	


#Final Conditions
	param xn:=2000;				#[m] final x position
	param zn:=10;			        #[m] final z position
	param vxn:=0;			        #[m/s] final velocity x
	param vzn:=-2.5;			#[m/s] final velocity z

#Other
	param tn:=(vz0 + sqrt(vz0^2 + 2*z0*gmoon))/gmoon; #[s] Guess for the final time
#-------------------------------------------------------------

#Sets---------------------------------------------------------
	set I := 1..n;					#This is the main grid [0 .. dt/2 .... 3/2dt .... (n-2)dt-dt/2 ..(n-2)dt]
	set J := 2..n-1;				#Here the impulses are collocated
	set K := 2..n-2;				#This is used for the AMPL for loop defining the dynamic constraints
	set L := 1..n-1;				#This is used for the initial guess
#-------------------------------------------------------------

#Variables of the OCP------------------------------------------
	var tf, >=0;					#final time
	var x{i in I};					#state 1
	var vx{i in I};					#state 2
	var z{i in I}; 					#state 3
	var vz{i in I};					#state 4
	var mass{i in I};               		#state 5
	var ux{i in J};    	   			#control 1
	var uz{i in J};					#control 2
#-------------------------------------------------------------

#Time definitions--------------------------------------------
	var dt = tf/(n-2);		
	param starttime := 0;
	var time{i in J} = starttime + dt/2 + (i-2)*dt;
#-------------------------------------------------------------

#Objective---------------------------------------------------
	#minimize tempo: tf;
	#maximize massa: mass[n];
        minimize fuel: (m0 - mass[n]);
#-------------------------------------------------------------

#Constraints--------------------------------------------------
	#Discontinuity created by the impulses
	var Dm{i in J}  =  -sqrt(ux[i]*ux[i] + uz[i]*uz[i])/(Isp*g0);
	var Dvz{i in J} =  (uz[i])/(mass[i] + Dm[i]/2);
	var Dvx{i in J} =  (ux[i])/(mass[i] + Dm[i]/2);

	#First propagation of dt/2
	subject to dynamic1s: x[2] = x[1] + vx[1] * dt/2;
	subject to dynamic2s: vx[2] = vx[1];
	subject to dynamic3s: z[2] = z[1] + vz[1] * dt/2 -0.5*gmoon*(dt^2/4);
	subject to dynamic4s: vz[2] = vz[1] - gmoon*dt/2;
	subject to dynamic5s: mass[2] = mass[1];
	
	#All other propagations
	subject to dynamic1{i in K}: x[i+1] = x[i] + (vx[i] + Dvx[i]) * dt;
	subject to dynamic2{i in K}: vx[i+1] = vx[i] + Dvx[i];
	subject to dynamic3{i in K}: z[i+1] = z[i] + (vz[i] + Dvz[i]) * dt -0.5*gmoon * (dt^2);
	subject to dynamic4{i in K}: vz[i+1] = vz[i] + Dvz[i] - gmoon*dt;
	subject to dynamic5{i in K}: mass[i+1] = mass[i] + Dm[i];
	
	#Last propagation of dt/2
	subject to dynamic1f: x[n] = x[n-1] + (vx[n-1] + Dvx[n-1]) * dt/2;
	subject to dynamic2f: vx[n] = vx[n-1] + Dvx[n-1];
	subject to dynamic3f: z[n] = z[n-1] + (vz[n-1] + Dvz[n-1]) * dt/2 -0.5*gmoon*(dt^2/4);
	subject to dynamic4f: vz[n] = vz[n-1] + Dvz[n-1] - gmoon*dt/2;
	subject to dynamic5f: mass[n] = mass[n-1] + Dm[n-1];

	#Boundary Conditions
	subject to InitialPositionx: x[1] = x0;
	subject to InitialPositionz: z[1] = z0;
	subject to InitialVelocityx: vx[1] = vx0;
	subject to InitialVelocityz: vz[1] = vz0;
	subject to InitialMass: mass[1] = m0;

#Un-comment the line below to activate pin-point landing
#	subject to FinalPositionx: x[n] = xn;
	subject to FinalPositionz: z[n] <= zn;       
	subject to FinalVelocityz: vz[n] >= vzn;
	subject to FinalVelocityx: vx[n] <= vxn;
	
#Un-comment the lines below to keep the optic flow constant
#	subject to ConstOF{i in K}: z[i]  = z0/vx0 * (vx[i] + Dvx[i]) - (vz[i] + Dvz[i]) * dt/2 - 0.5 * gmoon * ((dt)^2/4); 
#	subject to ConstOFend: z[n] = z0/vx0 * vx[n]; 

       
	#Conrol magnitude Constraint
	subject to controlmagnitude1{i in J}: sqrt(ux[i]*ux[i] + uz[i]*uz[i]) - maxthrust*dt <=0;  
#-------------------------------------------------------------

#-------------------------------------------------------------	
#Initial Conditions
	let tf := tn;
	let vz[1] := vz0;
	let {i in I} vx[i] := vx0;
	let {i in L} vz[i+1] := vz[i] - gmoon*dt;
	let {i in L} x[i+1] := x[i] + vx[i]*dt;
	let {i in L} z[i+1] := z[i] + vz[i]*dt - 0.5*gmoon*(dt^2);
	let {i in I} mass[i] := m0;
	let {i in J} ux[i] := 0.5;
#-------------------------------------------------------------

#Solver Options-----------------------------------------------
	option solver snopt;
	option substout 1;
	option show_stats 1;
	options snopt_options "outlev=2 Major_iterations=15000 Major_optimality_tolerance=1.0e-6 Major_feasibility_tolerance=1e-9";
#-------------------------------------------------------------

#Solve!!!-----------------------------------------------------
	solve;
#-------------------------------------------------------------

#Print the solution in the correct grid points-----------------
	var xplot{i in 2..n-1} = x[i] + (vx[i] + Dvx[i]) * dt/2;
	var vxplot{i in 2..n-1} = vx[i] + Dvx[i];
	var zplot{i in 2..n-1} = z[i] + (vz[i] + Dvz[i]) * dt/2 - 0.5 * gmoon * ((dt)^2/4);
	var vzplot{i in 2..n-1} = vz[i] + Dvz[i] - gmoon*dt/2;
	var massplot{i in 2..n-1} = mass[i] + Dm[i];
	
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", 0, mass[1], x[1], vx[1], z[1], vz[1], 0, 0 > out/sol.out;
	for {i in 2..n-1}
	{
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", (i-1)*dt, massplot[i], xplot[i], vxplot[i], zplot[i], vzplot[i], ux[i]/dt, uz[i]/dt> out/sol.out;
	}
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", tf, mass[n], x[n], vx[n], z[n], vz[n], 0, 0 > out/sol.out;
	close out/sol.out;
#-------------------------------------------------------------
	
	
	
