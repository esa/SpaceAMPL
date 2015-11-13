#####################################################################################
#
# Problem:  Optimal Control Problem (OCP)
# Dynamics: Variable mass point in 2D with one reaction wheel controlling the pitch rate 
# Transcription: Hermite-Simpson
#
# Author: Dario Izzo (Nov 2012)
#
#####################################################################################


#Parameters---------------------------------------------------
#Generic
	param n:=20;					#Numbers of nodes
	param g0:=9.81;					#[m/s^2] Earth gravity constant at sea level
	param gmoon:=1.6229;				#[m/s^2] Moon gravity constant

#Spacecarft
	param m0:=9472.06;         			#[kg] initial mass
	param Isp:=311;					#[s] Specific Impulse of the engine, Low Thrust = 2000s, Chemical = 300s 
	param maxthrust:=45760;				#[N] max Thrust
	param minthrust:=0;				
	param maxthetadot:=0.0698;			#[rad/s] max Pitch rate

#Initial Conditions
	param x0:=0;					#[m] Initial x
	param z0:=2300;					#[m] Initial z
	param vx0:=150;					#[m/s] initial velocity in x 
	param vz0:=-44;	 				#[m/s] initial Velocity in z
	param theta0:=-55*3.14159/180;			#[rad] initial pitch

#Final Conditions
	param xn:=2000;					#[m] final x position
	param zn:=10;				        #[m] final z position
	param vxn:=0;				        #[m/s] final velocity x
	param vzn:=-2.5;				#[m/s] final velocity z
	param thetan:= -16*3.14159/180;			#[rad] initial pitch

#Other
	param tn:=(vz0 + sqrt(vz0^2 + 2*z0*gmoon))/gmoon; #[s] Guess for the final time
#-------------------------------------------------------------	

#Sets---------------------------------------------------------
	set I := 1..n;
	set J := 1..n-1;
#-------------------------------------------------------------

#Variables---------------------------------------------------
	var x{i in I};
	var vx{i in I};
	var z{i in I};
	var vz{i in I};
	var theta{i in I}, <= 3.1416/2, >= -3.1416/2;   #state 5;
	var m{i in I},  >=0;
	var u1{i in I}, >=0;
	var u2{i in I};
	var u1m{i in J}, >=0;
	var u2m{i in J};
#-------------------------------------------------------------

#Time variables-----------------------------------------------
	var tf, >=0;
	var dt = tf/(n-1);
	var timegrid{i in I} = dt*(i-1);
#-------------------------------------------------------------

#Objective----------------------------------------------------
	#minimize tiempo: tf;
	minimize fuel: m0 - m[n];
#-------------------------------------------------------------

#Dynamic at the grid points-----------------------------------
	var f1{i in I} = vx[i];
	var f2{i in I} = u1[i] * sin(theta[i]) / m[i];
	var f3{i in I} = vz[i];
	var f4{i in I} = u1[i] * cos(theta[i]) / m[i] - gmoon;
	var f5{i in I} = u2[i];
	var f6{i in I} = -u1[i] / (Isp*g0);
#-----------------------------------------------------------------------

#State definition at mid-points via Simpson interpolation---------------
	var xm{i in J} 		= 	(x[i] + x[i+1])/2 + tf/(n-1)/8 * (f1[i] - f1[i+1]);
	var vxm{i in J} 	= 	(vx[i] + vx[i+1])/2 + tf/(n-1)/8 * (f2[i] - f2[i+1]);
	var zm{i in J} 		= 	(z[i] + z[i+1])/2 + tf/(n-1)/8 * (f3[i] - f3[i+1]);
	var vzm{i in J} 	= 	(vz[i] + vz[i+1])/2 + tf/(n-1)/8 * (f4[i] - f4[i+1]);
	var thetam{i in J}	= 	(theta[i] + theta[i+1])/2 + tf/(n-1)/8 * (f5[i] - f5[i+1]);
	var mm{i in J}	= 	(m[i] + m[i+1])/2 + tf/(n-1)/8 * (f6[i] - f6[i+1]);
#-----------------------------------------------------------------------

#Dynamic at the mid-points----------------------------------------------
	var f1m{i in J} = vxm[i];
	var f2m{i in J} = u1m[i]*sin(thetam[i])/mm[i];
	var f3m{i in J} = vzm[i];
	var f4m{i in J} = u1m[i]*cos(thetam[i])/mm[i] - gmoon;
	var f5m{i in J} = u2m[i];
	var f6m{i in J} = -u1m[i] / (Isp*g0);
#-----------------------------------------------------------------------

#Hermite Formula---------------------------------------------------------
subject to 
	dynamicx{i in J}:  x[i]  = x[i+1] - tf/(n-1)/6*(f1[i] + f1[i+1] + 4*f1m[i]);
	dynamicvx{i in J}:  vx[i]  = vx[i+1]  - tf/(n-1)/6*(f2[i] + f2[i+1] + 4*f2m[i]);
	dynamicz{i in J}:  z[i]  = z[i+1]  - tf/(n-1)/6*(f3[i] + f3[i+1] + 4*f3m[i]);
	dynamicvz{i in J}: vz[i] = vz[i+1] - tf/(n-1)/6*(f4[i] + f4[i+1] + 4*f4m[i]);
	dynamictheta{i in J}: theta[i] = theta[i+1] - tf/(n-1)/6*(f5[i] + f5[i+1] + 4*f5m[i]);
	dynamicmass{i in J}: m[i] = m[i+1] - tf/(n-1)/6*(f6[i] + f6[i+1] + 4*f6m[i]);
#--------------------------------------------------------------------------	

#Constraints------------------------------------------
	#Boundary Conditions
	subject to InitialPositionx: x[1] = x0;
	subject to InitialPositionz: z[1] = z0;
	subject to InitialVelocityx: vx[1] = vx0;
	subject to InitialVelocityz: vz[1] = vz0;
#	subject to InitialPitch: theta[1] = theta0;
	subject to InitialMass: m[1] = m0;

#Un-comment the line below to try a pin point landing
	#subject to FinalPositionx: x[n] = xn;
	subject to FinalPositionz: z[n] <= zn;       
	subject to FinalVelocityz: vz[n] >= vzn;
	subject to FinalVelocityx: vx[n] <= vxn;
#	subject to FinalPitch: theta[n] = thetan;
      
	
	#Control Constraint
	#Control Magnitude Constraint
	subject to controlmagnitude1{i in I}: u1[i] - maxthrust <=0;
#	subject to controlmagnitude2{i in I}: u1[i] - minthrust >=0;
	subject to controlmagnitude3{i in I}: u2[i] - maxthetadot <=0;
	subject to controlmagnitude4{i in I}: u2[i] + maxthetadot >=0;
	
	subject to controlmagnitude1m{i in J}: u1m[i] - maxthrust <=0;
#	subject to controlmagnitude2m{i in J}: u1m[i] - minthrust >=0;
	subject to controlmagnitude3m{i in J}: u2m[i] - maxthetadot <=0;
	subject to controlmagnitude4m{i in J}: u2m[i] + maxthetadot >=0;
#-------------------------------------------------------------

#Guess-------------------------------------------------------
	let tf := tn;
	let {i in I} m[i] := m0;
	let {i in I} vx[i] := vx0;
	let {i in I} x[i] := x0 + vx[i]*i*dt;
	let {i in I} u1[i] :=0.5;
	let {i in I} u2[i] :=0.5;
#-------------------------------------------------------------

#Solver Options-----------------------------------------------
	option solver ipopt;
	option substout 0;
	option show_stats 1;
	options ipopt_options "outlev=2";
	options snopt_options "outlev=2 Major_iterations=15000 Major_optimality_tolerance=1.0e-6 Major_feasibility_tolerance=1e-9";
#-------------------------------------------------------------

#Solve!!!-----------------------------------------------------
	solve;
#-------------------------------------------------------------

#Print the Solution with midpoints---------------------------
	for {i in J}
	{
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[i],  m[i], x[i], vx[i], z[i], vz[i], theta[i], u1[i], u2[i] > out/sol.out;
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[i] + dt/2, mm[i], xm[i], vxm[i], zm[i], vzm[i], thetam[i], u1m[i], u2m[i] > out/sol.out;
	}
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[n],  m[n], x[n], vx[n], z[n], vz[n], theta[n], u1[n], u2[n] > out/sol.out;
	close out/sol.out;
#-------------------------------------------------------------
