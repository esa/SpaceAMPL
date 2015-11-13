#####################################################################################
#
# Problem:  Optimal Control Problem (OCP)
# Dynamics: Variable mass point in 3D
# Transcription: Hermite-Simpson
#
# Author: Dario Izzo (Nov 2012)
#
#####################################################################################


#Parameters---------------------------------------------------
#Generic 
	param n:=20;				#Numbers of nodes
	param g0:=9.81;				#[m/s^2] Earth gravity constant at sea level
	param gmoon:=1.6229;			#[m/s^2] Moon gravity constant
	param pi:= 4*atan(1);			#definition of pi!!!!

#Spacecraft
	param Isp:=311;				#[s] Specific Impulse of the engine
	param maxthrust:=45760;			#[N] max Thrust

#Initial Conditions
	param x0:=0;				#[m] Initial x
	param y0:=0;				#[m] Initial y
	param z0:=2300;				#[m] Initial z (height)
	param vx0:=150;                 	#[m/s] Initial velocity in x 
	param vy0:=0;                 		#[m/s] Initial velocity in y
	param vz0:=-44;			 	#[m/s] Initial Velocity in z
	param m0:=9472.06;      		#[kg] initial total mass	

#Final Conditions
	param xn:=2000;				#[m] final x position
	param yn:=150;				#[m] final y position
	param zn:=10;			        #[m] final z position
	param vxn:=0;			        #[m/s] final velocity x
	param vyn:=0;			        #[m/s] final velocity y
	param vzn:=-2.5;			#[m/s] final velocity z

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
	var y{i in I};
	var vy{i in I};
	var z{i in I}, >=0;
	var vz{i in I}, <=0;
	var m{i in I},  >=0;
	var u{i in I}, >= 0;
	var phi{i in I}, <= pi, >= -pi;
	var theta{i in I}, <= pi/2, >= -pi/2;  
	var um{i in J}, >= 0;;
	var phim{i in J}, <= pi, >= -pi;
	var thetam{i in J}, <= pi/2, >= -pi/2;  
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
	var f2{i in I} = u[i] * sin(theta[i]) * cos(phi[i]) / m[i];
	var f3{i in I} = vy[i];
	var f4{i in I} = u[i] * sin(theta[i]) * sin(phi[i]) / m[i];
	var f5{i in I} = vz[i];
	var f6{i in I} = u[i] * cos(theta[i]) / m[i] - gmoon;
	var f7{i in I} = -u[i] / (Isp*g0);
#-----------------------------------------------------------------------

#State definition at mid-points via Simpson interpolation---------------
	var xm{i in J} 		= 	(x[i]  + x[i+1])/2  + tf/(n-1)/8 * (f1[i] - f1[i+1]);
	var vxm{i in J} 	= 	(vx[i] + vx[i+1])/2 + tf/(n-1)/8 * (f2[i] - f2[i+1]);
	var ym{i in J} 		= 	(y[i]  + y[i+1])/2  + tf/(n-1)/8 * (f3[i] - f3[i+1]);
	var vym{i in J} 	= 	(vy[i] + vy[i+1])/2 + tf/(n-1)/8 * (f4[i] - f4[i+1]);
	var zm{i in J} 		= 	(z[i]  + z[i+1])/2  + tf/(n-1)/8 * (f5[i] - f5[i+1]);
	var vzm{i in J} 	= 	(vz[i] + vz[i+1])/2 + tf/(n-1)/8 * (f6[i] - f6[i+1]);
	var mm{i in J}		= 	(m[i]  + m[i+1])/2  + tf/(n-1)/8 * (f7[i] - f7[i+1]);
#-----------------------------------------------------------------------

#Dynamic at the mid-points----------------------------------------------
	var f1m{i in J} = vxm[i];
	var f2m{i in J} = um[i]*sin(thetam[i])*cos(phim[i]) / mm[i];
	var f3m{i in J} = vym[i];
	var f4m{i in J} = um[i]*sin(thetam[i])*sin(phim[i]) / mm[i];
	var f5m{i in J} = vzm[i];
	var f6m{i in J} = um[i]*cos(thetam[i]) / mm[i] - gmoon;
	var f7m{i in J} = -um[i] / (Isp*g0);
#-----------------------------------------------------------------------

#Hermite Formula---------------------------------------------------------
subject to 
	dynamicx{i in J}:  x[i]  = x[i+1] -  tf/(n-1)/6*(f1[i] + f1[i+1] + 4*f1m[i]);
	dynamicvx{i in J}: vx[i] = vx[i+1] - tf/(n-1)/6*(f2[i] + f2[i+1] + 4*f2m[i]);
	dynamicy{i in J}:  y[i]  = y[i+1]  - tf/(n-1)/6*(f3[i] + f3[i+1] + 4*f3m[i]);
	dynamicvy{i in J}: vy[i] = vy[i+1] - tf/(n-1)/6*(f4[i] + f4[i+1] + 4*f4m[i]);
	dynamicz{i in J}:  z[i]  = z[i+1]  - tf/(n-1)/6*(f5[i] + f5[i+1] + 4*f5m[i]);
	dynamicvz{i in J}: vz[i] = vz[i+1] - tf/(n-1)/6*(f6[i] + f6[i+1] + 4*f6m[i]);
	dynamicmass{i in J}: m[i] = m[i+1] - tf/(n-1)/6*(f7[i] + f7[i+1] + 4*f7m[i]);
#--------------------------------------------------------------------------	

#Constraints------------------------------------------
#Boundary Conditions
	subject to InitialPositionx: x[1] = x0;
	subject to InitialPositiony: y[1] = y0;
	subject to InitialPositionz: z[1] = z0;
	subject to InitialVelocityx: vx[1] = vx0;
	subject to InitialVelocityy: vy[1] = vy0;
	subject to InitialVelocityz: vz[1] = vz0;
	subject to InitialMass: m[1] = m0;

#Un-comment the line below to activate pin-point landing
#	subject to FinalPositionx: x[n] = xn;
#	subject to FinalPositiony: y[n] = yn;
	subject to FinalPositionz: z[n] <= zn;       
	subject to FinalVelocityx: vx[n] <= vxn;
	subject to FinalVelocityy: vy[n] <= vyn;
	subject to FinalVelocityz: vz[n] >= vzn;
      
	
	#Control Constraint
	#Control Magnitude Constraint
	subject to controlmagnitude{i in I}: u[i] - maxthrust <=0;  
	subject to controlmagnitudem{i in J}: um[i] - maxthrust <=0; 
#-------------------------------------------------------------

#Guess-------------------------------------------------------
	let tf := tn;
	let {i in I} m[i] := m0;
	let {i in I} vx[i] := vx0;
	let {i in I} x[i] := x0;
	let {i in I} vy[i] := vy0;
	let {i in I} y[i] := y0;
	let {i in I} vz[i] := vz0;
	let {i in I} z[i] := z0;
	let {i in I} u[i] :=0.5;
	let {i in J} um[i] :=0.5;
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
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[i],  m[i], x[i], vx[i], y[i], vy[i], z[i], vz[i], u[i], theta[i], phi[i] > out/sol.out;
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[i] + dt/2, mm[i], xm[i], vxm[i], ym[i], vym[i], zm[i], vzm[i], um[i], thetam[i], phim[i] > out/sol.out;
	}
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[n],  m[n], x[n], vx[n], y[n], vy[n], z[n], vz[n], u[n], theta[n], phi[n] > out/sol.out;
	close out/sol.out;
#-------------------------------------------------------------
