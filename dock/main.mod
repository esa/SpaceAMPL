#Optimal Docking problem, impulsive transcription
#Equations of motion are, essentially, 2-D Hill-Cohessey-Wilthshire equations
#with one added degree of freedom for the attitude 
#Written by Nicolas Weiss and Dario Izzo (November 2009)


#Parameters---------------------------------------------------
#Generic 
	param n:=100;					#Numbers of nodes
	param eta:=0.08;				#[1/s] Tidal gravity constant
	param g0:=9.81;

#Spacecraft
	param m:=1.5;              	   		#[kg] initial mass
	param R:=0.5;					#[m] Radius of the spacecraft (disc)
	param Isp:=300;					#[s] Specific Impulse of the engine, Low Thrust = 2000s, Chemical = 300s 
	param maxthrustL:=0.1;				#[N] max lateral Thrust Left
	param maxthrustR:=0.1;	          		#[N] max lateral Thrust Right

#Initial Conditions
	param x0:=-2;					#[m] Initial x
	param z0:=0;					#[m] Initial z
	param vx0:=0;                  			#[m/s] initial velocity in x 
	param vz0:=0;	 				#[m/s] initial Velocity in z
	param theta0:=0;				#[rad] initial pitch
	param omega0:=0;				#[rad/s] initial pitch rate

#Final Conditions
	param xn:=0;					#[m] final x position
	param zn:=0;					#[m] final z position
	param vxn:=0;					#[m/s] final velocity x
	param vzn:=0;					#[m/s] final velocity z
	param thetan:=3.14;				#[rad] final pitch angle
	param omegan:=0;				#[rad/s] final pitch rate

#Other
	param tn:=10; 					#[s] Guess for the final time
#-------------------------------------------------------------

#Sets---------------------------------------------------------
	set I := 1..n;					#This is the main grid [0 .. dt/2 .... 3/2dt .... (n-2)dt-dt/2 ..(n-2)dt]
	set J := 2..n-1;				#Here the impulses are collocated
	set K := 2..n-2;				#This is used for the AMPL for loop defining the dynamic constraints
	set L := 1..n-1;				#This is used for the initial guess

#-------------------------------------------------------------

#Variables----------------------------------------------------
	var tf, >=0;					#final time
	var x{i in I};					#state 1
	var vx{i in I};					#state 2
	var z{i in I}; 					#state 3
	var vz{i in I};					#state 4
	var theta{i in I};			 	#state 5
	var omega{i in I};              		#state 6
	var uL{i in J}; 	      			#control 1
	var uR{i in J};					#control 2
#-------------------------------------------------------------
#Time variables

	var dt = tf/(n-2);		
	param starttime := 0;
	var time{i in J} = starttime + dt/2 + (i-2)*dt;

#Objective----------------------------------------------------
	
	minimize tiempo: tf;
	#maximize massa: m[n];
	#minimize power: exp ( - sum{h in J}  (abs(uL[h]) + abs(uR[h]))  /  (Isp*g0) ); 
#-------------------------------------------------------------

#Constraints--------------------------------------------------

	#Discontinuity created by the impulses
	var Dvx{i in J} =  (uL[i] + uR[i]) * cos(theta[i]);
	var Dvz{i in J} =  (uL[i] + uR[i]) * sin(theta[i]);
	var Domega{i in J} = (uL[i] - uR[i]) / (m*R);
       
	#First propagation of dt/2
	subject to dynamic1s: x[2] = x[1] + vx[1] * dt/2 + (2 * eta * vz[1] + 3 * (eta^2) * x[1]) * (dt^2)/4;
	subject to dynamic2s: vx[2] = vx[1] + (2 * eta * vz[1] + 3 * (eta^2) * x[1]) *dt/2;
	subject to dynamic3s: z[2] = z[1] + vz[1] * dt/2 - 2 * eta * vx[1] * (dt^2)/4;
	subject to dynamic4s: vz[2] = vz[1] - 2 * eta * vx[1] * dt/2;
        subject to dynamic5s: theta[2] = theta[1] + omega[1] * dt/2;
	subject to dynamic6s: omega[2] = omega[1];
	
	#Dynamic explicit solution constraints
	subject to dynamic1{i in K}: x[i+1] = x[i] + (vx[i] + Dvx[i]) * dt + (2 * eta * vz[i] + 3 * (eta^2) * x[i]) * dt^2;
	subject to dynamic2{i in K}: vx[i+1] = vx[i] + Dvx[i] + (2 * eta * vz[i] + 3 * (eta^2) * x[i]) * dt;
	subject to dynamic3{i in K}: z[i+1] = z[i] + (vz[i] + Dvz[i]) * dt - 2*eta*vx[i]*dt^2;
	subject to dynamic4{i in K}: vz[i+1] = vz[i] + Dvz[i] - 2*eta*vx[i] * dt;
        subject to dynamic5{i in K}: theta[i+1] = theta[i] + (omega[i] + Domega[i]) * dt;
	subject to dynamic6{i in K}: omega[i+1] = omega[i] + Domega[i];
	
	#Last propagation of dt/2
	subject to dynamic1f: x[n] = x[n-1] + (vx[n-1] + Dvx[n-1])* dt/2 + (2*eta*vz[n-1] + 3*(eta^2)*x[n-1]) * (dt^2)/4;
	subject to dynamic2f: vx[n] = vx[n-1] + Dvx[n-1] + (2*eta*vz[n-1] + 3*(eta^2)*x[n-1]) *dt/2;
	subject to dynamic3f: z[n] = z[n-1] + (vz[n-1] + Dvz[n-1]) * dt/2 - 2*eta*vx[n-1]*(dt^2)/4;
	subject to dynamic4f: vz[n] = vz[n-1] + Dvz[n-1] - 2*eta*vx[n-1] * dt/2;
        subject to dynamic5f: theta[n] = theta[n-1] + (omega[n-1] + Domega[n-1]) * dt/2;
	subject to dynamic6f: omega[n] = omega[n-1] + Domega[n-1];

	#Boundary Conditions
	subject to InitialPositionx: x[1] = x0;
	subject to InitialPositionz: z[1] = z0;
	subject to InitialVelocityx: vx[1] = vx0;
        subject to InitialVelocityz: vz[1] = vz0;
        subject to InitialPitch: theta[1] = theta0;
	subject to InitialPitchRate: omega[1] = omega0;

	subject to FinalPositionx: x[n]^2+z[n]^2 <=0.01;
	subject to FinalVelocityz: vz[n]^2 + vx[n]^2 <= 0.01;
        subject to FinalPitch1: theta[n] <= 3.1415/16;
	subject to FinalPitch2: theta[n] >= -3.1415/16;


	#subject to FinalPitchRate: omega[n] = omegan;


	#Conrol Constraint
	subject to controlmagnitude1{i in J}: uL[i] - maxthrustR*dt <=0;
	subject to controlmagnitude2{i in J}: uL[i] + maxthrustR*dt >=0;
	subject to controlmagnitude3{i in J}: uR[i] - maxthrustL*dt <=0;
        subject to controlmagnitude4{i in J}: uR[i] + maxthrustL*dt >=0;
       

#-------------------------------------------------------------

#-------------------------------------------------------------	
#Initial Conditions
	let tf := tn;
	let vz[1] := vz0;
	let theta[1] := theta0;
	let omega[1] := omega0;
        let {i in I} vx[i] := vx0;
        let {i in L} theta[i+1] := theta[i];
        let {i in L} omega[i+1] := 0;
        let {i in L} vz[i+1] :=0;
	let {i in L} x[i+1] := ((xn - x0)/tn) * dt + x[i];
        let {i in L} z[i+1] := ((zn - z0)/tn) * dt + z[i];


	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", starttime, x[1], vx[1], z[1], vz[1], theta[1], omega[1], 0, 0> out/guess.out;
	for {i in J}
	{printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", time[i], x[i], vx[i], z[i], vz[i], theta[i], omega[i], uL[i]/dt, uR[i]/dt> out/guess.out;
	}
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", tf, x[n], vx[n], z[n], vz[n], theta[n], omega[n], 0, 0 > out/guess.out;
	close guess.out;

#-------------------------------------------------------------

#-------------------------------------------------------------
#Solver Options
	option solver snopt;
	option substout 1;
	option show_stats 1;
	options snopt_options "outlev=2 Major_iterations=1500";
#-------------------------------------------------------------

#-------------------------------------------------------------
#Solve!!!
	solve;
	#drop massa; 
	#subject to heightpos {i in I}: z[i] >= 10;
	#solve;
	#drop time;
	#maximize rw: mass[n] - sum{i in J} u2[i]^2 - tf;
	#solve;
#-------------------------------------------------------------
#Print the AMPL raw solution
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", starttime, x[1], vx[1], z[1], vz[1], theta[1], omega[1], 0, 0> out/sol_raw.out;
	for {i in J}
	{
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", time[i], x[i], vx[i], z[i], vz[i], theta[i], omega[i], uL[i]/dt, uR[i]/dt > out/sol_raw.out;
	}
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n",tf, x[n], vx[n], z[n], vz[n], theta[n], omega[n], 0, 0  > out/sol_raw.out;
	close sol_raw.out;
	
#Print the AMPL solution in the good grid points
	var xplot{i in 2..n-1} = x[i] + (vx[i] + Dvx[i])* dt/2 + (2*eta*vz[i] + 3*(eta^2)*x[i]) * (dt^2)/4;
	var vxplot{i in 2..n-1} = vx[i] + Dvx[i] + (2*eta*vz[1] + 3*(eta^2)*x[i]) * dt/2;
	var zplot{i in 2..n-1} = z[i] + (vz[i] + Dvz[i]) * dt/2 - 2*eta*vx[i] * (dt^2)/4;
	var vzplot{i in 2..n-1} = vz[i] + Dvz[i] - 2*eta*vx[i] * dt/2;
	var thetaplot{i in 2..n-1} = theta[i] + (omega[i] + Domega[i]) * dt/2;
	var omegaplot{i in 2..n-1} = omega[i] + Domega[i];
	
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", 0, x[1], vx[1], z[1], vz[1], theta[1], omega[1], 0, 0> out/sol.out;
	for {i in 2..n-1}
	{
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", (i-1)*dt, xplot[i], vxplot[i], zplot[i], vzplot[i], thetaplot[i], omegaplot[i], uL[i]/dt, uR[i]/dt > out/sol.out;
	}
	close sol.out;
	
	
