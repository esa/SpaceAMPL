#####################################################################################
#
# Problem:  Optimal Control Problem (OCP)
# Dynamics: 6DOF spacecraft with one main engine and 6 lateral engines for attitude
# Transcription: Hermite-Simpson
#
# Author: Dario Izzo (Nov 2012)
#
#####################################################################################


#Parameters---------------------------------------------------
#Generic 
	param n:=11;					#Numbers of nodes
	param g0:=9.81;					#[m/s^2] Earth gravity constant at sea level
	param gmoon:=1.6229;				#[m/s^2] Earth gravity constant
	param pi:= 4*atan(1);				#definition of pi!!!!

#Spacecarft
	param Isp:=311;					#[s] Specific Impulse of the engine
	param maxthrust:=44000;				#[N] max Thrust
	param maxthrustL:=440;				#[N] max Thrust (lateral engines)
	param mass:=9472.06;                 		#[kg] wet mass
	param R:=3;					#[m] Lateral thrusters offset
	param A:= 30550 * 0.62;				#[kg m^2] Ixx
	param B:= 33014 * 0.62;				#[kg m^2] Iyy
	param C:= 33826 * 0.62;				#[kg m^2] Izz

#Initial Conditions
	param x0:=0;					#[m] Initial x
	param y0:=500;					#[m] Initial y
	param z0:=2300;					#[m] Initial z
	param vx0:=150;					#[m/s] initial velocity in x 
	param vy0:=0;					#[m/s] initial velocity in y
	param vz0:=-44;	 				#[m/s] initial velocity in z
	param m0:=mass;					#[kg] initial total mass
	param p0:=0;					#[kg] initial angular velocity x (rad/s)
	param q0:=0;					#[kg] initial angular velocity y (rad/s)
	param r0:=0;					#[kg] initial angular velocity z (rad/s)
	param e10:=0;					#initial euler parameter
	param e20:=0;					#initial euler parameter
	param e30:=0;					#initial euler parameter
	param eta0:=1;					#initial euler parameter
					
#Final Conditions
	param xn:=4500;					#[m] final x position
	param yn:=0;					#[m] final y position
	param zn:=10;					#[m] final z position
	param vxn:=0.1;					#[m/s] final velocity x (tolerance)
	param vyn:=0.1;					#[m/s] final velocity y (tolerance)
	param vzn:=-2.5;				#[m/s] final velocity z (absolute value)



#Other
	param tn:=(vz0 + sqrt(vz0^2 + 2*z0*gmoon))/gmoon; #[s] Guess for the final time
#-------------------------------------------------------------	

#Sets---------------------------------------------------------
	set I := 1..n;
	set J := 1..n-1;
#-------------------------------------------------------------

#Variables---------------------------------------------------
	var x{i in I}, >=0;
	var vx{i in I};
	var y{i in I};
	var vy{i in I};
	var z{i in I} >=0;
	var vz{i in I} <=0;
	var p{i in I}, >=-1, <=1;
	var q{i in I}, >=-1, <=1;
	var r{i in I}, >=-1, <=1;
	var e1{i in I}, >=-1, <=1;
	var e2{i in I}, >=-1, <=1;
	var e3{i in I}, >=-1, <=1;
	var eta{i in I}, >=-1, <=1;
	var m{i in I}, >=0;
	
	var uT{i in I}, >= 0.1, <= maxthrust;
	var uTm{i in J}, >= 0.1, <= maxthrust;
	
	var upL{i in I}, >= 0.1, <= maxthrustL;
	var upLm{i in J}, >= 0.1, <= maxthrustL;
	var upR{i in I}, >= 0.1, <= maxthrustL;
	var upRm{i in J}, >= 0.1, <= maxthrustL;

	var uqL{i in I}, >= 0.1, <= maxthrustL;
	var uqLm{i in J}, >= 0.1, <= maxthrustL;
	var uqR{i in I}, >= 0.1, <= maxthrustL;
	var uqRm{i in J}, >= 0.1, <= maxthrustL;
	
	var urL{i in I}, >= 0.1, <= maxthrustL;
	var urLm{i in J}, >= 0.1, <= maxthrustL;
	var urR{i in I}, >= 0.1, <= maxthrustL;
	var urRm{i in J}, >= 0.1, <= maxthrustL;
#-------------------------------------------------------------

#Time variables-----------------------------------------------
	var tf, >=0;
	var dt = tf/(n-1);
	var timegrid{i in I} = dt*(i-1);
#-------------------------------------------------------------

#Objective----------------------------------------------------
	#minimize tiempo: tf;
	minimize fuel: (m0 - m[n]);
#-------------------------------------------------------------

#Thrust projection to the inertial frame----------------------
	var uTx{i in I} = (uT[i]+upL[i]+upR[i]+uqL[i]+uqR[i]) * 2 * (e1[i]*e3[i] - e2[i]*eta[i]) + (urR[i] + urL[i]) * 2 * (e1[i]*e2[i] + e3[i]*eta[i]);
	var uTy{i in I} = (uT[i]+upL[i]+upR[i]+uqL[i]+uqR[i]) * 2 * (e2[i]*e3[i] + e1[i]*eta[i])  + (urR[i] + urL[i]) * (1 - 2*(e3[i]*e3[i]+e1[i]*e1[i]));
	var uTz{i in I} = (uT[i]+upL[i]+upR[i]+uqL[i]+uqR[i]) * (1 - 2 * (e2[i]*e2[i]+e1[i]*e1[i]))  + (urR[i] + urL[i]) * 2 * (e3[i]*e2[i] - e1[i]*eta[i]);
	
	var uTsum{i in I} = (uT[i]+upL[i]+upR[i]+uqL[i]+uqR[i]+urL[i]+urR[i]);
#-------------------------------------------------------------

#Dynamic at the grid points-----------------------------------
	var f1{i in I} = vx[i];
	var f2{i in I} = uTx[i] / m[i];
	var f3{i in I} = vy[i];
	var f4{i in I} = uTy[i] / m[i];
	var f5{i in I} = vz[i];
	var f6{i in I} = uTz[i] / m[i] - gmoon;
	var f7{i in I} = -uTsum[i] / (Isp*g0);
	var f8{i in I} = (B-C) / A * q[i] * r[i]  + (upR[i]-upL[i]) * R / A;
	var f9{i in I} = (C-A) / B * p[i] * r[i]  + (uqR[i]-uqL[i]) * R / B;
	var f10{i in I} = (A-B) / C * q[i] * p[i] + (urR[i]-urL[i]) * R / C;
	var f11{i in I} =   0.5 * ( e2[i]*r[i] - e3[i]*q[i] + eta[i]*p[i]);
	var f12{i in I} =   0.5 * (-e1[i]*r[i] + e3[i]*p[i] + eta[i]*q[i]);
	var f13{i in I} =   0.5 * ( e1[i]*q[i] - e2[i]*p[i] + eta[i]*r[i]);
	var f14{i in I} = - 0.5 * ( e1[i]*p[i] + e2[i]*q[i] +  e3[i]*r[i]);
#-----------------------------------------------------------------------

#State definition at mid-points via Simpson interpolation---------------
	var xm{i in J} 		= 	(x[i]  + x[i+1])/2  + tf/(n-1)/8     * (f1[i] - f1[i+1]);
	var vxm{i in J} 	= 	(vx[i] + vx[i+1])/2 + tf/(n-1)/8     * (f2[i] - f2[i+1]);
	var ym{i in J} 		= 	(y[i]  + y[i+1])/2  + tf/(n-1)/8     * (f3[i] - f3[i+1]);
	var vym{i in J} 	= 	(vy[i] + vy[i+1])/2 + tf/(n-1)/8     * (f4[i] - f4[i+1]);
	var zm{i in J} 		= 	(z[i]  + z[i+1])/2  + tf/(n-1)/8     * (f5[i] - f5[i+1]);
	var vzm{i in J} 	= 	(vz[i] + vz[i+1])/2 + tf/(n-1)/8     * (f6[i] - f6[i+1]);
	var mm{i in J}		= 	(m[i]  + m[i+1])/2  + tf/(n-1)/8     * (f7[i] - f7[i+1]);
	var pm{i in J}		= 	(p[i]  + p[i+1])/2  + tf/(n-1)/8     * (f8[i] - f8[i+1]);
	var qm{i in J}		= 	(q[i]  + q[i+1])/2  + tf/(n-1)/8     * (f9[i] - f9[i+1]);
	var rm{i in J}		= 	(r[i]  + r[i+1])/2  + tf/(n-1)/8     * (f10[i] - f10[i+1]);	
	var e1m{i in J}		= 	(e1[i]  + e1[i+1])/2  + tf/(n-1)/8   * (f11[i] - f11[i+1]);
	var e2m{i in J}		= 	(e2[i]  + e2[i+1])/2  + tf/(n-1)/8   * (f12[i] - f12[i+1]);
	var e3m{i in J}		= 	(e3[i]  + e3[i+1])/2  + tf/(n-1)/8   * (f13[i] - f13[i+1]);
	var etam{i in J}	= 	(eta[i]  + eta[i+1])/2  + tf/(n-1)/8 * (f14[i] - f14[i+1]);
	
#-----------------------------------------------------------------------

	var uTxm{i in J} = (uTm[i]+upLm[i]+upRm[i]+uqLm[i]+uqRm[i]) * 2 * (e1m[i]*e3m[i] - e2m[i]*etam[i]) + (urRm[i] + urLm[i]) * 2 * (e1m[i]*e2m[i] + e3m[i]*etam[i]);
	var uTym{i in J} = (uTm[i]+upLm[i]+upRm[i]+uqLm[i]+uqRm[i]) * 2 * (e2m[i]*e3m[i] + e1m[i]*etam[i])  + (urRm[i] + urLm[i]) * (1 - 2*(e3m[i]*e3m[i]+e1m[i]*e1m[i]));
	var uTzm{i in J} = (uTm[i]+upLm[i]+upRm[i]+uqLm[i]+uqRm[i]) * (1 - 2 * (e2m[i]*e2m[i]+e1m[i]*e1m[i]))  + (urRm[i] + urLm[i]) * 2 * (e3m[i]*e2m[i] - e1m[i]*etam[i]);
	var uTsumm{i in J} = (uTm[i]+upLm[i]+upRm[i]+uqLm[i]+uqRm[i]+urLm[i]+urRm[i]);

#Dynamic at the mid-points----------------------------------------------
	var f1m{i in J} = vxm[i];
	var f2m{i in J} = uTxm[i] / mm[i];
	var f3m{i in J} = vym[i];
	var f4m{i in J} = uTym[i] / mm[i];
	var f5m{i in J} = vzm[i];
	var f6m{i in J} = uTzm[i] / mm[i] - gmoon;
	var f7m{i in J} = -uTsumm[i] / (Isp*g0);
	var f8m{i in J} =  (B-C) / A * qm[i] * rm[i] + (upRm[i]-upLm[i]) * R / A;
	var f9m{i in J} =  (C-A) / B * pm[i] * rm[i] + (uqRm[i]-uqLm[i]) * R / A;
	var f10m{i in J} = (A-B) / C * qm[i] * pm[i] + (urRm[i]-urLm[i]) * R / A;
	var f11m{i in J} =   0.5 * ( e2m[i]*rm[i] - e3m[i]*qm[i] + etam[i]*pm[i]);
	var f12m{i in J} =   0.5 * (-e1m[i]*rm[i] + e3m[i]*pm[i] + etam[i]*qm[i]);
	var f13m{i in J} =   0.5 * ( e1m[i]*qm[i] - e2m[i]*pm[i] + etam[i]*rm[i]);
#	var f14m{i in J} = - 0.5 * ( e1m[i]*pm[i] + e2m[i]*qm[i] +  e3m[i]*rm[i]);
#-----------------------------------------------------------------------

#Hermite Formula---------------------------------------------------------
subject to 
	dynamicx{i in J}:    x[i]   = x[i+1]   - tf/(n-1)/6*(f1[i] + f1[i+1]   + 4*f1m[i]);
	dynamicvx{i in J}:   vx[i]  = vx[i+1]  - tf/(n-1)/6*(f2[i] + f2[i+1]   + 4*f2m[i]);
	dynamicy{i in J}:    y[i]   = y[i+1]   - tf/(n-1)/6*(f3[i] + f3[i+1]   + 4*f3m[i]);
	dynamicvy{i in J}:   vy[i]  = vy[i+1]  - tf/(n-1)/6*(f4[i] + f4[i+1]   + 4*f4m[i]);
	dynamicz{i in J}:    z[i]   = z[i+1]   - tf/(n-1)/6*(f5[i] + f5[i+1]   + 4*f5m[i]);
	dynamicvz{i in J}:   vz[i]  = vz[i+1]  - tf/(n-1)/6*(f6[i] + f6[i+1]   + 4*f6m[i]);
	dynamicmass{i in J}: m[i]   = m[i+1]   - tf/(n-1)/6*(f7[i] + f7[i+1]   + 4*f7m[i]);
	dynamicp{i in J}:    p[i]   = p[i+1]   - tf/(n-1)/6*(f8[i] + f8[i+1]   + 4*f8m[i]);
	dynamicq{i in J}:    q[i]   = q[i+1]   - tf/(n-1)/6*(f9[i] + f9[i+1]   + 4*f9m[i]);
	dynamicr{i in J}:    r[i]   = r[i+1]   - tf/(n-1)/6*(f10[i] + f10[i+1] + 4*f10m[i]);
	dynamice1{i in J}:   e1[i]  = e1[i+1]  - tf/(n-1)/6*(f11[i] + f11[i+1] + 4*f11m[i]);
	dynamice2{i in J}:   e2[i]  = e2[i+1]  - tf/(n-1)/6*(f12[i] + f12[i+1] + 4*f12m[i]);
	dynamice3{i in J}:   e3[i]  = e3[i+1]  - tf/(n-1)/6*(f13[i] + f13[i+1] + 4*f13m[i]);
#	dynamiceta{i in J}:  eta[i] = eta[i+1] - tf/(n-1)/6*(f14[i] + f14[i+1] + 4*f14m[i]);
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
	subject to Initialp: p[1] = p0;
	subject to Initialq: q[1] = q0;
	subject to Initialr: r[1] = r0;
	subject to Initiale1: e1[1] = e10;
	subject to Initiale2: e2[1] = e20;
	subject to Initiale3: e3[1] = e30;
	subject to Initialeta: eta[1] = eta0;

#Un-comment the two lines below to activate pin-point landing
#	subject to FinalPositionx: x[n] = xn;
#	subject to FinalPositiony: y[n] = yn;

	subject to FinalPositionz: z[n] <= zn;       
	subject to FinalVelocityx: vx[n] <= vxn;
	subject to FinalVelocityx2: vx[n] >= vxn;
	subject to FinalVelocityy: vy[n] <= vyn;
	subject to FinalVelocityy2: vy[n] >= vyn;
	subject to FinalVelocityz: vz[n] >= vzn;
#	subject to Finale1: e1[n] = e10;
#	subject to Finale2: e2[n] = e20;
#	subject to Finale3: e3[n] = e30;
#	subject to Finaleta: eta[n] = eta0;
      
	subject to quaternion{i in I}: e1[i]^2 + e2[i]^2 +e3[i]^2 +eta[i]^2 =1;
	
#-------------------------------------------------------------

#Guess-------------------------------------------------------
	let tf := tn;
	let {i in I} m[i] := m0;
	let {i in I} vx[i] := vx0;
	let {i in I} x[i] := x0 + vx0*(i-1)*dt;
	let {i in I} vy[i] := vy0;
	let {i in I} y[i] := y0 + vy0*(i-1)*dt;;
	let {i in I} vz[i] := vz0 - gmoon*(i-1)*dt;
	let {i in I} z[i] := z0 + vz0*(i-1)*dt - 0.5*gmoon*(i-1)^2*(dt^2);
	let {i in I} eta[i]:=1;
	let {i in I} uT[i] :=maxthrust/2;
	let {i in J} uTm[i] :=maxthrust/2;
	let {i in I} upL[i] :=maxthrustL/2;
	let {i in J} upLm[i] :=maxthrustL/2;
	let {i in I} uqL[i] :=maxthrustL/2;
	let {i in J} uqLm[i] :=maxthrustL/2;
	let {i in I} urL[i] :=maxthrustL/2;
	let {i in J} urLm[i] :=maxthrustL/2;
	let {i in I} upR[i] :=maxthrustL/2;
	let {i in J} upRm[i] :=maxthrustL/2;
	let {i in I} uqR[i] :=maxthrustL/2;
	let {i in J} uqRm[i] :=maxthrustL/2;
	let {i in I} urR[i] :=maxthrustL/2;
	let {i in J} urRm[i] :=maxthrustL/2;


#-------------------------------------------------------------

#Solver Options-----------------------------------------------
	option solver ipopt;
	option substout 1;
	option show_stats 1;
	options ipopt_options "outlev=5 max_iter=10000 tol=1e-5";
	options snopt_options "outlev=2 Major_iterations=10000 Major_optimality_tolerance=1.0e-6 Major_feasibility_tolerance=1e-10";
#-------------------------------------------------------------

#Solve!!!-----------------------------------------------------
	solve;
#-------------------------------------------------------------

#Print the Solution with midpoints---------------------------
	for {i in J}
	{
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[i],  m[i], x[i], vx[i],  y[i], vy[i], z[i], vz[i], p[i], q[i], r[i], e1[i], e2[i], e3[i], eta[i], upL[i],upR[i], uqL[i],uqR[i], urL[i], urR[i], uT[i]  > out/sol.out;
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[i] + dt/2, mm[i], xm[i], vxm[i],  ym[i], vym[i], zm[i], vzm[i], pm[i], qm[i], rm[i], e1m[i], e2m[i], e3m[i], etam[i], upLm[i],upRm[i], uqLm[i],uqRm[i], urLm[i], urRm[i], uTm[i] > out/sol.out;
	}
	printf "%17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e %17.16e\n", timegrid[n], m[n], x[n], vx[n],  y[n], vy[n], z[n], vz[n], p[n], q[n], r[n], e1[n], e2[n], e3[n], eta[n], upL[n], upR[n], uqL[n],uqR[n], urL[n], urR[n], uT[n] > out/sol.out;
	close out/sol.out;
#-------------------------------------------------------------
