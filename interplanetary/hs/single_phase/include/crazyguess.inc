#Single phase transfer (iniguess.inc)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and contains the calculations that AMPL
#needs to do in order to build an initial guess for the transfer by building
#up a trajectory with tangential thrust satisfying the boundary conditions



#--------------------------------------------------------------------------
#Problem definition
subject to 
	InitialMass:      m[1] =  1;
	InitialPositionx: x[1] = x0;
	InitialPositiony: y[1] = y0;
	InitialPositionz: z[1] = z0;
	InitialVelocityx: dx[1] = dx0;
	InitialVelocityy: dy[1] = dy0;
	InitialVelocityz: dz[1] = dz0;
	
#--------------------------------------------------------------------------
#Initial Guess for the Tangential Thrust traj
let {i in I} x[i] :=   (x0+(xf-x0)*(i-1)/(n-1))*cos(tf/(n-1)*(i-1)) - (y0+(yf-y0)*(i-1)/(n-1))*sin(tf/(n-1)*(i-1));
let {i in I} y[i] :=   (x0+(xf-x0)*(i-1)/(n-1))*sin(tf/(n-1)*(i-1)) + (y0+(yf-y0)*(i-1)/(n-1))*cos(tf/(n-1)*(i-1));
let {i in I} z[i] :=   (z0+(zf-z0)*(i-1)/(n-1));
let {i in I} dx[i] :=   (xf-x0)/tf*cos(tf/(n-1)*(i-1)) - (xf-x0)*(i-1)/(n-1)*sin(tf/(n-1)*(i-1))- ((yf-y0)/tf*sin(tf/(n-1)*(i-1)) + (yf-y0)*(i-1)/(n-1)*cos(tf/(n-1)*(i-1)));
let {i in I} dy[i] :=   (xf-x0)/tf*sin(tf/(n-1)*(i-1)) + (xf-x0)*(i-1)/(n-1)*cos(tf/(n-1)*(i-1))+ ((yf-y0)/tf*cos(tf/(n-1)*(i-1)) - (yf-y0)*(i-1)/(n-1)*sin(tf/(n-1)*(i-1)));
let {i in I} dz[i] :=   (zf-z0)/tf;
let {i in I} m[i] := 1-0.05*(i-1)/(n-1);

#--------------------------------------------------------------------------
#Print the Initial Guess
let {i in I} ux[i]  := Tmax * dx[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);
let {i in I} uy[i]  := Tmax * dy[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);
let {i in I} uz[i]  := Tmax * dz[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);

let {i in J} uxm[i]  := Tmax * dx[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);
let {i in J} uym[i]  := Tmax * dy[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);
let {i in J} uzm[i]  := Tmax * dz[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);

for {i in J} {
printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n",
x[i],y[i],z[i],dx[i],dy[i],dz[i],m[i],ux[i],uy[i],uz[i],x[i],y[i],z[i],dx[i],dy[i],dz[i],m[i],ux[i],uy[i],uz[i]> out/InitialGuess.out;
}
#------------------------------------------------------------------------

#Print the initial and final times
printf "%17.16e, %17.16e \n", ti/d2u , tF-ti/d2u > out/TimesGuess.out;

#--------------------------------------------------------------------------
#As to calculate the TT initial guess we made use of a simplified dynamic with 
#sipmplified objectives, we need to dispose some variables and constraints
#(note that AMPL will automatically dispose also those variables linked
# to the one listed here



purge InitialMass,InitialPositionx,InitialPositiony,InitialPositionz,InitialVelocityx,InitialVelocityy,InitialVelocityz;
#--------------------------------------------------------------------------

