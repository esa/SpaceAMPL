#Single phase transfer (iniguess.inc)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and contains the calculations that AMPL
#needs to do in order to build an initial guess for the transfer

#--------------------------------------------------------------------------
#Problem definition
subject to 
	InitialPositionx2: x2[1] = x02;
	InitialPositiony2: y2[1] = y02;
	InitialPositionz2: z2[1] = z02;
	FinalPositionx2  : x2[n2] = xf2;
	FinalPositiony2  : y2[n2] = yf2;
	FinalPositionz2  : z2[n2] = zf2;
#-----------------------------------------------------------------------		


#--------------------------------------------------------------------------
#Initial Guess for the Tangential Thrust traj
let {i in I2} x2[i] :=   (x02+(xf2-x02)*(i-1)/(n2-1))*cos(tf2/(n2-1)*(i-1)) - (y02+(yf2-y02)*(i-1)/(n2-1))*sin(tf2/(n2-1)*(i-1));
let {i in I2} y2[i] :=   (x02+(xf2-x02)*(i-1)/(n2-1))*sin(tf2/(n2-1)*(i-1)) + (y02+(yf2-y02)*(i-1)/(n2-1))*cos(tf2/(n2-1)*(i-1));
let {i in I2} z2[i] :=   (z02+(zf2-z02)*(i-1)/(n2-1));
let {i in I2} dx2[i] :=   (xf2-x02)/tf2*cos(tf2/(n2-1)*(i-1)) - (xf2-x02)*(i-1)/(n2-1)*sin(tf2/(n2-1)*(i-1))- ((yf2-y02)/tf2*sin(tf2/(n2-1)*(i-1)) + (yf2-y02)*(i-1)/(n2-1)*cos(tf2/(n2-1)*(i-1)));
let {i in I2} dy2[i] :=   (xf2-x02)/tf2*sin(tf2/(n2-1)*(i-1)) + (xf2-x02)*(i-1)/(n2-1)*cos(tf2/(n2-1)*(i-1))+ ((yf2-y02)/tf2*cos(tf2/(n2-1)*(i-1)) - (yf2-y02)*(i-1)/(n2-1)*sin(tf2/(n2-1)*(i-1)));
let {i in I2} dz2[i] :=   (zf2-z02)/tf2;


let {i in I2} m2[i] := 1-0.05*(i-1)/(n2-1);

#--------------------------------------------------------------------------
#Print the Initial Guess
let {i in I2} ux2[i]  := Tmax * dx2[i]/(dx2[i]**2+dy2[i]**2+dz2[i]**2)^(1/2);
let {i in I2} uy2[i]  := Tmax * dy2[i]/(dx2[i]**2+dy2[i]**2+dz2[i]**2)^(1/2);
let {i in I2} uz2[i]  := Tmax * dz2[i]/(dx2[i]**2+dy2[i]**2+dz2[i]**2)^(1/2);

let {i in J2} uxm2[i]  := Tmax * dx2[i]/(dx2[i]**2+dy2[i]**2+dz2[i]**2)^(1/2);
let {i in J2} uym2[i]  := Tmax * dy2[i]/(dx2[i]**2+dy2[i]**2+dz2[i]**2)^(1/2);
let {i in J2} uzm2[i]  := Tmax * dz2[i]/(dx2[i]**2+dy2[i]**2+dz2[i]**2)^(1/2);
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------


#--------------------------------------------------------------------------
#Print the Initial Guess
for {i in I2} {
printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e,%17.16e\n",x2[i],y2[i],z2[i],dx2[i],dy2[i],dz2[i],m2[i],ux2[i],uy2[i],uz2[i]>out/InitialGuess2.out;
}
#------------------------------------------------------------------------

#Print the initial and final times
printf "%17.16e, %17.16e \n",  ti2/d2u , tF2-ti2/d2u > out/TimesGuess2.out;



#--------------------------------------------------------------------------
#As to calculate the LP initial guess we made use of a simplified dynamic we 
#we need to dispose some variables (note that AMPL will automatically dispose
#also those variables linked to the one listed here


purge
InitialPositionx2,InitialPositiony2,InitialPositionz2,FinalPositionx2,FinalPositiony2,FinalPositionz2;
#--------------------------------------------------------------------------

