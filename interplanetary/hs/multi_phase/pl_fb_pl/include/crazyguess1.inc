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
	InitialPositionx1: x1[1] = x01;
	InitialPositiony1: y1[1] = y01;
	InitialPositionz1: z1[1] = z01;
	FinalPositionx1  : x1[n1] = xf1;
	FinalPositiony1  : y1[n1] = yf1;
	FinalPositionz1  : z1[n1] = zf1;
#-----------------------------------------------------------------------		


#--------------------------------------------------------------------------
#Initial Guess for the Tangential Thrust traj
let {i in I1} x1[i] :=   (x01+(xf1-x01)*(i-1)/(n1-1))*cos(tf1/(n1-1)*(i-1)) - (y01+(yf1-y01)*(i-1)/(n1-1))*sin(tf1/(n1-1)*(i-1));
let {i in I1} y1[i] :=   (x01+(xf1-x01)*(i-1)/(n1-1))*sin(tf1/(n1-1)*(i-1)) + (y01+(yf1-y01)*(i-1)/(n1-1))*cos(tf1/(n1-1)*(i-1));
let {i in I1} z1[i] :=   (z01+(zf1-z01)*(i-1)/(n1-1));
let {i in I1} dx1[i] :=   (xf1-x01)/tf1*cos(tf1/(n1-1)*(i-1)) - (xf1-x01)*(i-1)/(n1-1)*sin(tf1/(n1-1)*(i-1))- ((yf1-y01)/tf1*sin(tf1/(n1-1)*(i-1)) + (yf1-y01)*(i-1)/(n1-1)*cos(tf1/(n1-1)*(i-1)));
let {i in I1} dy1[i] :=   (xf1-x01)/tf1*sin(tf1/(n1-1)*(i-1)) + (xf1-x01)*(i-1)/(n1-1)*cos(tf1/(n1-1)*(i-1))+ ((yf1-y01)/tf1*cos(tf1/(n1-1)*(i-1)) - (yf1-y01)*(i-1)/(n1-1)*sin(tf1/(n1-1)*(i-1)));
let {i in I1} dz1[i] :=   (zf1-z01)/tf1;


let {i in I1} m1[i] := 1-0.05*(i-1)/(n1-1);

#--------------------------------------------------------------------------
#Print the Initial Guess
let {i in I1} ux1[i]  := Tmax * dx1[i]/(dx1[i]**2+dy1[i]**2+dz1[i]**2)^(1/2);
let {i in I1} uy1[i]  := Tmax * dy1[i]/(dx1[i]**2+dy1[i]**2+dz1[i]**2)^(1/2);
let {i in I1} uz1[i]  := Tmax * dz1[i]/(dx1[i]**2+dy1[i]**2+dz1[i]**2)^(1/2);

let {i in J1} uxm1[i]  := Tmax * dx1[i]/(dx1[i]**2+dy1[i]**2+dz1[i]**2)^(1/2);
let {i in J1} uym1[i]  := Tmax * dy1[i]/(dx1[i]**2+dy1[i]**2+dz1[i]**2)^(1/2);
let {i in J1} uzm1[i]  := Tmax * dz1[i]/(dx1[i]**2+dy1[i]**2+dz1[i]**2)^(1/2);
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------


#--------------------------------------------------------------------------
#Print the Initial Guess
for {i in I1} {
printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e,%17.16e\n",x1[i],y1[i],z1[i],dx1[i],dy1[i],dz1[i],m1[i],ux1[i],uy1[i],uz1[i]>out/InitialGuess1.out;
}
#------------------------------------------------------------------------

#Print the initial and final times
printf "%17.16e, %17.16e \n",  ti1/d2u , tF1-ti1/d2u > out/TimesGuess1.out;



#--------------------------------------------------------------------------
#As to calculate the LP initial guess we made use of a simplified dynamic we 
#we need to dispose some variables (note that AMPL will automatically dispose
#also those variables linked to the one listed here


purge
InitialPositionx1,InitialPositiony1,InitialPositionz1,FinalPositionx1,FinalPositiony1,FinalPositionz1;
#--------------------------------------------------------------------------

