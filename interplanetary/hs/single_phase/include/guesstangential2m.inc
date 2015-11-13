#Single phase transfer (iniguess.inc)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and contains the calculations that AMPL
#needs to do in order to build an initial guess for the transfer by building
#up a trajectory with tangential thrust satisfying the boundary conditions


var Th{i in I} <= Tmax, >= -Tmax, := -Tmax;
var Thm{i in J} <= Tmax, >= -Tmax, := -Tmax;

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
	
#-----------------------------------------------------------------------		

#Dynamic definition at grid points--------------------------------------
var f1{i in I} = dx[i];
var f2{i in I} = dy[i];
var f3{i in I} = dz[i];
var f4{i in I} = -x[i]/(x[i]**2+y[i]**2+z[i]**2)^(3/2) + Th[i]/m[i] * dx[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);
var f5{i in I} = -y[i]/(x[i]**2+y[i]**2+z[i]**2)^(3/2) + Th[i]/m[i] * dy[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);
var f6{i in I} = -z[i]/(x[i]**2+y[i]**2+z[i]**2)^(3/2) + Th[i]/m[i] * dz[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);
var f7{i in I} = -abs(Th[i])/Isp/g0;
#-----------------------------------------------------------------------

#State definition at mid-points via Simpson interpolation---------------
var xm{i in J}  = (x[i]+x[i+1])/2    + tf/(n-1)/8   * (f1[i] - f1[i+1]);
var ym{i in J}  = (y[i]+y[i+1])/2    + tf/(n-1)/8   * (f2[i] - f2[i+1]);
var zm{i in J}  = (z[i]+z[i+1])/2    + tf/(n-1)/8   * (f3[i] - f3[i+1]);
var dxm{i in J} = (dx[i]+dx[i+1])/2  + tf/(n-1)/8   * (f4[i] - f4[i+1]);
var dym{i in J} = (dy[i]+dy[i+1])/2  + tf/(n-1)/8   * (f5[i] - f5[i+1]);
var dzm{i in J} = (dz[i]+dz[i+1])/2  + tf/(n-1)/8   * (f6[i] - f6[i+1]);
var mm{i in J}  = (m[i]+m[i+1])/2    + tf/(n-1)/8   * (f7[i] - f7[i+1]);
#-----------------------------------------------------------------------

#Dynamic definition at mid-points---------------------------------------
var f1m{i in J} = dxm[i];
var f2m{i in J} = dym[i];
var f3m{i in J} = dzm[i];
var f4m{i in J} = -xm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2) + Thm[i]/mm[i] * dxm[i]/(dxm[i]**2+dym[i]**2+dzm[i]**2)^(1/2);
var f5m{i in J} = -ym[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2) + Thm[i]/mm[i] * dym[i]/(dxm[i]**2+dym[i]**2+dzm[i]**2)^(1/2);
var f6m{i in J} = -zm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2) + Thm[i]/mm[i] * dzm[i]/(dxm[i]**2+dym[i]**2+dzm[i]**2)^(1/2);
var f7m{i in J} = -abs(Thm[i])/Isp/g0;
#-----------------------------------------------------------------------

#Hermite Formula--------------------------------------------------------

	subject to 
		dynamicx{i in J}:  x[i]   = x[i+1]  - tf/(n-1)/6*(f1[i] + f1[i+1] + 4*f1m[i]);
		dynamicy{i in J}:  y[i]   = y[i+1]  - tf/(n-1)/6*(f2[i] + f2[i+1] + 4*f2m[i]);
		dynamicz{i in J}:  z[i]   = z[i+1]  - tf/(n-1)/6*(f3[i] + f3[i+1] + 4*f3m[i]);
		dynamicdx{i in J}: dx[i]  = dx[i+1] - tf/(n-1)/6*(f4[i] + f4[i+1] + 4*f4m[i]);
		dynamicdy{i in J}: dy[i]  = dy[i+1] - tf/(n-1)/6*(f5[i] + f5[i+1] + 4*f5m[i]);
		dynamicdz{i in J}: dz[i]  = dz[i+1] - tf/(n-1)/6*(f6[i] + f6[i+1] + 4*f6m[i]);
		dynamicm{i in J}:  m[i]   = m[i+1]  - tf/(n-1)/6*(f7[i] + f7[i+1] + 4*f7m[i]);
#--------------------------------------------------------------------------	

#--------------------------------------------------------------------------
#Initial Guess for the Tangential Thrust traj
let {i in I} m[i] := 1;
let {i in I} x[i] :=   x0*cos(tf/(n-1)*(i-1)) - y0*sin(tf/(n-1)*(i-1));
let {i in I} y[i] :=   x0*sin(tf/(n-1)*(i-1)) + y0*cos(tf/(n-1)*(i-1));
let {i in I} dx[i] :=  -x0*sin(tf/(n-1)*(i-1)) - y0*cos(tf/(n-1)*(i-1));
let {i in I} dy[i] :=  x0*cos(tf/(n-1)*(i-1)) - y0*sin(tf/(n-1)*(i-1));
#-----------------------------------------------------------------------


minimize Final: sqrt((dx[n]-dxf)**2 + (dy[n]-dyf)**2 + (dz[n]-dzf)**2) + sqrt((x[n]-xf)**2 + (y[n]-yf)**2 +(z[n]-zf)**2);
solve;

#--------------------------------------------------------------------------
#Print the Initial Guess
let {i in I} ux[i]  := Th[i] * dx[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);
let {i in I} uy[i]  := Th[i] * dy[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);
let {i in I} uz[i]  := Th[i] * dz[i]/(dx[i]**2+dy[i]**2+dz[i]**2)^(1/2);
let {i in J} uxm[i] := Thm[i] * dxm[i]/(dxm[i]**2+dym[i]**2+dzm[i]**2)^(1/2);
let {i in J} uym[i] := Thm[i] * dym[i]/(dxm[i]**2+dym[i]**2+dzm[i]**2)^(1/2);
let {i in J} uzm[i] := Thm[i] * dzm[i]/(dxm[i]**2+dym[i]**2+dzm[i]**2)^(1/2);


for {i in J} {
printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n",
x[i],y[i],z[i],dx[i],dy[i],dz[i],m[i],ux[i],uy[i],uz[i],xm[i],ym[i],zm[i],dxm[i],dym[i],dzm[i],mm[i],uxm[i],uym[i],uzm[i]> out/InitialGuess.out;
}
#------------------------------------------------------------------------

#Print the initial and final times
printf "%17.16e, %17.16e \n", ti/d2u , tF-ti/d2u > out/TimesGuess.out;

#--------------------------------------------------------------------------
#As to calculate the TT initial guess we made use of a simplified dynamic with 
#sipmplified objectives, we need to dispose some variables and constraints
#(note that AMPL will automatically dispose also those variables linked
# to the one listed here


purge Final;
purge f1,f2,f3,f4,f5,f6,f7,f1m,f2m,f3m,f4m,f5m,f6m,f7m,Th,Thm;
purge Th,Thm,InitialMass,InitialPositionx,InitialPositiony,InitialPositionz,InitialVelocityx,InitialVelocityy,InitialVelocityz;
#--------------------------------------------------------------------------

