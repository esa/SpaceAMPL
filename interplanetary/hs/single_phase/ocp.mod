#Single phase transfer (ocp.mod)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and contains the calculations that AMPL
#needs to do in order to solve the Optimal Control Problem. It maybe also
#called from ampl with model ocp.mod; after the objective function or other
#constraints are changed
#--------------------------------------------------------------------------

#Problem definition 
subject to 
	InitialMass: m[1] = m0;
	InitialPositionx: x[1] = x0;
	InitialPositiony: y[1] = y0;
	InitialPositionz: z[1] = z0;
	FinalPositionx  : x[n] = xf;
	FinalPositiony  : y[n] = yf;
	FinalPositionz  : z[n] = zf;
	
	ControlMagnitude{i in I} : sqrt(ux[i]**2  + uy[i]**2  + uz[i]**2)  <= Tmax;
	ControlMagnitudem{i in J}: sqrt(uxm[i]**2 + uym[i]**2 + uzm[i]**2) <= Tmax; 	 		
	
	FinalC3:	 sqrt((dx[n]-dxf)**2 + (dy[n]-dyf)**2 + (dz[n]-dzf)**2) <= sqrt(C3arr);
	InitialC3:	 sqrt((dx[1]-dx0)**2 + (dy[1]-dy0)**2 + (dz[1]-dz0)**2) <= sqrt(C3dep);
	
#--------------------------------------------------------------------------		

#--------------------------------------------------------------------------
#Dynamic at the grid points
var f1{i in I} = dx[i];
var f2{i in I} = dy[i];
var f3{i in I} = dz[i];
var f4{i in I} = -x[i]/(x[i]**2+y[i]**2+z[i]**2)^(3/2) + ux[i]/m[i];
var f5{i in I} = -y[i]/(x[i]**2+y[i]**2+z[i]**2)^(3/2) + uy[i]/m[i];
var f6{i in I} = -z[i]/(x[i]**2+y[i]**2+z[i]**2)^(3/2) + uz[i]/m[i];
var f7{i in I} = -sqrt(ux[i]**2+uy[i]**2+uz[i]**2)/Isp/g0;
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#State definition at mid-points via Simpson interpolation
var xm{i in J} = (x[i]+x[i+1])/2 + tf/(n-1)/8 * (f1[i] - f1[i+1]);
var ym{i in J} = (y[i]+y[i+1])/2 + tf/(n-1)/8 * (f2[i] - f2[i+1]);
var zm{i in J} = (z[i]+z[i+1])/2 + tf/(n-1)/8 * (f3[i] - f3[i+1]);
var dxm{i in J} = (dx[i]+dx[i+1])/2 + tf/(n-1)/8 * (f4[i] - f4[i+1]);
var dym{i in J} = (dy[i]+dy[i+1])/2 + tf/(n-1)/8 * (f5[i] - f5[i+1]);
var dzm{i in J} = (dz[i]+dz[i+1])/2 + tf/(n-1)/8 * (f6[i] - f6[i+1]);
var mm{i in J} = (m[i]+m[i+1])/2 + tf/(n-1)/8 * (f7[i] - f7[i+1]);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Dynamic at the mid-points
var f1m{i in J} = dxm[i];
var f2m{i in J} = dym[i];
var f3m{i in J} = dzm[i];
var f4m{i in J} = -xm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2) + uxm[i]/mm[i];
var f5m{i in J} = -ym[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2) + uym[i]/mm[i];
var f6m{i in J} = -zm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2) + uzm[i]/mm[i];
var f7m{i in J} = -sqrt(uxm[i]**2+uym[i]**2+uzm[i]**2)/Isp/g0;
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Hermite Formula
subject to 
	dynamicx{i in J}:  x[i]  = x[i+1]  - tf/(n-1)/6*(f1[i] + f1[i+1] + 4*f1m[i]);
	dynamicy{i in J}:  y[i]  = y[i+1]  - tf/(n-1)/6*(f2[i] + f2[i+1] + 4*f2m[i]);
	dynamicz{i in J}:  z[i]  = z[i+1]  - tf/(n-1)/6*(f3[i] + f3[i+1] + 4*f3m[i]);
	dynamicdx{i in J}: dx[i] = dx[i+1] - tf/(n-1)/6*(f4[i] + f4[i+1] + 4*f4m[i]);
	dynamicdy{i in J}: dy[i] = dy[i+1] - tf/(n-1)/6*(f5[i] + f5[i+1] + 4*f5m[i]);
	dynamicdz{i in J}: dz[i] = dz[i+1] - tf/(n-1)/6*(f6[i] + f6[i+1] + 4*f6m[i]);
	dynamicm{i in J}:  m[i]  = m[i+1]  - tf/(n-1)/6*(f7[i] + f7[i+1] + 4*f7m[i]);
#--------------------------------------------------------------------------	
	
		

var Thrust{i in I} = sqrt(ux[i]**2+uy[i]**2+uz[i]**2);
var Thrustm{i in J} = sqrt(uxm[i]**2+uym[i]**2+uzm[i]**2);

#minimize obj: tfmod;
maximize obj: m[n];
#minimize obj: 1000 * sum{i in J}(Thrust[i]**2+Thrust[i+1]**2+4*Thrustm[i]**2);


solve;


#--------------------------------------------------------------------------
#Print the Solution
for {i in J} {
printf "%17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e %17.16e\n%17.16e %17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e\n",
x[i], y[i],z[i],dx[i],dy[i],dz[i],m[i],ux[i],uy[i],uz[i],xm[i],ym[i],zm[i],dxm[i],dym[i],dzm[i],mm[i],uxm[i],uym[i],uzm[i]>out/Solution.out;
}
printf "%17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e  %17.16e\n",x[n],y[n],z[n],dx[n],dy[n],dz[n],m[n],ux[n],uy[n],uz[n]>out/Solution.out;
#------------------------------------------------------------------------
close out/Solution.out;

#Print the initial and final times
printf "%17.16e  %17.16e \n", ti/d2u , tF-ti/d2u > out/Times.out;
close out/Times.out;
#Print the final mass
printf "%17.16e \n", m[n] > out/mass.out;
close out/mass.out;


#--------------------------------------------------------------------------
#To be able to run ocp.mod again from a previously obtained solution we
#need to dispose some variables and constraints

purge f1,f2,f3,f4,f5,f6,f7,f1m,f2m,f3m,f4m,f5m,f6m,f7m, Thrust, Thrustm;
purge obj,InitialMass,ControlMagnitude,ControlMagnitudem,FinalC3,InitialC3;
purge FinalPositionx,FinalPositiony,FinalPositionz,InitialPositionx,InitialPositiony,InitialPositionz;
#--------------------------------------------------------------------------

