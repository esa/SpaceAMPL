#Single phase transfer (ocp.mod)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in ocp.mod and contains definition of the first phase
#--------------------------------------------------------------------------

#Problem definition 
subject to 
	InitialMass1: m1[1] = 1;
	InitialPositionx1: x1[1] = x01;
	InitialPositiony1: y1[1] = y01;
	InitialPositionz1: z1[1] = z01;
	FinalPositionx1  : x1[n1] = xf1;
	FinalPositiony1  : y1[n1] = yf1;
	FinalPositionz1  : z1[n1] = zf1;
	
	ControlMagnitude1{i in I1} : sqrt(ux1[i]**2  + uy1[i]**2  + uz1[i]**2)  <= Tmax;
	ControlMagnitudem1{i in J1}: sqrt(uxm1[i]**2 + uym1[i]**2 + uzm1[i]**2) <= Tmax; 	 		
	
	InitialC3:	 sqrt((dx1[1]-dx01)**2 + (dy1[1]-dy01)**2 + (dz1[1]-dz01)**2) <= sqrt(C3dep);
	
#--------------------------------------------------------------------------		

#--------------------------------------------------------------------------
#Dynamic at the grid points
var f11{i in I1} = dx1[i];
var f21{i in I1} = dy1[i];
var f31{i in I1} = dz1[i];
var f41{i in I1} = -x1[i]/(x1[i]**2+y1[i]**2+z1[i]**2)^(3/2) + ux1[i]/m1[i];
var f51{i in I1} = -y1[i]/(x1[i]**2+y1[i]**2+z1[i]**2)^(3/2) + uy1[i]/m1[i];
var f61{i in I1} = -z1[i]/(x1[i]**2+y1[i]**2+z1[i]**2)^(3/2) + uz1[i]/m1[i];
var f71{i in I1} = -sqrt(ux1[i]**2+uy1[i]**2+uz1[i]**2)/Isp/g0;
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#State definition at mid-points via Simpson interpolation
var xm1{i in J1} = (x1[i]+x1[i+1])/2 + tf1/(n1-1)/8 * (f11[i] - f11[i+1]);
var ym1{i in J1} = (y1[i]+y1[i+1])/2 + tf1/(n1-1)/8 * (f21[i] - f21[i+1]);
var zm1{i in J1} = (z1[i]+z1[i+1])/2 + tf1/(n1-1)/8 * (f31[i] - f31[i+1]);
var dxm1{i in J1} = (dx1[i]+dx1[i+1])/2 + tf1/(n1-1)/8 * (f41[i] - f41[i+1]);
var dym1{i in J1} = (dy1[i]+dy1[i+1])/2 + tf1/(n1-1)/8 * (f51[i] - f51[i+1]);
var dzm1{i in J1} = (dz1[i]+dz1[i+1])/2 + tf1/(n1-1)/8 * (f61[i] - f61[i+1]);
var mm1{i in J1} = (m1[i]+m1[i+1])/2 + tf1/(n1-1)/8 * (f71[i] - f71[i+1]);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Dynamic at the mid-points
var f1m1{i in J1} = dxm1[i];
var f2m1{i in J1} = dym1[i];
var f3m1{i in J1} = dzm1[i];
var f4m1{i in J1} = -xm1[i]/(xm1[i]**2+ym1[i]**2+zm1[i]**2)^(3/2) + uxm1[i]/mm1[i];
var f5m1{i in J1} = -ym1[i]/(xm1[i]**2+ym1[i]**2+zm1[i]**2)^(3/2) + uym1[i]/mm1[i];
var f6m1{i in J1} = -zm1[i]/(xm1[i]**2+ym1[i]**2+zm1[i]**2)^(3/2) + uzm1[i]/mm1[i];
var f7m1{i in J1} = -sqrt(uxm1[i]**2+uym1[i]**2+uzm1[i]**2)/Isp/g0;
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Hermite Formula
subject to 
	dynamicx1{i in J1}:  x1[i]  = x1[i+1]  - tf1/(n1-1)/6*(f11[i] + f11[i+1] + 4*f1m1[i]);
	dynamicy1{i in J1}:  y1[i]  = y1[i+1]  - tf1/(n1-1)/6*(f21[i] + f21[i+1] + 4*f2m1[i]);
	dynamicz1{i in J1}:  z1[i]  = z1[i+1]  - tf1/(n1-1)/6*(f31[i] + f31[i+1] + 4*f3m1[i]);
	dynamicdx1{i in J1}: dx1[i] = dx1[i+1] - tf1/(n1-1)/6*(f41[i] + f41[i+1] + 4*f4m1[i]);
	dynamicdy1{i in J1}: dy1[i] = dy1[i+1] - tf1/(n1-1)/6*(f51[i] + f51[i+1] + 4*f5m1[i]);
	dynamicdz1{i in J1}: dz1[i] = dz1[i+1] - tf1/(n1-1)/6*(f61[i] + f61[i+1] + 4*f6m1[i]);
	dynamicm1{i in J1}:  m1[i]  = m1[i+1]  - tf1/(n1-1)/6*(f71[i] + f71[i+1] + 4*f7m1[i]);
#--------------------------------------------------------------------------	
	
#--------------------------------------------------------------------------	
var Thrust1{i in I1} = sqrt(ux1[i]**2+uy1[i]**2+uz1[i]**2);
var Thrustm1{i in J1} = sqrt(uxm1[i]**2+uym1[i]**2+uzm1[i]**2);

