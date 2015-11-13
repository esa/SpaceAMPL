#Single phase transfer (ocp.mod)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in ocp.mod and contains definition of the second phase
#--------------------------------------------------------------------------

#Problem definition 
subject to 
	InitialMass2: m2[1] = 1;
	InitialPositionx2: x2[1] = x02;
	InitialPositiony2: y2[1] = y02;
	InitialPositionz2: z2[1] = z02;
	FinalPositionx2  : x2[n2] = xf2;
	FinalPositiony2  : y2[n2] = yf2;
	FinalPositionz2  : z2[n2] = zf2;
	
	ControlMagnitude2{i in I2} : sqrt(ux2[i]**2  + uy2[i]**2  + uz2[i]**2)  <= Tmax;
	ControlMagnitudem2{i in J2}: sqrt(uxm2[i]**2 + uym2[i]**2 + uzm2[i]**2) <= Tmax; 	 		
	
	FinalC3:	 sqrt((dx2[n2]-dxf2)**2 + (dy2[n2]-dyf2)**2 + (dz2[n2]-dzf2)**2) <= sqrt(C3arr);
	
#--------------------------------------------------------------------------		

#--------------------------------------------------------------------------
#Dynamic at the grid points
var f12{i in I2} = dx2[i];
var f22{i in I2} = dy2[i];
var f32{i in I2} = dz2[i];
var f42{i in I2} = -x2[i]/(x2[i]**2+y2[i]**2+z2[i]**2)^(3/2) + ux2[i]/m2[i];
var f52{i in I2} = -y2[i]/(x2[i]**2+y2[i]**2+z2[i]**2)^(3/2) + uy2[i]/m2[i];
var f62{i in I2} = -z2[i]/(x2[i]**2+y2[i]**2+z2[i]**2)^(3/2) + uz2[i]/m2[i];
var f72{i in I2} = -sqrt(ux2[i]**2+uy2[i]**2+uz2[i]**2)/Isp/g0;
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#State definition at mid-points via Simpson interpolation
var xm2{i in J2} =  (x2[i]+x2[i+1])/2   + tf2/(n2-1)/8 * (f12[i] - f12[i+1]);
var ym2{i in J2} =  (y2[i]+y2[i+1])/2   + tf2/(n2-1)/8 * (f22[i] - f22[i+1]);
var zm2{i in J2} =  (z2[i]+z2[i+1])/2   + tf2/(n2-1)/8 * (f32[i] - f32[i+1]);
var dxm2{i in J2} = (dx2[i]+dx2[i+1])/2 + tf2/(n2-1)/8 * (f42[i] - f42[i+1]);
var dym2{i in J2} = (dy2[i]+dy2[i+1])/2 + tf2/(n2-1)/8 * (f52[i] - f52[i+1]);
var dzm2{i in J2} = (dz2[i]+dz2[i+1])/2 + tf2/(n2-1)/8 * (f62[i] - f62[i+1]);
var mm2{i in J2} =  (m2[i]+m2[i+1])/2   + tf2/(n2-1)/8 * (f72[i] - f72[i+1]);
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Dynamic at the mid-points
var f1m2{i in J2} = dxm2[i];
var f2m2{i in J2} = dym2[i];
var f3m2{i in J2} = dzm2[i];
var f4m2{i in J2} = -xm2[i]/(xm2[i]**2+ym2[i]**2+zm2[i]**2)^(3/2) + uxm2[i]/mm2[i];
var f5m2{i in J2} = -ym2[i]/(xm2[i]**2+ym2[i]**2+zm2[i]**2)^(3/2) + uym2[i]/mm2[i];
var f6m2{i in J2} = -zm2[i]/(xm2[i]**2+ym2[i]**2+zm2[i]**2)^(3/2) + uzm2[i]/mm2[i];
var f7m2{i in J2} = -sqrt(uxm2[i]**2+uym2[i]**2+uzm2[i]**2)/Isp/g0;
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Hermite Formula
subject to 
	dynamicx2{i in J2}:  x2[i]  = x2[i+1]  - tf2/(n2-1)/6*(f12[i] + f12[i+1] + 4*f1m2[i]);
	dynamicy2{i in J2}:  y2[i]  = y2[i+1]  - tf2/(n2-1)/6*(f22[i] + f22[i+1] + 4*f2m2[i]);
	dynamicz2{i in J2}:  z2[i]  = z2[i+1]  - tf2/(n2-1)/6*(f32[i] + f32[i+1] + 4*f3m2[i]);
	dynamicdx2{i in J2}: dx2[i] = dx2[i+1] - tf2/(n2-1)/6*(f42[i] + f42[i+1] + 4*f4m2[i]);
	dynamicdy2{i in J2}: dy2[i] = dy2[i+1] - tf2/(n2-1)/6*(f52[i] + f52[i+1] + 4*f5m2[i]);
	dynamicdz2{i in J2}: dz2[i] = dz2[i+1] - tf2/(n2-1)/6*(f62[i] + f62[i+1] + 4*f6m2[i]);
	dynamicm2{i in J2}:  m2[i]  = m2[i+1]  - tf2/(n2-1)/6*(f72[i] + f72[i+1] + 4*f7m2[i]);
#--------------------------------------------------------------------------	
	
#-----------------------------------------------------------------------
var Thrust2{i in I2} = sqrt(ux2[i]**2+uy2[i]**2+uz2[i]**2);
var Thrustm2{i in J2} = sqrt(uxm2[i]**2+uym2[i]**2+uzm2[i]**2);


