#Two phases phase transfer (ocp.mod)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and contains the construction of the 
#optimisation problem with two phases
#--------------------------------------------------------------------------

include ocp1.mod;
include ocp2.mod;
purge InitialMass2;

#--------------------------------------------------------------------------
#Here we set the fly-by constraints:

var Vrelin2	=	(dx1[n1]-dxf1)**2 + (dy1[n1]-dyf1)**2 +(dz1[n1]-dzf1)**2;
var Vrelout2	=	(dx2[1] -dx02)**2 + (dy2[1] -dy02)**2 + (dz2[1]-dz02)**2;
var emax	=	1 + rmin/mupla*Vrelin2;
var alfa	=	acos(
			     ((dx1[n1]-dxf1)*(dx2[1]-dx02) + 
			      (dy1[n1]-dyf1)*(dy2[1]-dy02) + 
			      (dz1[n1]-dzf1)*(dz2[1]-dz02)) / sqrt(Vrelin2) / sqrt(Vrelout2)
			);

subject to 
	samemass: m2[1] = m1[n1];
	sametime: ti1 + tf1 = ti2;
	samerelvel: abs(Vrelin2 - Vrelout2) <= 0.00001;
	angularconstraint: alfa <= 2 * asin(1/emax);
#--------------------------------------------------------------------------


minimize QC: 1000 * sum{i in J1}(Thrust1[i]**2+Thrust1[i+1]**2+4*Thrustm1[i]**2) +
 	     1000 * sum{i in J2}(Thrust2[i]**2+Thrust2[i+1]**2+4*Thrustm2[i]**2);

solve;
purge QC;


maximize finalmass: m2[n2];
solve;


#--------------------------------------------------------------------------
#Print the Solution
for {i in J1} {
printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n",
x1[i],y1[i],z1[i],dx1[i],dy1[i],dz1[i],m1[i],ux1[i],uy1[i],uz1[i],xm1[i],ym1[i],zm1[i],dxm1[i],dym1[i],dzm1[i],mm1[i],uxm1[i],uym1[i],uzm1[i]>out/Solution1.out;
}
printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n",
x1[n1],y1[n1],z1[n1],dx1[n1],dy1[n1],dz1[n1],m1[n1],ux1[n1],uy1[n1],uz1[n1]>out/Solution1.out;

for {i in J2} {
printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n",
x2[i],y2[i],z2[i],dx2[i],dy2[i],dz2[i],m2[i],ux2[i],uy2[i],uz2[i],xm2[i],ym2[i],zm2[i],dxm2[i],dym2[i],dzm2[i],mm2[i],uxm2[i],uym2[i],uzm2[i]>out/Solution2.out;
}

printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n",
x2[n2],y2[n2],z2[n2],dx2[n2],dy2[n2],dz2[n2],m2[n2],ux2[n2],uy2[n2],uz2[n2]>out/Solution2.out;

#------------------------------------------------------------------------

#Print the initial and final times
printf "%17.16e, %17.16e \n", ti1/d2u , tF1-ti1/d2u > out/Times1.out;
printf "%17.16e, %17.16e \n", ti2/d2u , tF2-ti2/d2u > out/Times2.out;
#Print the final mass
printf "%17.16e \n", m1[n1] > out/mass1.out;
printf "%17.16e \n", m2[n2] > out/mass2.out;
#Print flyby radius and fly-by DV
printf "%17.16e \n %17.16e \n", mupla/Vrelin2*(1/sin(alfa/2)-1)*R, sqrt((dx1[n1]-dx2[1])**2 +(dy1[n1]-dy2[1])**2+(dz1[n1]-dz2[1])**2)*V  > out/flyby.out;

