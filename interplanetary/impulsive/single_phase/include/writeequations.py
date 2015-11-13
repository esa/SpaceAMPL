import sys;

file = open("equations.inc","w")
file2 = open("writeinitialguess.inc","w")
file3 = open("writesolution.inc","w")
file4 = open("guesstangential.inc","w")
n=int(sys.argv[1]);


file.write("#------------------------------------------------------------------------\n")
file.write("#Optimisation Variables\n\n")

file.write("#Impulsive DVs\n")
file.write("var ux{i in 2..n-1};\n")
file.write("var uy{i in 2..n-1};\n")
file.write("var uz{i in 2..n-1};\n")
file.write("var uT{i in 2..n-1} = sqrt(ux[i]**2+uy[i]**2+uz[i]**2);\n\n")

file.write("#Starting VINF\n")
file.write("var VINFx:=0.0001;\n")
file.write("var VINFy:=0.0001;\n")
file.write("var VINFz:=0.0001;\n")
file.write("var VINF = sqrt(VINFx^2+VINFy^2+VINFz^2);\n\n")

file.write("#Ending VINF\n")
file.write("var VINFxf:=0.0001;\n")
file.write("var VINFyf:=0.0001;\n")
file.write("var VINFzf:=0.0001;\n")
file.write("var VINFf = sqrt(VINFxf^2+VINFyf^2+VINFzf^2);\n\n")

file.write("#Eccentric Anomaly Differences between nodes\n")
file.write("var DE{i in J};\n\n")

file.write("#Initial time\n")
file.write("var timod    :=   tI * d2u * f, <= (tI+tbnd)*d2u*f, >= (tI-tbnd)*d2u*f;	\n")
file.write("#Time of flight \n")
file.write("var tfmod    :=   tT * d2u * f, <= (tT+tbnd)*d2u*f, >= (tT-tbnd)*d2u*f;	\n")
file.write("#--------------------------------------------------------------------------\n\n")

file.write("#--------------------------------------------------------------------------\n")
file.write("#We here introduce some time variables that simplifies the formulas	\n")
file.write("var ti	     =   timod /f;			#Initial time non dimensional\n")
file.write("var tf	     =   tfmod /f;			#Time of flight non dimensional\n")
file.write("var tF	     =   ti/d2u + tf/d2u;	#Arrival time (MJD2000)\n")
file.write("var dt = tf/(n-1);					#Inter-node temporal separation\n")
file.write("#--------------------------------------------------------------------------\n\n")

file.write("#--------------------------------------------------------------------------\n")
file.write("#Planet ephemerides are set and evaluated in tI, tI+tT\n")
file.write("include include/ephcalc.inc;\n")
file.write("fix timod;\n")
file.write("fix tfmod;\n")
file.write("solve;\n")
file.write("unfix timod;\n")
file.write("unfix tfmod;\n")
file.write("#--------------------------------------------------------------------------\n\n\n\n")


file.write("#--------------------------------------------------------------------------\n")
file.write("# Node 1: Departure Node\n")
file.write("var x1 =  x0;\n")
file.write("var y1 =  y0;\n")
file.write("var z1 =  z0;\n")
file.write("var dx1 =  dx0 + VINFx;\n")
file.write("var dy1 =  dy0 + VINFy;\n")
file.write("var dz1 =  dz0 + VINFz;\n\n")

file.write("#Basic definitions\n")
file.write("var r1 = sqrt(x1^2+y1^2+z1^2);\n")
file.write("var v1 = sqrt(dx1^2+dy1^2+dz1^2);\n")
file.write("var a1 = 1 / (2/r1 - v1^2);\n")
file.write("var sigma1 = x1*dx1+y1*dy1+z1*dz1;\n")
file.write("var meanmotion1 = sqrt(1/a1^3);\n")
file.write("var DM1 = meanmotion1 * dt/2;\n\n")

file.write("#Lagrange Coefficients\n")
file.write("var rvar1 = a1 + (r1-a1)*cos(DE[1]) + sigma1*sqrt(a1)*sin(DE[1]);\n")
file.write("var F1 = 1 - a1/r1 * (1-cos(DE[1]));\n")
file.write("var G1 = a1*sigma1*(1-cos(DE[1])) + r1*sqrt(a1)*sin(DE[1]);\n")
file.write("var Ft1 = -sqrt(a1)/(r1*rvar1) * sin(DE[1]);\n")
file.write("var Gt1 = 1 - a1/rvar1*(1-cos(DE[1]));\n\n")

file.write("subject to KeplerEquations1: \n")
file.write("	DM1 - DE[1] - sigma1/sqrt(a1) * (1 - cos(DE[1])) + (1 - r1/a1)*sin(DE[1]) = 0;\n")
file.write("#--------------------------------------------------------------------------\n\n")


for i in range(2,n-1):
	file.write("#--------------------------------------------------------------------------\n")
	file.write("# Node " +str(i)+"\n")
	file.write("var x" +str(i)+" = F" +str(i-1)+"*x" +str(i-1)+" + G" +str(i-1)+"*dx" +str(i-1)+";\n")
	file.write("var y" +str(i)+" = F" +str(i-1)+"*y" +str(i-1)+" + G" +str(i-1)+"*dy" +str(i-1)+";\n")
	file.write("var z" +str(i)+" = F" +str(i-1)+"*z" +str(i-1)+" + G" +str(i-1)+"*dz" +str(i-1)+";\n")
	file.write("var dx" +str(i)+" = Ft" +str(i-1)+"*x" +str(i-1)+" + Gt" +str(i-1)+"*dx" +str(i-1)+" + ux[" +str(i)+"];\n")
	file.write("var dy" +str(i)+" = Ft" +str(i-1)+"*y" +str(i-1)+" + Gt" +str(i-1)+"*dy" +str(i-1)+" + uy[" +str(i)+"];\n")
	file.write("var dz" +str(i)+" = Ft" +str(i-1)+"*z" +str(i-1)+" + Gt" +str(i-1)+"*dz" +str(i-1)+" + uz[" +str(i)+"];\n\n")

	file.write("#Basic definitions\n")
	file.write("var r" +str(i)+" = sqrt(x" +str(i)+"^2+y" +str(i)+"^2+z" +str(i)+"^2);\n")
	file.write("var v" +str(i)+" = sqrt(dx" +str(i)+"^2+dy" +str(i)+"^2+dz" +str(i)+"^2);\n")
	file.write("var a" +str(i)+" = 1 / (2/r" +str(i)+" - v" +str(i)+"^2);\n")
	file.write("var sigma" +str(i)+" = x" +str(i)+"*dx" +str(i)+"+y" +str(i)+"*dy" +str(i)+"+z" +str(i)+"*dz" +str(i)+";\n")
	file.write("var meanmotion" +str(i)+" = sqrt(1/a" +str(i)+"^3);\n")
	file.write("var DM" +str(i)+" = meanmotion" +str(i)+" * dt;\n\n")

	file.write("#Lagrange Coefficients\n")
	file.write("var rvar" +str(i)+" = a" +str(i)+" + (r" +str(i)+"-a" +str(i)+")*cos(DE[" +str(i)+"]) + sigma" +str(i)+"*sqrt(a" +str(i)+")*sin(DE[" +str(i)+"]);\n")
	file.write("var F" +str(i)+" = 1 - a" +str(i)+"/r" +str(i)+" * (1-cos(DE[" +str(i)+"]));\n")
	file.write("var G" +str(i)+" = a" +str(i)+"*sigma" +str(i)+"*(1-cos(DE[" +str(i)+"])) + r" +str(i)+"*sqrt(a" +str(i)+")*sin(DE[" +str(i)+"]);\n")
	file.write("var Ft" +str(i)+" = -sqrt(a" +str(i)+")/(r" +str(i)+"*rvar" +str(i)+") * sin(DE[" +str(i)+"]);\n")
	file.write("var Gt" +str(i)+" = 1 - a" +str(i)+"/rvar" +str(i)+"*(1-cos(DE[" +str(i)+"]));\n\n")

	
	file.write("subject to KeplerEquations" +str(i)+": \n")
	file.write("	DM" +str(i)+" - DE[" +str(i)+"] - sigma" +str(i)+"/sqrt(a" +str(i)+") * (1 - cos(DE[" +str(i)+"])) + (1 - r" +str(i)+"/a" +str(i)+")*sin(DE[" +str(i)+"]) = 0;\n")
	file.write("#--------------------------------------------------------------------------\n\n")

i=n-1
file.write("#--------------------------------------------------------------------------\n")
file.write("# Node " +str(i)+"\n")
file.write("var x" +str(i)+" = F" +str(i-1)+"*x" +str(i-1)+" + G" +str(i-1)+"*dx" +str(i-1)+";\n")
file.write("var y" +str(i)+" = F" +str(i-1)+"*y" +str(i-1)+" + G" +str(i-1)+"*dy" +str(i-1)+";\n")
file.write("var z" +str(i)+" = F" +str(i-1)+"*z" +str(i-1)+" + G" +str(i-1)+"*dz" +str(i-1)+";\n")
file.write("var dx" +str(i)+" = Ft" +str(i-1)+"*x" +str(i-1)+" + Gt" +str(i-1)+"*dx" +str(i-1)+" + ux[" +str(i)+"];\n")
file.write("var dy" +str(i)+" = Ft" +str(i-1)+"*y" +str(i-1)+" + Gt" +str(i-1)+"*dy" +str(i-1)+" + uy[" +str(i)+"];\n")
file.write("var dz" +str(i)+" = Ft" +str(i-1)+"*z" +str(i-1)+" + Gt" +str(i-1)+"*dz" +str(i-1)+" + uz[" +str(i)+"];\n\n")

file.write("#Basic definitions\n")
file.write("var r" +str(i)+" = sqrt(x" +str(i)+"^2+y" +str(i)+"^2+z" +str(i)+"^2);\n")
file.write("var v" +str(i)+" = sqrt(dx" +str(i)+"^2+dy" +str(i)+"^2+dz" +str(i)+"^2);\n")
file.write("var a" +str(i)+" = 1 / (2/r" +str(i)+" - v" +str(i)+"^2);\n")
file.write("var sigma" +str(i)+" = x" +str(i)+"*dx" +str(i)+"+y" +str(i)+"*dy" +str(i)+"+z" +str(i)+"*dz" +str(i)+";\n")
file.write("var meanmotion" +str(i)+" = sqrt(1/a" +str(i)+"^3);\n")
file.write("var DM" +str(i)+" = meanmotion" +str(i)+" * dt/2;\n\n")

file.write("#Lagrange Coefficients\n")
file.write("var rvar" +str(i)+" = a" +str(i)+" + (r" +str(i)+"-a" +str(i)+")*cos(DE[" +str(i)+"]) + sigma" +str(i)+"*sqrt(a" +str(i)+")*sin(DE[" +str(i)+"]);\n")
file.write("var F" +str(i)+" = 1 - a" +str(i)+"/r" +str(i)+" * (1-cos(DE[" +str(i)+"]));\n")
file.write("var G" +str(i)+" = a" +str(i)+"*sigma" +str(i)+"*(1-cos(DE[" +str(i)+"])) + r" +str(i)+"*sqrt(a" +str(i)+")*sin(DE[" +str(i)+"]);\n")
file.write("var Ft" +str(i)+" = -sqrt(a" +str(i)+")/(r" +str(i)+"*rvar" +str(i)+") * sin(DE[" +str(i)+"]);\n")
file.write("var Gt" +str(i)+" = 1 - a" +str(i)+"/rvar" +str(i)+"*(1-cos(DE[" +str(i)+"]));\n\n")

	
file.write("subject to KeplerEquations" +str(i)+": \n")
file.write("	DM" +str(i)+" - DE[" +str(i)+"] - sigma" +str(i)+"/sqrt(a" +str(i)+") * (1 - cos(DE[" +str(i)+"])) + (1 - r" +str(i)+"/a" +str(i)+")*sin(DE[" +str(i)+"]) = 0;\n")
file.write("#--------------------------------------------------------------------------\n\n")
		
	
	
	
file.write("#--------------------------------------------------------------------------\n")
file.write("# Node n: Arrival node\n")
file.write("var xn = F" +str(n-1)+"*x" +str(n-1)+" + G" +str(n-1)+"*dx" +str(n-1)+";\n")
file.write("var yn = F" +str(n-1)+"*y" +str(n-1)+" + G" +str(n-1)+"*dy" +str(n-1)+";\n")
file.write("var zn = F" +str(n-1)+"*z" +str(n-1)+" + G" +str(n-1)+"*dz" +str(n-1)+";\n")
file.write("var dxn = Ft" +str(n-1)+"*x" +str(n-1)+" + Gt" +str(n-1)+"*dx" +str(n-1)+"+ VINFxf;\n")
file.write("var dyn = Ft" +str(n-1)+"*y" +str(n-1)+" + Gt" +str(n-1)+"*dy" +str(n-1)+"+ VINFyf;\n")
file.write("var dzn = Ft" +str(n-1)+"*z" +str(n-1)+" + Gt" +str(n-1)+"*dz" +str(n-1)+"+ VINFzf;\n\n")

file.write("#Basic definitions\n")
file.write("var rn = sqrt(xn^2+yn^2+zn^2);\n")
file.write("var vn = sqrt(dxn^2+dyn^2+dzn^2);\n")
file.write("var an = 1 / (2/rn - vn^2);\n")
file.write("#--------------------------------------------------------------------------\n\n")

file.write("#--------------------------------------------------------------------------\n")
file.write("#Match Constraint\n")
file.write("subject to \n")
file.write("	FinalPositionx  : xn = xf;\n")
file.write("	FinalPositiony  : yn = yf;\n")
file.write("	FinalPositionz  : zn = zf;\n")
file.write("	FinalVelocityx  : dxn = dxf;\n")
file.write("	FinalVelocityy  : dyn = dyf;\n")
file.write("	FinalVelocityz  : dzn = dzf;\n")
file.write("#--------------------------------------------------------------------------\n")

#file2.write("printf \"%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\\n\",x1,y1,z1,dx1,dy1,dz1,1,VINFx,VINFy,VINFz>out/InitialGuess.out;\n")
for i in range(2,n):
	file2.write("printf \"%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\\n\\n\",x"+str(i)+",y"+str(i)+",z"+str(i)+",dx"+str(i)+",dy"+str(i)+",dz"+str(i)+",m["+str(i)+"],ux["+str(i)+"],uy["+str(i)+"],uz["+str(i)+"]>out/InitialGuess.out;\n")
#file2.write("printf \"%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\\n\",xn,yn,zn,dxn,dyn,dzn,m[n-1],VINFxf,VINFyf,VINFzf>out/InitialGuess.out;\n")
file2.write("close out/InitialGuess.out;")

#file3.write("printf \"%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\\n\",x1,y1,z1,dx1,dy1,dz1,1,VINFx,VINFy,VINFz>out/solution.out;\n")
for i in range(2,n):
	file3.write("printf \"%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\\n\\n\",x"+str(i)+",y"+str(i)+",z"+str(i)+",dx"+str(i)+",dy"+str(i)+",dz"+str(i)+",m["+str(i)+"],ux["+str(i)+"],uy["+str(i)+"],uz["+str(i)+"]>out/solution.out;\n")
#file3.write("printf \"%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\\n\",xn,yn,zn,dxn,dyn,dzn,m[n-1],VINFxf,VINFyf,VINFzf>out/solution.out;\n")
file3.write("close out/solution.out;")

file4.write("let {i in 2..n-1} ux[i]:=Tmax*0.0000001;\n")
file4.write("let {i in 2..n-1} uy[i]:=Tmax*0.0000001;\n")
file4.write("let {i in 2..n-1} uz[i]:=Tmax*0.0000001;\n\n")

#Tangentialguess
file4.write("#--------------------------------------------------------------------------\n")
file4.write("#Initial Guess for the DE variables\n")
file4.write("let  {i in J} DE[i] := DM1;\n")
file4.write("#-----------------------------------------------------------------------\n\n")

for i in range(2,n-1):
	file4.write("let ux["+str(i)+"]:=dx"+str(i)+"/v"+str(i)+"*Tmax/2* tf/(n-1);\n")
	file4.write("let uy["+str(i)+"]:=dy"+str(i)+"/v"+str(i)+"*Tmax/2* tf/(n-1);\n")
	file4.write("let uz["+str(i)+"]:=dz"+str(i)+"/v"+str(i)+"*Tmax/2* tf/(n-1);\n")



file4.write("subject to\n")
file4.write("	thrustON{i in 2..n-1}: uT[i] <= Tmax*tf/(n-1);\n\n")
	
file4.write("minimize\n")
file4.write("	position: (xf-xn)^2+(yf-yn)^2+(zf-zn)^2+(dxf-dxn)^2+(dyf-dyn)^2+(dzf-dzn)^2;\n\n")

file4.write("drop FinalPositionx;\n")
file4.write("drop FinalPositiony;\n")
file4.write("drop FinalPositionz;\n")
file4.write("drop FinalVelocityx;\n")
file4.write("drop FinalVelocityy;\n")
file4.write("drop FinalVelocityz;\n")

file4.write("#--------------------------------------------------------------------------\n")
file4.write("solve;\n")
file4.write("#-----------------------------------------------------------------------\n")

file4.write("#--------------------------------------------------------------------------\n")
file4.write("#Print The Initial Guess x,y,z,dx,dy,dz,m,ux,uy,uz variables\n\n")

file4.write("param m{i in I} := 1;\n")
file4.write("include include/writeinitialguess.inc;\n")
file4.write("purge m;\n\n")

file4.write("#Print the initial and final times\n")
file4.write("printf \"%17.16e, %17.16e \\n\", ti/d2u , tF-ti/d2u > out/TimesGuess.out;\n")
file4.write("close out/TimesGuess.out;\n")
file4.write("#------------------------------------------------------------------------\n\n")

file4.write("#--------------------------------------------------------------------------\n")
file4.write("#Clean up\n")
file4.write("unfix timod;\n")
file4.write("unfix tfmod;\n")
file4.write("restore FinalPositionx;\n")
file4.write("restore FinalPositiony;\n")
file4.write("restore FinalPositionz;\n")
file4.write("restore FinalVelocityx;\n")
file4.write("restore FinalVelocityy;\n")
file4.write("restore FinalVelocityz;\n")
file4.write("drop thrustON;\n")
file4.write("drop position;\n")
