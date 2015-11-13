#Single phase transfer (ephcalcfixed.inc)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file may be included in main.mod and contains the calculations that AMPL
#needs to do in order to start from a given fixed state and reach the final
#planet. It is used to refine the precision of the optimiser
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#This parameter is the scale factor for the time variables
param f		:=  1/100;
						 	
#Initial time  (this is a variable the optimiser plays with)
var timod;	

#Time of flight (this is a variable the optimiser plays with)
var tfmod    :=   tT * d2u * f;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#We here introduce some other time variables that simplifies the following formulae	
var ti	     =   timod /f;			#Initial time non dimensional
var tf	     =   tfmod /f;			#Time of flight non dimensional
var tF	     =   ti/d2u + tf/d2u;		#Arrival time (MJD2000)
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
#Starting point as given by checkprecision.m

param x0 :=0.5755104919256089;
param y0 :=-0.8365325972973815;
param z0 :=0.0000332069863071;
param dx0 :=0.8046751018101274;
param dy0 :=0.5594674551485674;
param dz0 :=-0.0000899621760164;
param m0 :=0.9780025444928832;
let timod :=10070.0170130085607525 * d2u * f;
let tfmod :=30 * d2u * f;

fix   timod;
#--------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Calculations to find the position of the target object at tF
	

var Mf		=	((el[target,'M']*d2r + sqrt(1/el[target,'a']^3) * ((tF + 51544) - el[target,'epoch'])*d2u)
mod (2*pi)+2*pi) mod (2*pi);

var Ef >=0, <=2*pi;
let Ef		:=	Mf + el[target,'e']*cos(Mf);  	#Initial guess for Kepler equation

var thetaf	=	2*atan(sqrt((1+el[target,'e'])/(1-el[target,'e']))*tan(Ef/2));

var gammaf	=	atan( el[target,'e'] * sin(thetaf) / (1 + el[target,'e']*cos(thetaf)) );

var rf		=	(el[target,'a'] * (1 - el[target,'e']^2)) / (1 + el[target,'e'] * cos(thetaf));

var vf		=	sqrt(2/rf - 1/el[target,'a']);
	
var xf		=	rf * ( cos(thetaf+el[target,'om']*d2r)*cos(el[target,'Om']*d2r) -
				sin(thetaf+el[target,'om']*d2r)*cos(el[target,'inc']*d2r)*sin(el[target,'Om']*d2r) );
				
var yf		=	rf * ( cos(thetaf+el[target,'om']*d2r)*sin(el[target,'Om']*d2r) +
				sin(thetaf+el[target,'om']*d2r)*cos(el[target,'inc']*d2r)*cos(el[target,'Om']*d2r) );
			
var zf		=	rf * ( sin(thetaf+el[target,'om']*d2r) * sin(el[target,'inc']*d2r) );
	
var dxf		=	vf * (-sin(thetaf+el[target,'om']*d2r-gammaf)*cos(el[target,'Om']*d2r)
				-cos(thetaf+el[target,'om']*d2r-gammaf)*cos(el[target,'inc']*d2r)*sin(el[target,'Om']*d2r));
				
var dyf		=	vf * (-sin(thetaf+el[target,'om']*d2r-gammaf)*sin(el[target,'Om']*d2r)
				+cos(thetaf+el[target,'om']*d2r-gammaf)*cos(el[target,'inc']*d2r)*cos(el[target,'Om']*d2r));
				
var dzf		=	vf * (cos(thetaf+el[target,'om']*d2r-gammaf)*sin(el[target,'inc']*d2r));
	
subject to
	KeplerEquationf: 1 * (Mf - Ef + el[target,'e']*sin(Ef)) = 0;
#--------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------
#Print the non dimensional units
printf "%17.16e %17.16e %17.16e \n ",R,V,M > out/units.out;

#------------------------------------------------------------------------ 


