#Single phase transfer (ephcalc.inc)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and contains the calculations that AMPL
#needs to do in order to locate the position of the departure and target 
#object at the initial and final time
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#This parameter is the scale factor for the time variables
param f		:=  1/100;

#This parameter allows to define intervals where tI and tT are constrained
param tbnd	:=  400;		
					 	
#Initial time  (this is a variable the optimiser plays with)
var timod    :=   tI * d2u * f, <= (tI+tbnd)*d2u*f, >= (tI-tbnd)*d2u*f;	

#Time of flight (this is a variable the optimiser plays with)
var tfmod    :=   tT * d2u * f, <= (tT+tbnd)*d2u*f, >= (tT-tbnd)*d2u*f;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#We here introduce some other time variables that simplifies the following formulae	
var ti	     =   timod /f;			#Initial time non dimensional
var tf	     =   tfmod /f;			#Time of flight non dimensional
var tF	     =   ti/d2u + tf/d2u;		#Arrival time (MJD2000)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Calculations to find the position of the departure object at tI


var M0		=	((el[depart,'M']*d2r + sqrt(1/el[depart,'a']^3) * ((ti/d2u + 51544) - el[depart,'epoch'])*d2u)
mod (2*pi)+2*pi) mod (2*pi);

var E0 >=0, <=2*pi;
let E0		:=	M0 + el[depart,'e']*cos(M0);  	#Initial guess for Kepler equation

var theta0	=	2*atan(sqrt((1+el[depart,'e'])/(1-el[depart,'e']))*tan(E0/2));

var gamma0	=	atan( el[depart,'e'] * sin(theta0) / (1 + el[depart,'e']*cos(theta0)) );

var r0		=	(el[depart,'a'] * (1 - el[depart,'e']^2)) / (1 + el[depart,'e'] * cos(theta0));

var v0		=	sqrt(2/r0 - 1/el[depart,'a']);
	
var x0		=	r0 * ( cos(theta0+el[depart,'om']*d2r)*cos(el[depart,'Om']*d2r) -
				sin(theta0+el[depart,'om']*d2r)*cos(el[depart,'inc']*d2r)*sin(el[depart,'Om']*d2r) );
				
var y0		=	r0 * ( cos(theta0+el[depart,'om']*d2r)*sin(el[depart,'Om']*d2r) +
				sin(theta0+el[depart,'om']*d2r)*cos(el[depart,'inc']*d2r)*cos(el[depart,'Om']*d2r) );
			
var z0		=	r0 * ( sin(theta0+el[depart,'om']*d2r) * sin(el[depart,'inc']*d2r) );
	
var dx0		=	v0 * (-sin(theta0+el[depart,'om']*d2r-gamma0)*cos(el[depart,'Om']*d2r)
				-cos(theta0+el[depart,'om']*d2r-gamma0)*cos(el[depart,'inc']*d2r)*sin(el[depart,'Om']*d2r));
				
var dy0		=	v0 * (-sin(theta0+el[depart,'om']*d2r-gamma0)*sin(el[depart,'Om']*d2r)
				+cos(theta0+el[depart,'om']*d2r-gamma0)*cos(el[depart,'inc']*d2r)*cos(el[depart,'Om']*d2r));
				
var dz0		=	v0 * (cos(theta0+el[depart,'om']*d2r-gamma0)*sin(el[depart,'inc']*d2r));
	
subject to
	KeplerEquation0: 1 * (M0 - E0 + el[depart,'e']*sin(E0)) = 0;
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
param m0 := 1;

#--------------------------------------------------------------------------
#Print the non dimensional units
printf "%17.16e %17.16e %17.16e \n ",R,V,M > out/units.out;

#------------------------------------------------------------------------ 


