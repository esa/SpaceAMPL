#Single phase transfer (eph.inc)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and contains the calculations that AMPL
#needs to do in order to locate the position of the flyby target at the initial and final time
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#This parameter is the scale factor for the time variables
param factor2	:=  1/100;

#This parameter allows to define intervals where tI and tT are constrained
param tbnd2	:=  250;		
					 	
#Initial time  (this is a variable the optimiser plays with)
var timod2    :=   tI2 * d2u * factor2, <= (tI2+tbnd2)*d2u*factor2, >= (tI2-tbnd2)*d2u*factor2;	

#Time of flight (this is a variable the optimiser plays with)
var tfmod2    :=   tT2 * d2u * factor2, <= (tT2+tbnd2)*d2u*factor2, >= (tT2-tbnd2)*d2u*factor2;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#We here introduce some other time variables that simplifies the following formulae	
var ti2	     =   timod2 /factor2;			#Initial time non dimensional
var tf2	     =   tfmod2 /factor2;			#Time of flight non dimensional
var tF2	     =   ti2/d2u + tf2/d2u;		#Arrival time (MJD2000)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Calculations to find the position of the flyby object at tI2

var M02		=	((el[flyby,'M']*d2r + sqrt(1/el[flyby,'a']^3) * ((ti2 / d2u + 51544) - el[flyby,'epoch'])*d2u)
mod (2*pi)+2*pi) mod (2*pi);

var E02 >=0, <=2*pi;
let E02		:=	M02 + el[flyby,'e']*cos(M02);  	#Initial guess for Kepler equation

var theta02	=	2*atan(sqrt((1+el[flyby,'e'])/(1-el[flyby,'e']))*tan(E02/2));

var gamma02	=	atan( el[flyby,'e'] * sin(theta02) / (1 + el[flyby,'e']*cos(theta02)) );

var r02		=	(el[flyby,'a'] * (1 - el[flyby,'e']^2)) / (1 + el[flyby,'e'] * cos(theta02));

var v02		=	sqrt(2/r02 - 1/el[flyby,'a']);
	
var x02		=	r02 * ( cos(theta02+el[flyby,'om']*d2r)*cos(el[flyby,'Om']*d2r) -
				sin(theta02+el[flyby,'om']*d2r)*cos(el[flyby,'inc']*d2r)*sin(el[flyby,'Om']*d2r) );
				
var y02		=	r02 * ( cos(theta02+el[flyby,'om']*d2r)*sin(el[flyby,'Om']*d2r) +
				sin(theta02+el[flyby,'om']*d2r)*cos(el[flyby,'inc']*d2r)*cos(el[flyby,'Om']*d2r) );
			
var z02		=	r02 * ( sin(theta02+el[flyby,'om']*d2r) * sin(el[flyby,'inc']*d2r) );
	
var dx02	=	v02 * (-sin(theta02+el[flyby,'om']*d2r-gamma02)*cos(el[flyby,'Om']*d2r)
				-cos(theta02+el[flyby,'om']*d2r-gamma02)*cos(el[flyby,'inc']*d2r)*sin(el[flyby,'Om']*d2r));
				
var dy02	=	v02 * (-sin(theta02+el[flyby,'om']*d2r-gamma02)*sin(el[flyby,'Om']*d2r)
				+cos(theta02+el[flyby,'om']*d2r-gamma02)*cos(el[flyby,'inc']*d2r)*cos(el[flyby,'Om']*d2r));
				
var dz02	=	v02 * (cos(theta02+el[flyby,'om']*d2r-gamma02)*sin(el[flyby,'inc']*d2r));
	
subject to
	KeplerEquation02: 1*(M02 - E02 + el[flyby,'e']*sin(E02)) = 0;
#--------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Calculations to find the position of the flyby object at tF1
	


var Mf2		=	((el[target,'M']*d2r + sqrt(1/el[target,'a']^3) * ((tF2 + 51544) - el[target,'epoch'])*d2u)
mod (2*pi)+2*pi) mod (2*pi);

var Ef2 >=0, <=2*pi;
let Ef2		:=	Mf2 + el[target,'e']*cos(Mf2);  	#Initial guess for Kepler equation

var thetaf2	=	2*atan(sqrt((1+el[target,'e'])/(1-el[target,'e']))*tan(Ef2/2));

var gammaf2	=	atan( el[target,'e'] * sin(thetaf2) / (1 + el[target,'e']*cos(thetaf2)) );

var rf2		=	(el[target,'a'] * (1 - el[target,'e']^2)) / (1 + el[target,'e'] * cos(thetaf2));

var vf2		=	sqrt(2/rf2 - 1/el[target,'a']);
	
var xf2		=	rf2 * ( cos(thetaf2+el[target,'om']*d2r)*cos(el[target,'Om']*d2r) -
				sin(thetaf2+el[target,'om']*d2r)*cos(el[target,'inc']*d2r)*sin(el[target,'Om']*d2r) );
				
var yf2		=	rf2 * ( cos(thetaf2+el[target,'om']*d2r)*sin(el[target,'Om']*d2r) +
				sin(thetaf2+el[target,'om']*d2r)*cos(el[target,'inc']*d2r)*cos(el[target,'Om']*d2r) );
			
var zf2		=	rf2 * ( sin(thetaf2+el[target,'om']*d2r) * sin(el[target,'inc']*d2r) );
	
var dxf2	=	vf2 * (-sin(thetaf2+el[target,'om']*d2r-gammaf2)*cos(el[target,'Om']*d2r)
				-cos(thetaf2+el[target,'om']*d2r-gammaf2)*cos(el[target,'inc']*d2r)*sin(el[target,'Om']*d2r));
				
var dyf2	=	vf2 * (-sin(thetaf2+el[target,'om']*d2r-gammaf2)*sin(el[target,'Om']*d2r)
				+cos(thetaf2+el[target,'om']*d2r-gammaf2)*cos(el[target,'inc']*d2r)*cos(el[target,'Om']*d2r));
				
var dzf2	=	vf2 * (cos(thetaf2+el[target,'om']*d2r-gammaf2)*sin(el[target,'inc']*d2r));
	
subject to
	KeplerEquationf2: 1*(Mf2 - Ef2 + el[target,'e']*sin(Ef2)) = 0;
#--------------------------------------------------------------------------------------------------------------------




	
