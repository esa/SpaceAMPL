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
param factor1	:=  1/100;

#This parameter allows to define intervals where tI and tT are constrained
param tbnd1	:=  250;		
					 	
#Initial time  (this is a variable the optimiser plays with)
var timod1    :=   tI1 * d2u * factor1, <= (tI1+tbnd1)*d2u*factor1, >= (tI1-tbnd1)*d2u*factor1;	

#Time of flight (this is a variable the optimiser plays with)
var tfmod1    :=   tT1 * d2u * factor1, <= (tT1+tbnd1)*d2u*factor1, >= (tT1-tbnd1)*d2u*factor1;
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#We here introduce some other time variables that simplifies the following formulae	
var ti1	     =   timod1 /factor1;			#Initial time non dimensional
var tf1	     =   tfmod1 /factor1;			#Time of flight non dimensional
var tF1	     =   ti1/d2u + tf1/d2u;			#Arrival time (MJD2000)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Calculations to find the position of the departure object at tI1

var M01		=	((el[depart,'M']*d2r + sqrt(1/el[depart,'a']^3) * ((ti1 / d2u+ 51544) - el[depart,'epoch'])*d2u)
mod (2*pi)+2*pi) mod (2*pi);

var E01 >=0, <=2*pi;
let E01		:=	M01 + el[depart,'e']*cos(M01);  	#Initial guess for Kepler equation

var theta01	=	2*atan(sqrt((1+el[depart,'e'])/(1-el[depart,'e']))*tan(E01/2));

var gamma01	=	atan( el[depart,'e'] * sin(theta01) / (1 + el[depart,'e']*cos(theta01)) );

var r01		=	(el[depart,'a'] * (1 - el[depart,'e']^2)) / (1 + el[depart,'e'] * cos(theta01));

var v01		=	sqrt(2/r01 - 1/el[depart,'a']);
	
var x01		=	r01 * ( cos(theta01+el[depart,'om']*d2r)*cos(el[depart,'Om']*d2r) -
				sin(theta01+el[depart,'om']*d2r)*cos(el[depart,'inc']*d2r)*sin(el[depart,'Om']*d2r) );
				
var y01		=	r01 * ( cos(theta01+el[depart,'om']*d2r)*sin(el[depart,'Om']*d2r) +
				sin(theta01+el[depart,'om']*d2r)*cos(el[depart,'inc']*d2r)*cos(el[depart,'Om']*d2r) );
			
var z01		=	r01 * ( sin(theta01+el[depart,'om']*d2r) * sin(el[depart,'inc']*d2r) );
	
var dx01	=	v01 * (-sin(theta01+el[depart,'om']*d2r-gamma01)*cos(el[depart,'Om']*d2r)
				-cos(theta01+el[depart,'om']*d2r-gamma01)*cos(el[depart,'inc']*d2r)*sin(el[depart,'Om']*d2r));
				
var dy01	=	v01 * (-sin(theta01+el[depart,'om']*d2r-gamma01)*sin(el[depart,'Om']*d2r)
				+cos(theta01+el[depart,'om']*d2r-gamma01)*cos(el[depart,'inc']*d2r)*cos(el[depart,'Om']*d2r));
				
var dz01	=	v01 * (cos(theta01+el[depart,'om']*d2r-gamma01)*sin(el[depart,'inc']*d2r));
	
subject to
	
	KeplerEquation01: 1*(M01 - E01 + el[depart,'e']*sin(E01)) = 0;
#--------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Calculations to find the position of the flyby object at tF1
	
var Mf1		=	((el[flyby,'M']*d2r + sqrt(1/el[flyby,'a']^3) * ((tF1 + 51544) - el[flyby,'epoch'])*d2u)
mod (2*pi)+2*pi) mod (2*pi);

var Ef1 >=0, <=2*pi;
let Ef1		:=	Mf1 + el[flyby,'e']*cos(Mf1);  	#Initial guess for Kepler equation

var thetaf1	=	2*atan(sqrt((1+el[flyby,'e'])/(1-el[flyby,'e']))*tan(Ef1/2));

var gammaf1	=	atan( el[flyby,'e'] * sin(thetaf1) / (1 + el[flyby,'e']*cos(thetaf1)) );

var rf1		=	(el[flyby,'a'] * (1 - el[flyby,'e']^2)) / (1 + el[flyby,'e'] * cos(thetaf1));

var vf1		=	sqrt(2/rf1 - 1/el[flyby,'a']);
	
var xf1		=	rf1 * ( cos(thetaf1+el[flyby,'om']*d2r)*cos(el[flyby,'Om']*d2r) -
				sin(thetaf1+el[flyby,'om']*d2r)*cos(el[flyby,'inc']*d2r)*sin(el[flyby,'Om']*d2r) );
				
var yf1		=	rf1 * ( cos(thetaf1+el[flyby,'om']*d2r)*sin(el[flyby,'Om']*d2r) +
				sin(thetaf1+el[flyby,'om']*d2r)*cos(el[flyby,'inc']*d2r)*cos(el[flyby,'Om']*d2r) );
			
var zf1		=	rf1 * ( sin(thetaf1+el[flyby,'om']*d2r) * sin(el[flyby,'inc']*d2r) );
	
var dxf1	=	vf1 * (-sin(thetaf1+el[flyby,'om']*d2r-gammaf1)*cos(el[flyby,'Om']*d2r)
				-cos(thetaf1+el[flyby,'om']*d2r-gammaf1)*cos(el[flyby,'inc']*d2r)*sin(el[flyby,'Om']*d2r));
				
var dyf1	=	vf1 * (-sin(thetaf1+el[flyby,'om']*d2r-gammaf1)*sin(el[flyby,'Om']*d2r)
				+cos(thetaf1+el[flyby,'om']*d2r-gammaf1)*cos(el[flyby,'inc']*d2r)*cos(el[flyby,'Om']*d2r));
				
var dzf1	=	vf1 * (cos(thetaf1+el[flyby,'om']*d2r-gammaf1)*sin(el[flyby,'inc']*d2r));
	
subject to
	KeplerEquationf1: 1*(Mf1 - Ef1 + el[flyby,'e']*sin(Ef1)) = 0;
#--------------------------------------------------------------------------------------------------------------------



	
