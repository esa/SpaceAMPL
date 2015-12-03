# Atmospheric Entry Problem - Non-Rotating Earth - Attitude Rate Controls
# Written by Alejandro Gonzalez Puerta (November 2015)

# AMPL Parameters
param n   			:= 20;			# Number of Nodes

# Constants
param pi			:= acos(-1);					# Definition of PI
param mu_E 			:= 3.986e14; 					# [m^3/s^2]		- Earth Gravitational Parameter

# Non-Dimentional Units
param length 		:= 6378136; 					# [m] 			- Earth Radius 
param time 			:= 2*pi*sqrt(length^3/mu_E); 	# [s] 			- Time
param velocity		:= sqrt(mu_E/length);			# [m/s]			- Velocity
param density		:= 1.225;						# [kg/m3]		- Sea Level Density	
param mass 			:= density*length^3;          	# [kg]          - Reference Mass
param acceleration	:= mu_E/length^2;				# [m/s^2]		- Gravitational Acceleration 
param force 		:= mass*acceleration;			# [N]

# Body Parameters (Orion Capsule)
param m 			:= (8913)/mass;					# [-] 			- Vehicle Mass
param d            	:= 5/length; 					# [-]			- Capsule Diameter
param S            	:= pi*d^2/4;           		 	# [-]           - Reference Area
param RN           	:= 6.035/length;                # [-]           - Nose Radius

# Entry Conditions
param he          	:= (100e3)/length;              # [-]           - Altitude
param long_e     	:= 0;               			# [rad]         - Longitude
param lat_e     	:= 0;               			# [rad]         - Latitude
param ve          	:= (8e3)/velocity;              # [-]         	- Velocity
param gamma_e     	:= - 5*pi/180;               	# [rad]         - Flight Path Angle
param psi_e     	:= -45*pi/180;               	# [rad]         - Heading Angle

# Control Constants
param bank_e 		:= 0; 							# [rad]         - Bank Angle
param max_alpharate := 10*pi/180*time;				# [rad/s]       - Maximum Angle of Attack Time Rate
param max_sigmarate := 30*pi/180*time;				# [rad/s]       - Maximum Bank Angle Time Rate
param alpha_zero 	:= 0; 							# [rad] 		- Entry Angle of Attack
param sigma_zero 	:= 0; 							# [rad] 		- Entry Bank Angle
param min_alpha 	:= -40*pi/180; 				    # [rad] 		- Minimum Angle of Attack

# Other Constants
param g0   			:= 1; 							# [-]			- Gravitational Acceleration
param R_E 			:= 1; 							# [-] 			- Earth Radius  
param SLdensity		:= 1; 							# [-] 			- Base Density
param H 			:= 287*250/9.81/length;			# [-] 			- Atmospheric Scale Height
param tfinal		:= 400/time; 					# [-] 			- Final Time Guess
param heatflux_C    := 1.7415e-4; 					# [-] 			- Heat Flux Constant
#---------------------------------------------------------------------------------------------------


# Sets
set I   			:= 1..n;
set J   			:= 1..n-1;
set G   			:= 2..n;
#---------------------------------------------------------------------------------------------------


# Variables
var     tf >= 0;
var     uu;

var    	 r{i in I};
var  theta{i in I};
var    phi{i in I};
var  	 V{i in I};
var	 gamma{i in I};
var	   psi{i in I};

var 	f1{i in I}; 
var 	f2{i in I}; 
var 	f3{i in I};
var 	f4{i in I};
var 	f5{i in I};
var 	f6{i in I};
var 	f7{i in I};
var 	f8{i in I};
 
var	    u1{i in I};
var	    u2{i in I};
var	    u3{i in I};

var  	rm{i in J};
var thetam{i in J};
var   phim{i in J};
var 	Vm{i in J};
var	gammam{i in J};
var	  psim{i in J};

var	   f1m{i in J}; 
var	   f2m{i in J}; 
var	   f3m{i in J};
var	   f4m{i in J};
var	   f5m{i in J};
var	   f6m{i in J};
var	   f7m{i in J};
var	   f8m{i in J};

var	   u1m{i in J};
var	   u2m{i in J};
var	   u3m{i in J};

# Supporting Variables
var   alpha{i in I};
var  alpham{i in J};

var   sigma{i in I};
var  sigmam{i in J};

var 	 h{i in I} 	=  r[i] - R_E;
var 	hm{i in J} 	= rm[i] - R_E;

var    rho{i in I}	= SLdensity*exp( -h[i]/H);
var   rhom{i in J} 	= SLdensity*exp(-hm[i]/H);

var     Cd{i in I} 	= - 0.0718 *  alpha[i]^5 - 0.6771 *  alpha[i]^4 - 2.1608 *  alpha[i]^3 - 2.4253 *  alpha[i]^2 - 0.1581 *  alpha[i] + 1.7072;
var    Cdm{i in J} 	= - 0.0718 * alpham[i]^5 - 0.6771 * alpham[i]^4 - 2.1608 * alpham[i]^3 - 2.4253 * alpham[i]^2 - 0.1581 * alpham[i] + 1.7072;

var     Cl{i in I} 	= 0.0777 *  alpha[i]^6  + 0.7629 *  alpha[i]^5 + 2.6522 *  alpha[i]^4 + 3.6395 *  alpha[i]^3 + 1.0828 *  alpha[i]^2 - 0.8108 *  alpha[i] + 0.0037;
var    Clm{i in J} 	= 0.0777 * alpham[i]^6  + 0.7629 * alpham[i]^5 + 2.6522 * alpham[i]^4 + 3.6395 * alpham[i]^3 + 1.0828 * alpham[i]^2 - 0.8108 * alpham[i] + 0.0037;

var    ClCd{i in I} = Cl[i]/Cd[i];
var   ClCdm{i in J} = Clm[i]/Cdm[i];

var   heatflux{i in I} = heatflux_C*sqrt( rho[i]/RN)* V[i]^3;
var  heatfluxm{i in J} = heatflux_C*sqrt(rhom[i]/RN)*Vm[i]^3;

var      q{i in I} 	= 0.5*rho[i]*V[i]^2;
var      L{i in I} 	= q[i]*Cl[i]*S;
var      D{i in I} 	= q[i]*Cd[i]*S;

var     qm{i in J} 	=  0.5*rhom[i]*Vm[i]^2;
var     Lm{i in J} 	= qm[i]*Clm[i]*S;
var     Dm{i in J} 	= qm[i]*Cdm[i]*S;

var   aeroload{i in I} = sqrt( L[i]^2 +  D[i]^2)/m/g0;
var  aeroloadm{i in J} = sqrt(Lm[i]^2 + Dm[i]^2)/m/g0;
#---------------------------------------------------------------------------------------------------


# Re-Entry Dynamics (Nodes)
subject to
	dynamics1{i in I}: f1[i] =  V[i]*sin(gamma[i]);
	dynamics2{i in I}: f2[i] =  V[i]*cos(gamma[i])*sin(psi[i])/r[i]/cos(phi[i]);
	dynamics3{i in I}: f3[i] =  V[i]*cos(gamma[i])*cos(psi[i])/r[i];
	dynamics4{i in I}: f4[i] = -D[i]/m - g0*sin(gamma[i]); 
	dynamics5{i in I}: f5[i] =  L[i]*cos(sigma[i])/V[i]/m - g0*cos(gamma[i])/V[i] + V[i]*cos(gamma[i])/r[i];
	dynamics6{i in I}: f6[i] =  L[i]*sin(sigma[i])/V[i]/cos(gamma[i])/m - V[i]*cos(gamma[i])*cos(psi[i])*tan(phi[i])/r[i];
	dynamics7{i in I}: f7[i] =  u1[i];
	dynamics8{i in I}: f8[i] =  u2[i];
#---------------------------------------------------------------------------------------------------


# Mid-Point State (Hermite)
subject to
	Hermite1{i in J}:
			rm[i] = (r[i]+r[i+1])/2 				+ tf/(n-1)/8 * (f1[i]-f1[i+1]);
	Hermite2{i in J}:
		thetam[i] = (theta[i]+theta[i+1])/2 		+ tf/(n-1)/8 * (f2[i]-f2[i+1]);
	Hermite3{i in J}:
		  phim[i] = (phi[i]+phi[i+1])/2 			+ tf/(n-1)/8 * (f3[i]-f3[i+1]);
	Hermite4{i in J}:
			Vm[i] = (V[i]+V[i+1])/2 				+ tf/(n-1)/8 * (f4[i]-f4[i+1]);
	Hermite5{i in J}:
		gammam[i] = (gamma[i]+gamma[i+1])/2 		+ tf/(n-1)/8 * (f5[i]-f5[i+1]);
	Hermite6{i in J}:
		  psim[i] = (psi[i]+psi[i+1])/2 			+ tf/(n-1)/8 * (f6[i]-f6[i+1]);
	Hermite7{i in J}:
	    alpham[i] = (alpha[i]+alpha[i+1])/2 		+ tf/(n-1)/8 * (f7[i]-f7[i+1]);
	Hermite8{i in J}:
	    sigmam[i] = (sigma[i]+sigma[i+1])/2 		+ tf/(n-1)/8 * (f8[i]-f8[i+1]);
#---------------------------------------------------------------------------------------------------


# Re-Entry Dynamics (Mid-Point)
subject to 
	dynamics1m{i in J}: f1m[i] =  Vm[i]*sin(gammam[i]);
	dynamics2m{i in J}: f2m[i] =  Vm[i]*cos(gammam[i])*sin(psim[i])/rm[i]/cos(phim[i]);
	dynamics3m{i in J}: f3m[i] =  Vm[i]*cos(gammam[i])*cos(psim[i])/rm[i];
	dynamics4m{i in J}: f4m[i] = -Dm[i]/m - g0*sin(gammam[i]); 
	dynamics5m{i in J}: f5m[i] =  Lm[i]*cos(sigmam[i])/Vm[i]/m - g0*cos(gammam[i])/Vm[i] + Vm[i]*cos(gammam[i])/rm[i];
	dynamics6m{i in J}: f6m[i] =  Lm[i]*sin(sigmam[i])/Vm[i]/cos(gammam[i])/m - Vm[i]*cos(gammam[i])*cos(psim[i])*tan(phim[i])/rm[i];
	dynamics7m{i in J}: f7m[i] =  u1m[i];
	dynamics8m{i in J}: f8m[i] =  u2m[i];
#---------------------------------------------------------------------------------------------------


# Simpsons Rule
subject to
	SimpsonRule1{i in J}:
			r[i] = r[i+1] 		- tf/(n-1)/6 * (f1[i] + f1[i+1] + 4*f1m[i]);
	SimpsonRule2{i in J}:
		theta[i] = theta[i+1] 	- tf/(n-1)/6 * (f2[i] + f2[i+1] + 4*f2m[i]);
	SimpsonRule3{i in J}:
		  phi[i] = phi[i+1] 	- tf/(n-1)/6 * (f3[i] + f3[i+1] + 4*f3m[i]);
	SimpsonRule4{i in J}:
			V[i] = V[i+1] 		- tf/(n-1)/6 * (f4[i] + f4[i+1] + 4*f4m[i]);
	SimpsonRule5{i in J}:
		gamma[i] = gamma[i+1] 	- tf/(n-1)/6 * (f5[i] + f5[i+1] + 4*f5m[i]);
	SimpsonRule6{i in J}:
		  psi[i] = psi[i+1] 	- tf/(n-1)/6 * (f6[i] + f6[i+1] + 4*f6m[i]);
	SimpsonRule7{i in J}:
	    alpha[i] = alpha[i+1] 	- tf/(n-1)/6 * (f7[i] + f7[i+1] + 4*f7m[i]);
	SimpsonRule8{i in J}:
	    sigma[i] = sigma[i+1] 	- tf/(n-1)/6 * (f8[i] + f8[i+1] + 4*f8m[i]);
#---------------------------------------------------------------------------------------------------


# Control Constraints
subject to   aoa_limits{i in I}: 	     min_alpha <=  alpha[i] <= 0;		
subject to  aoa_limitsm{i in J}: 	     min_alpha <= alpham[i] <= 0;		

subject to  controls_u1{i in I}:	-max_alpharate <=  	  u1[i] <= max_alpharate;
subject to controls_u1m{i in J}:	-max_alpharate <= 	 u1m[i] <= max_alpharate;

subject to  controls_u2{i in I}:	-max_sigmarate <=  	  u2[i] <= max_sigmarate;
subject to controls_u2m{i in J}:	-max_sigmarate <= 	 u2m[i] <= max_sigmarate;

subject to controls_uu:			  				uu 	 = gamma_e;
#---------------------------------------------------------------------------------------------------


# Initial Conditions
subject to initial_r:			  r[1] = he + R_E;
subject to initial_theta:	  theta[1] = long_e;
subject to initial_phi:		    phi[1] = lat_e;
subject to initial_V:			  V[1] = ve;
subject to initial_gamma:	  gamma[1] = uu;
subject to initial_psi:		  	psi[1] = psi_e;

subject to initial_alpha:     alpha[1] = alpha_zero;
subject to initial_alpham:   alpham[1] = alpha_zero;

subject to initial_sigma:     sigma[1] = sigma_zero;
subject to initial_sigmam:   sigmam[1] = sigma_zero;
#---------------------------------------------------------------------------------------------------


# Limits Conditions
subject to hlim{i in I}:	r[i]  >= R_E;
subject to Vlim{i in I}:	V[i]  >= 0;
subject to final_r:			r[n]   = R_E;
subject to final_V:			V[n]  <= 200/velocity;

#---------------------------------------------------------------------------------------------------

# Initial Guess (Linear decrease in Altitude and Velocity) 
let 	 				r[1]  := he + R_E;
let {i in G} 			r[i]  := (1 - i/n)*he + R_E;
let {i in I}	 	theta[i]  := 0;
let {i in I}	 	  phi[i]  := 0;
let						V[1]  := ve;
let {i in G} 			V[i]  := (1 - i/n)*ve + 10/velocity;
let {i in I}	 	gamma[i]  := uu;
let {i in I}	 	  psi[i]  := psi_e;
let {i in I} 		alpha[i]  := alpha_zero;
let {i in I} 		sigma[i]  := sigma_zero;

let {i in J} 	   	   rm[i]  := (1 - i/n)*he + R_E;
let {i in J}	   thetam[i]  := 0;
let {i in J}	 	 phim[i]  := 0;
let {i in J} 	   	   Vm[i]  := (1 - i/n)*ve + 10/velocity;
let {i in J}	   gammam[i]  := uu;
let {i in J}	 	 psim[i]  := psi_e;
let {i in J} 	   alpham[i]  := alpha_zero;
let {i in J} 	   sigmam[i]  := sigma_zero;


let 	   			  	  uu  := gamma_e;
let {i in I} 	   	   u1[i]  := 0;
let {i in J} 	  	  u1m[i]  := 0;
let {i in I} 	   	   u2[i]  := 0;
let {i in J} 	  	  u2m[i]  := 0;


let tf 			     	  := tfinal;
#---------------------------------------------------------------------------------------------------


# Objective Function
#minimize 			cum_acc: sum{i in I} f4[i];
minimize  aerodynamic_load: sum{i in I} aeroload[i];
#---------------------------------------------------------------------------------------------------

	
# Solver Options
option solver snopt;
option substout 1;
option show_stats 1;
options snopt_options "outlev=2 Major_iterations=7000 Major_optimality_tolerance=1.0e-9 Major_feasibility_tolerance=1e-11";
solve;
#---------------------------------------------------------------------------------------------------

display{i in I}		  					   length*h[i];
display{i in I}		  					   	  theta[i];
display{i in I}		  					        phi[i];
display{i in I}		 	 	  			 velocity*V[i];
display{i in I}	   		 		   			  gamma[i];
display{i in I}	   		 		   			  	psi[i];
display{i in I}	   			 		    density*rho[i];
display{i in I}			    		 180/pi/time*u1[i];
display{i in I}			    		 180/pi/time*u2[i];
display{i in I}							   sqrt(density)*velocity^3/sqrt(length)*heatflux[i];

display 180/pi*uu;
display time*tf;


# Print
for {i in I}{
printf "%e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e \n",length*r[i],theta[i],phi[i],velocity*V[i],gamma[i],psi[i],1/time*u1[i],1/time*u2[i],uu,density*rho[i],force*L[i],force*D[i],aeroload[i],alpha[i],sigma[i],n,time*tf > State_rate.out;
}