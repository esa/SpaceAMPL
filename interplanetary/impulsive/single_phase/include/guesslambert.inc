#--------------------------------------------------------------------------
#Single phase transfer (guesslambert.inc)
#Written by Dario Izzo (February 2009)
#
#This file may be included in main.mod and generates an initial guess solving a
#Lambert problem defined by the departure and arrival planets and tI and tT
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Problem Definition (only the impulse at the first node is allowed, the rest
#are set to be zero.

let {i in 2..n-1} ux[i]:=Tmax*0.0000001;
let {i in 2..n-1} uy[i]:=Tmax*0.0000001;
let {i in 2..n-1} uz[i]:=Tmax*0.0000001;

fix timod;
fix tfmod;
fix {i in 2..n-1} ux[i];
fix {i in 2..n-1} uy[i];
fix {i in 2..n-1} uz[i];

#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Initial Guess for the DE variables
let  {i in J} DE[i] := DM1;
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
solve;
#-----------------------------------------------------------------------

#--------------------------------------------------------------------------
#Print The Initial Guess x,y,z,dx,dy,dz,m,ux,uy,uz variables

param m{i in I} := 1;
include include/writeinitialguess.inc;
purge m;

#Print the initial and final times
printf "%17.16e, %17.16e \n", ti/d2u , tF-ti/d2u > out/TimesGuess.out;
close out/TimesGuess.out;
#------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Clean up
unfix timod;
unfix tfmod;
unfix ux;
unfix uy;
unfix uz;
#--------------------------------------------------------------------------

