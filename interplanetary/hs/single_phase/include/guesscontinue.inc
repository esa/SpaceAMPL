#Refinement of a single phase transfer (guesscontinue.inc)
#
#Written by Marco del Rey (November 2007)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and reads Solution.out and uses it
#as an initial guess.

set NUMBER := 1..10;
set LINES := 1..n*2-1;
param dumpmat{LINES,NUMBER};
read{lin in LINES,num in NUMBER} dumpmat[lin,num] < out/Solution.out;
close out/Solution.out;
let {j in J} x[j]:=dumpmat[2*j-1,1];
let {j in J} y[j]:=dumpmat[2*j-1,2];
let {j in J} z[j]:=dumpmat[2*j-1,3];
let {j in J} dx[j]:=dumpmat[2*j-1,4];
let {j in J} dy[j]:=dumpmat[2*j-1,5];
let {j in J} dz[j]:=dumpmat[2*j-1,6];
let {j in J} m[j]:=dumpmat[2*j-1,7];
let {j in J} ux[j]:=dumpmat[2*j-1,8];
let {j in J} uy[j]:=dumpmat[2*j-1,9];
let {j in J} uz[j]:=dumpmat[2*j-1,10];
let {j in J} uxm[j]:=dumpmat[2*j,8];
let {j in J} uym[j]:=dumpmat[2*j,9];
let {j in J} uzm[j]:=dumpmat[2*j,10];
let x[n]:=dumpmat[n*2-1,1];
let y[n]:=dumpmat[n*2-1,2];
let z[n]:=dumpmat[n*2-1,3];
let dx[n]:=dumpmat[n*2-1,4];
let dy[n]:=dumpmat[n*2-1,5];
let dz[n]:=dumpmat[n*2-1,6];
let m[n]:=dumpmat[n*2-1,7];
let ux[n]:=dumpmat[n*2-1,8];
let uy[n]:=dumpmat[n*2-1,9];
let uz[n]:=dumpmat[n*2-1,10];

purge dumpmat, LINES, NUMBER;
#--------------------------------------------------------------------------

