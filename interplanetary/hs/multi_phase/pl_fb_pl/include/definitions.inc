#Single phase transfer (definitions.inc)
#
#Written by Dario Izzo (October 2007)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and contains definitions of quantities
#that should not be visible to the user. The files also read the content of
#data.dat that should be formatted accordingly 




#Sets-------------------------------------------------------------------
set I1 := 1..n1;
set J1 := 1..n1-1;

set I2 := 1..n2;
set J2 := 1..n2-1;

set OscElem := {'a','e','inc','Om','om','M','epoch'};
set Planet := {'mercury','venus','earth','mars','jupiter','saturn','neptune','uranus','pluto', 
'ast1','ast2','ast3','ast4','ast5','ast6','ast7','ast8','ast9','ast10','ast11','ast12','ast13','ast14','ast15','ast16','ast17','ast18','ast19','ast20','ast21','ast22','ast23','ast24','ast25','ast26','ast27','ast28','ast29','ast30','ast31','ast32','ast33','ast34','ast35','ast36','ast37','ast38','ast39','ast40','ast41','ast42','ast43','ast44','ast45','ast46','ast47','ast48','ast49','ast50','ast51','ast52','ast53','ast54','ast55','ast56','ast57','ast58','ast59','ast60','ast61','ast62','ast63','ast64','ast65','ast66','ast67','ast68','ast69','ast70','ast71','ast72','ast73','ast74','ast75','ast76','ast77','ast78','ast79','ast80','ast81','ast82','ast83','ast84','ast85','ast86','ast87','ast88','ast89','ast90','ast91','ast92','ast93','ast94','ast95','ast96','ast97','ast98','ast99','ast100','ast101','ast102','ast103','ast104','ast105','ast106','ast107','ast108','ast109','ast110','ast111','ast112','ast113','ast114','ast115','ast116','ast117','ast118','ast119','ast120','ast121','ast122','ast123','ast124','ast125','ast126','ast127','ast128','ast129','ast130','ast131','ast132','ast133','ast134','ast135','ast136','ast137','ast138','ast139','ast140','ast141'};

#-----------------------------------------------------------------------

#Parameters-------------------------------------------------------------
param el{Planet,OscElem};		#This contains the data of the planets
read{i in Planet,j in OscElem} el[i,j] < data/data.dat;	#This reads the data of the planets/asteroids

param R:=1.49597870691E8; 		# Length unit (in km)
param MU:=1.32712440018E11;		# Grav. parameter unit(km^3/sec^2)
param V:=sqrt(MU/R);			# km/sec
param T:=R/V;				# sec
param A:=V^2/R;				# km/sec^2

param g0:=9.80665e-3/A;			#non-dimensional gravitational constant
param pi:=acos(-1);

param d2u:=60*60*24/T;			#conversion factors: days to unit 
param d2r:=pi/180;			#conversion factors: degrees to radians 
#------------------------------------------------------------------------

#------------------------------------------------------------------------
#Variables

#Phase 1
var x1{i in I1} <= 1.5, >= -1.5;
var y1{i in I1} <= 1.5, >= -1.5;
var z1{i in I1} <= 1.5, >= -1.5;
var dx1{i in I1};
var dy1{i in I1};
var dz1{i in I1};
var m1{i in I1} <=1;
var ux1{i in I1};
var uy1{i in I1};
var uz1{i in I1};
var uxm1{i in J1};
var uym1{i in J1};
var uzm1{i in J1};

#Phase 2
var x2{i in I2}<= 1.5, >= -1.5;
var y2{i in I2}<= 1.5, >= -1.5;
var z2{i in I2}<= 1.5, >= -1.5;
var dx2{i in I2};
var dy2{i in I2};
var dz2{i in I2};
var m2{i in I2} <=1;
var ux2{i in I2};
var uy2{i in I2};
var uz2{i in I2};
var uxm2{i in J2};
var uym2{i in J2};
var uzm2{i in J2};
#------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Print the orbital parameters of the departure/targets
printf "%8.7e, %8.7e, %8.7e, %8.7e, %8.7e, %8.7e, %8.7e\n%8.7e, %8.7e, %8.7e,%8.7e, %8.7e, %8.7e, %8.7e\n%8.7e, %8.7e, %8.7e, %8.7e, %8.7e, %8.7e, %8.7e\n",
el[depart,'a'],el[depart,'e'],el[depart,'inc'],el[depart,'Om'],el[depart,'om'],el[depart,'M'],el[depart,'epoch'], 
el[flyby,'a'],el[flyby,'e'],el[flyby,'inc'],el[flyby,'Om'],el[flyby,'om'],el[flyby,'M'],el[flyby,'epoch'],
el[target,'a'],el[target,'e'],el[target,'inc'],el[target,'Om'],el[target,'om'],el[target,'M'],el[target,'epoch']> out/PlaParam.out;

#------------------------------------------------------------------------ 

