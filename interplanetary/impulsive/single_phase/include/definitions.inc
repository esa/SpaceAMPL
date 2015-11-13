#Single phase transfer (definitions.inc)
#
#Written by Dario Izzo (February 2009)
#
#--------------------------------------------------------------------------
#This file is included in main.mod and contains definitions of quantities
#that should not be visible to the user. The files also read the content of
#data.dat that should be formatted accordingly 
#--------------------------------------------------------------------------



#Sets-------------------------------------------------------------------
set I := 1..n;
set J := 1..n-1;
set OscElem := {'a','e','inc','Om','om','M','epoch'};
set Planet := {'mercury','venus','earth','mars','jupiter','saturn','neptune','uranus','pluto', 
'ast1','ast2','ast3','ast4','ast5','ast6','ast7','ast8','ast9','ast10','ast11','ast12','ast13','ast14','ast15','ast16','ast17','ast18','ast19','ast20','ast21','ast22','ast23','ast24','ast25','ast26','ast27','ast28','ast29','ast30','ast31','ast32','ast33','ast34','ast35','ast36','ast37','ast38','ast39','ast40','ast41','ast42','ast43','ast44','ast45','ast46','ast47','ast48','ast49','ast50','ast51','ast52','ast53','ast54','ast55','ast56','ast57','ast58','ast59','ast60','ast61','ast62','ast63','ast64','ast65','ast66','ast67','ast68','ast69','ast70','ast71','ast72','ast73','ast74','ast75','ast76','ast77','ast78','ast79','ast80','ast81','ast82','ast83','ast84','ast85','ast86','ast87','ast88','ast89','ast90','ast91','ast92','ast93','ast94','ast95','ast96','ast97','ast98','ast99','ast100','ast101','ast102','ast103','ast104','ast105','ast106','ast107','ast108','ast109','ast110','ast111','ast112','ast113','ast114','ast115','ast116','ast117','ast118','ast119','ast120','ast121','ast122','ast123','ast124','ast125','ast126','ast127','ast128','ast129','ast130','ast131','ast132','ast133','ast134','ast135','ast136','ast137','ast138','ast139','ast140','ast141'};
#-----------------------------------------------------------------------

#Parameters-------------------------------------------------------------
param el{Planet,OscElem};		#This contains the data of the planets
read{i in Planet,j in OscElem} el[i,j] < data/data.dat;	#This reads the data of the planets
close data/data.dat;
param R:=1.49597870691E8; 		# Length unit (in km)
param MU:=1.32712440018E11;		# Grav. parameter unit(km^3/sec^2)
param V:=sqrt(MU/R);			# km/sec
param T:=R/V;					# sec
param A:=V^2/R;					# km/sec^2

param g0:=9.80665e-3/A;			#non-dimensional gravitational constant
param pi:=acos(-1);

param d2u:=60*60*24/T;			#conversion factors: days to unit 
param d2r:=pi/180;			#conversion factors: degrees to radians 


param f		:=  1/1;		#This parameter is the scale factor for the time variables
param tbnd	:=  400;		#This parameter allows to define intervals where tI and tT are constrained
#------------------------------------------------------------------------



#--------------------------------------------------------------------------
#Print the orbital parameters of the departure/targets
printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\n%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e",
el[depart,'a'],el[depart,'e'],el[depart,'inc'],el[depart,'Om'],el[depart,'om'],el[depart,'M'],el[depart,'epoch'], 
el[target,'a'],el[target,'e'],el[target,'inc'],el[target,'Om'],el[target,'om'],el[target,'M'],el[target,'epoch'] > out/PlaParam.out;
#------------------------------------------------------------------------ 

