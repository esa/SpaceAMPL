#Generator of a Multiphase AMPL solver
#Created by Marco del Rey Zapatero
#February 2008


#Parameters of the current problem
M=2000
Isp=3000
Tmax=0.15
nodes=20
C3lau=0.5
C3rnv=0.001
days=60
TOTALTIME=3652.6
#Function definitions
def maingen(sequence, flyby,times):


    #Opening main.mod file
    main=open('main.mod','w')
    
    #Writing the description of main.mod
    main.write('#Multiphase transfer (main.mod)\n')
    main.write('#\n')
    main.write('#(February 2008)\n')
    main.write('#\n')
    main.write('#--------------------------------------------------------------------------\n')
    main.write('#This AMPL model transcribe the NEP low-thrust OCP multiphase problem into an\n')
    main.write('#NLP using the Hermite-Simpson Formula. \n')
    main.write('#\n')
    main.write('\n')
    
    #Writing the arrival and departure bodies for each phase
    main.write('#--------------------------------------------------------------------------\n')
    main.write('#Departure and Arrival Choice\n')
    main.write('\n')
    for x in range(len(sequence)-1):
        main.write('param depart'+str(x)+':=\''+sequence[x]+'\' symbolic;\n')
        main.write('param target'+str(x)+':=\''+sequence[x+1]+'\' symbolic;\n')
    

    #Times initially set to 0 (Later will be written in the input.dat file
    main.write('#--------------------------------------------------------------------------\n')
    main.write('#Times (used as initial guesses)\n')
    main.write('\n')
    for x in range(len(sequence)-1):
        main.write('param tI'+str(x)+'	     :=   '+str(times[2*x])+';	        #Departure time (MJD2000)\n')
        main.write('param tT'+str(x)+'	     :=   '+str(times[2*x+1])+';		#Transfer time (MJD2000)\n')
    
    #Grid definition
    main.write('#--------------------------------------------------------------------------\n')
    main.write('#Grid options\n')
    main.write('\n')
    for x in range(len(sequence)-1):
        #Grid depending on timelength
        main.write('param n'+str(x)+':=floor(25+(tT'+str(x)+'-200)*40/1100);\n')
        #Fixed Grid
        #main.write('param n'+str(x)+':='+str(nodes)+';\n')
    
    #Loading of definitions file
    
    main.write('#--------------------------------------------------------------------------\n')
    main.write('#Loading definitions file\n')
    main.write('\n')
    main.write('include include/definitions.inc;\n')
    
    #Spacecraft Design Parameters
    main.write('#--------------------------------------------------------------------------\n')
    main.write('#Spacecraft Design\n')
    main.write('\n')
    main.write('param M      :=   ('+str(M)+');	\n')
    main.write('param Isp    :=   ('+str(Isp)+')   /T;		#Engine specific impulse (sec)\n')
    main.write('param Tmax   :=   ('+str(Tmax)+')   /A/M/1000;	#Maximum Thrust level 	 (N)\n')
    
    #Mission Design Parameters
    main.write('#--------------------------------------------------------------------------\n')
    main.write('#C3 conditions\n')
    main.write('param C3lau     :=   '+str(C3lau)+'**2 /V**2;	#Launch C3 (km^2/sec^2)\n')
    main.write('param C3rnv     :=   '+str(C3rnv)+'**2 /V**2;	#Rendezvous C3 (km^2/sec^2)\n')
  
    
    #Solver options for ephemerides
    main.write('#--------------------------------------------------------------------------\n')
    main.write('#Ephemerides Solver options\n')
    main.write('\n')
    main.write('options snopt_options "outlev=2 Major_iterations=100";\n')
    main.write('options conopt_options "outlev=3";\n')
    main.write('option log_file \'out/log.txt\';\n')
    main.write('options solver snopt;\n')

    #Loading ephemerides files
    main.write('#--------------------------------------------------------------------------\n')
    main.write('#Loading ephemerides files\n')
    main.write('\n')
    for x in range(len(sequence)-1):
        main.write('include include/ephcalc'+str(x)+'.inc;\n')
    
    #Printing the non dimensional units
    main.write('#--------------------------------------------------------------------------\n')
    main.write('#Printing the non dimensional units\n')
    main.write('\n')
    main.write('printf "%17.16e, %17.16e, %17.16e \\n ",R,V,M > out/units.out;\n')



    #Printing snopt options for total solver
    main.write('#------------------------------------------------------------------------ \n')
    main.write('#Regular Solver options\n')
    main.write('\n')
    main.write('options snopt_options "Superbasics_limit=3000 outlev=2 Major_iterations=2000";\n')
    
    #Printing sequence of guesses used
    for x in range(len(sequence)-1):
        main.write('include include/guess'+str(x)+'.inc;\n')
    
    #Final solve
    main.write('include ocp.mod;\n')
    main.close()
    return

def definitionsgen(sequence, flyby):
    definitions=open('./include/definitions.inc','w')
    definitions.write('\n')
    definitions.write('#Single phase transfer (definitions.inc)\n')
    definitions.write('#\n')
    definitions.write('#Written by Dario Izzo (October 2007)\n')
    definitions.write('#\n')
    definitions.write('#--------------------------------------------------------------------------\n')
    definitions.write('#This file is included in main.mod and contains definitions of quantities\n')
    definitions.write('#that should not be visible to the user. The files also read the content of\n')
    definitions.write('#data.dat that should be formatted accordingly \n')
    definitions.write('\n')

    
    definitions.write('#--------------------------------------------------------------------------\n')
    definitions.write('#Sets\n')
    definitions.write('\n')
    for x in range(len(sequence)-1):
        definitions.write('set I'+str(x)+' := 1..n'+str(x)+';\n')
        definitions.write('set J'+str(x)+' := 1..n'+str(x)+'-1;\n')
       
    definitions.write('#--------------------------------------------------------------------------\n')
    definitions.write('#Labels\n')
    definitions.write('\n')
    definitions.write('set OscElem := {\'a\',\'e\',\'inc\',\'Om\',\'om\',\'M\',\'epoch\'};\n')
    definitions.write('set Planet := {\'mercury\',\'venus\',\'earth\',\'mars\',\'jupiter\',\'saturn\',\'neptune\',\'uranus\',\'pluto\', \'ast1\',\'ast2\',\'ast3\',\'ast4\',\'ast5\',\'ast6\',\'ast7\',\'ast8\',\'ast9\',\'ast10\',\'ast11\',\'ast12\',\'ast13\',\'ast14\',\'ast15\',\'ast16\',\'ast17\',\'ast18\',\'ast19\',\'ast20\',\'ast21\',\'ast22\',\'ast23\',\'ast24\',\'ast25\',\'ast26\',\'ast27\',\'ast28\',\'ast29\',\'ast30\',\'ast31\',\'ast32\',\'ast33\',\'ast34\',\'ast35\',\'ast36\',\'ast37\',\'ast38\',\'ast39\',\'ast40\',\'ast41\',\'ast42\',\'ast43\',\'ast44\',\'ast45\',\'ast46\',\'ast47\',\'ast48\',\'ast49\',\'ast50\',\'ast51\',\'ast52\',\'ast53\',\'ast54\',\'ast55\',\'ast56\',\'ast57\',\'ast58\',\'ast59\',\'ast60\',\'ast61\',\'ast62\',\'ast63\',\'ast64\',\'ast65\',\'ast66\',\'ast67\',\'ast68\',\'ast69\',\'ast70\',\'ast71\',\'ast72\',\'ast73\',\'ast74\',\'ast75\',\'ast76\',\'ast77\',\'ast78\',\'ast79\',\'ast80\',\'ast81\',\'ast82\',\'ast83\',\'ast84\',\'ast85\',\'ast86\',\'ast87\',\'ast88\',\'ast89\',\'ast90\',\'ast91\',\'ast92\',\'ast93\',\'ast94\',\'ast95\',\'ast96\',\'ast97\',\'ast98\',\'ast99\',\'ast100\',\'ast101\',\'ast102\',\'ast103\',\'ast104\',\'ast105\',\'ast106\',\'ast107\',\'ast108\',\'ast109\',\'ast110\',\'ast111\',\'ast112\',\'ast113\',\'ast114\',\'ast115\',\'ast116\',\'ast117\',\'ast118\',\'ast119\',\'ast120\',\'ast121\',\'ast122\',\'ast123\',\'ast124\',\'ast125\',\'ast126\',\'ast127\',\'ast128\',\'ast129\',\'ast130\',\'ast131\',\'ast132\',\'ast133\',\'ast134\',\'ast135\',\'ast136\',\'ast137\',\'ast138\',\'ast139\',\'ast140\',\'ast141\'};\n')
                      
       
    definitions.write('#--------------------------------------------------------------------------\n')
    definitions.write('#Parameters\n')
    definitions.write('\n')
    definitions.write('param el{Planet,OscElem};		#This contains the data of the planets\n')
    definitions.write('read{i in Planet,j in OscElem} el[i,j] < data/data.dat;	#This reads the data of the planets/asteroids\n')
    definitions.write('\n')
    definitions.write('param R:=1.49597870691E8; 		# Length unit (in km)\n')
    definitions.write('param MU:=1.32712440018E11;		# Grav. parameter unit(km^3/sec^2)\n')
    definitions.write('param V:=sqrt(MU/R);			# km/sec\n')
    definitions.write('param T:=R/V;				# sec\n')
    definitions.write('param A:=V^2/R;				# km/sec^2\n')
    definitions.write('\n')
    definitions.write('param g0:=9.80665e-3/A;			#non-dimensional gravitational constant\n')
    definitions.write('param pi:=acos(-1);\n')
    definitions.write('\n')
    definitions.write('param d2u:=60*60*24/T;			#conversion factors: days to unit \n')
    definitions.write('param d2r:=pi/180;			#conversion factors: degrees to radians \n')


    definitions.write('#--------------------------------------------------------------------------\n')
    definitions.write('#Fly-by\n')
    definitions.write('param rmin   :=   (6871)   / R;		#Minimum fly-by altitude (km)\n')
    definitions.write('param mupla  :=   (398600) / MU;	#Gravitational parameter of the planet (km^3/sec^2)\n')
    definitions.write('#--------------------------------------------------------------------------\n')


    definitions.write('#------------------------------------------------------------------------\n')
    definitions.write('#Variables\n')
    for x in range(len(sequence)-1):
        definitions.write('#Phase '+str(x)+'\n')
        definitions.write('var x'+str(x)+'{i in I'+str(x)+'} <= 1.5, >= -1.5;\n')
        definitions.write('var y'+str(x)+'{i in I'+str(x)+'} <= 1.5, >= -1.5;\n')
        definitions.write('var z'+str(x)+'{i in I'+str(x)+'} <= 1.5, >= -1.5;\n')
        definitions.write('var dx'+str(x)+'{i in I'+str(x)+'};\n')
        definitions.write('var dy'+str(x)+'{i in I'+str(x)+'};\n')
        definitions.write('var dz'+str(x)+'{i in I'+str(x)+'};\n')
        definitions.write('var m'+str(x)+'{i in I'+str(x)+'} <=1;\n')
        definitions.write('var ux'+str(x)+'{i in I'+str(x)+'};\n')
        definitions.write('var uy'+str(x)+'{i in I'+str(x)+'};\n')
        definitions.write('var uz'+str(x)+'{i in I'+str(x)+'};\n')
        definitions.write('var uxm'+str(x)+'{i in J'+str(x)+'};\n')
        definitions.write('var uym'+str(x)+'{i in J'+str(x)+'};\n')
        definitions.write('var uzm'+str(x)+'{i in J'+str(x)+'};\n')
        definitions.write('\n')
    definitions.write('#------------------------------------------------------------------------\n')
    for x in range(len(sequence)-1):
        definitions.write('#--------------------------------------------------------------------------\n')
        definitions.write('#Print the orbital parameters of the departure/targets\n')
        definitions.write('printf "%8.7e, %8.7e, %8.7e,%8.7e, %8.7e, %8.7e, %8.7e\\n%8.7e, %8.7e, %8.7e, %8.7e, %8.7e, %8.7e, %8.7e\\n",\n')
        definitions.write('el[depart'+str(x)+',\'a\'],el[depart'+str(x)+',\'e\'],el[depart'+str(x)+',\'inc\'],el[depart'+str(x)+',\'Om\'],el[depart'+str(x)+',\'om\'],el[depart'+str(x)+',\'M\'],el[depart'+str(x)+',\'epoch\'], \n')
        definitions.write('el[target'+str(x)+',\'a\'],el[target'+str(x)+',\'e\'],el[target'+str(x)+',\'inc\'],el[target'+str(x)+',\'Om\'],el[target'+str(x)+',\'om\'],el[target'+str(x)+',\'M\'],el[target'+str(x)+',\'epoch\']> out/PlaParam'+str(x)+'.out;\n')
    definitions.close()
    return


def ephcalcgen(sequence, flyby):
    for x in range(len(sequence)-1):    
        ephcalc=open('./include/ephcalc'+str(x)+'.inc','w')
        ephcalc.write('#Single phase transfer (eph.inc)\n')
        ephcalc.write('#\n')
        ephcalc.write('#Written by Dario Izzo (October 2007)\n')
        ephcalc.write('#\n')
        ephcalc.write('#--------------------------------------------------------------------------\n')
        ephcalc.write('#This file is included in main.mod and contains the calculations that AMPL\n')
        ephcalc.write('#needs to do in order to locate the position of the initial and final bodies at the initial and final time\n')
        ephcalc.write('#--------------------------------------------------------------------------\n')
        ephcalc.write('\n')
        ephcalc.write('#--------------------------------------------------------------------------\n')
        ephcalc.write('#This parameter is the scale factor for the time variables\n')
        ephcalc.write('param factor'+str(x)+'	:=  1/100;\n')
        ephcalc.write('\n')
        ephcalc.write('#This parameter allows to define intervals where tI and tT are constrained\n')
        ephcalc.write('param tbnd'+str(x)+'	:=  20;		\n')
        ephcalc.write('					 	\n')
        ephcalc.write('#Initial time  (this is a variable the optimiser plays with)\n')
        ephcalc.write('var timod'+str(x)+'    :=   tI'+str(x)+' * d2u * factor'+str(x)+';	\n')
                      #,<= (tI'+str(x)+'+tbnd'+str(x)+')*d2u*factor'+str(x)+', >= (tI'+str(x)+'-tbnd'+str(x)+')*d2u*factor'+str(x)+';	\n')
        ephcalc.write('\n')
        ephcalc.write('#Time of flight (this is a variable the optimiser plays with)\n')
        ephcalc.write('var tfmod'+str(x)+'    :=   tT'+str(x)+' * d2u * factor'+str(x)+';	\n')
                      #, <= (tT'+str(x)+'+tbnd'+str(x)+')*d2u*factor'+str(x)+', >= (tT'+str(x)+'-tbnd'+str(x)+')*d2u*factor'+str(x)+';\n')
        ephcalc.write('#--------------------------------------------------------------------------\n')
        ephcalc.write('\n')
        ephcalc.write('#--------------------------------------------------------------------------\n')
        ephcalc.write('#We here introduce some other time variables that simplifies the following formulae	\n')
        ephcalc.write('var ti'+str(x)+'	     =   timod'+str(x)+' /factor'+str(x)+';			#Initial time non dimensional\n')
        ephcalc.write('var tf'+str(x)+'	     =   tfmod'+str(x)+' /factor'+str(x)+';			#Time of flight non dimensional\n')
        ephcalc.write('var tF'+str(x)+'	     =   ti'+str(x)+'/d2u + tf'+str(x)+'/d2u;			#Arrival time (MJD2000)\n')
        ephcalc.write('#--------------------------------------------------------------------------\n')
        ephcalc.write('\n')
        ephcalc.write('#--------------------------------------------------------------------------\n')
        ephcalc.write('#Calculations to find the position of the departure object at tI'+str(x)+'\n')
        ephcalc.write('\n')
        ephcalc.write('var M0'+str(x)+'		=	((el[depart'+str(x)+',\'M\']*d2r + sqrt(1/el[depart'+str(x)+',\'a\']^3) * ((ti'+str(x)+' / d2u+ 51544) - el[depart'+str(x)+',\'epoch\'])*d2u)\n')
        ephcalc.write('mod (2*pi)+2*pi) mod (2*pi);\n')
        ephcalc.write('\n')
        ephcalc.write('var E0'+str(x)+' >=0, <=2*pi;\n')
        ephcalc.write('let E0'+str(x)+'		:=	M0'+str(x)+' + el[depart'+str(x)+',\'e\']*cos(M0'+str(x)+');  	#Initial guess for Kepler equation\n')
        ephcalc.write('\n')
        ephcalc.write('var theta0'+str(x)+'	=	2*atan(sqrt((1+el[depart'+str(x)+',\'e\'])/(1-el[depart'+str(x)+',\'e\']))*tan(E0'+str(x)+'/2));\n')
        ephcalc.write('\n')
        ephcalc.write('var gamma0'+str(x)+'	=	atan( el[depart'+str(x)+',\'e\'] * sin(theta0'+str(x)+') / (1 + el[depart'+str(x)+',\'e\']*cos(theta0'+str(x)+')) );\n')
        ephcalc.write('\n')
        ephcalc.write('var r0'+str(x)+'		=	(el[depart'+str(x)+',\'a\'] * (1 - el[depart'+str(x)+',\'e\']^2)) / (1 + el[depart'+str(x)+',\'e\'] * cos(theta0'+str(x)+'));\n')
        ephcalc.write('\n')
        ephcalc.write('var v0'+str(x)+'		=	sqrt(2/r0'+str(x)+' - 1/el[depart'+str(x)+',\'a\']);\n')
        ephcalc.write('	\n')
        ephcalc.write('var x0'+str(x)+'		=	r0'+str(x)+' * ( cos(theta0'+str(x)+'+el[depart'+str(x)+',\'om\']*d2r)*cos(el[depart'+str(x)+',\'Om\']*d2r) -\n')
        ephcalc.write('				sin(theta0'+str(x)+'+el[depart'+str(x)+',\'om\']*d2r)*cos(el[depart'+str(x)+',\'inc\']*d2r)*sin(el[depart'+str(x)+',\'Om\']*d2r) );\n')
        ephcalc.write('				\n')
        ephcalc.write('var y0'+str(x)+'		=	r0'+str(x)+' * ( cos(theta0'+str(x)+'+el[depart'+str(x)+',\'om\']*d2r)*sin(el[depart'+str(x)+',\'Om\']*d2r) +\n')
        ephcalc.write('				sin(theta0'+str(x)+'+el[depart'+str(x)+',\'om\']*d2r)*cos(el[depart'+str(x)+',\'inc\']*d2r)*cos(el[depart'+str(x)+',\'Om\']*d2r) );\n')
        ephcalc.write('			\n')
        ephcalc.write('var z0'+str(x)+'		=	r0'+str(x)+' * ( sin(theta0'+str(x)+'+el[depart'+str(x)+',\'om\']*d2r) * sin(el[depart'+str(x)+',\'inc\']*d2r) );\n')
        ephcalc.write('	\n')
        ephcalc.write('var dx0'+str(x)+'	=	v0'+str(x)+' * (-sin(theta0'+str(x)+'+el[depart'+str(x)+',\'om\']*d2r-gamma0'+str(x)+')*cos(el[depart'+str(x)+',\'Om\']*d2r)\n')
        ephcalc.write('				-cos(theta0'+str(x)+'+el[depart'+str(x)+',\'om\']*d2r-gamma0'+str(x)+')*cos(el[depart'+str(x)+',\'inc\']*d2r)*sin(el[depart'+str(x)+',\'Om\']*d2r));\n')
        ephcalc.write('				\n')
        ephcalc.write('var dy0'+str(x)+'	=	v0'+str(x)+' * (-sin(theta0'+str(x)+'+el[depart'+str(x)+',\'om\']*d2r-gamma0'+str(x)+')*sin(el[depart'+str(x)+',\'Om\']*d2r)\n')
        ephcalc.write('				+cos(theta0'+str(x)+'+el[depart'+str(x)+',\'om\']*d2r-gamma0'+str(x)+')*cos(el[depart'+str(x)+',\'inc\']*d2r)*cos(el[depart'+str(x)+',\'Om\']*d2r));\n')
        ephcalc.write('				\n')
        ephcalc.write('var dz0'+str(x)+'	=	v0'+str(x)+' * (cos(theta0'+str(x)+'+el[depart'+str(x)+',\'om\']*d2r-gamma0'+str(x)+')*sin(el[depart'+str(x)+',\'inc\']*d2r));\n')
        ephcalc.write('	\n')
        ephcalc.write('subject to\n')
        ephcalc.write('	\n')
        ephcalc.write('	KeplerEquation0'+str(x)+': 1*(M0'+str(x)+' - E0'+str(x)+' + el[depart'+str(x)+',\'e\']*sin(E0'+str(x)+')) = 0;\n')
        ephcalc.write('#--------------------------------------------------------------------------------------------------------------------\n')
        ephcalc.write('\n')
        ephcalc.write('#--------------------------------------------------------------------------\n')
        ephcalc.write('#Calculations to find the position of the target'+str(x)+' object at tF'+str(x)+'\n')
        ephcalc.write('	\n')
        ephcalc.write('var Mf'+str(x)+'		=	((el[target'+str(x)+',\'M\']*d2r + sqrt(1/el[target'+str(x)+',\'a\']^3) * ((tF'+str(x)+' + 51544) - el[target'+str(x)+',\'epoch\'])*d2u)\n')
        ephcalc.write('mod (2*pi)+2*pi) mod (2*pi);\n')
        ephcalc.write('\n')
        ephcalc.write('var Ef'+str(x)+' >=0, <=2*pi;\n')
        ephcalc.write('let Ef'+str(x)+'		:=	Mf'+str(x)+' + el[target'+str(x)+',\'e\']*cos(Mf'+str(x)+');  	#Initial guess for Kepler equation\n')
        ephcalc.write('\n')
        ephcalc.write('var thetaf'+str(x)+'	=	2*atan(sqrt((1+el[target'+str(x)+',\'e\'])/(1-el[target'+str(x)+',\'e\']))*tan(Ef'+str(x)+'/2));\n')
        ephcalc.write('\n')
        ephcalc.write('var gammaf'+str(x)+'	=	atan( el[target'+str(x)+',\'e\'] * sin(thetaf'+str(x)+') / (1 + el[target'+str(x)+',\'e\']*cos(thetaf'+str(x)+')) );\n')
        ephcalc.write('\n')
        ephcalc.write('var rf'+str(x)+'		=	(el[target'+str(x)+',\'a\'] * (1 - el[target'+str(x)+',\'e\']^2)) / (1 + el[target'+str(x)+',\'e\'] * cos(thetaf'+str(x)+'));\n')
        ephcalc.write('\n')
        ephcalc.write('var vf'+str(x)+'		=	sqrt(2/rf'+str(x)+' - 1/el[target'+str(x)+',\'a\']);\n')
        ephcalc.write('	\n')
        ephcalc.write('var xf'+str(x)+'		=	rf'+str(x)+' * ( cos(thetaf'+str(x)+'+el[target'+str(x)+',\'om\']*d2r)*cos(el[target'+str(x)+',\'Om\']*d2r) -\n')
        ephcalc.write('				sin(thetaf'+str(x)+'+el[target'+str(x)+',\'om\']*d2r)*cos(el[target'+str(x)+',\'inc\']*d2r)*sin(el[target'+str(x)+',\'Om\']*d2r) );\n')
        ephcalc.write('\n')
        ephcalc.write('var yf'+str(x)+'		=	rf'+str(x)+' * ( cos(thetaf'+str(x)+'+el[target'+str(x)+',\'om\']*d2r)*sin(el[target'+str(x)+',\'Om\']*d2r) +\n')
        ephcalc.write('				sin(thetaf'+str(x)+'+el[target'+str(x)+',\'om\']*d2r)*cos(el[target'+str(x)+',\'inc\']*d2r)*cos(el[target'+str(x)+',\'Om\']*d2r) );\n')
        ephcalc.write('			\n')
        ephcalc.write('var zf'+str(x)+'		=	rf'+str(x)+' * ( sin(thetaf'+str(x)+'+el[target'+str(x)+',\'om\']*d2r) * sin(el[target'+str(x)+',\'inc\']*d2r) );\n')
        ephcalc.write('	\n')
        ephcalc.write('var dxf'+str(x)+'	=	vf'+str(x)+' * (-sin(thetaf'+str(x)+'+el[target'+str(x)+',\'om\']*d2r-gammaf'+str(x)+')*cos(el[target'+str(x)+',\'Om\']*d2r)\n')
        ephcalc.write('				-cos(thetaf'+str(x)+'+el[target'+str(x)+',\'om\']*d2r-gammaf'+str(x)+')*cos(el[target'+str(x)+',\'inc\']*d2r)*sin(el[target'+str(x)+',\'Om\']*d2r));\n')
        ephcalc.write('				\n')
        ephcalc.write('var dyf'+str(x)+'	=	vf'+str(x)+' * (-sin(thetaf'+str(x)+'+el[target'+str(x)+',\'om\']*d2r-gammaf'+str(x)+')*sin(el[target'+str(x)+',\'Om\']*d2r)\n')
        ephcalc.write('				+cos(thetaf'+str(x)+'+el[target'+str(x)+',\'om\']*d2r-gammaf'+str(x)+')*cos(el[target'+str(x)+',\'inc\']*d2r)*cos(el[target'+str(x)+',\'Om\']*d2r));\n')
        ephcalc.write('				\n')
        ephcalc.write('var dzf'+str(x)+'	=	vf'+str(x)+' * (cos(thetaf'+str(x)+'+el[target'+str(x)+',\'om\']*d2r-gammaf'+str(x)+')*sin(el[target'+str(x)+',\'inc\']*d2r));\n')
        ephcalc.write('	\n')
        ephcalc.write('subject to\n')
        ephcalc.write('	KeplerEquationf'+str(x)+': 1*(Mf'+str(x)+' - Ef'+str(x)+' + el[target'+str(x)+',\'e\']*sin(Ef'+str(x)+')) = 0;\n')
        ephcalc.write('#--------------------------------------------------------------------------------------------------------------------\n')
        ephcalc.close()
    return


def ocpsgen(sequence, flyby):
    for x in range(len(sequence)-1): 
        ocps=open('ocp'+str(x)+'.mod','w')
        ocps.write('#Single phase transfer (ocp.mod)\n')
        ocps.write('#\n')
        ocps.write('#--------------------------------------------------------------------------\n')
        ocps.write('#This file is included in ocp.mod and contains definition of the first phase\n')
        ocps.write('#--------------------------------------------------------------------------\n')
	
   	ocps.write('\n')
	ocps.write('#--------------------------------------------------------------------------\n')
        ocps.write('#Problem definition\n')
        ocps.write('subject to \n')
        ocps.write('	InitialMass'+str(x)+': m'+str(x)+'[1] = 1;\n')
        ocps.write('	InitialPositionx'+str(x)+': x'+str(x)+'[1] = x0'+str(x)+';\n')
        ocps.write('	InitialPositiony'+str(x)+': y'+str(x)+'[1] = y0'+str(x)+';\n')
        ocps.write('	InitialPositionz'+str(x)+': z'+str(x)+'[1] = z0'+str(x)+';\n')
        ocps.write('	FinalPositionx'+str(x)+'  : x'+str(x)+'[n'+str(x)+'] = xf'+str(x)+';\n')
        ocps.write('	FinalPositiony'+str(x)+'  : y'+str(x)+'[n'+str(x)+'] = yf'+str(x)+';\n')
        ocps.write('	FinalPositionz'+str(x)+'  : z'+str(x)+'[n'+str(x)+'] = zf'+str(x)+';\n')
        ocps.write('	\n')
        ocps.write('	ControlMagnitude'+str(x)+'{i in I'+str(x)+'} : sqrt(ux'+str(x)+'[i]**2  + uy'+str(x)+'[i]**2  + uz'+str(x)+'[i]**2)  <= Tmax;\n')
        ocps.write('	ControlMagnitudem'+str(x)+'{i in J'+str(x)+'}: sqrt(uxm'+str(x)+'[i]**2 + uym'+str(x)+'[i]**2 + uzm'+str(x)+'[i]**2) <= Tmax; 	 		\n')
        ocps.write('	\n')
        
	ocps.write('#C3allowed-----------------------------------------------------------------\n')
	
	if (x == 0):
		ocps.write('InitialC3'+str(x)+':	 sqrt((dx'+str(x)+'[1]-dx0'+str(x)+')**2 + (dy'+str(x)+'[1]-dy0'+str(x)+')**2 + (dz'+str(x)+'[1]-dz0'+str(x)+')**2) <= sqrt(C3lau);\n')
	elif (flyby[x-1]!=str(1)):
		ocps.write('InitialC3'+str(x)+':	 sqrt((dx'+str(x)+'[1]-dx0'+str(x)+')**2 + (dy'+str(x)+'[1]-dy0'+str(x)+')**2 + (dz'+str(x)+'[1]-dz0'+str(x)+')**2) <= sqrt(C3rnv);\n')

	if (x == len(sequence)-2):
		ocps.write('FinalC3'+str(x)+':	 sqrt((dx'+str(x)+'[n'+str(x)+']-dxf'+str(x)+')**2 + (dy'+str(x)+'[n'+str(x)+']-dyf'+str(x)+')**2 + (dz'+str(x)+'[n'+str(x)+']-dzf'+str(x)+')**2) <= sqrt(C3rnv);\n')
	elif (flyby[x]!=str(1)):
		ocps.write('FinalC3'+str(x)+':	 sqrt((dx'+str(x)+'[n'+str(x)+']-dxf'+str(x)+')**2 + (dy'+str(x)+'[n'+str(x)+']-dyf'+str(x)+')**2 + (dz'+str(x)+'[n'+str(x)+']-dzf'+str(x)+')**2) <= sqrt(C3rnv);\n')	

        ocps.write('\n')
        ocps.write('#--------------------------------------------------------------------------\n')
        ocps.write('#Dynamic at the grid points\n')
        ocps.write('var f1'+str(x)+'{i in I'+str(x)+'} = dx'+str(x)+'[i];\n')
        ocps.write('var f2'+str(x)+'{i in I'+str(x)+'} = dy'+str(x)+'[i];\n')
        ocps.write('var f3'+str(x)+'{i in I'+str(x)+'} = dz'+str(x)+'[i];\n')
        ocps.write('var f4'+str(x)+'{i in I'+str(x)+'} = -x'+str(x)+'[i]/(x'+str(x)+'[i]**2+y'+str(x)+'[i]**2+z'+str(x)+'[i]**2)^(3/2) + ux'+str(x)+'[i]/m'+str(x)+'[i];\n')
        ocps.write('var f5'+str(x)+'{i in I'+str(x)+'} = -y'+str(x)+'[i]/(x'+str(x)+'[i]**2+y'+str(x)+'[i]**2+z'+str(x)+'[i]**2)^(3/2) + uy'+str(x)+'[i]/m'+str(x)+'[i];\n')
        ocps.write('var f6'+str(x)+'{i in I'+str(x)+'} = -z'+str(x)+'[i]/(x'+str(x)+'[i]**2+y'+str(x)+'[i]**2+z'+str(x)+'[i]**2)^(3/2) + uz'+str(x)+'[i]/m'+str(x)+'[i];\n')
        ocps.write('var f7'+str(x)+'{i in I'+str(x)+'} = -sqrt(ux'+str(x)+'[i]**2+uy'+str(x)+'[i]**2+uz'+str(x)+'[i]**2)/Isp/g0;\n')
        ocps.write('#-----------------------------------------------------------------------\n')
        ocps.write('\n')
        ocps.write('#--------------------------------------------------------------------------\n')
        ocps.write('#State definition at mid-points via Simpson interpolation\n')
        ocps.write('var xm'+str(x)+'{i in J'+str(x)+'} = (x'+str(x)+'[i]+x'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f1'+str(x)+'[i] - f1'+str(x)+'[i+1]);\n')
        ocps.write('var ym'+str(x)+'{i in J'+str(x)+'} = (y'+str(x)+'[i]+y'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f2'+str(x)+'[i] - f2'+str(x)+'[i+1]);\n')
        ocps.write('var zm'+str(x)+'{i in J'+str(x)+'} = (z'+str(x)+'[i]+z'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f3'+str(x)+'[i] - f3'+str(x)+'[i+1]);\n')
        ocps.write('var dxm'+str(x)+'{i in J'+str(x)+'} = (dx'+str(x)+'[i]+dx'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f4'+str(x)+'[i] - f4'+str(x)+'[i+1]);\n')
        ocps.write('var dym'+str(x)+'{i in J'+str(x)+'} = (dy'+str(x)+'[i]+dy'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f5'+str(x)+'[i] - f5'+str(x)+'[i+1]);\n')
        ocps.write('var dzm'+str(x)+'{i in J'+str(x)+'} = (dz'+str(x)+'[i]+dz'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f6'+str(x)+'[i] - f6'+str(x)+'[i+1]);\n')
        ocps.write('var mm'+str(x)+'{i in J'+str(x)+'} = (m'+str(x)+'[i]+m'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f7'+str(x)+'[i] - f7'+str(x)+'[i+1]);\n')
        ocps.write('#-----------------------------------------------------------------------\n')
        ocps.write('\n')
        ocps.write('#--------------------------------------------------------------------------\n')
        ocps.write('#Dynamic at the mid-points\n')
        ocps.write('var f1m'+str(x)+'{i in J'+str(x)+'} = dxm'+str(x)+'[i];\n')
        ocps.write('var f2m'+str(x)+'{i in J'+str(x)+'} = dym'+str(x)+'[i];\n')
        ocps.write('var f3m'+str(x)+'{i in J'+str(x)+'} = dzm'+str(x)+'[i];\n')
        ocps.write('var f4m'+str(x)+'{i in J'+str(x)+'} = -xm'+str(x)+'[i]/(xm'+str(x)+'[i]**2+ym'+str(x)+'[i]**2+zm'+str(x)+'[i]**2)^(3/2) + uxm'+str(x)+'[i]/mm'+str(x)+'[i];\n')
        ocps.write('var f5m'+str(x)+'{i in J'+str(x)+'} = -ym'+str(x)+'[i]/(xm'+str(x)+'[i]**2+ym'+str(x)+'[i]**2+zm'+str(x)+'[i]**2)^(3/2) + uym'+str(x)+'[i]/mm'+str(x)+'[i];\n')
        ocps.write('var f6m'+str(x)+'{i in J'+str(x)+'} = -zm'+str(x)+'[i]/(xm'+str(x)+'[i]**2+ym'+str(x)+'[i]**2+zm'+str(x)+'[i]**2)^(3/2) + uzm'+str(x)+'[i]/mm'+str(x)+'[i];\n')
        ocps.write('var f7m'+str(x)+'{i in J'+str(x)+'} = -sqrt(uxm'+str(x)+'[i]**2+uym'+str(x)+'[i]**2+uzm'+str(x)+'[i]**2)/Isp/g0;\n')
        ocps.write('#-----------------------------------------------------------------------\n')
        ocps.write('\n')
        ocps.write('#--------------------------------------------------------------------------\n')
        ocps.write('#Hermite Formula\n')
        ocps.write('subject to \n')
        ocps.write('	dynamicx'+str(x)+'{i in J'+str(x)+'}:  x'+str(x)+'[i]  = x'+str(x)+'[i+1]  - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f1'+str(x)+'[i] + f1'+str(x)+'[i+1] + 4*f1m'+str(x)+'[i]);\n')
        ocps.write('	dynamicy'+str(x)+'{i in J'+str(x)+'}:  y'+str(x)+'[i]  = y'+str(x)+'[i+1]  - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f2'+str(x)+'[i] + f2'+str(x)+'[i+1] + 4*f2m'+str(x)+'[i]);\n')
        ocps.write('	dynamicz'+str(x)+'{i in J'+str(x)+'}:  z'+str(x)+'[i]  = z'+str(x)+'[i+1]  - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f3'+str(x)+'[i] + f3'+str(x)+'[i+1] + 4*f3m'+str(x)+'[i]);\n')
        ocps.write('	dynamicdx'+str(x)+'{i in J'+str(x)+'}: dx'+str(x)+'[i] = dx'+str(x)+'[i+1] - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f4'+str(x)+'[i] + f4'+str(x)+'[i+1] + 4*f4m'+str(x)+'[i]);\n')
        ocps.write('	dynamicdy'+str(x)+'{i in J'+str(x)+'}: dy'+str(x)+'[i] = dy'+str(x)+'[i+1] - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f5'+str(x)+'[i] + f5'+str(x)+'[i+1] + 4*f5m'+str(x)+'[i]);\n')
        ocps.write('	dynamicdz'+str(x)+'{i in J'+str(x)+'}: dz'+str(x)+'[i] = dz'+str(x)+'[i+1] - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f6'+str(x)+'[i] + f6'+str(x)+'[i+1] + 4*f6m'+str(x)+'[i]);\n')
        ocps.write('	dynamicm'+str(x)+'{i in J'+str(x)+'}:  m'+str(x)+'[i]  = m'+str(x)+'[i+1]  - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f7'+str(x)+'[i] + f7'+str(x)+'[i+1] + 4*f7m'+str(x)+'[i]);\n')
        ocps.write('#--------------------------------------------------------------------------	\n')
        ocps.write('	\n')
        ocps.write('#--------------------------------------------------------------------------	\n')
        ocps.write('var Thrust'+str(x)+'{i in I'+str(x)+'} = sqrt(ux'+str(x)+'[i]**2+uy'+str(x)+'[i]**2+uz'+str(x)+'[i]**2);\n')
        ocps.write('var Thrustm'+str(x)+'{i in J'+str(x)+'} = sqrt(uxm'+str(x)+'[i]**2+uym'+str(x)+'[i]**2+uzm'+str(x)+'[i]**2);      \n')
        ocps.write('                   \n')
        ocps.write('#--------------------------------------------------------------------------------------------------------------------\n')
        ocps.close()
    return

def ocpgen(sequence, flyby):
    ocp=open('ocp.mod','w')
    ocp.write('#Multiphase transfer (ocp.mod)\n')
    ocp.write('#\n')
    ocp.write('#--------------------------------------------------------------------------\n')
    ocp.write('#This file is included in main.mod and contains the construction of the \n')
    ocp.write('#optimisation problem with two phases\n')
    ocp.write('#--------------------------------------------------------------------------\n')
    ocp.write('\n')

    for x in range(len(sequence)-1):
        ocp.write('include ocp'+str(x)+'.mod;\n')
	ocp.write('minimize QCpart'+str(x)+': 1000 * sum{i in J'+str(x)+'}(Thrust'+str(x)+'[i]**2+Thrust'+str(x)+'[i+1]**2+4*Thrustm'+str(x)+'[i]**2);\n')
	ocp.write('solve;\n')
	ocp.write('purge QCpart'+str(x)+';\n')


    #Purge Initial Mass
    for x in range(1,len(sequence)-1):
    	ocp.write('purge InitialMass'+str(x)+';\n')
    ocp.write('\n')

    #Patching Conditions
    for x in range(len(sequence)-2):
    
        if (flyby[x]==str(1)):
	    	ocp.write('#--------------------------------------------------------------------------\n')
    		ocp.write('#Flyby Constraints patching phases '+str(x)+' and '+str(x+1)+';\n')
	    	ocp.write('\n')
    		ocp.write('var Vrelin'+str(x+1)+'	=	(dx'+str(x)+'[n'+str(x)+']-dxf'+str(x)+')**2 + (dy'+str(x)+'[n'+str(x)+']-dyf'+str(x)+')**2 +(dz'+str(x)+'[n'+str(x)+']-dzf'+str(x)+')**2;\n')
    		ocp.write('var Vrelout'+str(x+1)+'	=	(dx'+str(x+1)+'[1] -dx0'+str(x+1)+')**2 + (dy'+str(x+1)+'[1] -dy0'+str(x+1)+')**2 + (dz'+str(x+1)+'[1]-dz0'+str(x+1)+')**2;\n')
    		ocp.write('var emax'+str(x)+'	=	1 + rmin/mupla*Vrelin'+str(x+1)+';\n')
    		ocp.write('var alfa'+str(x)+'	=	acos(\n')
    		ocp.write('			     ((dx'+str(x)+'[n'+str(x)+']-dxf'+str(x)+')*(dx'+str(x+1)+'[1]-dx0'+str(x+1)+') + \n')
    		ocp.write('			      (dy'+str(x)+'[n'+str(x)+']-dyf'+str(x)+')*(dy'+str(x+1)+'[1]-dy0'+str(x+1)+') + \n')
    		ocp.write('			      (dz'+str(x)+'[n'+str(x)+']-dzf'+str(x)+')*(dz'+str(x+1)+'[1]-dz0'+str(x+1)+')) / sqrt(Vrelin'+str(x+1)+') / sqrt(Vrelout'+str(x+1)+')\n')
    		ocp.write('			);\n')
    		ocp.write('\n')
    		ocp.write('subject to \n')
    		ocp.write('	samemass'+str(x)+': m'+str(x+1)+'[1] = m'+str(x)+'[n'+str(x)+'];\n')
    		ocp.write('	sametime'+str(x)+': ti'+str(x+1)+' = ti'+str(x)+' + tf'+str(x)+';\n')
    		ocp.write('	samerelvel'+str(x)+': abs(Vrelin'+str(x+1)+' - Vrelout'+str(x+1)+') <= 0.00001;\n')
    		ocp.write('	angularconstraint'+str(x)+': alfa'+str(x)+' <= 2 * asin(1/emax'+str(x)+');\n')
    		ocp.write('#--------------------------------------------------------------------------\n')

	else :
	
	    	ocp.write('#--------------------------------------------------------------------------\n')
    		ocp.write('#Waiting time Constraints patching phases '+str(x)+' and '+str(x+1)+';\n')
		ocp.write('\n')
    		ocp.write('subject to \n')
    		ocp.write('	samemass'+str(x)+': m'+str(x+1)+'[1] = m'+str(x)+'[n'+str(x)+'];\n')
    		ocp.write('	waittime'+str(x)+': ti'+str(x+1)+' >= '+str(days)+' * d2u + ti'+str(x)+' + tf'+str(x)+' ;\n')
    		ocp.write('#--------------------------------------------------------------------------\n')
##  for x in range(1,len(sequence)-2):
##     	ocp.write('minimize QC'+str(x)+': 1000 * sum{i in J'+str(x)+'}(Thrust'+str(x)+'[i]**2+Thrust'+str(x)+'[i+1]**2+4*Thrustm'+str(x)+'[i]**2);\n')
##	ocp.write('solve;\n')
##	ocp.write('purge QC'+str(x)+';\n')
      	ocp.write('maximize mass'+str(x)+': m'+str(x)+'[n'+str(x)+'];\n')
	ocp.write('solve;\n')
  	ocp.write('purge mass'+str(x)+';\n')
#  	ocp.write('purge QC'+str(x)+';\n')
#	if (flyby[x]==str(1)):
#    		ocp.write('drop	samemass'+str(x)+';\n')
#    		ocp.write('drop	sametime'+str(x)+';\n')
#    		ocp.write('drop	samerelvel'+str(x)+';\n')
#    		ocp.write('drop	angularconstraint'+str(x)+';\n')
#	else :
#    		ocp.write('drop	samemass'+str(x)+';\n')
#    		ocp.write('drop	waittime'+str(x)+';\n')
 
#   for x in range(len(sequence)-2):
#    
#	if (flyby[x]==str(1)):
#    		ocp.write('subject to \n')
#    		ocp.write('	samemass'+str(x)+': m'+str(x+1)+'[1] = m'+str(x)+'[n'+str(x)+'];\n')
#    		ocp.write('	sametime'+str(x)+': ti'+str(x+1)+' = ti'+str(x)+' + tf'+str(x)+';\n')
#    		ocp.write('	samerelvel'+str(x)+': abs(Vrelin'+str(x+1)+' - Vrelout'+str(x+1)+') <= 0.00001;\n')
#    		ocp.write('	angularconstraint'+str(x)+': alfa'+str(x)+' <= 2 * asin(1/emax'+str(x)+');\n')
#	else :
#    		ocp.write('subject to \n')
#    		ocp.write('	samemass'+str(x)+': m'+str(x+1)+'[1] = m'+str(x)+'[n'+str(x)+'];\n')
#    		ocp.write('	waittime'+str(x)+': ti'+str(x+1)+' >= '+str(days)+' * d2u + ti'+str(x)+' + tf'+str(x)+' ;\n')
    ocp.write('#--------------------------------------------------------------------------\n')

    ocp.write('#Total time constraint  \n')
    ocp.write('subject to  \n')
    ocp.write('    totaltime: ti0 + TOTALTIME * d2u >= tf'+str(len(sequence)-2)+';\n')
##  ocp.write('#--------------------------------------------------------------------------\n')
##  ocp.write('#FINAL QC Optimization \n')
##  ocp.write('\n')
##  ocp.write('minimize QC: \n')
##  for x in range(len(sequence)-2):
##      ocp.write('		1000 * sum{i in J'+str(x)+'}(Thrust'+str(x)+'[i]**2+Thrust'+str(x)+'[i+1]**2+4*Thrustm'+str(x)+'[i]**2) +')
##  ocp.write('		1000 * sum{i in J'+str(len(sequence)-2)+'}(Thrust'+str(len(sequence)-2)+'[i]**2+Thrust'+str(len(sequence)-2)+'[i+1]**2+4*Thrustm'+str(len(sequence)-2)+'[i]**2); \n')
##  ocp.write('solve;\n')
##  ocp.write('purge QC;\n')
    ocp.write('#--------------------------------------------------------------------------\n')
    ocp.write('#FINAL OBJFUN Optimization \n')
    ocp.write('\n')
    ocp.write('maximize finalmass: m'+str(len(sequence)-2)+'[n'+str(len(sequence)-2)+'];\n')
    ocp.write('solve;\n')
    ocp.write('\n')
    ocp.write('#--------------------------------------------------------------------------\n')
    ocp.write('#PRINTING THE GLOBAL SOLUTION\n')
    for x in range(len(sequence)-1):
        ocp.write('#--------------------------------------------------------------------------\n')
        ocp.write('#SOLUTION PHASE '+str(x)+'\n')
        ocp.write('\n')       
        ocp.write('for {i in J'+str(x)+'} {\n')
        ocp.write('printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\\n%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e\\n",\n')
        ocp.write('x'+str(x)+'[i],y'+str(x)+'[i],z'+str(x)+'[i],dx'+str(x)+'[i],dy'+str(x)+'[i],dz'+str(x)+'[i],m'+str(x)+'[i],ux'+str(x)+'[i],uy'+str(x)+'[i],uz'+str(x)+'[i],xm'+str(x)+'[i],ym'+str(x)+'[i],zm'+str(x)+'[i],dxm'+str(x)+'[i],dym'+str(x)+'[i],dzm'+str(x)+'[i],mm'+str(x)+'[i],uxm'+str(x)+'[i],uym'+str(x)+'[i],uzm'+str(x)+'[i]>out/Solution'+str(x)+'.out;\n')
        ocp.write('}\n')
        ocp.write('printf "%17.16e, %17.16e \\n", ti'+str(x)+'/d2u , tF'+str(x)+'-ti'+str(x)+'/d2u > out/Times'+str(x)+'.out;\n')
        ocp.write('printf "%17.16e \\n", m'+str(x)+'[n'+str(x)+'] > out/mass'+str(x)+'.out;\n')
        if (x<=len(sequence)-3):
            if (flyby[x]==str(1)):
                ocp.write('printf "%17.16e \\n %17.16e \\n", mupla/Vrelin'+str(x+1)+'*(1/sin(alfa'+str(x)+'/2)-1)*R, sqrt((dx'+str(x)+'[n'+str(x)+']-dx'+str(x+1)+'[1])**2 +(dy'+str(x)+'[n'+str(x)+']-dy'+str(x+1)+'[1])**2+(dz'+str(x)+'[n'+str(x)+']-dz'+str(x+1)+'[1])**2)*V  > out/flyby'+str(x)+'.out;\n')
    ocp.write('#--------------------------------------------------------------------------\n')
    ocp.close()
    return




def guessgen(sequence, flyby):
    for x in range(len(sequence)-1):    
        guess=open('./include/guess'+str(x)+'.inc','w')
        guess.write('#Single phase transfer (iniguess.inc)\n')
        guess.write('#\n')
        guess.write('#--------------------------------------------------------------------------\n')
        guess.write('#This file is included in main.mod and contains the calculations that AMPL\n')
        guess.write('#needs to do in order to build an initial guess for the transfer\n')
        guess.write('\n')
        guess.write('#--------------------------------------------------------------------------\n')
        guess.write('#Problem definition\n')
        guess.write('subject to \n')
        guess.write('	InitialPositionx'+str(x)+': x'+str(x)+'[1] = x0'+str(x)+';\n')
        guess.write('	InitialPositiony'+str(x)+': y'+str(x)+'[1] = y0'+str(x)+';\n')
        guess.write('	InitialPositionz'+str(x)+': z'+str(x)+'[1] = z0'+str(x)+';\n')
        guess.write('	FinalPositionx'+str(x)+'  : x'+str(x)+'[n'+str(x)+'] = xf'+str(x)+';\n')
        guess.write('	FinalPositiony'+str(x)+'  : y'+str(x)+'[n'+str(x)+'] = yf'+str(x)+';\n')
        guess.write('	FinalPositionz'+str(x)+'  : z'+str(x)+'[n'+str(x)+'] = zf'+str(x)+';\n')
        guess.write('#-----------------------------------------------------------------------		\n')
        guess.write('\n')
        guess.write('#--------------------------------------------------------------------------\n')
        guess.write('#Dynamic at the grid points\n')
        guess.write('var f1{i in I'+str(x)+'} = dx'+str(x)+'[i];\n')
        guess.write('var f2{i in I'+str(x)+'} = dy'+str(x)+'[i];\n')
        guess.write('var f3{i in I'+str(x)+'} = dz'+str(x)+'[i];\n')
        guess.write('var f4{i in I'+str(x)+'} = -x'+str(x)+'[i]/(x'+str(x)+'[i]**2+y'+str(x)+'[i]**2+z'+str(x)+'[i]**2)^(3/2);\n')
        guess.write('var f5{i in I'+str(x)+'} = -y'+str(x)+'[i]/(x'+str(x)+'[i]**2+y'+str(x)+'[i]**2+z'+str(x)+'[i]**2)^(3/2);\n')
        guess.write('var f6{i in I'+str(x)+'} = -z'+str(x)+'[i]/(x'+str(x)+'[i]**2+y'+str(x)+'[i]**2+z'+str(x)+'[i]**2)^(3/2);\n')
        guess.write('#-----------------------------------------------------------------------\n')
        guess.write('\n')
        guess.write('#--------------------------------------------------------------------------\n')
        guess.write('#State definition at mid-points via Simpson interpolation\n')
        guess.write('var xm{i in J'+str(x)+'} = (x'+str(x)+'[i]+x'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f1[i] - f1[i+1]);\n')
        guess.write('var ym{i in J'+str(x)+'} = (y'+str(x)+'[i]+y'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f2[i] - f2[i+1]);\n')
        guess.write('var zm{i in J'+str(x)+'} = (z'+str(x)+'[i]+z'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f3[i] - f3[i+1]);\n')
        guess.write('var dxm{i in J'+str(x)+'} = (dx'+str(x)+'[i]+dx'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f4[i] - f4[i+1]);\n')
        guess.write('var dym{i in J'+str(x)+'} = (dy'+str(x)+'[i]+dy'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f5[i] - f5[i+1]);\n')
        guess.write('var dzm{i in J'+str(x)+'} = (dz'+str(x)+'[i]+dz'+str(x)+'[i+1])/2 + tf'+str(x)+'/(n'+str(x)+'-1)/8 * (f6[i] - f6[i+1]);\n')
        guess.write('#-----------------------------------------------------------------------\n')
        guess.write('\n')
        guess.write('#--------------------------------------------------------------------------\n')
        guess.write('#Dynamic at the mid-points\n')
        guess.write('var f1m{i in J'+str(x)+'} = dxm[i];\n')
        guess.write('var f2m{i in J'+str(x)+'} = dym[i];\n')
        guess.write('var f3m{i in J'+str(x)+'} = dzm[i];\n')
        guess.write('var f4m{i in J'+str(x)+'} = -xm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);\n')
        guess.write('var f5m{i in J'+str(x)+'} = -ym[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);\n')
        guess.write('var f6m{i in J'+str(x)+'} = -zm[i]/(xm[i]**2+ym[i]**2+zm[i]**2)^(3/2);\n')
        guess.write('#-----------------------------------------------------------------------\n')
        guess.write('\n')
        guess.write('#--------------------------------------------------------------------------\n')
        guess.write('#Hermite Formula\n')
        guess.write('subject to \n')
        guess.write('	dynamicx'+str(x)+'{i in J'+str(x)+'}: x'+str(x)+'[i] = x'+str(x)+'[i+1] - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f1[i] + f1[i+1] + 4*f1m[i]);\n')
        guess.write('	dynamicy'+str(x)+'{i in J'+str(x)+'}: y'+str(x)+'[i] = y'+str(x)+'[i+1] - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f2[i] + f2[i+1] + 4*f2m[i]);\n')
        guess.write('	dynamicz'+str(x)+'{i in J'+str(x)+'}: z'+str(x)+'[i] = z'+str(x)+'[i+1] - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f3[i] + f3[i+1] + 4*f3m[i]);\n')
        guess.write('	dynamicdx'+str(x)+'{i in J'+str(x)+'}: dx'+str(x)+'[i] = dx'+str(x)+'[i+1] - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f4[i] + f4[i+1] + 4*f4m[i]);\n')
        guess.write('	dynamicdy'+str(x)+'{i in J'+str(x)+'}: dy'+str(x)+'[i] = dy'+str(x)+'[i+1] - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f5[i] + f5[i+1] + 4*f5m[i]);\n')
        guess.write('	dynamicdz'+str(x)+'{i in J'+str(x)+'}: dz'+str(x)+'[i] = dz'+str(x)+'[i+1] - tf'+str(x)+'/(n'+str(x)+'-1)/6*(f6[i] + f6[i+1] + 4*f6m[i]);\n')
        guess.write('		\n')
        guess.write('#--------------------------------------------------------------------------\n')
        guess.write('#Initial Guess for the LP\n')
        guess.write('let {i in I'+str(x)+'} x'+str(x)+'[i] :=    x0'+str(x)+'*cos(tf'+str(x)+'/(n'+str(x)+'-1)*(i-1)) - y0'+str(x)+'*sin(tf'+str(x)+'/(n'+str(x)+'-1)*(i-1));\n')
        guess.write('let {i in I'+str(x)+'} y'+str(x)+'[i] :=    x0'+str(x)+'*sin(tf'+str(x)+'/(n'+str(x)+'-1)*(i-1)) + y0'+str(x)+'*cos(tf'+str(x)+'/(n'+str(x)+'-1)*(i-1));\n')
        guess.write('let {i in I'+str(x)+'} dx'+str(x)+'[i] :=  -x0'+str(x)+'*sin(tf'+str(x)+'/(n'+str(x)+'-1)*(i-1)) - y0'+str(x)+'*cos(tf'+str(x)+'/(n'+str(x)+'-1)*(i-1));\n')
        guess.write('let {i in I'+str(x)+'} dy'+str(x)+'[i] :=   x0'+str(x)+'*cos(tf'+str(x)+'/(n'+str(x)+'-1)*(i-1)) - y0'+str(x)+'*sin(tf'+str(x)+'/(n'+str(x)+'-1)*(i-1));\n')
        guess.write('#-----------------------------------------------------------------------\n')
        guess.write('\n')
        guess.write('solve;\n')
        guess.write('\n')
        guess.write('#--------------------------------------------------------------------------\n')
        guess.write('#Initial Guess for the remaining variables\n')
        guess.write('let {i in I'+str(x)+'} m'+str(x)+'[i] := 1;\n')
        guess.write('\n')
        guess.write('#The following guess is needed for snopt as otherwise snopt starts with an\n')
        guess.write('#arithmetic error\n')
        guess.write('let {i in I'+str(x)+'} ux'+str(x)+'[i]:=Tmax*0.0001;\n')
        guess.write('let {i in I'+str(x)+'} uy'+str(x)+'[i]:=Tmax*0.0001;\n')
        guess.write('let {i in I'+str(x)+'} uz'+str(x)+'[i]:=Tmax*0.0001;\n')
        guess.write('let {i in J'+str(x)+'} uxm'+str(x)+'[i]:=Tmax*0.0001;\n')
        guess.write('let {i in J'+str(x)+'} uym'+str(x)+'[i]:=Tmax*0.0001;\n')
        guess.write('let {i in J'+str(x)+'} uzm'+str(x)+'[i]:=Tmax*0.0001;\n')
        guess.write('\n')
        guess.write('#-----------------------------------------------------------------------\n')
        guess.write('\n')
        guess.write('\n')
        guess.write('#--------------------------------------------------------------------------\n')
        guess.write('#Print the Initial Guess\n')
        guess.write('for {i in I'+str(x)+'} {\n')
        guess.write('printf "%17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e, %17.16e,%17.16e\\n",x'+str(x)+'[i],y'+str(x)+'[i],z'+str(x)+'[i],dx'+str(x)+'[i],dy'+str(x)+'[i],dz'+str(x)+'[i],m'+str(x)+'[i],ux'+str(x)+'[i],uy'+str(x)+'[i],uz'+str(x)+'[i]>out/InitialGuess'+str(x)+'.out;\n')
        guess.write('}\n')
        guess.write('#------------------------------------------------------------------------\n')
        guess.write('\n')
        guess.write('#Print the initial and final times\n')
        guess.write('printf "%17.16e, %17.16e \\n", ti'+str(x)+'/d2u , tF'+str(x)+'-ti'+str(x)+'/d2u > out/TimesGuess'+str(x)+'.out;\n')
        guess.write('\n')
        guess.write('\n')
        guess.write('\n')
        guess.write('#--------------------------------------------------------------------------\n')
        guess.write('#As to calculate the LP initial guess we made use of a simplified dynamic we \n')
        guess.write('#we need to dispose some variables (note that AMPL will automatically dispose\n')
        guess.write('#also those variables linked to the ones listed here\n')
        guess.write('\n')
        guess.write('\n')
        guess.write('purge\n')
        guess.write('f1,f2,f3,f4,f5,f6,f1m,f2m,f3m,f4m,f5m,f6m,InitialPositionx'+str(x)+',InitialPositiony'+str(x)+',InitialPositionz'+str(x)+',FinalPositionx'+str(x)+',FinalPositiony'+str(x)+',FinalPositionz'+str(x)+';\n')
        guess.write('#--------------------------------------------------------------------------\n')
    return





#------------------------------------------------------------------------
#MAIN PROGRAM
#------------------------------------------------------------------------


#Importing modules needed: re(regular expressions), os (operative system)
import re, os
#Opening file input.dat
inp=open('input.dat','r')
#Reading sequence of bodies from input.dat 
sequence= re.findall('(\w+)',inp.readline())
#Reading flyby,times conditions from input.dat and closing file
#BE careful, with the next two lines, the variables become strings
#REMEMBER TO MODIFY THE REGULAR EXPRESSION TO MATCH ALL POSSIBLE NUMBER NOTATIONS
flyby=re.findall('([\de\-\.]+)',inp.readline())
times=re.findall('([\de\-\.]+)',inp.readline())

inp.close()
#Generating main.mod
maingen(sequence,flyby,times)
#Creating folders include, out
os.system('mkdir include')
os.system('mkdir out')
#Generating ephcalcfiles and definitions
definitionsgen(sequence,flyby)
ephcalcgen(sequence,flyby)
ocpsgen(sequence,flyby)
ocpgen(sequence,flyby)
guessgen(sequence,flyby)
