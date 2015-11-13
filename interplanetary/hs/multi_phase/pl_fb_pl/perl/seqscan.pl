#!/usr/bin/perl

#Including Libraries

use File::Copy;
use Time::HiRes qw(gettimeofday);

#Getting the inputs

open(INPUT,"./inputfly.dat");
@inp=<INPUT>;
$_=$inp[0];
@asteroids=/(\d+)/g;
$_=$inp[1];
@timedata=/([+-\d.Ee]+)/g;
$_=$inp[2];
@mass_arr1_dep2=/([+-\d.Ee]+)/g;
close(INPUT);

#Modifying variables to complete the sequence
$asteroids[2]=$asteroids[1];
$asteroids[1]=141;

#Initialising mass, minimum waiting time and other variables
@guesses=("crazyguess.inc","guesslambert.inc","guesslambert2.inc","guesstangential.inc","guesstangential2.inc","guesstangential2m.inc","guesstangential3.inc","guesstangential3m.inc");
@corrections=(0,-50,+50);
#Creating folder structure
system("mkdir ./out");
system("mkdir ./out/".$asteroids[0]."_".$asteroids[2]);
#Opening OUTPUT file
open(OUTPUT,">./out/".$asteroids[0]."_".$asteroids[2]."/output.dat");

#Creating a temporary file as copy
copy("./main.mod","./main.bak");

#Writing Sequence in the Output file
print OUTPUT "Sequence ".$asteroids[0]." - Flyby Earth - ".$asteroids[2]."\n\n";
@best=(0, 0, 0);
$max=0;
for($l=0;$l<=$#corrections;$l++) {
  for($m=0;$m<=$#corrections;$m++) {
   for($n=0;$n<=$#corrections;$n++) {
	if ($timedata[0]+$timedata[1]+$timedata[2]+$corrections[$l]+$corrections[$m]+$corrections[$n]>$mass_arr1_dep2[2]) { 
		last	;}
	if ($timedata[0]+$corrections[$l]<$mass_arr1_dep2[1]+60) { 
		last	;}
	open(MAINOLD,"main.bak");
	open(MAINNEW,">main.mod");
	foreach (<MAINOLD>){
	    if (/param depart/) { 
	      $_="param depart:='ast".$asteroids[0]."' symbolic;\n";
	    }
	    elsif (/param target/) {
	      $_="param target:='ast".$asteroids[2]."' symbolic;\n";
	    }
	    elsif (/param M/) {
	     $_="param M	     :=   (".($mass_arr1_dep2[0]).");	\n";
	    }
	    elsif (/param C3dep/) {
	      if ($asteroids[$0]==141) {	
	      $_="param C3dep     :=   (0.5**2)  /V**2;	#Initial Launch C3 (km^2/sec^2)\n";
	      }
	      else {
	      $_="param C3dep     :=   (0.001**2)  /V**2;	#Initial Launch C3 (km^2/sec^2)\n";
	      }
	    }
	    elsif (/param C3arr/) {
	      $_="param C3arr     :=   (0.001**2)  /V**2;	#Arrival C3 allowed (km^2/sec^2)\n";
	    }
	    elsif (/var tI1/) {
	      $_="var tI1	     :=   ".($timedata[0]+$corrections[$l]).";		#Departure time (MJD2000)\n";
	    } 
	    elsif (/var tT1/) {
	      $_="var tT1	     :=   ".($timedata[1]+$corrections[$m]).";		#Arrival time (MJD2000)\n";
	    }
	    elsif (/var tI2/) {
	      $_="var tI2	     :=   ".($timedata[0]+$timedata[1]+$corrections[$l]+$corrections[$m]).";		#Departure time (MJD2000)\n";
	    } 

	    elsif (/var tT2/) {
	      $_="var tT2	     :=   ".($timedata[2]+$corrections[$n]).";		#Arrival time (MJD2000)\n";

	    }

	    print MAINNEW $_;
	}
	close(MAINOLD);
	close(MAINNEW);
	system("ampl main.mod");
	system("cp main.mod ./out/".$asteroids[0]."_".$asteroids[2]);
	system("cp out/*  ./out/".$asteroids[0]."_".$asteroids[2]);
	
	open(TIME,"out/Times1.out");
	foreach (<TIME>){
	    @time1=/([+-\d.Ee]+)/g;
	    last;
	}
	close(TIME);
	
	open(TIME,"out/Times2.out");
	foreach (<TIME>){
	    @time2=/([+-\d.Ee]+)/g;
	    last;
	}
	close(TIME);
	
		open(MASS,"out/mass2.out");
		foreach (<MASS>){
		    /([+-\d.Ee]+)/;
		    $objfunguess=$1;
		    last;
		}
		close(MASS);
		$value=$objfunguess;

		if (&infeasible) {
		print OUTPUT "ERROR: INFEASIBLE\n";
		$objfunguess=$objfunguess-0.4;
		}

		if ($objfunguess>$max) {
			$max=$objfunguess;
			@best=($l, $m, $n);
		}
		print OUTPUT "DEPARTURE ". ($time1[0])." \n";
		print OUTPUT "PERIOD1 ". ($time1[1])." \n";
		print OUTPUT "FLYBY ". ($time2[0])." \n";
		print OUTPUT "PERIOD2 ". ($time2[1])." \n";
		print OUTPUT "MASS ".($mass_arr1_dep2[0]*$value)."\n";
		

		print OUTPUT "OBJFUN ".($value)."\n\n";
		&imprimearray(\@best);
	}
	       
    }

}

print " BEST CORRECTIONS GIVE A VALUE OF : ".$max."\n";
print OUTPUT " BEST CORRECTIONS GIVE A VALUE OF : ".$max."\n";
&imprimearray(\@best);

print "============================".$objfun."======================================\n";
close(OUTPUT);

open(VALUE,">sol");
print VALUE $objfun;
close(VALUE);

sub min {
    local ($x, $y) = @_;
    return $x < $y ? $x : $y;
}



sub infeasible {
    open(LOGFILE,"out/log.txt");
    @LOG=<LOGFILE>;
    close(LOGFILE);
    if (($LOG[$#LOG-30] =~ /i/)  || !($LOG[$#LOG-30] =~ /\([+-\d.Ee]+\)/) ){
    return 1;
    }
    else {
    return 0;
    }
}



sub max {
    local ($x, $y) = @_;
    return $x > $y ? $x : $y;
}

sub imprimearray {
	local (@array) = @{$_[0]};
	local ($variable);
	local ($i);
	for ($i=0; $i <= $#array; $i++) {
		$variable= $array[$i];
		if ($i == $#array) {
				print OUTPUT $corrections[$variable]."\n";
		}
		else {
				print OUTPUT $corrections[$variable]."  ";
		}
	}
}

sub inarray {
	local ($valor) =0;
	local ($elemento) = $_[0];
	local ($posicion) = -1;
	local (@array) = @{$_[1]};
	for ($j=0; $j <= $#array; $j++) {
		if ($array[$j]==($elemento/1)) {
			$valor = 1;
			$posicion = $j;
		}
	}
	return ($valor, $posicion) ;
}

sub maxarray {
	local ($i) =0;
	local ($pos) =0;
	local (@array) = @{$_[0]};
 	local ($maxim) = $array[0];
	for ($i=0; $i <= $#array; $i++) {
		if ($maxim<=$array[$i]) {
			$maxim = $array[$i];
			$pos = $i;
		}

	}
	return ($pos, $array[$pos]) ;
}

