#!/usr/bin/perl

#Including Libraries

use File::Copy;
use Time::HiRes qw(gettimeofday);

#Guesses, Corrections, Locks

@guesses=("guesslambert.inc","guesstangential2m.inc");
@corrections=(0,-50,+50);
$lock1=-1;
$lock2=-1;
$lock3=-1;

#Getting the inputs

open(INPUT,"./input.dat");
@inp=<INPUT>;
$_=$inp[0];
@asteroids=/(\d+)/g;
$_=$inp[1];
@data=/([+-\d.Ee]+)/g;
close(INPUT);


#Modifying variables to complete the sequence
$asteroids[3]=$asteroids[2];
$asteroids[2]=$asteroids[1];
$asteroids[1]=$asteroids[0];
$asteroids[4]=141;
$asteroids[0]=141;

#Initialising mass, minimum waiting time and other variables
$mass=1;
$min=3652.5;

#Creating folder structure
system("mkdir ./out");
system("mkdir ./out/".$asteroids[1]."_".$asteroids[2]."_".$asteroids[3]);
system("mkdir ./out/".$asteroids[1]."_".$asteroids[2]."_".$asteroids[3]."/A0_A1");
system("mkdir ./out/".$asteroids[1]."_".$asteroids[2]."_".$asteroids[3]."/A1_A2");
system("mkdir ./out/".$asteroids[1]."_".$asteroids[2]."_".$asteroids[3]."/A2_A3");
system("mkdir ./out/".$asteroids[1]."_".$asteroids[2]."_".$asteroids[3]."/A3_A4");

#Opening OUTPUT files
open(OUTPUT,">./out/".$asteroids[1]."_".$asteroids[2]."_".$asteroids[3]."/outputTOTAL.dat");
open(OUTPUT2,">./out/".$asteroids[1]."_".$asteroids[2]."_".$asteroids[3]."/prune.dat");

#Creating a temporary file as copy
copy("./main.mod","./main.bak");

#Writing Sequence in the Output file
print OUTPUT "Sequence ".$asteroids[1]." ".$asteroids[2]." ".$asteroids[3]."\n\n";

#Loop! Each call modifies main.mod, and launches a full execution of ampl
for($j=0;$j<=3;$j++) {
   print OUTPUT "PHASE A".$j."_A".($j+1)."\n";
   for($k=0;$k<=$#guesses;$k++) { 
	
	@realbest=($l, $m);
	$realmax=0;
	for($l=0;$l<=$#corrections;$l++) {
		for($m=0;$m<=$#corrections;$m++) {
		if (($lock1==$j)||($lock2==$j)||($lock3==$j)) {
			if (($m!=0)||($l!=0)){
				last;
			}
		}
	
		#The variables mink, totalk and $objfunguess (mass) are evaluated for each of the guesses, they are kept to be stored when the best guess is selected
		$mink[$k]=$min;
		$totalk=$total;
		open(MAINOLD,"main.bak");
		open(MAINNEW,">main.mod");
		foreach (<MAINOLD>){
		    if (/param depart/) { 
		      $_="param depart:='ast".$asteroids[$j]."' symbolic;\n";
		    }
		    elsif (/param target/) {
		      $_="param target:='ast".$asteroids[$j+1]."' symbolic;\n";
		    }
		    elsif (/param M/) {
		     $_="param M	     :=   (".($mass*2000).");	\n";
		    }
		    elsif (/param C3dep/) {
		      if ($j==0) {	
		      $_="param C3dep     :=   (0.5**2)  /V**2;	#Initial Launch C3 (km^2/sec^2)\n";
		      }
		      else {
		      $_="param C3dep     :=   (0.001**2)  /V**2;	#Initial Launch C3 (km^2/sec^2)\n";
		      }
		    }
		    elsif (/param C3arr/) {
		      $_="param C3arr     :=   (0.001**2)  /V**2;	#Arrival C3 allowed (km^2/sec^2)\n";
		    }
		    elsif (/param tI/) {
		      $_="param tI	     :=   ".($data[$j*2]+$corrections[$l]).";		#Departure time (MJD2000)\n";
		    } 
		    elsif (/param tT/) {
		      $_="param tT	     :=   ".($data[$j*2+1]+$corrections[$m]).";		#Arrival time (MJD2000)\n";
		    }
		    if ((/include/) & (/guess/)) {
		      $_="include include/".$guesses[$k].";\n";
		    }
		    print MAINNEW $_;
		}
		close(MAINOLD);
		close(MAINNEW);
		system("ampl main.mod");
		system("cp main.mod ./out/".$asteroids[1]."_".$asteroids[2]."_".$asteroids[3]."/A".$j."_A".($j+1));
		system("cp out/*  ./out/".$asteroids[1]."_".$asteroids[2]."_".$asteroids[3]."/A".$j."_A".($j+1));
		
		open(TIME,"out/Times.out");
		foreach (<TIME>){
		    @time=/([+-\d.Ee]+)/g;
		    last;
		}
	
		close(TIME);
	
		$waiting=$time[0]-$total;
		$totalk[$k]=$time[0]+$time[1];
	
		open(MASS,"out/mass.out");
		foreach (<MASS>){
		    /([+-\d.Ee]+)/;
		    $objfunguess[$k]=$1;
		    last;
		}
		close(MASS);
	
		print OUTPUT "\n* GUESS  ".$guesses[$k]." \n";
		print OUTPUT "DEPARTURE ". $time[0]." \n";
		print OUTPUT "PERIOD ". $time[1]." \n";
		print OUTPUT "ARRIVAL ". $totalk[$k]." \n";
		print OUTPUT "INITIAL MASS ".(($mass*2000)/1)."\n";
		if ($min>=$waiting) {
			$mink[$k]=$waiting;
		}	
		
		print OUTPUT "MIN ". $mink[$k]." \n";
		print OUTPUT "MASS ".(($mass*$objfunguess[$k]*2000)/1)."\n";
		
	
		$massk[$k]=$objfunguess[$k];

		$objfunguess[$k]=$mass*$objfunguess[$k]+0.2*$mink[$k]/3652.5/(4-$j);
		
		if ($waiting<60) {
			print OUTPUT "ERROR: DOES NOT WAIT ENOUGH\n";
			$objfunguess[$k]=$objfunguess[$k]-0.5;
		}	
	
		if (&infeasible) {
		print OUTPUT "ERROR: INFEASIBLE\n";
		$objfunguess[$k]=$objfunguess[$k]-0.5;
		print OUTPUT2 "$k $corrections[$l] $corrections[$m] $time[0] $totalk[$k] 0 \n";
		}
		else{		
		print OUTPUT2 "$k $corrections[$l] $corrections[$m] $time[0] $totalk[$k] $massk[$k] \n";
		}

	 	if ($j==3) {
			$cond= (($data[0] + 3625.5) > $totalk[$k]);
          	 }
 	  	else {
			$cond= (($data[$j*2+2])>($totalk[$k]+60));
	  	 }

		
		if (!($cond)) {
			$objfunguess[$k]=$objfunguess[$k]-0.05;
			print OUTPUT "ERROR: POSSIBLE MISMATCHING WITH NEXT PHASE\n";
		}

		print OUTPUT "RELATIVE OBJFUN ".($objfunguess[$k])."\n";
		print OUTPUT "PARTIAL MASS OBJFUN ".($massk[$k])."\n";


		if ($objfunguess[$k]>$realmax) {
			$realmax=$objfunguess[$k];
			@realbest=($l, $m);
			$realmin= $mink[$k];
			$realmass=$massk[$k];
			$realtotal=$totalk[$k];
		}
		print OUTPUT "$corrections[$l], $corrections[$m] \n";
       	   }
	}
	
	$totalk[$k]=$realtotal;
	$massk[$k]=$realmass;
	$mink[$k]=$realmin;
	$objfunguess[$k]=$realmax;
	print OUTPUT "$corrections[$realbest[0]], $corrections[$realbest[1]] \n";
	print OUTPUT "=========================================================\n\n";
        }

@x=&maxarray(\@objfunguess);
$tempmass=$x[1];
$bestguess[$j]=$x[0];
$total=$totalk[$bestguess[$j]];
$min=$mink[$bestguess[$j]];
$mass=$mass*$massk[$bestguess[$j]];

print OUTPUT "\nBEST GUESS ".$guesses[$bestguess[$j]]."\n\n";

}

print OUTPUT "GUESS SEQUENCE : ";
&imprimearray(\@bestguess);


$objfun=$mass+0.2*$min/3652.5;
print OUTPUT "OBJECTIVE FUNCTION ".$objfun."\n";
print "=============================================================================\n";
print "=========================OBJECTIVE FUNCTION VALUE============================\n";
print "============================".$objfun."======================================\n";
close(OUTPUT);
close(OUTPUT2);

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
    for($p=0;$p<=40;$p++) {
	if ($LOG[$#LOG-20-$p] =~ /\A\s*\d+\s/) {
		last;
	}
    }
    if (($LOG[$#LOG-20-$p] =~ /i/) ||(!($LOG[$#LOG-20-$p] =~ /\([+-\d.Ee]+\)/))){
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
				print OUTPUT $guesses[$variable]."\n";
		}
		else {
				print OUTPUT $guesses[$variable]."  ";
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

