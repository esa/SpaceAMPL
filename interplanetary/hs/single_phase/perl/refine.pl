#!/usr/bin/perl

#Including Libraries

use File::Copy;
use Time::HiRes qw(gettimeofday);

#Creating a temporary file as copy
copy("./out/Solution.out","./out/Solution.old");
open(OLD,"out/Solution.old");
open(NEW,">out/Solution.out");
@VOLD=<OLD>;
for($i=0;$i<$#VOLD; $i++){
	$_=$VOLD[$i];
	@array1=/([+-\d.Ee]+)/g;
	$_=$VOLD[$i+1];
	@array2=/([+-\d.Ee]+)/g;
	for($j=0;$j<$#array1;$j++) {
		print NEW sprintf("%8.7e ",($array1[$j])); 
	}
	print NEW sprintf("%8.7e\n",($array1[$j])); 
	
	for($j=0;$j<$#array2;$j++) {
		print NEW sprintf("%8.7e ",(($array1[$j]+$array2[$j])/2)); 
	}
	print NEW sprintf("%8.7e\n",(($array1[$j]+$array2[$j])/2)); 
}
$_=$VOLD[$i];
@array1=/([+-\d.Ee]+)/g;
for($j=0;$j<$#array1;$j++) {
	print NEW sprintf("%8.7e ",($array1[$j])); 
}
print NEW sprintf("%8.7e\n",($array1[$j])); 

close(OLD);
close(NEW);

#Creating a temporary file as copy
copy("./main.mod","./main.bak");

open(MAINOLD,"main.bak");
open(MAINNEW,">main.mod");
foreach (<MAINOLD>){
		    if (/param n/) { 
		     /(\d+)/g;
		     $_="param n:=".($1*2-1).";\n";
		    }
		    elsif ((/include/) & (/guess/)) {
		      $_="include include/guesscontinue.inc;\n";
		    }
		    print MAINNEW $_;
}

close(MAINOLD);
close(MAINNEW);
system("ampl main.mod");



sub imprimearray {
	local (@array) = @{$_[0]};
	local ($variable);
	local ($i);
	for ($i=0; $i <= $#array; $i++) {
		$variable= $array[$i];
		if ($i == $#array) {
				print  $guesses[$variable]."\n";
		}
		else {
				print  $guesses[$variable]."  ";
		}
	}
}
