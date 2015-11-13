#!/usr/bin/perl

#Including Libraries

use File::Copy;
use Time::HiRes qw(gettimeofday);

#Smoothing
copy("./out/Solution.out","./out/Solution.old");
open(OLD,"out/Solution.old");
open(NEW,">out/Solution.out");
@VOLD=<OLD>;
	$i=0;
	$_=$VOLD[$i];
	@array1=/([+-\d.Ee]+)/g;
	$_=$VOLD[$i+1];
	@array2=/([+-\d.Ee]+)/g;
	for($j=0;$j<$#array2;$j++) {
		print NEW sprintf("%8.7e ",($array1[$j]+$array2[$j])/2); 
	}
	print NEW sprintf("%8.7e\n",($array1[$j])); 
for($i=1;$i<$#VOLD; $i++){
	$_=$VOLD[$i-1];
	@array1=/([+-\d.Ee]+)/g;
	$_=$VOLD[$i];
	@array2=/([+-\d.Ee]+)/g;
	$_=$VOLD[$i+1];
	@array3=/([+-\d.Ee]+)/g;
	for($j=0;$j<$#array2;$j++) {
		print NEW sprintf("%8.7e ",(($array1[$j]+$array2[$j]+$array3[$j])/3)); 
	}
	print NEW sprintf("%8.7e\n",(($array1[$j]+$array2[$j]+$array3[$j])/3)); 
}
$_=$VOLD[$i-1];
@array1=/([+-\d.Ee]+)/g;
$_=$VOLD[$i];
@array2=/([+-\d.Ee]+)/g;
for($j=0;$j<$#array1;$j++) {
	print NEW sprintf("%8.7e ",($array1[$j]+$array2[$j])/2); 
}
print NEW sprintf("%8.7e\n",($array1[$j]+$array2[$j])/2); 

close(OLD);
close(NEW);



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
