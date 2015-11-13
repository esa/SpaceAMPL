#!/usr/bin/perl

	open(TIME,"out/Times.out");
	foreach (<TIME>){
	    @time=/([+-\d.Ee]+)/g;
	    last;
	}

	close(TIME);

	open(MASS,"out/mass.out");
	foreach (<MASS>){
	    /([+-\d.Ee]+)/;
	    $objfunguess=$1;
	    last;
	}
	close(MASS);

	print "DEPARTURE ". $time[0]." \n";
	print "PERIOD ". $time[1]." \n";
	print "ARRIVAL ". ($time[0]+$time[1])."\n";
	print "OBJFUN ".$objfunguess."\n";
