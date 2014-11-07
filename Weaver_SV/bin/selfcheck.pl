#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	$F = 0;
	if(exists $hash{"$m[0]\t$m[2]\t$m[3]\t$m[5]"}){
		for $i (0 .. $#{$hash{"$m[0]\t$m[2]\t$m[3]\t$m[5]"}}){
			if(abs($m[1] - $hash{"$m[0]\t$m[2]\t$m[3]\t$m[5]"}[$i][0]) < 100 && abs($m[4] - $hash{"$m[0]\t$m[2]\t$m[3]\t$m[5]"}[$i][1]) < 100){
				$F = 1;
				last;
			}
		}
	}
	if($F == 1){
		next;
	}
	if(!exists  $KK{"$m[0]:$m[1]:$m[2]"}){
		$KK{"$m[0]:$m[1]:$m[2]"}=1;
		$T1 = $m[1]+1;
		$S1 = $m[1]-1;
		$KK{"$m[0]:$T1:$m[2]"}=1;
		$KK{"$m[0]:$S1:$m[2]"}=1;
	}
	else{
		$BLACK{"$m[0]:$m[1]:$m[2]"}=1;
		$T1 = $m[1]+1;
		$S1 = $m[1]-1;
		$BLACK{"$m[0]:$T1:$m[2]"}=1;
		$BLACK{"$m[0]:$S1:$m[2]"}=1;
	}
	if(!exists  $KK{"$m[3]:$m[4]:$m[5]"}){
		$KK{"$m[3]:$m[4]:$m[5]"}=1;
		$T1 = $m[4]+1;
		$S1 = $m[4]-1;
		$KK{"$m[3]:$T1:$m[5]"}=1;
		$KK{"$m[3]:$S1:$m[5]"}=1;
	}
	else{
		$BLACK{"$m[3]:$m[4]:$m[5]"}=1;
		$T1 = $m[4]+1;
		$S1 = $m[4]-1;
		$BLACK{"$m[3]:$T1:$m[5]"}=1;
		$BLACK{"$m[3]:$S1:$m[5]"}=1;

	}

	######INVERSION?
	#
	#
	#
	#
	#
	################
	###

	if(!exists  $MB{"$m[0]:$m[1]"}){
		$MB{"$m[0]:$m[1]"}=1;
		$T1 = $m[1]+1;
		$S1 = $m[1]-1;
		$MB{"$m[0]:$T1"}=1;
		$MB{"$m[0]:$S1"}=1;
	}
	else{
		$MBLACK{"$m[0]:$m[1]"}=1;
		$T1 = $m[1]+1;
		$S1 = $m[1]-1;
		$MBLACK{"$m[0]:$T1"}=1;
		$MBLACK{"$m[0]:$S1"}=1;
	}
	if(!exists  $MB{"$m[3]:$m[4]"}){
		$MB{"$m[3]:$m[4]"}=1;
		$T1 = $m[4]+1;
		$S1 = $m[4]-1;
		$MB{"$m[3]:$T1"}=1;
		$MB{"$m[3]:$S1"}=1;
	}
	else{
		$MBLACK{"$m[3]:$m[4]"}=1;
		$T1 = $m[4]+1;
		$S1 = $m[4]-1;
		$MBLACK{"$m[3]:$T1"}=1;
		$MBLACK{"$m[3]:$S1"}=1;
	}


	push @{$hash{"$m[0]\t$m[2]\t$m[3]\t$m[5]"}}, [$m[1], $m[4]];
	push @FF, $_;
	#print $_,"\n";
}
for $i (0 .. $#FF){
	@m = split(/\s+/, $FF[$i]);
	if((exists $BLACK{"$m[3]:$m[4]:$m[5]"} || exists $BLACK{"$m[0]:$m[1]:$m[2]"} || exists $MBLACK{"$m[3]:$m[4]"} || exists $MBLACK{"$m[0]:$m[1]"})){# && $m[6] + $m[7] < 20){ // never tolerate rep
		#print $FF[$i],"\n";
	}
	else{
		#if((exists $MBLACK{"$m[3]:$m[4]"} || exists $MBLACK{"$m[0]:$m[1]"})){
			#print $FF[$i],"\n";

			#}
			#else{
				print $FF[$i],"\n";
			#}
	}
}




