#!/usr/bin/perl
#
$F = shift@ARGV;
open(I,"<$F");
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($m[0] eq $m[3] && $m[2] eq "+" && $m[5] eq "-" && abs($m[1] - $m[4]) < 7000){
		push @{$C{$m[0]}}, $_;
		$ALL{$_}=1;
	}
	$S{$_}=1;


}

open(I,"<$F");
while(<I>){
	chomp;
	@m = split(/\s+/);
	for $i (0 .. $#{$C{$m[0]}}){
		@t = split(/\s+/, $C{$m[0]}[$i]);
		if($t[1] < $m[1] && $t[4] > $m[1]){
			$B{$C{$m[0]}[$i]}=1;
		}
	}
	for $i (0 .. $#{$C{$m[3]}}){
		@t = split(/\s+/, $C{$m[3]}[$i]);
		if($t[1] < $m[4] && $t[4] > $m[4]){
			$B{$C{$m[0]}[$i]}=1;
		}
	}
}
$L = shift@ARGV;
if($L == 1){
	foreach $key (keys %S){
		if(!exists $ALL{$key}){
			print $key,"\n";
		}
		else{
			if(exists $B{$key}){
				print $key,"\n";
			}
		}
	}
}
else{
	foreach $key (keys %ALL){
		if(!exists $B{$key}){
			print $key,"\n";
		}
	}
}



