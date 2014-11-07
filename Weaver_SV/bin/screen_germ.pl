#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
$FLAG = shift@ARGV; # 0 for number 1 for chr
$K = 0;
open(I,"<",shift@ARGV);
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($K == 0){
		$K = 1;
		if(substr($m[0],0,1) ne "c"){
			$FLAG = 1;
		}
	}
	if($FLAG == 0){
		$m[0] = substr($m[0],3);
		$m[3] = substr($m[3],3);
	}
	push @{$hash{"$m[0]\t$m[2]\t$m[3]\t$m[5]"}}, [$m[1], $m[4]];
	if($m[$#m-1] == 0 || $m[$#m] == 0){
		push @{$BLACK{$m[0]}}, $m[1];
		push @{$BLACK{$m[3]}}, $m[4];
	}
}
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
	if($m[$#m-1] == 0 || $m[$#m] == 0){
		if(exists $BLACK{$m[0]}){
			for $i (0 .. $#{$BLACK{$m[0]}}){
				if(abs($m[1] - $BLACK{$m[0]}[$i]) < 100){
					$F = 1;
					last;
				}
			}
		}
		if($F == 1){
			next;
		}
		if(exists $BLACK{$m[3]}){
			for $i (0 .. $#{$BLACK{$m[3]}}){
				if(abs($m[4] - $BLACK{$m[3]}[$i]) < 100){
					$F = 1;
					last;
				}
			}
		}
		if($F == 1){
			next;
		}
	}
	print $_,"\n";
}



