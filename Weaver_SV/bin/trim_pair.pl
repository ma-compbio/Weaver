#!/usr/bin/perl
#
#open(I,"<",shift@ARGV);
$CUT = shift@ARGV;
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	$F = 1;
	if(exists $hash{"$m[0]"}){
		for $i (0 .. $#{$hash{"$m[0]"}}){
			if(abs($m[1] - $hash{"$m[0]"}[$i]) < 1000){
				$F = 0;
				push @{$hash{"$m[0]"}}, $m[1];
				push @{$hash{"$m[3]"}}, $m[4];
				last;
			}
		}
	}
	if(exists $hash{"$m[3]"}){
		for $i (0 .. $#{$hash{"$m[3]"}}){
			if(abs($m[4] - $hash{"$m[3]"}[$i]) < 1000){
				push @{$hash{"$m[0]"}}, $m[1];
				push @{$hash{"$m[3]"}}, $m[4];
				$F = 0;
				last;
			}
		}
	}
	if($F == 0){
		next;
	}

	push @{$hash{"$m[0]"}}, $m[1];
	push @{$hash{"$m[3]"}}, $m[4];
	if($m[$#m-1] > $CUT){ # 10?
		print $_,"\n";
	}
}



