#!/usr/bin/perl
#7       108442443       +       7       108455069       +       8       3
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
		push @ALL,$_;
		if(exists $hash{"$m[0]\t$m[3]\t$m[2]\t$m[5]"}){
			for $i (0 .. $#{$hash{"$m[0]\t$m[3]\t$m[2]\t$m[5]"}}){
				if(abs($hash{"$m[0]\t$m[3]\t$m[2]\t$m[5]"}[$i][0] - $m[1]) < 20 && abs($hash{"$m[0]\t$m[3]\t$m[2]\t$m[5]"}[$i][1] - $m[4]) < 20){
					$black{"$m[0]\t$m[3]\t$m[1]\t$m[4]"}=1;
					$b = $hash{"$m[0]\t$m[3]\t$m[2]\t$m[5]"}[$i][0];
					$e = $hash{"$m[0]\t$m[3]\t$m[2]\t$m[5]"}[$i][1];
					$N++;
					last;
				}
			}
		}
		push @{$hash{"$m[0]\t$m[3]\t$m[2]\t$m[5]"}}, [$m[1],$m[4]];
		$BB{"$m[0]\t$m[3]\t$m[1]\t$m[4]"} = $_;
}

foreach $L (@ALL){
	chomp;
	@m = split(/\s+/,$L);
	if(exists $black{"$m[0]\t$m[3]\t$m[1]\t$m[4]"}){
		next;
	}
	print $L,"\n";
}




