#!/usr/bin/perl
open(O,">CROSS.SV");
$N  = 0;
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	push @ALL,$_;
	if(exists $hash{"$m[0]\t$m[3]"}){
		for $i (0 .. $#{$hash{"$m[0]\t$m[3]"}}){
			if(abs($hash{"$m[0]\t$m[3]"}[$i][0] - $m[1]) < 500 && abs($hash{"$m[0]\t$m[3]"}[$i][1] - $m[4]) < 500){
				$black{"$m[0]\t$m[3]\t$m[1]\t$m[4]"}=1;
				$b = $hash{"$m[0]\t$m[3]"}[$i][0];
				$e = $hash{"$m[0]\t$m[3]"}[$i][1];
				$black{"$m[0]\t$m[3]\t$b\t$e"}=1;
				$N++;
				print O $_,"\t$N\n";
				print O $BB{"$m[0]\t$m[3]\t$b\t$e"},"\t$N\n";
				last;
			}
		}
	}
	push @{$hash{"$m[0]\t$m[3]"}}, [$m[1],$m[4]];
	$BB{"$m[0]\t$m[3]\t$m[1]\t$m[4]"} = $_;
}

foreach $L (@ALL){
	chomp;
	@m = split(/\s+/,$L);
	if(exists $black{"$m[0]\t$m[3]\t$m[1]\t$m[4]"}){
		#print O $L,"\n";
		next;
	}
	print $L,"\n";
}




