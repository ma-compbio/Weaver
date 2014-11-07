#!/usr/bin/perl
#chrM%16%21%+%SALLY:276:D1AVYACXX:5:2316:12402:89904     -       chr1    1223645
#chr4%143193514%25%-%AMELIA:245:D1GH8ACXX:6:2310:14940:6745      +       chr4    143193489
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	@p = split(/%/,$m[0]);
	#substr($p[0],0,1) = 'c';
	if($p[3] eq "+"){
		if($m[1] eq "+"){
			$temp = $m[3];
			$L2 = "$m[2]\t$temp\t+";
		}
		else{
			$temp = $m[3]+$p[2]-1;
			$L2 = "$m[2]\t$temp\t-";
		}
	}
	else{
		if($m[1] eq "+"){
			$temp = $m[3]+$p[2]-1;
			$L2 = "$m[2]\t$temp\t-";
		}
		else{
			$temp = $m[3];
			$L2 = "$m[2]\t$temp\t+";
		}
	}
	if($p[0] eq $m[2]){
		if(abs($p[1] - $temp) < 5){
			next;
		}
	}
	if($p[3] eq "+"){
		$L1 = "$p[0]\t$p[1]\t-";## - means left
	}
	else{
		$L1 = "$p[0]\t$p[1]\t+";## extend to right
	}
	if($p[0] eq $m[2]){
		if($p[1] < $temp){
			$hash{"$L1\t$L2"} ++;
		}
		else{
			$hash{"$L2\t$L1"} ++;
		}
	}
	elsif($p[0] gt $m[2]){
		$hash{"$L1\t$L2"} ++;
	}
	else{
		$hash{"$L2\t$L1"} ++;
	}
}

foreach $key (keys %hash){
	#if($hash{$key} > 2){
		print $key,"\t",$hash{$key},"\n";
		#}
}


