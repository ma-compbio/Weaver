#!/usr/bin/perl
#
#open(I,"<",STDIN);
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	if($chr_temp_1 eq $m[0] && $chr_temp_2 eq $m[3] && abs($pos_1 - $m[1]) < 25 && abs($pos_2 - $m[4]) < 25 && $ori_1 eq $m[2] && $ori_2 eq $m[5]){
		if($ori_1 ne $ori_2  && abs($pos_2 - $pos_1 - $m[4]+ $m[1])<15 || $ori_1 eq $ori_2 && abs($pos_2 - $m[4] + $pos_1 - $m[1])<15 ){
			$hash{$temp}+=$m[$#m];
			next;
		}
	}
	$chr_temp_1 = $m[0];
	$chr_temp_2 = $m[3];
	$pos_1 = $m[1];
	$pos_2 = $m[4];
	$ori_1 = $m[2];
	$ori_2 = $m[5];
	$num = $m[$#m];
	$#m = 5;
	$temp = join("\t",@m);
	$hash{$temp}=$num;
}

foreach $key (keys %hash){
	if($hash{$key} > 2){
		print $key,"\t",$hash{$key},"\n";
	}
}
