#!/usr/bin/perl
#
$FLAG = shift@ARGV; # 1 for chr 0 for number
$Chr = shift@ARGV;

open(I,"<",shift@ARGV);# read in SNP file from cancer genome sequencing data
if($FLAG == 0){
	$C_temp = substr($Chr, 3);
}
else{
	$C_temp = $Chr;
}
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($m[0] ne $C_temp){
		next;
	}
	$hash{$m[0]}{$m[1]} = "$m[4]";
	if($m[0] eq $temp_chr){
		$LINK{$m[0]}{"$temp_pos\t$m[1]"}=1;
	}
	$temp_chr = $m[0];
	$temp_pos = $m[1];
}

$temp_chr = "";
$temp_pos = "";
open(I,"<",shift@ARGV);# read in 1000G file
while(<I>){
	chomp;
	@m = split(/\s+/);
	if($FLAG == 1){
		$chr = "chr".$m[0];
	}
	else{
		$chr = $m[0];
	}
	if(exists $hash{$chr}{$m[1]}){
		if($hash{$chr}{$m[1]} eq $m[3]){
			if(exists $LINK{$chr}{"$temp_pos\t$m[1]"}){
				@s = split(/\s+/,$SAVE);
				$AB1 = 0;
				$AB2 = 0;
				$AB3 = 0;
				$AB4 = 0;
				for $j (4 .. $#s){
					@k = split(/:/, $s[$j]);
					@p = split(/:/, $m[$j]);
					if($k[0] == 0 && $p[0] == 0 || $k[1] == 0 && $p[1] == 0){
						$AB1 ++;
					}
					if($k[0] == 1 && $p[0] == 1 || $k[1] == 1 && $p[1] == 1){
						$AB2 ++;
					}
					if($k[0] == 0 && $p[0] == 1 || $k[1] == 0 && $p[1] == 1){
						$AB3 ++;
					}
					if($k[0] == 1 && $p[0] == 0 || $k[1] == 1 && $p[1] == 0){
						$AB4 ++;
					}
				}
				$n1 = $AB1 * $AB2;
				$n2 = $AB3 * $AB4;
				#print $n1,"\t",$n2,"\n";
				if($n1 == 0 && $n2 == 0){
					$R = ($AB1+$AB2)/($AB1+$AB2+$AB3+$AB4);
				}
				else{
					$R = $n1/($n1+$n2);
				}
				print $chr,"\t",$temp_pos,"\t",$m[1],"\t",$R,"\n";
			}
			$temp_pos = $m[1];
			$SAVE = $_;
		}
	}
}

