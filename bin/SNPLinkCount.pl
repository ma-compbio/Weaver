#!/usr/bin/perl
#pull out SNP link information
#$CHR = shift@ARGV;
# $Bin/SNPLinkCount.pl SNP LINK > SNPLINK
open(I,"<",shift@ARGV);
$temp = 0;
$CHR_temp = "";
while(<I>){
	chomp;
	@m = split(/\s+/);
	#print $m[0],"\n";

	if($CHR_temp eq ""){
		$temp = 0;
		$CHR_temp = $m[0];
	}
	if($temp == 0){
		$temp = $m[1];
		next;
	}
#	print $CHR_temp,"\t",$m[0],"\n";
	if($CHR_temp eq $m[0]){
		$ALL{$m[0]}{"$temp\t$m[1]"}=1;
		$temp = $m[1];
	}
	else{
		$CHR_temp = $m[0];
		$temp = $m[1];
	}

}
open(I,"<",shift@ARGV);
while(<I>){
	chomp;
	@m = split(/\s+/);
	$chr = $m[1];

	if(!exists $ALL{$chr}){
		next;
	}
	for $i (2 .. $#m-1){
		@s = split(/:/,$m[$i]);
		@t = split(/:/,$m[$i+1]);
		if(exists $ALL{$chr}{"$s[0]\t$t[0]"}){

			if($s[1] == 0){
				$hash{$chr}{"$m[$i]\t$m[$i+1]"}++;
				$B{$chr}{"$s[0]\t$t[0]"}{"$m[$i]\t$m[$i+1]"}++;

			}
			else{
				#@t = split(/:/,$m[$i+1]);
				$k = 1-$t[1];
				$hash{$chr}{"$s[0]:0\t$t[0]:$k"}++;
				$B{$chr}{"$s[0]\t$t[0]"}{"$s[0]:0\t$t[0]:$k"}++;
			}
		}
	}
}
foreach $chr (keys %B){
	%A = %{$B{$chr}};
	foreach $key (keys %A){
		$F = 1;
		if(scalar keys %{$A{$key}} > 1){
			$sum=0;
			foreach $k (keys %{$A{$key}}){
				$sum+=$A{$key}{$k};
			}
			foreach $k (keys %{$A{$key}}){
				if($A{$key}{$k} < 3 || $A{$key}{$k}/$sum < 0.2){
					$F = 0;
					$A{$key}{$k}=0;
				}
			}
			if($F > 0){
				#print $key,"\n";#join("\t",keys %{$ALL{$key}}),"\n";	
			}
			else{
				foreach $k (keys %{$A{$key}}){
					if($A{$key}{$k} > 1){
						@sb = split(/:/, $k);
						print $chr,"\t";
						if($sb[2] == 0){
							print $sb[0],"\t+\n";
						}
						else{
							print $sb[0],"\t-\n";
						}

						#print $k,"\n";
					}
				}
			}
		}
		else{
			foreach $k (keys %{$A{$key}}){
				if($A{$key}{$k} > 1){
					@sb = split(/:/, $k);
					print $chr,"\t";
					if($sb[2] == 0){
						print $sb[0],"\t+\n";
					}
					else{
						print $sb[0],"\t-\n";
					}
				}
				#print $k,"\n";
			}
		}
	}
}
=cut
foreach $key (keys %hash){
	if($hash{$key} > 1){
		print $key,"\t",$hash{$key},"\n";
	}
}
