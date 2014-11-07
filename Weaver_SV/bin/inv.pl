#!/usr/bin/perl
$REF = shift@ARGV;
open(I,"<",$REF);# reference sequence
while(<I>){
        chomp;
        if(substr($_,0,1) eq ">"){
                if($S ne ""){
                        $hash{$chr} = $S;
                        $S = "";
                }
                $chr = substr($_,1);
                next;
        }
        $S .= $_;
}
$hash{$chr} = $S;

open(O,">INV.SV.vcf");
$N  = 0;
#7       108442443       +       7       108455069       +       8       3
while(<STDIN>){
	chomp;
	@m = split(/\s+/);
	if($m[0] eq $m[3] && $m[2] eq $m[5]){
		push @ALL,$_;
		if(exists $hash{"$m[0]\t$m[3]"}){
			for $i (0 .. $#{$hash{"$m[0]\t$m[3]"}}){
				if(abs($hash{"$m[0]\t$m[3]"}[$i][0] - $m[1]) < 500 && abs($hash{"$m[0]\t$m[3]"}[$i][1] - $m[4]) < 500){
					$black{"$m[0]\t$m[3]\t$m[1]\t$m[4]"}=1;
					$b = $hash{"$m[0]\t$m[3]"}[$i][0];
					$e = $hash{"$m[0]\t$m[3]"}[$i][1];
					$black{"$m[0]\t$m[3]\t$b\t$e"}=1;
					$N++;
					$BASE1 = substr($hash{$m[0]}, $m[1], 1);
					print O "$m[0]\t$m[1]\t.\t$BASE1\t<INV>\t100\tPASS\t";
					$SVLEN = $m[4] - $m[1];
					print O "IMPRECISE;CIPOS=-100,100;CIEND=-100,100;SOMATIC;END=$m[4];SVLEN=$SVLEN\n";
					#print O $_,"\t$N\n";
					#print O $BB{"$m[0]\t$m[3]\t$b\t$e"},"\t$N\n";
					last;
				}
			}
		}
		push @{$hash{"$m[0]\t$m[3]"}}, [$m[1],$m[4]];
		$BB{"$m[0]\t$m[3]\t$m[1]\t$m[4]"} = $_;
	}
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




