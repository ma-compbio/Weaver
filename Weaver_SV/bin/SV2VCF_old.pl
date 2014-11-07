#!/usr/bin/perl
#
#1       109494847       -       3       110413400       +       120     0
$REF = shift@ARGV;
open(I,"<",$REF);# reference sequence
while(<I>){
	chomp;
	if(substr($_,0,1) eq ">"){
		if($S ne ""){
			$Shash{$chr} = $S;
			$S = "";
		}
		$chr = substr($_,1);
		next;
	}
	$S .= $_;
}
$Shash{$chr} = $S;
$VERSION = 1;
print "##fileformat=VCFv4.1\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
print "##source=AEGIS_V$VERSION\n##reference=file:$REF\n";
print "##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n";
print "##FILTER=<ID=LowQual,Description=\"PE support below 3 or mapping quality below 20.\">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">
##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">
##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">
##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">
##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SPIKEIN
#1       1146618 .       C       <DUP>   100     PASS    IMPRECISE;CIPOS=-100,100;CIEND=-100,100;SOMATIC;SVTYPE=DUP;END=1159763;SVLEN=-47780     GT      ./.
open(I,"<",shift@ARGV);## SV file
while(<I>){
	chomp;
	@m = split(/\s+/);
	$N++;
	$BASE1 = substr($Shash{$m[0]}, $m[1], 1);
	$BASE2 = substr($Shash{$m[3]}, $m[4], 1);
	$SVLEN = $m[4] - $m[1];
	if($m[2] eq "+" && $m[5] eq "-" && $m[0] eq $m[3]){
		print "$m[0]\t$m[1]\t.\t$BASE1\t<DEL>\t100\tPASS\t";
	#	if($m[7] > 0){
	#		print "PRECISE;SOMATIC;END=$m[4];SVLEN=$SVLEN\n"
	#	}
	#	else{
			print "IMPRECISE;CIPOS=-100,100;CIEND=-100,100;SOMATIC;SVTYPE=DEL;END=$m[4];SVLEN=$SVLEN\n"
	#	}
	}
	if($m[2] eq "-" && $m[5] eq "+" && $m[0] eq $m[3]){
		print "$m[0]\t$m[1]\t.\t$BASE1\t<DUP>\t100\tPASS\t";
	#	if($m[7] > 0){
	#		print "PRECISE;SOMATIC;END=$m[4];SVLEN=-$SVLEN\n"
	#	}
	#	else{
			print "IMPRECISE;CIPOS=-100,100;CIEND=-100,100;SOMATIC;SVTYPE=DUP;END=$m[4];SVLEN=-$SVLEN\n"
	#	}
	}
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
                                        #$BASE1 = substr($Shash{$m[0]}, $m[1], 1);
                                        print "$m[0]\t$m[1]\t.\t$BASE1\t<INV>\t100\tPASS\t";
                                        $SVLEN = $m[4] - $m[1];
                                        print "IMPRECISE;CIPOS=-100,100;CIEND=-100,100;SOMATIC;SVTYPE=INV;END=$m[4];SVLEN=0\n";
                                        #print O $_,"\t$N\n";
                                        #print O $BB{"$m[0]\t$m[3]\t$b\t$e"},"\t$N\n";
                                        last;
                                }
                        }
                }
                push @{$hash{"$m[0]\t$m[3]"}}, [$m[1],$m[4]];
                $BB{"$m[0]\t$m[3]\t$m[1]\t$m[4]"} = $_;
        }
=cut
	if($m[2] eq "+" && $m[5] eq "+"){
		$M = $N+1;
		print "$m[0]\t$m[1]\tbnd_$N\t$BASE1\t$BASE1\]$m[3]:$m[4]\]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$M\n";
		print "$m[3]\t$m[4]\tbnd_$N\t$BASE2\t$BASE2\]$m[0]:$m[1]\]\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$N\n";
	}
	if($m[2] eq "+" && $m[5] eq "-"){
		$M = $N+1;
		print "$m[0]\t$m[1]\tbnd_$N\t$BASE1\t$BASE1\[$m[3]:$m[4]\[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$M\n";
		print "$m[3]\t$m[4]\tbnd_$N\t$BASE2\t\]$m[0]:$m[1]\]$BASE2\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$N\n";
	}
	if($m[2] eq "-" && $m[5] eq "+"){
		$M = $N+1;
		print "$m[3]\t$m[4]\tbnd_$N\t$BASE2\t$BASE2\[$m[0]:$m[1]\[\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$M\n";
		print "$m[0]\t$m[1]\tbnd_$N\t$BASE1\t\]$m[3]:$m[4]\]$BASE1\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$N\n";
	}
	if($m[2] eq "-" && $m[5] eq "-"){
		$M = $N+1;
		print "$m[0]\t$m[1]\tbnd_$N\t$BASE1\t\[$m[3]:$m[4]\[$BASE1\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$M\n";
		print "$m[3]\t$m[4]\tbnd_$N\t$BASE2\t\[$m[0]:$m[1]\[$BASE2\t6\tPASS\tSVTYPE=BND;MATEID=bnd_$N\n";
	}
	$N = $M;
=cut
}

