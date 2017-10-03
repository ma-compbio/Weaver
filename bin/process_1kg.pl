#!/usr/bin/perl
#
open(I,"<",shift@ARGV);
#1       10583   rs58108140      G       A       100     PASS    AVGPOST=0.7707;RSQ=0.4319;LDAF=0.2327;ERATE=0.0161;AN=2184;VT=SNP;AA=.;THETA=0.0046;AC=314;SNPSOURCE=LOWCOV;AF=0.14;ASN_AF=0.13;AMR_AF=0.17;AFR_AF=0.04;EUR_AF=0.21        GT:DS:GL        0|0:0.200:-0.18,-0.47,-2.42
while(<I>){
    chomp;
    if(substr($_,0,1) eq "#"){
        next;
    }
    @m = split(/\s+/);
    if(length($m[3]) != 1 || length($m[4]) != 1){
        next;
    }
    print $m[0],"\t",$m[1],"\t",$m[3],"\t",$m[4];
    for $i (9 .. $#m){
        #print $#m,"\n";
        @s = split(/\||:/, $m[$i]);
        print "\t",$s[0],":",$s[1];
    }
    print "\n";

}
