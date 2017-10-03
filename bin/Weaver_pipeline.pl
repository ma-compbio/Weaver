#!/usr/bin/perl
#
#
######
use Getopt::Long ;
use Getopt::Long qw(:config no_ignore_case);
use POSIX;
use FindBin qw($Bin);
use Parallel::ForkManager;
use Cwd 'abs_path';

$RUN_TYPE = shift@ARGV;
$MODEFLAG = shift@ARGV; # SNP:SV:WIG

print $RUN_TYPE,"\n";
GetOptions(
        #MANDATORY
        #OPTIONAL
        'p|thread=i'=>\$P,
        'g|gap=s'=>\$GAP, # with chr [MANDATORY]
	'b|bam=s'=>\$BAM, # [MANDATORY]
	'f|fa=s'=>\$FA, # no .fa [MANDATORY]
	'F|FullFa=s'=>\$FULLFA, # .fasta or .fa
        'h|help' =>\$help,
        'o|output=s'=>\$OUT_DIR,
	'k|onekg=s'=>\$ONEKG, # dir [MANDATORY]
	't=s'=>\$TorN,# Tumor or Normal
	's|sex=s'=>\$SEX, # M or F
        'C=i'=>\$cov);

#################################
$VERSION = 1.0;
#################################

if($BAM eq "" || $help || $FA eq "" || $ONEKG eq ""){
        if(defined $FA && ! -e "$FA.fa"){
                print "$FA.fa does not exist\n";
        }
        print "Weaver v$VERSION\nUsage:\n
        -p/--thread             number of cores
        -f/--fa			[MANDATORY] bowtie and bwa reference dir/name
        -g/--gap		[MANDATORY] bowtie reference dir/name, thus bowtie index should be dir/name.*.ebwt, reference genome should also be located here, e.g. dir/name.fa; chromosome line for dir/name.fa should be clean and there is no space within it, such as \">XXXXX\\n\"
        -b/--bam                bam file
        -o/--output             output dir
	-k/--onekg		1000 Gemomes Project data dir
	-s/--sex		Female (F) or Male (M). Y chromosome will not be used if the bam is from female tissue.
        -h/--help
        \n\n";

        exit(0);
}

print "Weaver\t",$RUN_TYPE,"\n";


if(!(-e $FULLFA)){
	if(-e "$FA.fa"){
		$FULLFA = "$FA.fa";
	}
	elsif(-e "$FA.fasta"){
		$FULLFA = "$FA.fasta";
	}
	else{
		die "fasta reference not found! Provide by -F \n";
	}
}


$P = $P || 50; # default thread
$TorN = $TorN || "T"; # default tumor sample
$SEX = $SEX || "F"; # default female Y discard
$GAP = $GAP || "$Bin/../data/GAP_20140416";

#die "$Bin/../../external_bin directory not working!" unless -x "$Bin/../../external_bin/bowtie";

die "GAP file not found!" unless -e $GAP;

die "bam file not exist!" unless -e $BAM;



open(RUNLOG,">>RUNLOG");
$NOW = time;
$now_string = localtime;
print RUNLOG "$now_string:\tStart\n";
$bwt = abs_path($bwt);
## SNP SNPLINK 1000_LINK
#
if ($MODEFLAG  =~ /SNP/){
$now_string = localtime;
print RUNLOG "$now_string:\tSNP\n";

#-------------------
#
#Germline point mutations and number of reads mapped on
#
#-------------------
#if($RUN_TYPE == "lite"){
#	system("$Bin/pipe.pl $RUN_TYPE $GAP $BAM $FULLFA $P $SEX");
#}
#else{

system("$Bin/pipe.pl $RUN_TYPE $GAP $BAM $FULLFA $P $ONEKG $SEX $TorN");


#}
$now_string = localtime;
print RUNLOG "$now_string:\tSNP done\n";
}
## SV finding
if ($MODEFLAG  =~ /SV/){
$now_string = localtime;
print RUNLOG "$now_string:\tSV\n";
system("$Bin/../Weaver_SV/bin/Weaver_SV.pl $BAM $FA $FULLFA $P");
$now_string = localtime;
print RUNLOG "$now_string:\tSV done\n";
}
## get wig/bw file
$now_string = localtime;
#-------------------
#
#Generate wiggle file which will be used in Weaver core program as read depth
#
#
#-------------------
if ($MODEFLAG  =~ /WIG/){
	print RUNLOG "$now_string:\twig\n";
	system("$Bin/bam2bw.pl $RUN_TYPE $BAM $P $SEX");
	$now_string = localtime;
	print RUNLOG "$now_string:\twig done\n";
}

