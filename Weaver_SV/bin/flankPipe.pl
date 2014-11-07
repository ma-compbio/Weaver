#!/usr/bin/perl
#
system("cat $BARCODE.ALL_SV | $Bin/flankSV.pl | sort -k 1,1 -k 2,2n");
