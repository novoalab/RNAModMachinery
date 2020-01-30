use strict;
use warnings;
my %mem=();
open A, $ARGV[0];
while (<A>) {
	chomp;
	my @a = split /\s+/,$_;
	$mem{$a[0]} = $a[1]."_".$a[2];
}

open B, $ARGV[1];

while (<B>) {
	if (/^>/) {
		my @a = split /\|/,$_;
		$a[1] = $mem{$a[1]} if exists $mem{$a[1]};
		my $l = join "|",@a;
		print ">$a[1]\n" ;
	}
	else {
		print $_;
	}
}
