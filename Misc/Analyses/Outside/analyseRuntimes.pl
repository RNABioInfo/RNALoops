#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my ($dir) = @ARGV;

my %results = ();
opendir(DIR, $dir) || die "can't read directory: $!";
	while (my $file = readdir(DIR)) {
		if ($file ne '.' && $file ne '..') {
			open (FILE, $dir.'/'.$file) || die "can't read file: $!";
				my $length = -1;
				my $resType = undef;
				my $exitCode = 'fail';
				my %times = ();
				my %ss = ();
				while (my $line = <FILE>) {
					#~ print $line;
					if ($line =~ m/^sequence: (\w+)$/) {
						$length = length($1);
					} elsif ($line =~ m/^#(.+?)#$/) {
						$resType = $1;
					} elsif ($line =~ m/Exit \[0\]/) {
						$exitCode = 'ok';
					} elsif ($line =~ m/^(.+?) user, (.+?) system, /) {
						if ($exitCode eq 'ok') {
							$times{$resType} = $1+$2;
						}
						$exitCode = 'fail';
					}
				}
				$results{$length} = \%times;
			close (FILE);
			#~ print Dumper \%times;
		}
	}
closedir(DIR);
#~ print Dumper \%results; die;
my $datafile = "tmpRuntimes.data";
open (OUT, "> $datafile") || die "can't write to '$datafile': $!\n";
	my @lengths = sort {$a <=> $b} keys(%results);
	my @programs = keys(%{$results{$lengths[0]}});
	print OUT "length";
	foreach my $program (@programs) {
		print OUT "\t".$program;
	}
	print OUT "\n";
	
	foreach my $length (@lengths) {
		my $missingValue = 'false';
		foreach my $program (@programs) {
			if (not (exists $results{$length}->{$program})) {
				$missingValue = 'true';
				last;
			}
		}
		if ($missingValue eq 'false') {
			print OUT $length;
			foreach my $program (@programs) {
				print OUT "\t".$results{$length}->{$program};
			}
			print OUT "\n";
		}
	}
close (OUT);

my $pdffile = "runtimes.pdf";
my @colors = ('"red"','"orange"','"green"','"darkgreen"','"black"','"cyan"','"gray"');
open (R, " | R --vanilla");
	print R 'require(gplots)'."\n";
	print R 'pdf("'.$pdffile.'", width=10, height=7)'."\n";
	print R 'require(gplots)'."\n";
	#~ print R 'par( mfrow = c( 1, 2 ));'."\n";
	print R 'data <- read.csv("'.$datafile.'", header=TRUE, sep="\t")'."\n";
	#~ print R 'tmp <- subset(data, (data$fake > 0) & (data$mfe > 0) & (data$pfunc > 0) & (data$rnafold > 0) & (data$newmfe > 0));'."\n";
	print R 'tmp = data;'."\n";
	print R 'par(mar=c(5.1, 4.1, 0.5, 0.5));'."\n";
	print R 'plot(log="y", tmp$oa_o_nodangle_pfunc ~ tmp$length, xlab="sequence length", ylab="run-time in sec.", cex=.0)'."\n";
	#~ print R 'plot(log="y", tmp$fake ~ tmp$length, ylim=c(min(tmp$mfe, tmp$fake, tmp$pfunc, tmp$rnafold, tmp$newmfe),max(tmp$mfe, tmp$fake, tmp$pfunc, tmp$rnafold, tmp$newmfe)), xlab="sequence length", ylab="run-time in sec.", cex=.0)'."\n";
	for (my $i = 0; $i < @programs; $i++) {
		print R 'lines(tmp$'.getRname($programs[$i]).' ~ tmp$length, col='.$colors[$i].',lwd=2)'."\n";
	}
	print R 'smartlegend(x="left",y="top", inset = 0.05, c('.join(',', (map {translateNames($_)} @programs)).'), fill=c('.join(',', @colors).'), bg="white");'."\n";
	print R 'dev.off()'."\n";
close (R);

sub translateNames {
	my ($name) = @_;
	
	if ($name eq 'RNAfold') {
		return '"RNAfold -d0"';
	} elsif ($name eq 'RNAfold-p') {
		return '"RNAfold -d0 -p"';
	} elsif ($name eq 'oa_i_nodangle_mfepp') {
		return 'expression(paste("G"[inside],"(I"[mfe], " * I"[db], ")"))';
	} elsif ($name eq 'oa_o_nodangle_pfunc') {
		return 'expression(paste("G"[outside],"(I"[o_bwe], ")"))';
	} elsif ($name eq 'oa_i_nodangle_pfunc') {
		return 'expression(paste("G"[inside],"(I"[bwe], ")"))';
	} else {
		return '"'.$name.'"';
	}
}

sub getRname {
	my ($name) = @_;
	my $progname = $name;
	$progname =~ s/-/\./g;
	return $progname;
}