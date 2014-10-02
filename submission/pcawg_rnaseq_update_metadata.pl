#!/usr/bin/perl
#modifies TCGA rnaseq metadata to be PCAWG 2.0 rnaseq metadata v0.1
use strict;

my $date = `date +%Y-%m-%dT%T`;
chomp($date);

my $TEMPLATE_ANALYSIS='analysis.pcawg_rnaseq.STAR_template.xml';
my $RUN_LABELS='<RUN_LABELS>';
my $RUN_LABELS_STOP='<\/RUN_LABELS>';
my $TARGETS='<TARGETS>';
my $TARGETS_STOP='<\/TARGETS>';
my $FAILED_READS='(<alignment_includes_failed_reads>.+<\/alignment_includes_failed_reads>)';

my $analysis_xml = shift;
my $filename = shift;
my $md5 = shift;

main();

sub main()
{
	#`rsync -av $analysis_xml $analysis_xml.old`;
	my $md_lines = extract_old_metadata_elements($analysis_xml);
	synthesize_new_analysis($md_lines,$filename,$md5,$TEMPLATE_ANALYSIS,$analysis_xml);
}

sub synthesize_new_analysis()
{
	my ($md_lines,$filename,$md5,$templateF,$analysisF) = @_;
	`rsync -av $templateF $analysisF.temp`;
	
	open(IN,"<$analysisF.temp");
	open(OUT,">$analysisF");
	while(my $line = <IN>)
	{
		chomp($line);
		$line =~ s/alias="[^"]*"/alias="$filename"/;
		$line =~ s/analysis_date="[^"]*"/analysis_date="$date"/;
		$line =~ s/data_block_name="[^"]*"/data_block_name="$filename"/;
		$line =~ s/DATA_BLOCK name="[^"]*"/DATA_BLOCK name="$filename"/;
		if($line =~ /<FILE/)
		{
			$line =~ s/checksum="[^"]*"/checksum="$md5"/;
			$line =~ s/filename="[^"]*"/filename="$filename"/;
			$line =~ s/filetype="[^"]*"/filetype="bam"/;
		}
		if($line =~ /$FAILED_READS/)
		{
			my $cur_val = $1;
			my $val = $md_lines->{$FAILED_READS};
			$line =~ s/$cur_val/$val/;
		}
		print OUT "$line\n";
		if($line =~ /<\/ASSEMBLY>/)
		{
			my $lines = $md_lines->{$RUN_LABELS};
			foreach my $i (@$lines)
			{
				print OUT "$i\n";
			}
		}
		elsif($line =~ /<\/ANALYSIS_TYPE>/)
		{
			my $lines = $md_lines->{$TARGETS};
			foreach my $i (@$lines)
			{
				print OUT "$i\n";
			}
		}
	}
	close(IN);
	close(OUT);
}		
		
	


sub extract_old_metadata_elements($)
{
	my $f = shift;
	my %lines;
	open(IN,"<$f");
	my $SAVE=undef;
	while(my $line = <IN>)
	{
		chomp($line);
		if($line =~ /$FAILED_READS/)
		{
			$lines{$FAILED_READS}=$1;
		}	
		if($line =~ /$RUN_LABELS/ || $line =~ /$TARGETS/)
		{
			$SAVE=$RUN_LABELS;
		}
		if($line =~ /$TARGETS/)
		{
			$SAVE=$TARGETS;
		}
		next unless($SAVE);
		push(@{$lines{$SAVE}},$line);
		if($line =~ /$RUN_LABELS_STOP/ || $line =~ /$TARGETS_STOP/)
		{
			$SAVE=undef;
		}
	}
	close(IN);
	return \%lines;	
}
