#!/usr/bin/perl
#modifies TCGA rnaseq metadata to be PCAWG 2.0 rnaseq metadata v0.1
#and does the submissiona and upload of the new realignments to CGHub
use strict;
use warnings;

#programs
my $CGSUBMIT='cgsubmit';
my $GTUPLOAD='gtupload';

#metadata constants
my $RUN_LABELS='<RUN_LABELS>';
my $RUN_LABELS_STOP='<\/RUN_LABELS>';
my $TARGETS='<TARGETS>';
my $TARGETS_STOP='<\/TARGETS>';
my $FAILED_READS='(<alignment_includes_failed_reads>.+<\/alignment_includes_failed_reads>)';

my $date = `date +%Y-%m-%dT%T`;
chomp($date);


#args
my $template_analysis_xml = shift;
my $file_list = shift;
die "must submit an analysis.xml template file AND a list of uuids,filenames, and md5s!\n" if(!$file_list || !$template_analysis_xml);
#optional, if not passed in, we don't submit
my $submit_key = shift;

#temporary output files, will be deleted if script finishes successfully
my $STDOUT_FILE="$file_list.stdout.$date";
my $STDERR_FILE="$file_list.stderr.$date";

sub main()
{
	open(IN,"<$file_list");
	my @id_map;
	while(my $line = <IN>)
	{
		chomp($line);
		my ($original_analysis_id,$new_filepath,$new_md5,$run_cmd,$aligner) = split(/\t/,$line);
		#dump original metadata
		run_command("dump_all_metadata.py $original_analysis_id",$STDOUT_FILE,$STDERR_FILE);
		#exract original metadata still relevant to PCAWG metadata
		my $md_lines = extract_old_metadata_elements($original_analysis_id);
		#create new metadata package from old metadata bits and template
		my $new_analysis_id = synthesize_new_analysis($md_lines,$new_filepath,$new_md5,$run_cmd,$template_analysis_xml,$original_analysis_id,$aligner);
		#do the actual validation->submission->upload
		if(validate_new_metadata($new_analysis_id) && $submit_key)
		{
			if(submit_new_metadata($new_analysis_id,$submit_key))
			{
				upload_data($new_analysis_id,$submit_key);
			}
		}
		print "successfully processed $original_analysis_id as $new_analysis_id\n";
		push(@id_map,[$original_analysis_id,$new_analysis_id]);
		#clean up extraneous files
		`rm -rf $original_analysis_id`;
		`rm $STDERR_FILE`;
		`rm $STDOUT_FILE`;
	}
	close(IN);
	open(IDS,">$file_list.id_map.tsv");
	print IDS "old_analysis_id\tnew_analysis_id\n";
	foreach my $a (@id_map)
	{
		my ($old,$new)=@$a;
		print IDS "$old\t$new\n";
	}
	close(IDS);	
}


sub validate_new_metadata()
{	
	my $analysis_id = shift;
	my $run = run_command("$CGSUBMIT --validate-only -u $analysis_id");
	print "$run\n";
}

sub submit_new_metadata()
{	
	my ($analysis_id,$submit_key) = @_;
	my $run = run_command("$CGSUBMIT -c $submit_key -u $analysis_id");
	print "$run\n";
}

sub upload_data()
{
	my ($analysis_id,$submit_key) = @_;
	my $run = run_command("$GTUPLOAD -c $submit_key -u $analysis_id/manifest.xml -vv",$STDOUT_FILE,"$analysis_id.gtupload");
	print "$run\n";
}

sub synthesize_new_analysis()
{
	my ($md_lines,$filepath,$md5,$run_cmd,$templateF,$original_analysis_id,$aligner) = @_;
    my @run_cmds = split(/\$/,$run_cmd);

    my $filename = "PCAWG.$original_analysis_id.$aligner.v1.bam";

	my $new_analysis_id = run_command('uuidgen');
	chomp($new_analysis_id);

	#create the new analysis id directory
	run_command("mkdir $new_analysis_id");
	#copy the original metadata into the new directory
	run_command("rsync -av $original_analysis_id/*.xml $new_analysis_id/");
	#link in the new realigned file into the new directory
	run_command("ln -s $filepath $new_analysis_id/$filename");
	
	open(TEMPLATE,"<$templateF");
	open(OUT,">$new_analysis_id/analysis.xml");
	while(my $line = <TEMPLATE>)
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
        my $tmp_cmd = join("\n", @run_cmds);
        $line =~ s/<NOTES>.*<\/NOTES>/<NOTES>$tmp_cmd<\/NOTES>/;
		if($line =~ /$FAILED_READS/)
		{
			my $cur_val = $1;
			my $val = $md_lines->{$FAILED_READS};
			$line =~ s/$cur_val/$val/;
		}
		print OUT "$line\n";
		if($line =~ /<\/ASSEMBLY>/)
		{
			my $labels = $md_lines->{$RUN_LABELS};
			my $length = scalar @$labels;
			for(my $i=0; $i < $length; $i++)
			{
				my $label = $labels->[$i];
				$label =~ s/^\s+/        / if($i==0 || $i == ($length - 1));
				$label =~ s/^\s+/            / if($i > 0 && $i < ($length - 1));
				print OUT "$label\n";
			}
		}
		elsif($line =~ /<\/ANALYSIS_TYPE>/)
		{
			my $targets = $md_lines->{$TARGETS};
			my $length = scalar @$targets;
			for(my $i=0; $i < $length; $i++)
			{
				my $target = $targets->[$i];
				$target =~ s/^\s+/    / if($i==0 || $i == ($length - 1));
				$target =~ s/^\s+/      / if($i > 0 && $i < ($length - 1));
				print OUT "$target\n";
			}
		}
	}
	close(TEMPLATE);
	close(OUT);
	return $new_analysis_id;
}		


sub extract_old_metadata_elements()
{
	my $uuid = shift;
	my $f = "$uuid/analysis.xml";
	my %lines;
	open(OLD,"<$f");
	my $SAVE=undef;
	while(my $line = <OLD>)
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
	close(OLD);
	return \%lines;	
}


sub run_command
{
	my ($command,$stdout_file,$stderr_file) = @_;

	$stdout_file = $STDOUT_FILE if(!$stdout_file);
	$stderr_file = $STDERR_FILE if(!$stderr_file);

	$command = "$command > $stdout_file 2>$stderr_file";
	print "running: $command\n";
	my $run = system($command);
	if($run != 0)
	{
		open(ERR,"<$stderr_file");
		my @err = <ERR>;
		close(ERR);
		chomp(@err);
		die "Command $command failed:\n".join("\n",@err)."\n";
	}
	open(INFILE,"<$stdout_file");
	my @output = <INFILE>;
	close(INFILE);
	chomp(@output);
	
	return join("\n",@output);	
}	

main();
