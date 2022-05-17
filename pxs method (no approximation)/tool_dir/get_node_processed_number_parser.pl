#!/usr/bin/perl

my $main_output_dir = $ARGV[0]; # like "./output_dir/demo_results"
my $node_processed_number_filename_in = $ARGV[1];

my $node_processed_number_filename = "${main_output_dir}/$node_processed_number_filename_in";

open(outfile, ">$node_processed_number_filename") or die "Cannot open output file $node_processed_number_filename\n";

my $bm_count = 0;


for(my $group_id = 1; $group_id <= 4; $group_id++){
	for(my $bm_id = 1; $bm_id <= 12; $bm_id++){
	
		my $bm_name = "bm${bm_id}.${group_id}";
		#if($bm_name eq "bm7.1" or $bm_name eq "bm7.2"){
		#	next;
		#}
		
		$bm_count++;

		my $output_subdir = "${main_output_dir}/${bm_name}";
		
		my $input_node_processed_number_filename = "${output_subdir}/node_processed_number_report.txt";

		open(infile, "<$input_node_processed_number_filename") or die "Can't open ${input_node_processed_number_filename} for ${bm_name}\n";

		my $nodeCount;
		while(<infile>){
			chomp;
			if(/processed node number = (\d+)/){
				print "$_\n";
				$nodeCount = $1;
				print "$nodeCount\n";
				print outfile "$bm_name $nodeCount\n";
			}
		}
		close(infile);
	}
}

close(outfile);







