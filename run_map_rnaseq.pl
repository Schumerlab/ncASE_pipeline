#perl! -w

if(@ARGV<6){
    print "perl run_map_rnaseq.pl current_job_file genome1.fa genome2.fa PE_or_SE bwa_or_STAR mapper_path gtf\n"; exit;
}#print usage

my $infile1=shift(@ARGV); chomp $infile1;
open IN1, $infile1 or die "cannot open list of reads\n";

my $genome1=shift(@ARGV); chomp $genome1;
my $genome2=shift(@ARGV); chomp $genome2;

my $read_type=shift(@ARGV); chomp $read_type;

my $mapper=shift(@ARGV); chomp $mapper;

my $path=shift(@ARGV); chomp $path;

my $gtf=shift(@ARGV); chomp $gtf;

while(my $line1 =<IN1>){

    chomp $line1;

    if(($read_type eq 'SE') && ($mapper eq 'bwa')){
	my $sam1 = "$line1".".par1.sam";
	my $sam2 = "$line1".".par2.sam";

	my $RG1="'"."\@RG"."\\t"."ID:hyb"."\\t"."SM:tn5"."\\t"."PL:illumina"."\\t"."LB:hyblib1"."\\t"."PU:LSIslowmode"."'";
        #print "$RG1\n";
	system("$path mem -M -R $RG1 $genome1 -t 3 $line1 > $sam1");
	system("$path mem -M -R $RG1 $genome2 -t 3 $line1 > $sam2");

    }#SE reads

    if(($read_type eq 'PE') && ($mapper eq 'bwa')){

        my @read_array=split(/\t/,$line1);

        my $read1=$read_array[0]; chomp $read1;
        my $read2=$read_array[1]; chomp $read2;

        my $sam1 = "$read1".".par1.sam";
        my $sam2 = "$read1".".par2.sam";

	my $RG1="'"."\@RG"."\\t"."ID:hyb"."\\t"."SM:tn5"."\\t"."PL:illumina"."\\t"."LB:hyblib1"."\\t"."PU:LSIslowmode"."'";
	#print "$RG1\n";
	system("$path mem -M -R $RG1 $genome1 -t 3 $read1 $read2 > $sam1");
	system("$path mem -M -R $RG1 $genome2 -t 3 $read1 $read2 > $sam2");

    }#PE reads

    if(($read_type eq 'SE') && ($mapper eq 'star')){
	my $prefix1="$line1"."_par1";
        system("$path --genomeDir ./par1 --readFilesIn $line1 --readFilesCommand zcat --outFileNamePrefix $prefix1 --sjdbGTFfile $gtf --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM");

	my $prefix2="$line1"."_par2";
	system("$path --genomeDir ./par2 --readFilesIn $line1 --readFilesCommand zcat --outFileNamePrefix $prefix2 --sjdbGTFfile $gtf --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM");

    }#SE reads STAR

    if(($read_type eq 'PE') && ($mapper eq 'star')){

	my @read_array=split(/\t/,$line1);

	my $read1=$read_array[0]; chomp $read1;
        my $read2=$read_array[1]; chomp $read2;

	if(length($gtf)>0){
	my $prefix1="$read1"."_par1";
	system("$path --genomeDir ./par1 --readFilesIn $read1 $read2 --readFilesCommand zcat --outFileNamePrefix $prefix1 --sjdbGTFfile $gtf --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM");

	my $prefix2="$read1"."_par2";
	system("$path --genomeDir ./par2 --readFilesIn $read1 $read2 --readFilesCommand zcat --outFileNamePrefix $prefix2 --sjdbGTFfile $gtf --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM");
	} else{
	    my $prefix1="$read1"."_par1";
	    system("$path --genomeDir ./par1 --readFilesIn $read1 $read2 --readFilesCommand zcat --outFileNamePrefix $prefix1 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM");
	    
	    my $prefix2="$read1"."_par2";
	    system("$path --genomeDir ./par2 --readFilesIn $read1 $read2 --readFilesCommand zcat --outFileNamePrefix $prefix2 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM");

	}#is a gtf file provided for STAR?
	
    }#PE reads STAR 

}#map each read in list
