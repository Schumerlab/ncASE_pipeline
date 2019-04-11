#perl! -w

my $config=shift(@ARGV); chomp $config;
open CONFIG, $config or die "cannot open configuration file\n";

my $genome1=""; my $genome2=""; my $read_type=""; my $read_list=""; my $read_length=""; my $number_indiv_per_job=""; my $gtf=""; my $mapper=""; my $gtf_status=""; my $aims_status=""; my $map_path=""; my $sam_path=""; my $allow_zero_counts=""; my $bias_threshold=""; my $skip=0; my $bedtools_path=""; my @jobs=();

my $num_jobs=""; my $job1_submit=""; my $job2_submit=""; my $job3_submit="";

my $provide_AIMs=""; my $provide_counts="";

my $save_files=0;

while (my $line=<CONFIG>){

    chomp $line;
    my @elements=split(/\=/,$line);
    
    if($line =~ /genome1=/g){
        $genome1=$elements[1]; chomp $genome1;
        print "parent genome 1 is $genome1\n";
    }#define genome1
    if($line =~ /genome2=/g){
        $genome2=$elements[1]; chomp $genome2;
        print "parent genome 2 is $genome2\n";
    }#define genome2
    if($line =~ /skip_genome_index/g){
	$skip=$elements[1]; chomp $skip;
    }#record whether to index genome
    if($line =~ /mapping_program/g){
	$mapper=$elements[1]; chomp $mapper;
	if(($mapper eq 'bwa') or ($mapper eq 'BWA')){
	    print "mapping RNAseq reads with bwa\n";
	    $mapper="bwa";
	} elsif(($mapper eq 'star') or ($mapper eq 'STAR')){
	    print "mapping RNAseq reads with star\n";
	    $mapper="star";
	} else{
	    die "mapper must be bwa or star\n";
	}#mapper must be one of the approved mapper
    }#define mapper
    if($line =~ /path_to_mapper/g){
	$map_path=$elements[1]; chomp $map_path;
	if(length($map_path) eq 0){
	    print "no mapper path provided, assuming global install\n";
	}#no path provided
    }#path to mapper
    if($line =~ /path_to_samtools/g){
	$sam_path=$elements[1];chomp $sam_path;
	if(length($sam_path) eq 0){
            print "no samtools path provided, assuming global install\n";
        }#no path provided 
    }#path to samtools
    if($line =~/path_to_bedtools/g){
	$bedtools_path=$elements[1]; chomp $bedtools_path;
	if(length($bedtools_path) eq 0){
	    print "no bedtools path provided, assuming global install\n";
	}#no path for bedtools
    }#path to bedtools
    if($line =~/gtf_file/g){
	$gtf=$elements[1]; chomp $gtf;
	if(length($gtf) eq 0){
	    $gtf_status=0; print "gtf file not provided\n";
	} elsif(length($gtf) > 0){
	    $gtf_status=1; print "gtf file is $gtf\n";
	}#check gft status
    }#define gtf file
    if($line =~ /allow_zero_counts/g){
	$allow_zero_counts=$elements[1]; chomp $allow_zero_counts;
    }#define counts limit
    if($line =~ /allelic_bias_threshold/g){
	$bias_threshold=$elements[1]; chomp $bias_threshold;
    }#define bias threshold;
    if($line =~ /read_type/g){
        $read_type=$elements[1]; chomp $read_type;
        if(($read_type ne 'SE') && ($read_type ne 'PE')){
            die "read type must be SE or PE\n";
        }
        print "read type is $read_type\n";
    }#define read type
    if($line =~ /read_length/g){
        $read_length=$elements[1]; chomp $read_length;
        print "expected read length is $read_length\n";
    }#save read length
    if($line =~ /read_list/g){
        $read_list=$elements[1]; chomp $read_list;
        print "read list is $read_list\n";
    }#define read  list
    if($line=~ /retain_intermediate_files/g){
        $save_files=$elements[1]; chomp $save_files;
    }#define save files
    if($line =~ /number_indiv_per_job/g){
        $number_indiv_per_job=$elements[1]; chomp $number_indiv_per_job;
       
        print "task being split into $number_indiv_per_job per job\n";
	$prefix="$read_list.";

	open SPLIT, ">split_jobs_list";

	my $all_read_list=qx(cat $read_list); chomp $all_read_list;
	my @read_list_array=split(/\n/,$all_read_list);

	my $total=qx(wc -l $read_list | perl -p -e 's/ +/\t/g' | cut -f 1); chomp $read_list;
	#!$total=$total-1;
	#print "total number of files is $total\n";

	my $current_indiv_rl=0; my $counter_rl=0; my $file_num=0;
	while($current_indiv_rl < $total){
	    #print "$current_indiv_rl\t$total\n";
	    if($counter_rl == 0){
		my $prefix_file_curr="$prefix"."part"."$file_num";
		print SPLIT "$prefix_file_curr\n";
		push(@jobs,$prefix_file_curr);
		open RL, ">$prefix_file_curr";
		$file_num=$file_num+1;
	    }#open new prefix file                                                                                                           
	    print RL "$read_list_array[$current_indiv_rl]\n";

	    $current_indiv_rl++; $counter_rl++;

	    if($counter_rl == $number_indiv_per_job){
		$counter_rl=0;
	    }#reset counter                                                                                                                  
	}#count jobs                                                                                                                     

#	system("ls $prefix* > split_jobs_list");
    
    }#batch parameters
    if($line =~ /slurm_command_map/g){
        my @job1=split(/\#/,$line);
        $job1_submit="#"."$job1[1]"."\n"."#"."$job1[2]"."\n"."#"."$job1[3]"."\n"."#"."$job1[4]"."\n"."#"."$job1[5]"."\n"."$job1[6]"."\n"."$job1[7]"."\n"; chomp $job1_submit;
        print "cluster command for job1 is: $job1_submit\n";
    }#mapping command
    if($line =~ /slurm_command_variant_call/g){
        my @job2=split(/\#/,$line);
        $job2_submit="#"."$job2[1]"."\n"."#"."$job2[2]"."\n"."#"."$job2[3]"."\n"."#"."$job2[4]"."\n"."#"."$job2[5]"."\n"."$job2[6]"."\n"."#"."$job2[7]"."\n"; chomp $job2_submit;
        print "cluster command for job2 is: $job2_submit\n";
    }#variant calling command 
    if($line =~ /slurm_command_ncASE/g){
        my @job3=split(/\#/,$line);
        $job3_submit="#"."$job3[1]"."\n"."#"."$job3[2]"."\n"."#"."$job3[3]"."\n"."#"."$job3[4]"."\n"."#"."$job3[5]"."\n"."$job2[6]"."\n"."#"."$job2[7]"."\n"; chomp $job3_submit;
        print "cluster command for job3 is: $job3_submit\n";
    }#mapping command 
    if($line =~ /provide_AIMs/){
        $provide_AIMs=$elements[1]; chomp $provide_AIMs;
	if(length($provide_AIMs) eq 0){
	    $aims_status=0;
	} elsif(length($provide_AIMs) > 0){
	    $aims_status=1; print "pre-defined aims file was provided: $provide_AIMs\n";
	}#aims are provided
    }#aims list if provided
   

}#read in the configuration file

####replace mapper and variant caller with full path variables
my $path_bcf="bcftools"; my $path_sam="samtools"; my $path_bedtools="bedtools intersect"; 
my $program=$mapper;
if(length($map_path)>0){
    if(($mapper eq 'bwa') or ($mapper eq 'BWA')){
    $mapper="$map_path"."/"."bwa";
    } elsif(($mapper eq 'star') or ($mapper eq 'STAR')){
    $mapper="$map_path"."/"."STAR";
    }#set correct path for mapping
}#replace
if(length($sam_path)>0){
    $path_bcf="$sam_path"."/"."$bcftools";
    $path_sam="$sam_path"."/"."$samtools";
}#replace
if(length($bedtools_path)>0){
    $path_bedtools="$bedtools_path"."/bin/"."intersectBed";
}#replace

####check if files exist
if(! -f $genome1){
    die "cannot find file $genome1\n";
}#genome2
if(! -f $genome2){
    die "cannot find file $genome2\n";
}#genome2
if(! -f $read_list){
    die "cannot find file $read_list\n";
}#read list
if(($gtf_status eq 1) && (! -f $gtf)){
    die "cannot find file $gtf\n";
}#gtf if defined
if(($aims_status eq 1) && (! -f $provide_AIMs)){
    die "cannot find file $provide_AIMs\n";
}#aims if defined

#####first index reference if needed

if($skip eq 0){

if($program eq 'bwa'){

    system("$mapper index $genome1");
    system("$mapper index $genome2");

} elsif($program eq 'star'){

    print "$mapper  --runMode genomeGenerate --genomeDir par1 --genomeFastaFiles $genome1 --sjdbGTFfile $gtf --outFileNamePrefix par1_trans --sjdbOverhang 99","\n";
    if((! -e 'par1') or (! -e 'par2')){
	system("mkdir par1"); system("mkdir par2");
    }#make genome output directories

    #change so that --runThreadN equals number specified in map slurm command
    system("$mapper  --runMode genomeGenerate --genomeDir par1 --genomeFastaFiles $genome1 --sjdbGTFfile $gtf --outFileNamePrefix par1_trans --sjdbOverhang 99");
    system("$mapper  --runMode genomeGenerate --genomeDir par2 --genomeFastaFiles $genome2 --sjdbGTFfile $gtf --outFileNamePrefix par2_trans  --sjdbOverhang 99");

}#check mapper and index appropriately

}#only index if requested

#####generate AIMs list or identify use-defined AIMs list
my $aims="all_AIMs_"."$genome1"."_"."$genome2";
$aims=~ s/\.\///g;
$aims=~ s/\//_/g;

if(length($provide_AIMs)==0){

    if((! -f $aims) or (! -f "current_aims_file")){
	system("perl identify_AIMs_two_genomes.pl $genome1 $genome2 > $aims");
	open AIMSFILE, ">current_aims_file";
	print AIMSFILE "$aims\n";
    }#if aims file and key do not exist, write them
    else{
print "aims files $aims and current_aims_file exist, not overwriting\n"
    }#warn about the re-use of these files

}#no aims are provided, generate
elsif(-f $provide_AIMs){

    system("cp $provide_AIMs $aims");
    system("perl -pi -e 's/\t/_/g' $aims"); #reformat for downstream compatibility
    open AIMSFILE, ">current_aims_file";
    print AIMSFILE "$aims\n";

}#aims file defined and exists
elsif(! -f $provide_AIMs){

    die "cannot open the user defined AIMs file $provide_AIMs\n";

}#aims file defined but does not exist


##two types of approaches, one for SE and one for PE
my @slurm_ids_map=();
my $slurm_sam_string="";
my $hyb_string=""; my $par_string="";


for my $j (0..scalar(@jobs)-1){
    my $current_job=$jobs[$j];
    #print "$current_job\n";
    my $mapscript="map_batch"."$j".".sh";
    open MAPSCRIPT, ">$mapscript";
    print MAPSCRIPT "$job1_submit\n";

####put appropriate mapping commands for STAR or BWA###
if($program eq 'bwa'){
    print MAPSCRIPT "perl run_map_rnaseq.pl $current_job $genome1 $genome2 $read_type bwa $mapper $gtf\n";
} elsif($program eq 'star'){
    print MAPSCRIPT "perl run_map_rnaseq.pl $current_job $genome1 $genome2 $read_type star $mapper $gtf\n";
}#map using appropriate commands

    my $id_temp=qx(sbatch $mapscript); chomp $id_temp;
    my @idarray=split(/\n/,$id_temp);
    $id=$idarray[0]; chomp $id;
    $id=~ s/\D//g;
    push(@slurm_ids_map,$id);
    print "submitting mapping batch id $id\n";
}#all mapping

for my $m (0..scalar(@jobs)-1){
    my $current_job=$jobs[$m];

    my $samscript="samtools_batch"."$m".".sh";
    open VARSCRIPT, ">$samscript";
    print VARSCRIPT "$job2_submit\n";

####PUT IN APPROPRIATE JOB FOR SAMTOOLS MPILEUP AND OTHER
####SAVE FILENAMES for final ncASE string
#!    print VARSCRIPT "perl run_samtools_to_hmm_v8.pl $current_job $genome1 $genome2 $read_length $save_files $max_align $focal_chrom $rec_M_per_bp\n";
    if($program eq 'bwa'){
    print VARSCRIPT "perl run_samtools_rnaseq.pl $current_job $genome1 $genome2 $read_length bwa $path_sam $path_bcf $path_bedtools $gtf\n";
    } elsif($program eq 'star'){
    print VARSCRIPT "perl run_samtools_rnaseq.pl $current_job $genome1 $genome2 $read_length star $path_sam $path_bcf $path_bedtools $gtf\n";
    }#use appropriate commands  

    my $map_depend=$slurm_ids_map[$m];

    my $id_temp=qx(sbatch --dependency=afterok:$map_depend $samscript); chomp $id_temp;
    my @idarray=split(/\n/,$id_temp);
    $id=$idarray[0]; chomp $id;
    $id=~ s/\D//g;
    if($m==0){
	$slurm_sam_string="$id";
    } else{
        $slurm_sam_string="$slurm_sam_string".","."$id";
    }
    print "submitting variant batch id $id\n";

}#all variants

open READS, $read_list or die "cannot open $read_list\n";
my $ase_string="";
while(my $reads=<READS>){
    chomp $reads;
    my @read_split=split(/\t/,$reads);
    my $file1_tmp=$read_split[0]; chomp $file1_tmp;
    my @readinfo=split(/\//,$file1_tmp);
    my $file1=$readinfo[-1];
    my $file2=$file1;

    $file1="$file1"."_par1_ASE_counts";
    $file2="$file2"."_par2_ASE_counts";

    $ase_string="$ase_string"." "."$file1"." "."$file2";

}#read in read data for summary

open ASESCRIPT, ">ncase_batch.sh";
print ASESCRIPT "$job3_submit\n";
##print ncASE_v3.R command
my $summary_file="$read_list"."_counts_summarized";
print ASESCRIPT "Rscript ncASE_pipeline_cmd.R $ase_string $read_length $allow_zero_counts $bias_threshold $summary_file\n";

system("sbatch --dependency=afterok:$slurm_sam_string ncase_batch.sh");
