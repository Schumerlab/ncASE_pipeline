#perl! -w
# perl $path/samtools_vcf_to_ASE_counts_v3.pl $overlap1
# input files are: .par1/.par2

my $infile=shift(@ARGV); chomp $infile;
open IN, $infile or die "cannot open infile\n";

open OUT, ">"."$infile"."_ASE_counts";

while(my $line=<IN>){
  if($line !~ '#'){
    chomp $line;
#    print "$line\n";
    my @fields=split(/\t/,$line);
    my $info=$fields[11]; chomp $info;
#    print "$info\n";
    my @infoarray=split(/;/,$info);
    my @dp4 = grep /DP4/, @infoarray; chomp $dp4;
    my $alldepths=$dp4[0]; chomp $alldepths;
    $alldepths=~ s/DP4=//g;
    my @indiv_depths=split(/,/,$alldepths);
    # in DP4 field, index 0 and 1 are ref allele counts
    #               index 2 and 3 are alt allele counts
    my $a1=$indiv_depths[0]+$indiv_depths[1];
    my $a2=$indiv_depths[2]+$indiv_depths[3];
    
    # ref base in pos 2, alt base in pos 3
    my $base1=$fields[2]; chomp $base1;
    my $base2=$fields[3]; chomp $base2;
    my $ref=$fields[7]; chomp $ref;
    my $alt=$fields[8]; chomp $alt;
#    print "$base1\t$base2\t$alt\t$fields[0]\t$a1\t$a2\n";
   
    # check that the bases are as expected from AIMs file,
    # organize according to ref/alt
    if((($alt eq '.') or ($alt eq $base2)) and ($ref eq $base1)){
	    print OUT $fields[0],"\t",$fields[1],"\t",$base1,"\t",$base2,"\t$a1\t$a2\n";
    }
    elsif((($alt eq '.') or ($alt eq $base1)) and ($ref eq $base2)){
      # counts switched to reflect flipped ref/alt alleles
	    print OUT $fields[0],"\t",$fields[1],"\t",$base1,"\t",$base2,"\t$a2\t$a1\n";
	  }

  }#not header line
}#for all lines
