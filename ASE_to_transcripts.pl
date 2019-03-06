#perl! -w

if(@ARGV<2){
    print "perl ASE_to_transcripts.pl intersect_file tag\n"; exit;
}#print usage

my $file=shift(@ARGV); chomp $file;
open IN, $file or die "cannot open intersected vcf\n";

my $tag=shift(@ARGV); chomp $tag;

while(my $line=<IN>){

    chomp $line;
    my @elements=split(/\t/,$line);
    my $info_field=$elements[18]; #print "$info_field\n";

    my $pos=$elements[1];

    my @info_elements=split(/;/,$info_field);
    my $found_tag=$info_elements[0]; chomp $found_tag; #print "$found_tag\n";
    $found_tag=~ s/"//g;

    my $class=$elements[12];

    if(($found_tag =~ /gene_id/g) && ($tag eq $class)){

	my @gene_name_info=split(/ /,$found_tag);
	my $gene_name=$gene_name_info[-1]; chomp $gene_name;

	print "$gene_name\t$pos\t$elements[2]\t$elements[3]\t$elements[4]\t$elements[5]\t$elements[6]\t$elements[7]\t$elements[8]\t$elements[9]\t$elements[10]\n";

    }#check for right annotation flag

}#all lines
