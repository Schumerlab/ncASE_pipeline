#perl! -w

if(@ARGV<1){
    print "perl cleanup_run.pl runs_read_list\n"; exit;
}#print usage

my $infile=shift(@ARGV); chomp $infile;
open IN, $infile or die "cannot open read list\n";

my $mod="$infile".".part"."*";
system("rm $mod");

while(my $line=<IN>){

    chomp $line;
    my @elements=split(/\t/,$line);
    my $line_mod="$elements[0]"."_"."*";
    system("rm -r $line_mod");

}
