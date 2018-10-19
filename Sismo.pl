#!usr/bin/perl 
use strict;
use Data::Dumper;
use List::Util 'sum';
use List::Util 'shuffle';
#use lib '/users/so/lzapata/perl5/lib/perl5';
use Statistics::R;
use Getopt::Std;

my %opt=();
getopts("i:f:n:",\%opt);

my $usage="$0
-i FASTA transcript files
-f mutational profile
-n number of simulated mutations
";

my $fasta = $opt{i} or die $usage;
my $profile= $opt{f} or die $usage;
my $nsim = $opt{n} or die $usage;

my %freq=getFreq($profile);

open (FASTA_FILE, $fasta) or die;  # we get the cds sequences. They MUST be inline (>ENST\nATGGATGCTAGCTAC blabla)
my %seq;
#open (OUTFILE, ">$ARGV[0].count") or die;

####print STDERR "Step 1, read fasta file\n";
while (<FASTA_FILE>){
	my $read = $_; 
	chomp $read;
	if ($read =~ m/>/){		
		my ($id) = ($read =~ m/>(\S+)/);
		$read = <FASTA_FILE>;
		chomp $read;
		my $sec = $read;
		$seq{$id} = $sec;
	}
}

##Create hash for storing transcript, triplet, and counts, and lengths per transcript

my %tripletcount;
my %genelength;
my %posingene;
my %value_count;

###print STDERR "Step 2, split and count triplets in fasta\n";
##go to each transcripr
foreach my $trans(keys %seq){
	#print $trans."\n";
	#$seq{$trans}."\t".length($seq{$trans})."\n";
	###Store genelength per transcript
	$genelength{$trans}=length($seq{$trans});
	
	#Define array of triplets per transcript
	my @trip;
	my $totalcount=0;
	#split all sequence in base pairs
	my @sequence_trans = split ('', $seq{$trans});
	
	#store triplets per gene with index from 0 to the total length of the gene and store the position of triplet in each gene
	for my $i (0 .. $#sequence_trans-2){	
		$trip[$i]=substr($seq{$trans},$i,3,);
		push (@{$posingene{$trans}{$trip[$i]}{pos}},$i);
		#print $i."\n";
	}
	
	##count number of triplet ocurrences per gene and store in hash count
	my %count;
	foreach my $element( @trip ) {
		++$count{$element};
	}
	
	##for each triplet in count store the value in to the hash associated to the transcript and the
	#triplet and calculate the frequency of the triplet based on the length of the gene
	foreach my $element( sort keys %count ) {
		$tripletcount{$trans}{$element}{count}=$count{$element};
		$tripletcount{$trans}{$element}{freq}=$count{$element}/(length($seq{$trans})-2);
	}
}

#print Dumper(%posingene);
#print Dumper(%tripletcount);

##Count total ocurrences of triplet in the genome and get the probabilti per gene to fall into a gene
my %counting;
foreach my $transcript(sort keys %tripletcount){
	#print Dumper($tripletcount{$transcript});
	foreach my $triplet(keys %{ $tripletcount{$transcript} }){
		$counting{$triplet} += $tripletcount{$transcript}{$triplet}{count}."\n";			
	}
	
}

my %prob;
foreach my $transcript(sort keys %tripletcount){
	#print Dumper($tripletcount{$transcript});
	foreach my $triplet(keys %{ $tripletcount{$transcript} }){
		$tripletcount{$transcript}{$triplet}{prob}= $tripletcount{$transcript}{$triplet}{count}/$counting{$triplet};
		push (@{$prob{$triplet}{name}},$transcript);
		push (@{$prob{$triplet}{prob}},$tripletcount{$transcript}{$triplet}{prob});
	}
}


my @triplets;
my @trip_probs;

###print STDERR "Step 3, get triplet frquencies per transcripts\n";

foreach my $el( keys %freq){
	foreach my $el2(keys %{ $freq{$el} }){
		push @triplets, $el."_".$el2;
		push @trip_probs, $freq{$el}{$el2}{freq};		
		}
}

#print Dumper(@trip_probs);

my $R = Statistics::R->new();
$R->set( 'x', \@triplets );
$R->set('p', \@trip_probs);
$R->set('n', $nsim);
$R->run( q'sample_for_perl = sample(x, n, replace = T, prob= p)' );
my $Rsample = $R->get('sample_for_perl');

$R->stop();

my @changes=split(" ","@$Rsample");

#open(VEP, ">results.vep");

##main simulation part once a triplet is obtained

##print STDERR "Step 4, Get one transcript based on R and get a random position\n";
foreach my $runs(@changes){
	my ($trip,$sub)=split("_",$runs);
	
	my @names = split(" ","@{$prob{$trip}{name}}");
	my @probs2 = split(" ","@{$prob{$trip}{prob}}" );
	
	my $R = Statistics::R->new();
	
	#print "@{$prob{$trip}{name}}"."\n";
	
	$R->set( 'x1', \@names );
	$R->set('p1', \@probs2 );
	$R->run( q'sample_for_perl2 = sample(x1, 1, replace = T, prob= p1)' );
	my $Rsample2 = $R->get('sample_for_perl2');

	$R->stop();
	my $position = getPosInGene($Rsample2,$trip);
	#print join ("\t", $Rsample2, ($position+1), $trip , $sub)."\n";
	my (@ch)=split("/",$sub);
	print $Rsample2.":c.".($position+1).$ch[0].">".$ch[1]."\n";
	#print Dumper($position);
	
}

#close VEP;

sub getPosInGene{
	my $transcript = shift;
	my $context = shift;	
	my @position = split(" ", "@{$posingene{$transcript}{$context}{pos}}");
	
	#print Dumper($position[10]);
	
	my $superpos= $position[rand @position] + 1;
	
	#print Dumper(int(rand scalar(@position)));
	return $superpos;
	
	
}



sub getFreq{
	open FH, shift;
	my %info;
	my $tripcount=0;
	
	while (my $line = <FH>) {
		chomp $line;
		my ($triplet,$change,$number)=split("\t",$line);
		$info{$triplet}{$change}{count}=$number;
		$tripcount+=$number;
	}
	
	foreach my $el( keys %info){
		foreach my $el2(keys %{ $info{$el} }){
			$info{$el}{$el2}{freq}=$info{$el}{$el2}{count}/$tripcount;			
		}
	}
	#print Dumper(%info);
	return %info;
}











