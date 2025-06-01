#!/usr/bin/perl -w
use strict;
use List::Util qw(sum);
use Getopt::Std;
my %opt;
&getopts("hv:w:p:d:h",\%opt);
if(!exists $opt{w} || !exists $opt{p} || !exists $opt{d})
{
	print "perl xxx.pl -p path_ind_vcf -d average_depth -w window_size\n";
	print "\t-p /strage_151/zhangwp_data/mapping_pipelines/Pipeline_Chu-NE/snp_vcf\n";
	print "\t-w 50000\n";
	exit;
}
my %depth;
open (IN,$opt{d})or die;
while(<IN>)
{
	/(\S+)\s+(\S+)/;
	$depth{$1}=$2;
}
close(IN);
printf "%-7s%-9s"."%-15s"x(keys %depth)."\n","CHROM","RANGE",(sort keys %depth);

my %Het;
my %Dep;
my @mapping = ('A'..'Z', 'a'..'z', '0'..'9', '!', '@', '#', '$', '%', '^', '&', '*', '(', ')', '-', '_', '+', '=', '[', ']', '{', '}', '|', '\\', ':', ';', '<', '>', ',', '.', '?', '/', '~');
my @files=split/\n/,`ls $opt{p}/*gz`;
foreach my $f (@files)
{
	if($f=~/([^\/]+).g.vcf.gz$/ && exists $depth{$1})
	{
		my $ind=$1;
		my ($min,$max)=($depth{$ind}/3,$depth{$ind}*2);
		if($depth{$ind}<=18){$min=6;}
		#my ($min,$max)=(3,30);
		open (IN, "gzip -dc $f |") or die;
		while (<IN>){
			if(/^Chr/){
				my @t=split/\t/,$_;
				my ($chr,$pos)=($t[0],$t[1]);
				my ($len,$label)=(1,"N");
				if(exists $Het{$chr}{$ind})
				{
					my $d=length($Het{$chr}{$ind})-$pos+1;
					$Het{$chr}{$ind}=~s/.{$d}$//;
					$Dep{$chr}{$ind}=~s/.{$d}$//;
				}
				$t[9]=~/^.*?:(\S+?):/;
				my $text=$1;
				my $dp=&sum(split/,/,$text);
				if($dp>90){$dp=90;}
				#depth 1/3 ~ 2
				if($dp<=$min || $dp>=$max){
					$label="X";
					if(/END=(\d+)/)	{$len=$1-$pos+1;}
					else{$len=length($t[3]);}
				# <1/2 or >2
				}else{
					##nonbiallele && indel
					if(length($t[3])!=1 || ($t[4]=~/,/ && length($t[4])!=11)){
						$label="X";
						$len=length($t[3]);
					#Hom
					}elsif($t[9]=~/^(0\/0|1\/1)/){
						if(/END=(\d+)/){$len=$1-$pos+1;}
					#Het
					}else{
						$label="#";
						$t[9]=~/^0\/1:(\d+),/;
						my ($limit_1,$limit_2) = (0.1 * $dp,0.9 * $dp);
                   		if($dp >= 20){($limit_1,$limit_2) = (0.2 * $dp,0.8 * $dp);}   
                     	if ( $1 >= $limit_2 ||  $1 <= $limit_1 ){$label="N";}
					}
				}
				$Het{$chr}{$ind}.=$label x $len ;
				$Dep{$chr}{$ind}.=$mapping[$dp] x $len;
				#print STDERR "$chr,$pos,$len,"."$label\n";
			}
		}
		close(IN);
	}
}

foreach my $chr(sort keys %Het)
{
	my %a=%{$Het{$chr}};
	my %b=%{$Dep{$chr}};
	my $len=length((values%a)[0]);
	#grep{print STDERR "$chr-$_:".length($a{$_})."\n";}keys %a;
	for (my $i=0; $i+$opt{w}<$len; $i=$i+$opt{w})
	{
		my @pro_het=();
		my @ave_dep=();
		foreach my $ind(sort keys %depth)
		{
			my $subseq=substr($a{$ind},$i,$opt{w});
			my $num_het=$subseq=~tr/#//;
			my $num_missing=$subseq=~tr/X//;
			if($opt{w}/2 < $num_missing){ push(@pro_het,"NA");}
			#if($opt{w}*0.8 < $num_missing){ push(@pro_het,"NA");}
			else{push(@pro_het,sprintf "%.6f",$num_het/($opt{w}-$num_missing));}
			my $subseq_dep=substr($b{$ind},$i,$opt{w});
			my $sum_dep=0;
			grep{$sum_dep+=index(join("",@mapping),$_);}split("",$subseq_dep);
			push(@ave_dep,sprintf "%.2f",$sum_dep/$opt{w});
		}
		my @comb=();
		grep{$comb[$_]=$pro_het[$_].":".$ave_dep[$_];} 0 .. $#pro_het;
		printf "%-7s%-9s"."%-15s"x (keys %depth) ."\n",$chr,$i+1,@comb;	
	}
}
