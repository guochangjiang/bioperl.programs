#!/usr/bin/perl -w

#############################################################################
##程序功能：通过改变氨基酸编码的同义密码子来查看错位终止密码子的分布
#############################################################################

#################################################################
##用法：perl xx.pl -in in.fasta -rand.time [n] -out outfile
#################################################################

my $infile=$ARGV[1];
my $rand_time=$ARGV[3];
my $outfile=$ARGV[5];

open IN,"<$infile";
my @data=<IN>;
chomp @data;
close IN;

open OUT,">>$outfile";
print OUT "cds name\trandom method\t2nd stop codon number\t3rd stop codon number\t2nd random stop codon number\t3rd random stop codon number\n";
close OUT;

my $i=0;
while($i<$#data){
        if($data[$i]=~/^>/) {
                my $cds_name=$data[$i];
                $cds_name=~s/^>//;
				my $prin=$i/2+1;
				print "$cds_name................($prin)";
                $data[$i+1]=~tr/agtc/AGTC/;
                my $cds=&DelSE($data[$i+1]);
                my $ori2nd=&Count2ndStopCodon($cds);
                my $ori3rd=&Count3rdStopCodon($cds);
				my $aaseq= &Translate($cds);
				#print "\n$aaseq\n";
				my $time=1;
				my $ran2nd=0;
				my $ran3rd=0;
				while($time<=$rand_time){
					my $new_cds=&SynCodonSeq($aaseq);
					#print "\n$new_cds\n";
                	$ran2nd+=&Count2ndStopCodon($new_cds);
                	$ran3rd+=&Count3rdStopCodon($new_cds);
					$time++;
				}
				$ran2nd=$ran2nd/$rand_time;
				$ran3rd=$ran3rd/$rand_time;
				open OUT,">>$outfile";
                print OUT "$cds_name\tsynonym codon\t$ori2nd\t$ori3rd\t$ran2nd\t$ran3rd\n";
				close OUT;
				print "\n";
        }
        $i++;
}

##删除起始密码子和终止密码子
sub DelSE {
        my($seq)=@_;
        $seq=~s/^ATG//;
        $seq=~s/TAA$//;
        $seq=~s/TAG$//;
        $seq=~s/TGA$//;
        return $seq;
}



sub Count2ndStopCodon {
        my $seq=$_[0];
        my $num=(length $seq)/3-1;
        my $count=0;
        my $i=1;
        while($i<=$num){
                my $codon=substr($seq,$i,3);
                if($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA") {$count++;}
                $i+=3;
        }
        return $count;
}

sub Count3rdStopCodon {
        my $seq=$_[0];
        my $num=(length $seq)/3-1;
        my $count=0;
        my $i=1;
        while($i<=$num){
                my $codon=substr($seq,$i+1,3);
                if($codon eq "TAA" || $codon eq "TAG" || $codon eq "TGA") {$count++;}
                $i+=3;
        }
        return $count;
}

sub Translate{
               my $cds=$_[0];
               my $tri=0;
               my $pro='';
               my $length=length $cds;
               while($tri<=$length-3){my $codon=substr($cds,$tri,3);
                                      my $aa=&Codon2aa($codon);
                                      $pro.=$aa;
                                      $tri+=3;
                                     }
               return $pro;
}

sub Codon2aa{
              my $codon=$_[0];
              if($codon=~/TTT|TTC/)  {return 'F';}
              if($codon=~/TTA|TTG|CTT|CTC|CTA|CTG/)  {return 'L';}
              if($codon=~/ATT|ATC|ATA/)  {return 'I';}
              if($codon=~/ATG/)  {return 'M';}
              if($codon=~/GTT|GTC|GTA|GTG/)  {return 'V';}
              if($codon=~/TCT|TCC|TCA|TCG/)  {return 'S';}
              if($codon=~/CCT|CCC|CCA|CCG/)  {return 'P';}
              if($codon=~/ACT|ACC|ACA|ACG/)  {return 'T';}
              if($codon=~/GCT|GCC|GCA|GCG/)  {return 'A';}
              if($codon=~/TAT|TAC/)  {return 'Y';}
              if($codon=~/TAA|TAG|TGA/)  {return '*';}
              if($codon=~/CAT|CAC/)  {return 'H';}
              if($codon=~/CAA|CAG/)  {return 'Q';}
              if($codon=~/AAT|AAC/)  {return 'N';}
              if($codon=~/AAA|AAG/)  {return 'K';}
              if($codon=~/GAT|GAC/)  {return 'D';}
              if($codon=~/GAA|GAG/)  {return 'E';}
              if($codon=~/TGT|TGC/)  {return 'C';}
              if($codon=~/TGG/)  {return 'W';}
              if($codon=~/CGT|CGC|CGA|CGG|AGA|AGG/)  {return 'R';}
              if($codon=~/AGT|AGC/)  {return 'S';}
              if($codon=~/GGT|GGC|GGA|GGG/)  {return 'G';}
}

sub SynCodonSeq{
	my $protein=$_[0];
	my $dnaseq="";
	my $i=0;
	while($i<length $protein){
		$dnaseq.=&SynCodon(substr($protein,$i,1));
		$i++;
	}
	return $dnaseq;
}

sub SynCodon {
	my $aa=$_[0];
	#print "$aa\t";
	my $codon="";
	if($aa eq "F") { my @codon=("TTT","TTC"); $codon=&RandGet(@codon);}
	elsif($aa eq "L") { my @codon=("TTA","TTG","CTT","CTC","CTA","CTG"); $codon=&RandGet(@codon);}
	elsif($aa eq "I") { my @codon=("ATT","ATC","ATA");$codon=&RandGet(@codon);}
	elsif($aa eq "M") { $codon="ATG";}
	elsif($aa eq "V") { my @codon=("GTT","GTC","GTA","GTG"); $codon=&RandGet(@codon);}
	elsif($aa eq "S") { my @codon=("TCT","TCC","TCA","TCG","AGT","AGC"); $codon=&RandGet(@codon);}
	elsif($aa eq "P") { my @codon=("CCT","CCC","CCA","CCG"); $codon=&RandGet(@codon);}
	elsif($aa eq "T") { my @codon=("ACT","ACC","ACA","ACG"); $codon=&RandGet(@codon);}
	elsif($aa eq "A") { my @codon=("GCT","GCC","GCA","GCG"); $codon=&RandGet(@codon);}
	elsif($aa eq "Y") { my @codon=("TAT","TAC"); $codon=&RandGet(@codon);}
	elsif($aa eq "H") { my @codon=("CAT","CAC"); $codon=&RandGet(@codon);}
	elsif($aa eq "Q") { my @codon=("CAA","CAG"); $codon=&RandGet(@codon);}
	elsif($aa eq "N") { my @codon=("AAT","AAC"); $codon=&RandGet(@codon);}
	elsif($aa eq "K") { my @codon=("AAA","AAG"); $codon=&RandGet(@codon);}
	elsif($aa eq "D") { my @codon=("GAT","GAC"); $codon=&RandGet(@codon);}
	elsif($aa eq "E") { my @codon=("GAA","GAG"); $codon=&RandGet(@codon);}
	elsif($aa eq "C") { my @codon=("TGT","TGC"); $codon=&RandGet(@codon);}
	elsif($aa eq "W") { $codon="TGG";}
	elsif($aa eq "R") { my @codon=("CGT","CGC","CGA","CGG","AGA","AGG"); $codon=&RandGet(@codon);}
	elsif($aa eq "G") { my @codon=("GGT","GGC","GGA","GGG"); $codon=&RandGet(@codon);}
	#print "$codon\t";
	return $codon;
}

sub RandGet {
	my @syncodon=@_;
	srand( time() ^ ($$ + ($$ << 15)) );
	my $i=int (rand (@syncodon));
	#print "$syncodon[$i]\t";
	return $syncodon[$i];
}
