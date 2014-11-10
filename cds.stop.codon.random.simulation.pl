#!/usr/bin/perl -w

####################################################################################
##程序功能：将cds序列的2位和3位密码子的终止密码子信息进行统计与随机codon重排序列进行
##          比较
####################################################################################

##################################################
##用法：perl xx.pl -in cds.fasta -random.time [n] -out outfile
##################################################

my $infile=$ARGV[1];
my $rand_time=$ARGV[3];
my $outfile=$ARGV[5];

open IN,"<$infile";
my @data=<IN>;
chomp @data;
close IN;

open OUT,">$outfile";
print OUT "cds name\t2nd stop codon number\t3rd stop codon number\t2nd random stop codon number\t3rd random stop codon number\n";
my $i=0;
while($i<$#data){
        if($data[$i]=~/^>/) {
                my $cds_name=$data[$i];
                $cds_name=~s/^>//;
				print "$cds_name................";
                $data[$i+1]=~tr/agtc/AGTC/;
                my $cds=&DelSE($data[$i+1]);
                my $ori2nd=&Count2ndStopCodon($cds);
                my $ori3rd=&Count3rdStopCodon($cds);
				my $time=1;
				my $ran2nd=0;
				my $ran3rd=0;
				while($time<=$rand_time){
					my $rand_cds=&CodonRandom($cds);
                	$ran2nd+=&Count2ndStopCodon($rand_cds);
                	$ran3rd+=&Count3rdStopCodon($rand_cds);
					$time++;
				}
				$ran2nd=$ran2nd/$rand_time;
				$ran3rd=$ran3rd/$rand_time;
                print OUT "$cds_name\t$ori2nd\t$ori3rd\t$ran2nd\t$ran3rd\n";
				print "\n";
        }
        $i++;
}
close OUT;

##删除起始密码子和终止密码子
sub DelSE {
        my($seq)=@_;
        $seq=~s/^ATG//;
        $seq=~s/TAA$//;
        $seq=~s/TAG$//;
        $seq=~s/TGA$//;
        return $seq;
}

sub CodonRandom{
        my $sequence=$_[0];
        ##codon提取
        my $i=0;
		my $codon_num=(length $sequence)/3;
		my @ori_num=(1..$codon_num);
        my @codons=();
        while($i<=$codon_num){
        push @codons,substr($sequence,$i*3,3);
        $i++;
        }

        my @random=();
        srand( time() ^ ($$ + ($$ << 15)) );
        my $r1=int (rand(@ori_num));
        push @random,$ori_num[$r1];
		$ori_num[$r1]="NA";
		@ori_num=grep(!(/NA/),@ori_num);
        while($#ori_num>=0){
                my $rand=int (rand(@ori_num));
                push @random,$ori_num[$rand];
				$ori_num[$rand]="NA";
				@ori_num=grep(!(/NA/),@ori_num);
        }
        my $randseq="";
        foreach my $r (@random){
                $randseq.=$codons[$r-1];
        }
        return $randseq;
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