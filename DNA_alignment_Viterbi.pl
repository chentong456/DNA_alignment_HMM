###############################################
#尝试自己写，写不出来再说                     #
#Finished！！！噢耶~~~~~~~~~~                 #                 
#DNA sequence alignment Viterbi HMM           #   
#Author:xiaofanpen                            #
#Research unit：Beijing Institute of Genomics #
#Chinese Academy of Sciences                  #
#Major：Bioengineering                        #
#Date:2017年6月23日14:02:20                   #
###############################################
#!perl -w
$seq1="ACGCTA";
$seq2="ACGCTA";
@seq1=split(//,$seq1);
@seq2=split(//,$seq2);
#print"@seq1 and @seq2\n";
$delta=0.2; $eta=0.5;

#### Transition probabilities ####
#—— M:match/mismatch
#—— X：gap in seq1
#—— Y：gap in seq2
#—— B：begin
$BX=$delta; $BY=$delta; $BM=1-2*$delta;
$MX=$delta; $MY=$delta; $MM=1-2*$delta;
$XX=$eta; $XY=0; $XM=1-$eta;
$YY=$eta; $YX=0; $YM=1-$eta;

#### Emission probabilities ####
#—— pgap：probability of gap by X or Y. (equal to 1)
#—— pm:probabilitity of match by M. (1-pm)
#—— pmm：probability of mismatch by M. (pm=0.6)
#$pgap=1; $pm=0.6; $pmm=1-$pm; 我觉着pgap=1,貌似不合理
$pgap=0.4; $pm=0.6; $pmm=1-$pm;
#### Initialize B ####
$B[0][0]=1; $B[0][1]=0; $B[1][0]=0;
for(my $i=1;$i<=length($seq2);$i++){
	for(my $j=1;$j<=length($seq1);$j++){
		$B[$i][$j]=0;
	}
}

#### Initialize M ####
$M[0][0]=0;
for(my $i=1;$i<=length($seq2);$i++){
	$M[$i][0]=0;
}
for(my $i=1;$i<=length($seq1);$i++){
	$M[0][$i]=0;
}

#### Initialize X ####
$X[0][0]=0;
$X[1][0]=$pgap*$BX*$B[0][0];
for(my $i=2;$i<length($seq2);$i++){
	$X[$i][0]=$pgap*$XX*$X[$i-1][0];
}
for(my $i=1;$i<length($seq1);$i++){
	$X[0][$i]=0;
}

#### Initialize Y ####
$Y[0][0]=0;
$Y[0][1]=$pgap*$BY*$B[0][0];
for(my $i=1;$i<length($seq2);$i++){
	$Y[$i][0]=0;
}
for(my $i=2;$i<length($seq1);$i++){
	$Y[0][$i]=$pgap*$YY*$Y[0][$i-1];
}

#### Initialize T ####
#感觉初始化有问题啊，改了结果输出不变，简直不为所动？？？
#OK,have solved!

for(my $i=1;$i<=length($seq2);$i++){
	$T[$i][0]='U';
}
for(my $i=1;$i<=length($seq1);$i++){
	$T[0][$i]='L';
}

#### Alignment HMM Viterbi ####
#for($i=1;$i<=length($seq2);$i++){
#	for($j=1;$j<=length($seq1);$j++){
		#calculate the value of M[$i][$j], Y[$i][$j], X[$i][$j]
		
#	}
#}
#### Recurrence ####
for(my $i=1;$i<=@seq2;$i++){
	for(my $j=1;$j<=@seq1;$j++){
		my $max_1;
		if($seq2[$i] eq $seq1[$j]){
			#M[i][j]=$pm*
			#for(my $k=0;$k<3;$k++)
			my @aa=($MM*$M[$i-1][$j-1],$XM*$X[$i-1][$j-1],$YM*$Y[$i-1][$j-1],$BM*$B[$i-1][$j-1]);
			$max_1=&Max(@aa);	
			$M[$i][$j]=$pm*$max_1;		
		}
		else{
			$M[$i][$j]=$pmm*$max_1;
		}
		my @bb=($MY*$M[$i][$j-1],$YY*$Y[$i][$j-1],$BY*$B[$i][$j-1]);
		my $max_2=&Max(@bb);
		$Y[$i][$j]=$pgap*$max_2;
		my @cc=($MX*$M[$i-1][$j],$XX*$X[$i-1][$j],$BX*$B[$i-1][$j]);
		my $max_3=&Max(@cc);
		$X[$i][$j]=$pgap*$max_3;
		if($M[$i][$j]>=$X[$i][$j] && $M[$i][$j]>=$Y[$i][$j]){
			$T[$i][$j]='D';
		}
		elsif($X[$i][$j]>=$M[$i][$j] && $X[$i][$j]>=$Y[$i][$j]){
			$T[$i][$j]='U';
		}
		elsif($Y[$i][$j]>=$M[$i][$j] && $Y[$i][$j]>=$X[$i][$j]){
			$T[$i][$j]='L';
		}
	}
}
print "---------M---------\n";
for(my $i=0;$i<length($seq2);$i++){
	for(my $j=0;$j<length($seq1);$j++){
		print "$M[$i][$j]  ";
	}
	print "\n";
}
print "\n\n\n---------X---------\n";

for(my $i=0;$i<length($seq2);$i++){
	for(my $j=0;$j<length($seq1);$j++){
		print "$X[$i][$j]  ";
	}
	print "\n";
}
print "\n\n\n---------Y---------\n";
for(my $i=0;$i<length($seq2);$i++){
	for(my $j=0;$j<length($seq1);$j++){
		print "$Y[$i][$j]  ";
	}
	print "\n";
}
print "\n\n\n---------T---------\n";
for(my $i=0;$i<=length($seq2);$i++){
	for(my $j=0;$j<=length($seq1);$j++){
		print "$T[$i][$j]  ";
	}
	print "\n";
	
}


#### Backtracking ####
$j=@seq1;
$i=@seq2;
$k=0;
while(1)
{
	if($i>0 && $j>0)
    {
		if($T[$i][$j] eq "D")
		{
			$align1[$k]=$seq1[$j-1];
            $align2[$k]=$seq2[$i-1];
        	$hiddenState[$k]="M";
            $k++;
            $j--;
            $i--;
            next;
        }
        elsif($T[$i][$j] eq "L")
		{
            $align1[$k]=$seq1[$j-1];
            $align2[$k]="-";
       	    $hiddenState[$k]="X";
            $k++;
            $j--;
            next;
        }
        elsif($T[$i][$j] eq "U")
		{
            $align1[$k]="-";
            $align2[$k]=$seq2[$i-1];
		    $hiddenState[$k]="Y";
            $k++;
            $i--;
            next;
        }
    }
    elsif($i==0&&$j!=0)
    {
		$align2[$k]=$seq1[$j-1];
		$align1[$k]= "-";
		$hiddenState[$k]="Y";
		$k++;
		$j--;
    }
    elsif($j==0&&$i!=0)
    {
		$align1[$k]=$seq2[$i-1];
		$align2[$k]="-";
		$hiddenState[$k]="X";
		$k++;
		$i--;
    }
    elsif($j==0&&$i==0){last;}
}

#print "@align1\n\n\n\n";
#print "@align2\n\n\n\n";

@align1=reverse(@align1);
@align2=reverse(@align2);
@hiddenState=reverse(@hiddenState);
print"The alignment result：\n";
print"Sequence 1：\n@align1\n";
print"Sequence 2：\n@align2\n";
print "Hidden state sequence：\n@hiddenState\n";

sub Max{
	my $maxValue=shift @_;
	for(@_){
		if($_ > $maxValue){
			$maxValue=$_;
		}
	}
	return $maxValue;
}










 