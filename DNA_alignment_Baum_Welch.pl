############################################################
#Description��Given two DNA sequences, find the most       #
#likely set of state transition and output probabilities.  #
# The problem can be solved by the Baum-Welch algorithm.   #
#Requirement��Find the best final alignment result.        #
#Author��xiaofanpen(С����)                                #
#Date��2017��6��21��23:24:06                               #
############################################################
#!perl  -w
$seq1="TGTTCC";
$seq2="AGTTTA";

@seq1=split(//,$seq1);
@seq2=split(//,$seq2);

#���巣�ֹ���
$match=1;
$mismatch=-1;
$gap=-1;

#��־����һ�к͵�һ�и�ֵ
#$score[0][0]={0};
$len1=@seq1;
$len2=@seq2;

#���Ÿ�ֵ
for $i(0..$len1)
{
	$score[0][$i]=-$i;

}
#���Ÿ�ֵ
for $j(0..$len2)
{
    $score[$j][0]=-$j;
}
#���÷־��������λ��
for $i(1..$len1)
{
   $letter2=$seq2[$i-1];
   for $j(1..$len2)
   {
        $letter1=$seq1[$j-1];
		if($letter1 eq $letter2) #·��б���ߣ�Ҫôƥ�䣬Ҫô���̣�
       	{
            $diagonal=$score[$i-1][$j-1]+$match;
       	}
		else
		{
	       $diagonal=$score[$i-1][$j-1]+$mismatch;
		}
            $up=$score[$i-1][$j]+$gap;#·��������
	        $left=$score[$i][$j-1]+$gap;#·�ߺ�����

		#�Ƚ�������������ŵ÷�
		if($diagonal>=$up)
		{
	       	if($diagonal>=$left)
   			{
                $score[$i][$j]=$diagonal;
 			    $trace[$i][$j]=1;
	       	}
			else
			{
			    $score[$i][$j]=$left;
                $trace[$i][$j]=2;
			}
		}
        else
		{
			if($up>=$left)
	       	{
	       		$score[$i][$j]=$up;
	      		$trace[$i][$j]=3;
	      	}
	      	else
	      	{
	      		$score[$i][$j]=$left;
	      		$trace[$i][$j]=2;
	     	}
		}
   }
}
####@score����ÿ�ж��ٸ�������ȶ�һ�¡�
for $j(0..$len2){
	for $i(0..$len1){
		print "$score[$j][$i]  ";
	}
	print "\n";
}
print "The highest score is $score[$len1][$len2]\n";

#�ڵ÷־����У��ҵ�����·���ı�ǣ��Ӷ��ҵ�����·��
$j = @seq1;
$i = @seq2;
$k = 0;
#@align1=();
#@align2=();
#@hiddenState=();

while(1)
{
	if($i>0&&$j>0)
    {
		if($trace[$i][$j]==1)
		{
			$align1[$k]=$seq1[$j-1];
            $align2[$k]=$seq2[$i-1];
        	$hiddenState[$k]="M";
            $k++;
            $j--;
            $i--;
            next;
        }
        elsif($trace[$i][$j]==2)
		{
            $align1[$k]=$seq1[$j-1];
            $align2[$k]="-";
       	    $hiddenState[$k]="X";
            $k++;
            $j--;
            next;
        }
        elsif($trace[$i][$j]==3)
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
@align1=reverse(@align1);
@align2=reverse(@align2);
@hiddenState=reverse(@hiddenState);
#for(int $i=0;$i<@align1;$i++){
#	if($align1[$i]!=)
#}
print"The alignment result��\n";
print"Sequence 1��\n@align1\n";
print"Sequence 2��\n@align2\n";
print "Hidden state sequence��\n@hiddenState\n";
