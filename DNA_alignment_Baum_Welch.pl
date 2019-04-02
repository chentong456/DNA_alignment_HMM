############################################################
#Description：Given two DNA sequences, find the most       #
#likely set of state transition and output probabilities.  #
# The problem can be solved by the Baum-Welch algorithm.   #
#Requirement：Find the best final alignment result.        #
#Author：xiaofanpen(小饭盆)                                #
#Date：2017年6月21日23:24:06                               #
############################################################
#!perl  -w
$seq1="TGTTCC";
$seq2="AGTTTA";

@seq1=split(//,$seq1);
@seq2=split(//,$seq2);

#定义罚分规则
$match=1;
$mismatch=-1;
$gap=-1;

#打分矩阵第一列和第一行赋值
#$score[0][0]={0};
$len1=@seq1;
$len2=@seq2;

#横着赋值
for $i(0..$len1)
{
	$score[0][$i]=-$i;

}
#竖着赋值
for $j(0..$len2)
{
    $score[$j][0]=-$j;
}
#填充得分矩阵的其他位置
for $i(1..$len1)
{
   $letter2=$seq2[$i-1];
   for $j(1..$len2)
   {
        $letter1=$seq1[$j-1];
		if($letter1 eq $letter2) #路线斜着走（要么匹配，要么容忍）
       	{
            $diagonal=$score[$i-1][$j-1]+$match;
       	}
		else
		{
	       $diagonal=$score[$i-1][$j-1]+$mismatch;
		}
            $up=$score[$i-1][$j]+$gap;#路线竖着走
	        $left=$score[$i][$j-1]+$gap;#路线横着走

		#比较四种情况的最优得分
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
####@score到底每行多少个，输出比对一下。
for $j(0..$len2){
	for $i(0..$len1){
		print "$score[$j][$i]  ";
	}
	print "\n";
}
print "The highest score is $score[$len1][$len2]\n";

#在得分矩阵中，找到最优路径的标记，从而找到最优路径
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
print"The alignment result：\n";
print"Sequence 1：\n@align1\n";
print"Sequence 2：\n@align2\n";
print "Hidden state sequence：\n@hiddenState\n";
