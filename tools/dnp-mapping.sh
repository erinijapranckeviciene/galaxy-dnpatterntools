#!/bin/sh
if test "$#" -ne 6; then

echo ""
echo " CALL  "
echo "   sh dnp-mapping.sh input.fasta input.pattern input.trimstart input.length output.file1 output.file2"
echo ""
echo " INPUT" 
echo "   input.fasta     - input fasta file "
echo "   input.pattern   - 'one or more columns with dinucleotide frequency pattern'"
echo "   input.trimstart - 'number of positions to trim from the start of the sequence'"
echo "   input.length    - 'sequence length to retain past trimming from start'"
echo ""
echo " OUTPUT"
echo "   output.file1   - tabular file with correlations  "
echo "   output.file2   - file to store the max correlation position  "

echo ""
echo " DESCRIPTION"
echo "   Each sequence in the fasta file is reduced by trimming  "
echo "   and retaining a given number of positions, but no less than 147." 
echo "   Correlation of the nucleosome's sequence with the patterns" 
echo "   is computed within the sliding window. Correlation coefficients "
echo "   of the patterns with the sequence starting at a position 73 - dyad " 
echo "   are computed and saved in output.file1. The maximum correlation position"
echo "   is saved in output.file2."
echo ""
echo " REQUIREMENT"
echo "   dnp-mapping installed"
echo "   conda install -c bioconda dnp-mapping"
echo ""
  exit 1
fi

faseqfile=$1
patternfile=$2
seqstart=$3
seqlength=$4

outfile1=$5
outfile2=$6

call=dnp-mapping


awk_program=$( cat << 'EOF' 
###################################################################
# position of maximum
# parameters: window=W (minimal distance between two peaks)
#             buffer=N (size of buffer)
###################################################################
function max_pos_funct(min_pos, max_pos)
{
  sum=0;
    start_position=min_pos;
    for(i=min_pos+window;i<=max_pos&&sum<1000;)
    {
      sum++;
      max=arr[i];
      pos=i;
      for(j=i-window;j<=i+window&&j<=max_pos;j++)
      {
        if(arr[j]>max)
        {
          max=arr[j]
          pos=j
        }
      }
      if(arr[pos]>=arr[pos-1]&&arr[pos]>=arr[pos+1]&&arr[pos]>0)
      {
        if(pos==i)
        {
          start_position=pos+window+1
          printf("%d %f\n", pos+buffer*num_buf, arr[pos]);
          i=pos+window*2+1;
        }
        else
        {
          if(pos>=start_position&&pos>min_pos)
          {
            i=pos;
          }
          else
          {
            i+=window*2+1;
          }
        }
      }
      else
      {
        i+=window+1;
      }
      if(sum==999)
      {
        i+=window*3;
        sum=1;
      }
    }
#  printf("\n");
}
{
  if(FNR==1)
    num_buf=0;
  pos_buf=int($1/buffer);
  if(pos_buf>num_buf)
  {
    max_pos_funct(1,buffer);
    num_buf=pos_buf;
  }
  arr[$1-num_buf*buffer]=$2;
}
END{
  max_pos_funct(1,$1-num_buf*buffer);
}

EOF
)

> ${outfile2}
> ${outfile1}

for seq in `cat ${faseqfile} | tr "\t" "="`; 
  do 
    echo $seq | sed 's/=.*$//' > id; 
    echo $seq | sed 's/^.*=//' > dseq; 
    dseq=`cat dseq`;
    echo ${dseq:${seqstart}:${seqlength}} > dseq 
    id=`cat id`;
    echo ${id}
    cat dseq
    ${call} -m ${patternfile} -s dseq | awk -v id=${id} '{print $0 "\t" id}'  >>  ${outfile1}

    #compute average correlation    
    cat ${outfile1} | gawk '{sum=0; for(i=2;i<=NF;i++) sum+=$i; print $1, sum/(NF-1);}' > avgc
    
    # compute most likely position of the nucleosome 
    cat avgc | awk "$awk_program" window=73 buffer=10000 | awk -v num=$patternfile -v id=$id '{print id "\t" num "\t" $0}' >> ${outfile2} ; 
 done 

rm id dseq avgc 
exit 0

