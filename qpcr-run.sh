#!bin/sh


###############  Reading user defined parameters  ###################
##Read in
while getopts 'd:s:m:xe:i:o:h' optname
  do
    case $optname in
      d)
        data=$OPTARG;;
      s)
        sheetname=$OPTARG;;
      m)
        mode=$OPTARG;;
      x)
        stat="stat";;
      e)
        ec=$OPTARG;;
      i)
        ic=$OPTARG;;
      o)
        out=$OPTARG;;
      h)
        echo "See help information: python qpcrCalculate.py -h"
        exit 1;;
      ?)
        echo "Unknow argument!!!"
        exit 1;;
    esac
  done




# # Read ABI machine output
python qpcrRead.py -d ${data} -o ${out}_extracted

#
echo "Calculating Delta Ct, DDelta Ct, Fold Changes."
python qpcrCalculate.py -d ${out}_extracted \
                        -i ${ic} \
                        -e ${ec} \
                        -m ${mode}
                        -o ${out}_delta
# extract delta ct
if [ $stat == "stat" ]
then
echo "Statistical testing"
python qpcrCalculate.py -d ${out}_delta.csv \
                        -i ${ic} \
                        -e ${ec} \
                        -m stat
                        -o ${out}_ttest
else
    echo "Skip statistical testing"
fi

echo "Finish pipeline."
