FNR == NR {
  arr[$1][$2]["min"]=$3
  arr[$1][$2]["max"]=$4
  arr[$1][$5]["min"]=$6
  arr[$1][$5]["max"]=$7
  next
}
#$2 in arr 
{
  split($1,type,"-")
  ( $4 > arr[$2][type[1]]["min"] && $4 < arr[$2][type[1]]["max"]) ? decision="WithinLimits" : decision="OutOfLimits"
  printf"%s\t%s-%s-%d-%d\n",$0,decision,type[1],arr[$2][type[1]]["min"],arr[$2][type[1]]["max"]
}
