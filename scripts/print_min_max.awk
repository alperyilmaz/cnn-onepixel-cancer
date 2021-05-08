{ 
  # thank you anubhava @stackoverflow for the help: https://stackoverflow.com/a/67442465/4579196
  split($2,arr,"-")
  if (max[$1][arr[1]] < $3) max[$1][arr[1]] = $3
  if (!min[$1][arr[1]] || min[$1][arr[1]] > $3) min[$1][arr[1]] = $3 
}
END{
  for (i in max) {
    printf "%s", i 
    for (j in max[i]) printf " %s %d %d", j, min[i][j], max[i][j]; print ""
  }
}
