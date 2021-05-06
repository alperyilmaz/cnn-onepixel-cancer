export PROCESSED="processed-data"

function _filtergenes() { awk 'FNR==NR {genes[$2]++; next} NR > 2 && FNR==1 {print $0;next}  {split($1,parse,".");$1=parse[1];if(parse[1] in genes){print $0}}' $PROCESSED/genes_high_deg -; }

export -f _filtergenes

function _filtersamples() { awk 'FNR==NR {samples[$1]++; next} ($2 in samples)' $PROCESSED/sample_labels - ; }

export -f _filtersamples

function _pivotlonger() { awk 'NR==1{split($0,header,"\t"); next} {for (j=2;j<NF;j++){printf"%s\t%s\t%s\n",$1,header[j],2^$j}}'; }

export -f _pivotlonger

function _pivotwider() { awk '{thedata[$2][$1]=$3}END{asorti(thedata,sindex); asorti(thedata[sindex[1]],gindex); printf"Sample\t"; for (g in gindex){printf"%s\t",gindex[g]};print"";for (s in sindex){ printf"%s\t",sindex[s]; for (g in gindex){printf"%s\t",thedata[sindex[s]][gindex[g]]};print ""; }}' | sed -e 's/\t$//'; }

export -f _pivotwider

function _max_rgb() { awk '{if ($2 > dict[$1]["r"]){ dict[$1]["r"]=$2 }; if ($3 > dict[$1]["g"]){ dict[$1]["g"]=$3 }; if ($4 > dict[$1]["b"]){ dict[$1]["b"]=$4 } }END{for (gene in dict){printf"%s\t%s\t%s\t%s\n",gene, dict[gene]["r"], dict[gene]["g"], dict[gene]["b"]}}' |sed -e 's/\t\t/\t0\t/' -e 's/\t\t/\t0\t/' ; }

export _max_rgb

function _convert_rgb() { awk '{b=and($3,255); g=and(rshift($3,8),255); r=and(rshift($3,16),255); print $1, r, g, b}' ; }

export _convert_rgb

function _get_max_tcga_gtex() { awk '{printf"%s\t%s\t%.0f\n",$1,$2,$3}'  | awk '{ split($2,sample,"-"); ($3 > count[$1][sample[1]]) ? count[$1][sample[1]]=$3 : 1 }END{for(gene in count){ printf"%s\t%s\t%s\n", gene, count[gene]["GTEX"], count[gene]["TCGA"]}}'; }

export -f _get_max_tcga_gtex

function _count_up_down_cases() { awk 'NR>3 {sum=0; for(j=4;j<=NF;j++){sum+=$j^2}; printf"%s\t%s\n",$1,sum}'; }

export -f _count_up_down_cases

function _get_ensembl() { awk 'FNR==NR {conv[$1]=$2; next} ($1 in conv) {printf"%s\t%s\t%s\n",$1,conv[$1],$2 }' $PROCESSED/entrez2ensembl - ; }

export -f _get_ensembl

function _has_deseq() { awk 'FNR==NR {split($1,parse,"."); deseq[parse[1]]++;next} ($2 in deseq)' $PROCESSED/genes-deseq2 - ; }

export -f _has_deseq
