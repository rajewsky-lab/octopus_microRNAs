BEGIN{FS="\t";OFS="\t"}{if ($6=="+") {$6="-"} else {$6="+"}} {print $0}
