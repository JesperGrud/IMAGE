seq 1 $(((lines=$(wc -l <STRINGMOTIFS))/PROCESSORS+1)) $lines | sed 'N;s|\(.*\)\(\n\)\(.*\)|\1d;\1,\3w STRINGDIR/\3\2\3|;P;$d;D' | sed -ne :nl -ne '/\n$/!{N;bnl}' -nf - STRINGMOTIFS
