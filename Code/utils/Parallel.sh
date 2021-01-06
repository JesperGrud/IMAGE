Processors=$1
Dir=$2
String=$3
cd $Dir
sed 's/>/\n>/g' utils/Collection.motif > tmp/$String.motif
cp utils/SplitMotifLibrary.sh tmp/$String.sh
cd tmp
stringdir=$String"_motifs"
stringfile=$String".motif"
mkdir $stringdir
sed "s/STRINGMOTIFS/$stringfile/g" $String.sh > $String.tmp
mv $String.tmp $String.sh 
sed "s/PROCESSORS/$Processors/g" $String.sh > $String.tmp
mv $String.tmp $String.sh 
sed "s/STRINGDIR/$stringdir/g" $String.sh > $String.tmp
mv $String.tmp $String.sh 
sh $String.sh
rm $String.sh
VAR=$(ls $stringdir/)
for i in $VAR; do echo $i >> $String.SplitFiles; done
cat $String.SplitFiles | sort -k1n,1 > $String.SplitSort
ProcessorsCount=$(($Processors-1))
for i in $(seq 1 $ProcessorsCount); do File=$(awk -v VAR=$i 'NR == VAR { print $0 }' $String.SplitSort); cat $stringdir/$File >> $String.Reassemble; done
rm $String.SplitFiles
rm $String.SplitSort
diff $String.motif $String.Reassemble | sed '1d' | tr -d '< ' > $stringdir/End
rm $String.Reassemble
rm $String.motif
cd $stringdir
VAR=$(ls)
for i in $VAR; do awk -v RS="" '{gsub (/\n>/,"")}1' $i > Tmp; mv Tmp $i; done
