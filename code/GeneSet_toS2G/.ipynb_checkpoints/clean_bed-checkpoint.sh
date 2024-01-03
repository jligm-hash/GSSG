bedops_cell=[BEDOPS_DIRECTORY]/BEDOPS/bin
bedtools_cell=[BEDTOOLS_DIRECTORY]/BEDTOOLS/bedtools2/bin
input_cell=$1  #sclinker_beds/{Celltype} [directory where bed files are located]
names=`ls $input_cell | cut -f 1 -d '.'`

#! sort and merge the bed files, enhancer chr, start, end -> target gene score
for name in $names
do
$bedtools_cell/bedtools sort -i $input_cell/$name.bed > $input_cell/$name.2.bed
$bedtools_cell/bedtools merge -i $input_cell/$name.2.bed -c 4 -o max > $input_cell/$name.3.bed
mv $input_cell/$name.3.bed $input_cell/$name.bed
rm $input_cell/$name.2.bed
done
