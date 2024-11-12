#!/usr/bin/bash

Usage()
{
    echo -e "本程序用于提取单个per_seq.txt文件的结果\nUsage: $0 -f contig_file -a protein_file -i per_seq_file -c cut_off=0"
}
cut_off=0
lines1='【------------------------------'
lines2='------------------------------】'

while getopts ':f:a:i:c:r:' OPT; do
    case $OPT in
        f) contig="$OPTARG";;
        a) protein="$OPTARG";;
        i) per_seq="$OPTARG";;
        c) cut_off="$OPTARG";;
        *) Usage; exit 1;;
    esac
done
if [ -z $contig ];then Usage; exit 1; fi
if [ -z $protein ];then Usage; exit 1; fi
if [ -z $per_seq ];then Usage; exit 1; fi

outdir=${per_seq%/*}
b=recombinases.txt

grep -P "^[^#]" $per_seq > ${outdir}/tmp
if [ `wc -l ${outdir}/tmp|awk '{print $1}'` == 0 ];then rm ${outdir}/tmp;echo "${lines1}No recombinase, exit.${lines2}";exit 0;fi
ind1=`awk '{print $3}' ${outdir}/tmp|grep Recombinase -c`
ind2=`awk '{print $3}' ${outdir}/tmp|grep Resolvase -c`
if [ $ind1 == 0 -a $ind2 == 0 ];then rm ${outdir}/tmp;echo "${lines1}No recombinase, exit.${lines2}";exit 0;fi

echo -n > ${outdir}/${b}
echo -n > ${outdir}/hitted_contigs.fa
echo "${lines1}开始写入文件：${outdir}/${b}，${outdir}/hitted_contigs.fa${lines2}"
while read line
do
    l=${line%%\#*}
    r=${line#*\#}
    hmm_start=`echo ${r}|awk -v FS="[ #]*" '{print $1}'`
    hmm_end=`echo ${r}|awk -v FS="[ #]*" '{print $2}'`
    hmm_strand=`echo ${r}|awk -v FS="[ #]*" '{print $3}'`
    line1="${l}\"#${r}\""
    i=${line1%% *}
    j=${line1#* }
    ii=${i%_*}

    contig_seq=`seqkit grep -I -w 0 -p $ii ${contig}|tail -1`
    protein_seq=`seqkit grep -I -w 0 -p $i ${protein}|tail -1`
    echo "$i $hmm_start $hmm_end $hmm_strand ${#contig_seq} $j $contig_seq $protein_seq" >> ${outdir}/${b}
    echo ">${i}" >> ${outdir}/hitted_contigs.fa
    echo $contig_seq >> ${outdir}/hitted_contigs.fa
done < ${outdir}/tmp
seqkit rmdup -n ${outdir}/hitted_contigs.fa -o ${outdir}/hitted_contigs.uniq.fa -w 0
rm ${outdir}/hitted_contigs.fa
echo "完成。"
rm ${outdir}/tmp