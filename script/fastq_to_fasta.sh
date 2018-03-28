if [ $# -ne 3 ];then
        echo "usage: [xx.fastq.gz] [xx.fa]"
        exit 1
fi
fastq_gz=$1
fasta=$2
awk '{if ((NR - 1)% 4 <= 1) {print $0}}' <(gzip -dc $fastq_gz) | sed 's/@/>/g' > $fasta
