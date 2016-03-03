# Bash code

**Sync with NCBI's ftp server**

```bash
# add time stamp to log file
lftp -c 'open -e "mirror -c -p --no-empty-dirs -I *assembly*.txt -P=5 --log=lftp_log.txt /genomes/genbank/bacteria/ ~/MGGen/ncbi_bacteria_mirror" ftp.ncbi.nlm.nih.gov'
```

**Get assembly summary for specific organism**

```bash
wget -O Acinetobacter_nosocomialis_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Acinetobacter_nosocomialis/assembly_summary.txt
awk -F "\t" '/Acinetobacter/ {print $20}' Acinetobacter_nosocomialis_summary.txt | \
sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/)(GCA_.+)|\1\2/\2_genomic.fna.gz|'>genome_urls.txt
wget -P ./genome_files/ --input genome_urls.txt
gunzip ./genome_files/*.gz
```

**Get assembly summary for all bacteria**

Since it downloads every summary as `assembly_summary.txt` in a separate directory named after the bacteria, we will have to copy and rename the files.  The files must remain in the structure that `lftp` puts them in in order to stay synced with NCBI.

```bash
mkdir bacteria_summaries && cd bacteria_summaries
wget -rl 0 -A "*assembly*.txt" -nH --cut-dirs=3 --no-remove-listing ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria
for subdir in *; do cp ./$subdir/assembly_summary.txt ./renamed/$subdir.txt; done;
```

*It has come to my attention that most organisms (~8K out of 12K) do NOT have a `assembly_summary.txt` file in the expected directory.  Running another test to find out where/if the equivalent file is located.*
