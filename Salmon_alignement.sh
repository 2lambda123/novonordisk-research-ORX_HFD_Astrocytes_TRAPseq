###############################
###############################

# Salmon indexes

# # https://github.com/COMBINE-lab/salmon/issues/603
# # similar to here https://github.com/COMBINE-lab/salmon/issues/214
# 
# mkdir /opt/genomes/mouse_mm10/Salmon_homemade/
# cd /opt/genomes/mouse_mm10/Salmon_homemade/
# 
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.transcripts.fa.gz
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
# 
# gunzip gencode.vM23.transcripts.fa.gz
# gunzip GRCm38.primary_assembly.genome.fa.gz
# 
# grep "^>" GRCm38.primary_assembly.genome.fa | cut -d " " -f 1 > GRCm38.decoys.txt
# sed -i 's/>//g' GRCm38.decoys.txt
# cat gencode.vM23.transcripts.fa GRCm38.primary_assembly.genome.fa | gzip > GRCm38.gentrome.fa.gz
# salmon index -t GRCm38.gentrome.fa.gz -d GRCm38.decoys.txt -p 12 -i salmon_index --gencode
# 
#
# 
mkdir ./Salmon_homemade_mm10_bootstraped_update/

ls ./Output/trim_fastq/*_trimmed.fq.gz >  ./sample_list_trimmedfastq.txt
 
while read p;
do
echo "$p"
salmon quant -i /opt/genomes/mouse_mm10/Salmon_homemade/salmon_index -l SF --noLengthCorrection --validateMappings --numBootstraps 100 -r <(zcat "$p") -o ./Salmon_homemade_mm10_bootstraped_update/"$p"
done < ./sample_list_trimmedfastq.txt

multiqc ./Salmon_homemade_mm10_bootstraped_update/  -d -f -s -v -p

