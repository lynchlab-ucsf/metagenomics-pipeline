#!/bin/bash
#$ -cwd
#$ -l mem_free=1G
#$ -l h_rt=12:00:00
## #$ -m e
## #$ -M kathryn.mccauley@ucsf.edu


software_location=/wynton/group/lynch/kmccauley/mySoftware/
## Note: first time this was run, three samples were still outstanding
FASTA_DIRECTORY=/wynton/group/lynch/kmccauley/Metagenomics_CREW/dyad_mapping_v2/


####
cd $FASTA_DIRECTORY

mkdir -p humann3_combined
for i in */humann3_results/*genefamilies.tsv; do
id=$(echo $i | cut -d '/' -f1)
cp $i humann3_combined/${id}_genefamilies.tsv
done

for i in */humann3_results/*pathabundance.tsv; do
id=$(echo $i | cut -d '/' -f1)
cp $i humann3_combined/${id}_pathabundance.tsv
done

singularity exec -B $PWD:$PWD $software_location/humann.img humann_join_tables -i humann3_combined -o combined_genefamilies.tsv --file_name genefamilies
singularity exec -B $PWD:$PWD $software_location/humann.img humann_renorm_table -i combined_genefamilies.tsv -o combined_genefamilies-cpm.tsv --units cpm --update-snames
singularity exec -B $PWD:$PWD,/wynton/group/lynch/Shared/humann_db/utility_mapping/:/mnt $software_location/humann.img humann_regroup_table -i combined_genefamilies-cpm.tsv -o level4ec-cpm.tsv -c /mnt/map_level4ec_uniref90.txt.gz
singularity exec -B $PWD:$PWD $software_location/humann.img humann_rename_table -i level4ec-cpm.tsv -o level4ec-cpm-named.tsv --names ec

singularity exec -B $PWD:$PWD $software_location/humann.img humann_join_tables -i humann3_combined -o combined_pathabundance.tsv --file_name pathabundance
singularity exec -B $PWD:$PWD $software_location/humann.img humann_renorm_table -i combined_pathabundance.tsv -o combined_pathabundance-cpm.tsv --units cpm --update-snames

