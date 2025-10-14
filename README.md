# OrthoFinder_GO

This analysis is mainly to asign GO terms for orthogroups and summarize based on UP/DOWN-exression.

Please note, in-house python scripts can be found in <[02.GO_annotation](https://github.com/chongjing/OrthoFinder_GO/tree/main/02.GO_annotation/1.interproscan)> and <[03.Summary_UP_DOWN_for_each_orthogroup](https://github.com/chongjing/OrthoFinder_GO/tree/main/03.Summary_UP_DOWN_for_each_orthogroup)> directories.

## 1.Orthogroups
```bash
cd /data/pathology/cxia/projects/Sebastian/Victor/1.orthofinder/
python /data/pathology/cxia/program/OrthoFinder_source/orthofinder.py -f ./1.compliantFasta --threads 32 --algthreads 32
```

## 2.GO annotation
```bash
cd /data/pathology/cxia/projects/Sebastian/Victor/1.orthofinder/1.compliantFasta/01.arabidopsis_introproscan
export JAVA_HOME=/data/pathology/program/java/jdk-21.0.2/
export PATH=/data/pathology/program/InterProScan/interproscan-5.66-98.0:/data/pathology/program/InterProScan/interproscan-5.66-98.0/bin:$JAVA_HOME/bin:$PATH
for file in ../*.faa; do
    echo "Processing interproscan for ${file}";
    /data/pathology/program/InterProScan/interproscan-5.66-98.0/interproscan.sh -i ${file} -f tsv -appl Pfam --goterms -pa --iprlookup --cpu 16;
done
```

## 3. Summarize UP and DOWN genes for each orthogroup
```bash
#using a set of inhouse script to format and append orthogroup and function to DE results: 7.DEGs_arab_log2FC_down.withOG_Function; 7.DEGs_arab_log2FC_up.withOG_Function
cd /data/pathology/cxia/projects/Sebastian/Victor/1.orthofinder/02.GO_annotation/1.interproscan
python 02.remove.transcript.versions.from.orthofinder.py ../../1.compliantFasta/OrthoFinder/Results_Feb05/Orthogroups/Orthogroups.tsv 5.Orthogroups.noTranscripts.no.space.tsv
python 03.append.OrthologGroup.py 5.Orthogroups.noTranscripts.no.space.tsv DEGs_arab_log2FC_down 6.DEGs_arab_log2FC_down.withOG.removed.space
# finally, add description to DEG results
python 04.extractGeneFunction.to.DEresults.py 6.DEGs_arab_log2FC_down.withOG.removed.space 4.Arabidopsis_thaliana.TAIR10.annotation.reformat.list 8.DEGs_arab_log2FC_down.withOG_Function

# To summarize UP and DOWN genes for each OrthologGroup
# arabidopsis
awk -F"\t" '{print $1"\t"$2}' ../Orthogroups.ID_Changed.lettuce.tsv | sed -e 's/,.........\.[2-9]//g' -e 's/,.........\.1[0-9]//g' -e 's/\..//g' -e 's/\...//g' > 01.Orthogroups.Arab.tsv
awk -F"\t" '{print $1"\t"$4}' Orthogroups.ID_Changed.lettuce.tsv | sed -e 's/,Zm.*._P00[2-9]//g' -e 's/,Zm.*._P01[0-9]//g' -e 's/,Zm.*._P02[0-9]//g' | sed 's/_P...//g' > 03.Maz/01.Orthogroups.Maz.tsv
# to label UP/DOWN for each gene
cd /data/pathology/cxia/projects/Sebastian/Victor/1.orthofinder/03.Summary_UP_DOWN_for_each_orthogroup
python 00.ortho_express_summery.py 07.Saff/01.Orthogroups.Saff.tsv DEGs_saff_log2FC_up DEGs_saff_log2FC_down 07.Saff/04.Orthogroups.Saff.Express_label.tsv
# to summarize numbers
python 00.ortho_express_summery.further.py 07.Saff/02.Orthogroups.Saff.Express_label.tsv 07.Saff/05.Orthogroups.Saff.Express_number.tsv
paste -d$'\t' <(tr -d '\r' < 01.Arab/02.Orthogroups.Arab.Express_label.tsv) <(tr -d '\r' < 02.Let/02.Orthogroups.Lettuce.Express_label.tsv) <(tr -d '\r' < 03.Maz/02.Orthogroups.Maz.Express_label.tsv) <(tr -d '\r' < 04.Med/02.Orthogroups.Med.Express_label.tsv) <(tr -d '\r' < 05.Mel/02.Orthogroups.Mel.Express_label.tsv) <(tr -d '\r' < 06.Rice/02.Orthogroups.Rice.Express_label.tsv) <(tr -d '\r' < 07.Saff/04.Orthogroups.Saff.Express_label.tsv) <(tr -d '\r' < 08.Tobac/02.Orthogroups.Tobac.Express_label.tsv) <(tr -d '\r' < 09.Tom/02.Orthogroups.Tom.Express_label.tsv) <(tr -d '\r' < 10.Wmel/02.Orthogroups.Wmel.Express_label.tsv) | awk -F"\t" '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12"\t"$14"\t"$16"\t"$18"\t"$20}' > Orthogroups.ExpressionLabel.1.tsv
paste -d$'\t' <(tr -d '\r' < 01.Arab/03.Orthogroups.Arab.Express_number.tsv) <(tr -d '\r' < 02.Let/03.Orthogroups.Lettuce.Express_number.tsv) <(tr -d '\r' < 03.Maz/02.Orthogroups.Maz.Express_number.tsv) <(tr -d '\r' < 04.Med/03.Orthogroups.Med.Express_number.tsv) <(tr -d '\r' < 05.Mel/03.Orthogroups.Mel.Express_number.tsv) <(tr -d '\r' < 06.Rice/03.Orthogroups.Rice.Express_number.tsv) <(tr -d '\r' < 07.Saff/05.Orthogroups.Saff.Express_number.tsv) <(tr -d '\r' < 08.Tobac/03.Orthogroups.Tobac.Express_number.tsv) <(tr -d '\r' < 09.Tom/03.Orthogroups.Tom.Express_number.tsv) <(tr -d '\r' < 10.Wmel/03.Orthogroups.Wmel.Express_number.tsv) | awk -F"\t" '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12"\t"$14"\t"$16"\t"$18"\t"$20}' > Orthogroups.ExpressionNumber.1.tsv
```

Expected result:
<img src="https://github.com/chongjing/OrthoFinder_GO/blob/main/03.Summary_UP_DOWN_for_each_orthogroup/Orthogroups.ExpressionNumber.1.jpg" alt="Image 1" width="600"/>

## 4. add GO terms for UPs and DOWNs
```bash
cd /data/pathology/cxia/projects/Sebastian/Victor/1.orthofinder/03.Summary_UP_DOWN_for_each_orthogroup/002.summary_with_GO
01: get only UP genes for each species
python 01.get_UP.py ../01.Arab/02.Orthogroups.Arab.Express_label.tsv 001.Arab.UP.label.tsv
python 01.get_DOWN.py ../01.Arab/02.Orthogroups.Arab.Express_label.tsv 001.Arab.DOWN.label.tsv
python 01.get_DOWN.py ../06.Rice/02.Orthogroups.Rice.Express_label.3.tsv 006.Rice.DOWN.label.tsv
python 01.get_DOWN.py ../07.Saff/04.Orthogroups.Saff.Express_label.tsv 007.Saff.DOWN.label.tsv

### protein name not match gene name for lettuce
egrep 'protein_id' edited_not_question_mark_GCF_002870075.4_Lsat_Salinas_v11_annotation.gtf | awk -F"\t" '{print $9}' | awk -F";" '{print $5"\t"$7}' | sed -e 's/gene //g' -e 's/protein_id //g' -e 's/"//g' | sort -k2,2 | uniq > edited_not_question_mark_GCF_002870075.4_Lsat_Salinas_v11_annotation.geneID_proteinID.list
awk 'NR==FNR{a[$2]=$1; next} $1 in a{$1=a[$1]}1' /home/vhm24/files_hosts/reads/multimapping/lettuce/genome_annotation/edited_not_question_mark_GCF_002870075.4_Lsat_Salinas_v11_annotation.geneID_proteinID.list lettuce.GO.list | uniq > lettuce.GO.transcriptID.list
### protein name not match gene name for lettuce
awk -F"\t" '$3 == "transcript" {print $9}' /home/vhm24/files_hosts/reads/multimapping/medicago/genome_annotation/medicago_truncatula.medtat17_4.0.54.gtf | awk -F";" '{print $1" "$2}' | awk -F" " '{print $2"\t"$4}' | sed 's/"//g' > /home/vhm24/files_hosts/reads/multimapping/medicago/genome_annotation/medicago_truncatula.medtat17_4.0.54_GeneID-TranscriptID.list
awk 'NR==FNR{a[$2]=$1; next} $1 in a{$1=a[$1]}1' /home/vhm24/files_hosts/reads/multimapping/medicago/genome_annotation/medicago_truncatula.medtat17_4.0.54_GeneID-TranscriptID.list ../../1.compliantFasta/01.medicago_introproscan/medicago.GO.list | sort -k1,1 | uniq > ../../1.compliantFasta/01.medicago_introproscan/medicago.GO.transcriptID.list

# add GO terms for orthogroups
python 02.add_GO_to_Orthogroup.py 001.Arab.DOWN.label.tsv ../../1.compliantFasta/01.arabidopsis_introproscan/arabidopsis.GO.list 101.Arab.DOWN.GO-label.tsv
python 02.add_GO_to_Orthogroup.py 001.Arab.UP.label.tsv ../../1.compliantFasta/01.arabidopsis_introproscan/arabidopsis.GO.list 101.Arab.UP.GO-label.tsv
python 02.add_GO_to_Orthogroup.py 002.Lettuce.DOWN.label.tsv ../../1.compliantFasta/01.lettuce_introproscan/lettuce.GO.transcriptID.list 102.Lettuce.DOWN.GO-label.tsv
python 02.add_GO_to_Orthogroup.py 002.Lettuce.UP.label.tsv ../../1.compliantFasta/01.lettuce_introproscan/lettuce.GO.transcriptID.list 102.Lettuce.DOWN.UP-label.tsv
python 02.add_GO_to_Orthogroup.py 004.Med.DOWN.label.tsv ../../1.compliantFasta/01.medicago_introproscan/medicago.GO.transcriptID.list 104.Med.DOWN.GO-label.tsv
python 02.add_GO_to_Orthogroup.py 004.Med.UP.label.tsv ../../1.compliantFasta/01.medicago_introproscan/medicago.GO.transcriptID.list 104.Med.UP.GO-label.tsv

#concate all GO-label for orthogroups
paste -d'\t' <(awk '{print $1"\t"$2}' 101.Arab.UP.GO-label.tsv) <(awk '{print $2}' 102.Lettuce.UP.GO-label.tsv) <(awk '{print $2}' 103.Maz.UP.GO-label.tsv) <(awk '{print $2}' 104.Med.UP.GO-label.tsv) <(awk '{print $2}' 105.Mel.UP.GO-label.tsv) <(awk '{print $2}' 106.Rice.UP.GO-label.tsv) <(awk '{print $2}' 107.Saff.UP.GO-label.tsv) <(awk '{print $2}' 108.Tobac.UP.GO-label.tsv) <(awk '{print $2}' 109.Tom.UP.GO-label.tsv) <(awk '{print $2}' 110.Wmel.UP.GO-label.tsv) > 20.UP.GO-label.tsv
```

Expected output:
<img src="https://github.com/chongjing/OrthoFinder_GO/blob/main/03.Summary_UP_DOWN_for_each_orthogroup/002.summary_with_GO/20.UP.GO-label.jpg" alt="Image 1" width="600"/>

## 5. ENA EMBL submission
5.1 install webin-cli https://ena-docs.readthedocs.io/en/latest/submit/general-guide/webin-cli.html
```bash
cd /home/cx264/rds/rds-csc_programmes-FTKWLWDeHys/programs/webin-cli
wget https://github.com/enasequence/webin-cli/releases/download/9.0.1/webin-cli-9.0.1.jar
``
5.2 login Webin submission Portal to register study and samples, via https://www.ebi.ac.uk/ena/submit/webin/login

| Type       | Accession    | Unique name (alias)                  |
|------------|--------------|--------------------------------------|
| Project    | PRJEB100529  | 26ca3a6e-6baa-4f04-a327-e29411fc3df5 |
| Submission | ERA35056725  | SUBMISSION-11-10-2025-07:43:28:740   |

5.3 Submission
After get sample accessions, prepare metadata table with header like

```bash
cd ~/rds/rds-csc_programmes-FTKWLWDeHys/chongjing/01.Sebastian/01.Victor
#prepare a sample file fastq2_template_1760166804427.tsv:
| sample_id   | study      | sample_acc  | instrument          | library_name | source         | selection   | strategy | layout | fq1             | fq1_md5                          | fq2             | fq2_md5                          |
|-------------|------------|-------------|---------------------|--------------|----------------|-------------|----------|--------|-----------------|----------------------------------|-----------------|----------------------------------|
| ERS27017987 | PRJEB100529 | ERS27017987 | Illumina HiSeq 2500 | arab_c1     | TRANSCRIPTOMIC | RANDOM PCR | RNA-Seq | PAIRED | arab_c1_1.fq.gz | a2dcae8f99493c8e84036a5b121a82d3 | arab_c1_2.fq.gz | 67b7a234abe3a1bcd07e6a4db967e3e8 |
| ERS27017988 | PRJEB100529 | ERS27017988 | Illumina HiSeq 2500 | arab_c2     | TRANSCRIPTOMIC | RANDOM PCR | RNA-Seq | PAIRED | arab_c2_1.fq.gz | d8f51409fd28daa3c7c33161e88fe2c8 | arab_c2_2.fq.gz | 749f619d33ab0a5e95b659157b36034a |
```
prepare `gene_manifests.py`:
```python
import csv, os

# Input TSV
INPUT = "fastq2_template_1760166804427.tsv"
# Output dirs
MANIFEST_DIR = "manifests"
CHECKSUM_FILE = "checksums.md5"

os.makedirs(MANIFEST_DIR, exist_ok=True)

checksums = []

with open(INPUT) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        manifest = f"""STUDY\t{row['study']}
SAMPLE\t{row['sample_acc']}
NAME\t{row['library_name']}
INSTRUMENT\t{row['instrument']}
LIBRARY_NAME\t{row['library_name']}
LIBRARY_SOURCE\t{row['source']}
LIBRARY_SELECTION\t{row['selection']}
LIBRARY_STRATEGY\t{row['strategy']}
FASTQ\t{row['fq1']}
FASTQ\t{row['fq2']}
"""
        fname = os.path.join(MANIFEST_DIR, f"{row['library_name']}.manifest.txt")
        with open(fname, "w") as out:
            out.write(manifest)
        print(f"Wrote {fname}")

        # Collect checksums
        checksums.append(f"{row['fq1_md5']}  {row['fq1']}")
        checksums.append(f"{row['fq2_md5']}  {row['fq2']}")

# Write combined checksums file
with open(CHECKSUM_FILE, "w") as out:
    out.write("\n".join(checksums) + "\n")

print(f"Wrote {CHECKSUM_FILE}")
```

prepare submission script:

```bash
#!/usr/bin/env bash
set -euo pipefail

WEBIN_USER="Webin-69760"
WEBIN_PASS="" ###password here
JAR="$HOME/rds/rds-csc_programmes-FTKWLWDeHys/programs/webin-cli/webin-cli-9.0.1.jar"
OUTDIR="./webin_out"

mkdir -p "$OUTDIR"

for m in manifests/*.manifest.txt; do
  echo "Validating $m"
  /rds/project/rds-FTKWLWDeHys/programs/java/jdk-17.0.10/bin/java -jar "$JAR" -context reads -manifest "$m" \
       -userName "$WEBIN_USER" -password "$WEBIN_PASS" \
       -outputDir "$OUTDIR" -validate

  echo "Submitting $m"
  /rds/project/rds-FTKWLWDeHys/programs/java/jdk-17.0.10/bin/java -jar "$JAR" -context reads -manifest "$m" \
       -userName "$WEBIN_USER" -password "$WEBIN_PASS" \
       -outputDir "$OUTDIR" -submit
done
```

Submission.
```bash
python3 generate_manifests.py
./submit_reads.sh
```
Successful submission should look like this:
```bash
% ./submit_reads.sh
Validating manifests/arab_c1.manifest.txt
INFO : Your application version is 9.0.1
INFO : Submission(s) validated successfully.
INFO : Creating report file: /rds/project/rds-FTKWLWDeHys/chongjing/01.Sebastian/01.Victor/./webin_out/./webin-cli.report
Submitting manifests/arab_c1.manifest.txt
INFO : Your application version is 9.0.1
INFO : Connecting to FTP server : webin2.ebi.ac.uk
INFO : Creating report file: /rds/project/rds-FTKWLWDeHys/chongjing/01.Sebastian/01.Victor/./webin_out/./webin-cli.report
INFO : Uploading file: /rds/project/rds-FTKWLWDeHys/chongjing/01.Sebastian/01.Victor/arab_c1_1.fq.gz
INFO : Uploading file: /rds/project/rds-FTKWLWDeHys/chongjing/01.Sebastian/01.Victor/arab_c1_2.fq.gz
INFO : Files have been uploaded to webin2.ebi.ac.uk.
INFO : The submission has been completed successfully. The following run accession was assigned to the submission: ERR15696041
INFO : The submission has been completed successfully. The following experiment accession was assigned to the submission: ERX15100159
```

