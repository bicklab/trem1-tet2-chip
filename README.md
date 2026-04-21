#  TREM1 signaling associated with increased cardiovascular risk in TET2 clonal hematopoiesis 

## Study Cohort
The UK Biobank is a large-scale, prospective cohort study comprising over 500,000 men and women aged 40 to 70 years at the time of recruitment (2006–2010) from across the United Kingdom. For this analysis, we included a subset of 487,409 participants with available whole-exome sequencing data (40x). The mean age at the time of DNA sample collection was 56.5 years (±8.1), and 222,905 participants (46%) were male. BioVU is a biobank at Vanderbilt University Medical Center, containing de-identified electronic health records (EHRs) linked to genetic data, spanning from 2006 to 2021. Whole-genome sequencing (40x) was available for 250,391 participants. The mean age at the time of DNA sample collection was 51.1 years (±22.1), with 103,866 (42%) male participants. 

## CHIP variant calls
Putative somatic SNPs and short indels were called with GATK Mutect2 (https://software.broadinstitute.org/gatk). Mutect2 first identifies candidate sites with evidence of variation and then performs local reassembly to refine variant calling. It employs an external reference, known as a “panel of normal samples,” to filter out recurrent sequencing artifacts and calls variants only at sites with evidence of somatic variation. The panel of normal samples comprised 100 randomly selected individuals under 40 years of age, with the absence of hotspot CHIP mutations confirmed prior to their inclusion. 

Variants included in the preliminary dataset met the following criteria: presence in a pre-established list of candidate CHIP variants, total sequencing depth ≥ 20, alternate allele read depth count ≥ 5, and representation in both sequencing directions (i.e., F1R2 ≥ 1 and F2R1 ≥ 1). CHIP mutations were defined as those with a variant allele fraction (VAF) ≥ 0.02.


## Calculation for genetically predicted TREM1 levels
We calculated genetically predicted scores of TREM1 expression using 49 variants from a model trained with Somalogic protein measurements from the INTERVAL study (OPGS000019)

Code:
  

    file_df = pd.read_table('filetable.txt')
    batch_df = pd.read_table('batch_file.txt') 

    
    for chrom in range(1,23): 
      bgen = batch_df.iloc[chrom][1]
      bgen_bgi = batch_df.iloc[chrom][2]
      sample = batch_df.iloc[chrom][3]
      for index, row in df.iterrows():
        name = row['file_name']
        file_id = row['file_id']
        command = f"dx run swiss-army-knife \
      --priority high \
      --instance-type "mem3_ssd1_v2_x16" \
      -iin=project-GbjJ860J2p5ZJqZKZZVY883Z:file-FxY60F8JkF63k38X9V5BFY55 \
      -iin=project-GbjJ860J2p5ZJqZKZZVY883Z:file-FxZ2byjJkF65g2vX9Vx8vgb3 \
      -iin=project-GbjJ860J2p5ZJqZKZZVY883Z:file-GbxY6YQJ5YqBqx9qX6Yyv4Bk \
      -iin=project-GbjJ860J2p5ZJqZKZZVY883Z:file-Gx33yX8J2p5bpFKJQ356G20G \
      -icmd="plink2 --bgen ukb22828_c22_b0_v3.bgen ref-first \
                  --sample ukb22828_c22_b0_v3.sample \
                  --score TREM1.txt 1 4 6 header list-variants cols=scoresums \
                  --out TREM1" \
      --destination "project-GbjJ860J2p5ZJqZKZZVY883Z:/omics_pred/score_files/" \
      --tag plink2_score \
      -y"



      %%sh
      dx run swiss-army-knife \
         --priority low \
         --instance-type "mem3_ssd1_v2_x16" \
         -iin=project-GbjJ860J2p5ZJqZKZZVY883Z:file-FxY60F8JkF63k38X9V5BFY55 \
         -iin=project-GbjJ860J2p5ZJqZKZZVY883Z:file-FxZ2byjJkF65g2vX9Vx8vgb3 \
         -iin=project-GbjJ860J2p5ZJqZKZZVY883Z:file-GbxY6YQJ5YqBqx9qX6Yyv4Bk \
         -iin=project-GbjJ860J2p5ZJqZKZZVY883Z:file-Gx33yX8J2p5bpFKJQ356G20G \
         -icmd="plink2 --bgen ukb22828_c22_b0_v3.bgen ref-first \
                  --sample ukb22828_c22_b0_v3.sample \
                  --score TREM1.txt 1 4 6 header list-variants cols=scoresums \
                  --out TREM1" \
         --destination "project-GbjJ860J2p5ZJqZKZZVY883Z:/omics_pred/score_files/" \
         --tag plink2_score \
         -y


## Association between gTREM1 and CVD incidence
BioVU Code available in: [BioVU code](https://github.com/bicklab/trem1-tet2-chip/blob/main/KZ_TREM1_CVD_BioVU.r)

UKB Code available in: [UKB code]([https://github.com/bicklab/trem1-tet2-chip/KZ_TREM1_CVD_UKB.r](https://github.com/bicklab/trem1-tet2-chip/blob/main/KZ_TREM1_CVD_UKB.r))

## Data
This analysis was performed on the [UK Biobank DNA Nexus Research Analysis Platform](https://ukbiobank.dnanexus.com) and BioVU Terra.bio environment.

## Acknowledgements
Individual-level sequence data and CHIP calls have been deposited with UK Biobank and are available to approved researchers by application (https://www.ukbiobank.ac.uk/register-apply/). Vanderbilt BioVU data are available through an application to the Vanderbilt Institute for Clinical and Translational Research (VICTR) BioVU Review Committee.

## Contact
Kun Zhao, kun.zhao@vumc.org; Yash Pershad, yash.pershad@vanderbilt.edu
