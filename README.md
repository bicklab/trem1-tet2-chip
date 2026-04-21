#  TREM1 signaling associated with increased cardiovascular risk in TET2 clonal hematopoiesis 

## Study Cohort
The UK Biobank is a large-scale, prospective cohort study comprising over 500,000 men and women aged 40 to 70 years at the time of recruitment (2006–2010) from across the United Kingdom. For this analysis, we included a subset of 487,409 participants with available whole-exome sequencing data (40x). The mean age at the time of DNA sample collection was 56.5 years (±8.1), and 222,905 participants (46%) were male. A total of 192 individuals were identified as carriers of JAK2-CHIP mutations based on our detection pipeline (detailed below). These carriers had a mean age of 61.2 years (±6.2), and 118 (61.5%) were male.

BioVU is a biobank at Vanderbilt University Medical Center, containing de-identified electronic health records (EHRs) linked to genetic data, spanning from 2006 to 2021. Whole-genome sequencing (40x) was available for 250,391 participants. The mean age at the time of DNA sample collection was 51.1 years (±22.1), with 103,866 (42%) male participants. Given the wide age range in BioVU, we restricted our CHIP detection to individuals aged 40 to 70 years. We identified 179 individuals carrying JAK2-CHIP mutations. The mean age of these carriers was 59.2 years (±8.8), and 84 (47%) were male. 


## CHIP variant calls
Putative somatic SNPs and short indels were called with GATK Mutect2 (https://software.broadinstitute.org/gatk). Mutect2 first identifies candidate sites with evidence of variation and then performs local reassembly to refine variant calling. It employs an external reference, known as a “panel of normal samples,” to filter out recurrent sequencing artifacts and calls variants only at sites with evidence of somatic variation. The panel of normal samples comprised 100 randomly selected individuals under 40 years of age, with the absence of hotspot CHIP mutations confirmed prior to their inclusion. 

Variants included in the preliminary dataset met the following criteria: presence in a pre-established list of candidate CHIP variants, total sequencing depth ≥ 20, alternate allele read depth count ≥ 5, and representation in both sequencing directions (i.e., F1R2 ≥ 1 and F2R1 ≥ 1). CHIP mutations were defined as those with a variant allele fraction (VAF) ≥ 0.02.


## Calculation for genetically predicted TREM1 levels
We calculated genetically predicted scores of IL-17RA expression using 76 variants from a model trained with Somalogic protein measurements from the INTERVAL study (OPGS000019)

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


## Association between gIL17RA and CVD
We tested the association between gIL-17RA and CVD prevalence with logistic regression and the association between JAK2-CHIP and CVD incidence with time-to-event regression. We used LASSO regression to identify gIL-17RA variants that most significantly abrogate the association between JAK2-CHIP and CVD. We tested the association between gIL-17RA and CVD prevalence with logistic regression. Given the relatively small number of events within each subtype, we applied Firth’s penalized logistic regression (via the logistf package in R). We adjusted for age, age-squared, genetic sex, smoking status, genetic ancestry, body mass index, low-density lipoprotein cholesterol level, and baseline status of diabetes, hypertension, and autoimmune diseases (composite of rheumatoid arthritis, psoriasis, psoriatic arthritis, ankylosing spondylitis, and systemic lupus erythematosus). To ensure numerical stability, we increased the maximum number of iterations for model convergence (maxit = 1000). As a sensitivity analysis, the patients with autoimmune diseases were removed.

Code available in: [Association code](https://github.com/bicklab/il17ra/blob/main/KZ_IL17RA_CVD_UKB.ipynb)

### gIL17RA and CVD among JAK2
We stratified JAK2-CHIP carriers into three groups of equal size based on variant allele fraction (VAF) and assessed the association between gIL-17RA and CVD within each group. The low VAF group (bottom 33 percentile) were defined as VAF < 0.27, the intermediate VAF group (middle 33 percentile) were defined as 0.27 ≤ VAF < 0.47, and the high VAF group (top 33 percentile) defined as VAF ≥ 0.47. 

<img width="370" alt="image" src="https://github.com/user-attachments/assets/efd93d5c-5fb5-4152-8906-8f7ca1f0cc85" />

The association between genetically predicted IL-17RA levels and cardiovascular disease in different participant groups in the UK Biobank stratified by VAF.

Code available in: [Association code](https://github.com/bicklab/il17ra/blob/main/KZ_IL17RA_CVD_UKB.ipynb)

### gIL-17RA modifies the hematologic manifestations of JAK2-CHIP
We analyzed five available hematologic parameters at baseline in the UK Biobank: white blood cell count (WBC), hemoglobin (HGB), mean corpuscular volume (MCV), platelet count (PLT), and red cell distribution width (RDW). As expected, we first confirmed strong positive associations between JAK2-CHIP and each of these blood traits (all FDR-adjusted P < 0.01), consistent with prior literature. 

<img width="582" alt="Screenshot 2025-07-02 at 14 36 09" src="https://github.com/user-attachments/assets/4923d7c9-6956-4874-b5cb-a6c88d2fd38f" />

Next, we tested for statistical interaction between JAK2-CHIP and gIL-17RA on each blood trait. Among the five parameters, we observed a significant interaction only for RDW (interaction β = -0.27, 95% CI: -0.46 to -0.08, FDR-adjusted P = 0.022), suggesting that higher gIL-17RA may attenuate the JAK2-CHIP–associated increase in RDW. This result may reflect a dampened systemic inflammatory state, given RDW’s established association with chronic inflammation and cardiovascular risk like coronary artery disease (Zalawadiya et al. Am J Cardiol. 2010), acute myocardial infarction (Wang et al. J Atheroscler Thromb. 2015), and heart failure (Borne et al. Eur J Heart Fail. 2011).

<img width="576" alt="Screenshot 2025-07-02 at 14 37 07" src="https://github.com/user-attachments/assets/d140f989-2d19-497e-b4d4-50994196742f" />

Code available in: [Association code](https://github.com/bicklab/il17ra/blob/main/KZ_IL17RA_CVD_UKB.ipynb)

## Kaplan–Meier survival curve
Kaplan–Meier survival curves were generated across tertiles (via the survminer package in R), where the red line represents the lowest 33 percentile of gIL-17RA score (Low), the blue line represents the middle 33 percentile of gIL-17RA score (Intermediate), and the green line represents highest 33 percentile of gIL-17RA score (High).

Code available in: [Kaplan–Meier survival code](https://github.com/bicklab/il17ra/blob/main/KZ_KMplot.R)


## Single-cell RNA sequencing (scRNAseq) analysis in JAK2 V617F mice with and without Ninjurin-1 deficiency
We generated atheroprone JAK2-V617F mice with (N=17) and without (N=18) hematopoietic-specific Ninjurin-1 (NINJ1) deficiency by transplanting bone marrow from Mx1-Cre JAK2-V617F donors (with or without NINJ1-knockout) into female LDL-receptor-deficient recipients. Chimeric female LDL-receptor-deficient mice with 80% wildtype and 20% Mx1-Cre JAK2-V617F (with or without NINJ1 knockout) bone marrow were allowed to recover for 4 weeks after bone marrow transplant. After recovery, mice received 2 intraperitoneal (i.p.) injections with 200 μg/mouse polyinosinic:polycytidylic acid (pIpC) to induce JAK2-V617F and knockout NINJ1. Mice were placed on a Western-type diet (WTD) for 12 weeks to promote atherosclerotic lesion formation. Atherosclerotic lesion characteristics, such as lesion size (micrometer-squared) and percentage necrotic core (%) were quantified in aortic root sections.  (scRNAseq) was conducted on a pool of 5 digested aortas per group. The adventitial layer of the aorta was removed prior to digestion, and the cells were then flow-sorted for live CD45+ cells. After sequencing and quality control, a total of 1,020 intimal-medial immune cells were identified and followed by clustering and cell-type annotation based on existing literature, enabling the quantification of cell type abundance of (1) Trem2-expressing foamy-like macrophages, (2) Th17-like T cells, and (3) neutrophils relative to lesional immune cells. Cell type annotation gene expression values are available [here](https://github.com/bicklab/il17ra/blob/main/il7ra_jak2_scrnaseq_supplement.xlsx). Atherosclerotic lesion characteristics were quantified in aortic root sections. scRNAseq code available [here](https://github.com/bicklab/il17ra/blob/main/The%20atheroprotective%20role%20of%20interleukin-17-receptor-A%20signaling%20in%20JAK2%20clonal%20hematopoiesis.Rmd).


Cell-cell communication with LIANA and Tensor-cell2cell evaluated signaling interactions between IL-17A-producing Th17 cells and myeloid cells, with code available [here](https://github.com/bicklab/il17ra/blob/main/EHJ_CCC.py).

## Data
This analysis was performed on the [UK Biobank DNA Nexus Research Analysis Platform](https://ukbiobank.dnanexus.com) and BioVU Terra.bio environment.

## Acknowledgements
Individual-level sequence data and CHIP calls have been deposited with UK Biobank and are available to approved researchers by application (https://www.ukbiobank.ac.uk/register-apply/). Vanderbilt BioVU data are available through an application to the Vanderbilt Institute for Clinical and Translational Research (VICTR) BioVU Review Committee.

## Contact
Kun Zhao, kun.zhao@vumc.org; Yash Pershad, yash.pershad@vanderbilt.edu
