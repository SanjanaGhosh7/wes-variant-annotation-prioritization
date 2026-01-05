# 🧬 ***Whole-Exome Variant Prioritization and Annotation (chr22, Healthy Sample)***



## 

## 📌 Project Overview



This project demonstrates an end-to-end variant annotation, prioritization, and interpretation workflow using whole-exome sequencing (WES) data from a healthy reference individual (NA12878), focusing on chromosome 22.



***The primary goal is not disease discovery, but rather to:***

-Build a technically sound and biologically interpretable variant prioritization pipeline

-Understand baseline human genetic variation

-Establish a reference framework that can later be extended to case–control, disease, or multi-omics studies



This mirrors how real-world bioinformatics pipelines are validated on healthy or reference samples before being applied to clinical or disease cohorts.



## 

## 📚 Dataset Description



Data type: ***Whole Exome Sequencing (WES)***



Sample: ***Healthy sample (NA12878)***



Reference genome: ***Human (GRCh38)***



Scope: ***Chromosome 22 only (selected as a reduced yet biologically relevant test case)***



⚠️ Raw FASTQ files, reference genome data, and other intermediate files are intentionally excluded from this repository due to size constraints. The pipeline is fully reproducible using the scripts and instructions provided.



## 

## 🧪 Biological Context



Whole-Exome Sequencing enables the detection of coding and splice-region variants that may affect protein function. WES data lets us find the genetic causes of rare diseases, identify cancer mutations for targeted treatments, understand complex disorders, and predict drug responses by focusing on the protein-coding parts (exons) of the genome, which harbor most disease-causing mutations. It helps diagnose developmental delays, epilepsy, heart issues, and other inherited conditions, improving patient outcomes and guiding personalized medicine.



###### Why a Healthy Sample is used for here?

Using a healthy individual is intentional and scientifically meaningful:

-A healthy genome still contains thousands of variants, including:

    -Rare variants

    -Missense changes

    -Even occasional predicted loss-of-function variants



###### This allows us to:

    -Validate annotation and filtering logic

    -Study background genetic variation

    -Avoid bias from disease-driven assumptions



It provides a baseline expectation against which disease cohorts can later be compared. Importantly, the presence of HIGH or MODERATE impact variants in healthy individuals is biologically normal, reinforcing why variant impact ≠ pathogenicity.



## 

## 🔬 Pipeline Overview



Raw WES FASTQ files
        │
        ▼
Quality Control (FastQC)
        │
        ▼
Alignment to Reference Genome (GRCh38)
        │
        ▼
Variant Calling (VCF generation)
        │
        ▼
Variant Annotation (Ensembl VEP)
        │
        ▼
Variant Prioritization
  ├─ Canonical transcripts
  ├─ Protein-coding variants
  ├─ HIGH / MODERATE impact
  └─ Rare alleles (gnomAD AF < 0.01)
        │
        ▼
Exploratory Analysis & Visualization
  (Impact vs AF, ClinVar overlay, summaries)



## 

## 🗂 Project Structure



wes_chr22_project/
│
├── scripts/
│   ├── wes_variant_calling_pipeline.sh   # Alignment, post-processing & variant calling
│   ├── vep_annotation.sh                 # Ensembl VEP-based variant annotation
│   └── wes_variant_analysis.ipynb        # Downstream analysis & visualization
│
├── results/
│   ├── prioritized_variants.tsv          # Final filtered & annotated variant table
│   └── *.png                             # Generated plots & figures
│
├── README.md                             # Project documentation
└── .gitignore                            # Excluded intermediate & system files



## 

## ⚙️ Methodology





##### 1\. Variant Calling(GATK Best Practices) and Initial Filtering



**Variant calling was performed using standard best-practice steps:**



   -Alignment to GRCh38



   -Duplicate marking



   -Base Quality Score Recalibration (BQSR)



   -Variant calling



   -Hard filtering to retain only PASS variants



**Hard filtering ensures:**



   -Removal of low-confidence calls



   -High technical reliability for downstream interpretation



This step was performed using bash scripts (`/scripts/wes_variant_calling_pipeline.sh`) and standard command-line tools, without Conda, as system-installed tools were sufficient and stable.







##### 2\. Variant Annotation with Ensembl VEP



**Variants were annotated using Ensembl Variant Effect Predictor (VEP):**



   -Reference genome: GRCh38



   -Offline cache



   -Plugins and population frequency annotations enabled

&nbsp; 

&nbsp;  -Produced vep.vcf and vep.summary files



**Conda was used for VEP. Why?**



VEP has:



   -Complex Perl dependencies



   -Plugin compatibility constraints



   -Cache version coupling



Using Conda ensured:



  -Reproducibility



  -Dependency isolation



  -Clean environment management


This step was performed using bash script (`/scripts/vep_annotation.sh`).




##### 3\. To get Prioritized Variants



***Key Insight* 🔑**



VEP stores detailed annotations inside the CSQ field, which must be structured before filtering. All these steps were performed in CLI.



###### 

###### Step 1: Split CSQ into Structured Fields



**Command(bash):**



*bcftools +split-vep NA12878.chr22.vep.vcf \\*

*-c Allele,Consequence,IMPACT,SYMBOL,Gene,Feature,BIOTYPE,CANONICAL,EXON,INTRON,HGVSc,HGVSp,gnomADe\_AF,SIFT,PolyPhen,CLIN\_SIG \\*

*-o vep.split.vcf*



What this does:



  -Extracts CSQ subfields into usable INFO tags



  -Enables biologically meaningful filtering







###### Step 2: Impact + Canonical Transcript Filtering



**Command(bash):**



*bcftools view -i \\*

*'INFO/CANONICAL="YES" \&\& INFO/BIOTYPE="protein\_coding" \&\& (INFO/IMPACT="HIGH" || INFO/IMPACT="MODERATE")' \\*

*vep.split.vcf -o vep.prioritized.vcf*



Why this step?



  -Focuses on interpretable protein-coding effects



  -Retains canonical transcript consequences



  -Removes regulatory and low-confidence noise



*⚠️ **Note:** Variants may still carry multiple transcript annotations.*

*Final impact prioritization is handled downstream.*







###### Step 3: Population Frequency Filtering



***We apply an allele frequency threshold to enrich for rare variants, because common variants are unlikely to have strong functional or clinical impact. An AF cutoff of 1% (0.01) is a widely used research standard to separate common population polymorphisms from potentially interesting variants.***



**Command(bash):**



*bcftools view -i \\*

*'INFO/gnomADe\_AF<0.01 || INFO/gnomADe\_AF="."' \\*

*vep.prioritized.vcf -o vep.prioritized.rare.vcf*



**Why is AF filtering?**



Without it, our dataset contains:



  -Millions of background polymorphisms



  -Variants that everyone has



  -Noise that obscures meaningful biology



Why allow "." (missing AF)?



  -Because Missing frequency ≠ rare nor common



  -Prevents premature exclusion of novel variants







###### Step 4: Create the Prioritized Variant Table



**Command(bash):**



*bcftools query -f \\*

*'%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/SYMBOL\\t%INFO/Consequence\\t%INFO/IMPACT\\t%INFO/HGVSc\\t%INFO/HGVSp\\t%INFO/gnomADe\_AF\\t%INFO/SIFT\\t%INFO/PolyPhen\\t%INFO/CLIN\_SIG\\n' \\*

*vep.prioritized.rare.vcf > prioritized\_variants.tsv*



This table serves as the input for all downstream analysis.







## 📊 Downstream Analysis (Jupyter Notebook)



The notebook *wes\_variant\_analysis.ipynb* performs exploratory analysis using the prioritized variant table:







###### Analyses include:



* Variant impact distribution



* Variant burden per gene (top 10)



* SIFT vs PolyPhen functional predictions



* Allele frequency spectrum (gnomAD)



* Allele frequency vs impact assessment



* ClinVar significance overlay on impact vs allele frequency (Extracted the most severe impact for each variant)



* Gene-level deep dive (e.g., PPP6R2)





###### 

###### Key Visualization:



**-ClinVar significance overlaid on Impact vs Allele Frequency plot:** This visualization highlights how most ClinVar-annotated variants cluster around Moderate-impact, with almost no rare variants in High-impact, which is expected for a healthy sample and demonstrates realistic biological patterns.



*⚠️ **Note:** Refer to the Notebook for clear understanding!*







* **All plots generated by the notebook are saved to:** `/results*`*



Also, using Conda for the Notebook analysis ensured:



  -Matching Python, pandas, seaborn versions



  -Reproducibility



  -Clean separation from system Python



## 

## 🛠 Technologies \& Tools Used



* Linux / Bash scripting



* BWA / SAMtools (alignment \& processing)



* Variant calling tools



* Conda



* Ensembl VEP



* Python (pandas, matplotlib, seaborn)



* JupyterLab



* Git \& GitHub





## 🧠 Key Skills Demonstrated



* Linux \& shell scripting



* Variant calling workflows



* Variant annotation using VEP



* Interpretation of functional consequences



* Variant prioritization strategies



* Jupyter-based exploratory data analysis



* Reproducible project structuring

## 

## 🚫 N.B.





&nbsp; -This project is research-focused, not clinical



&nbsp; -No diagnostic claims are made

  

&nbsp; -No phenotype-based interpretation



&nbsp; -ACMG classification was intentionally excluded



  -No multi-sample cohort analysis



These were intentionally excluded to keep the project focused on variant annotation and prioritization fundamentals.





## 

## 📝 Project Conclusion



This project successfully demonstrates a research-grade variant prioritization workflow starting from WES data of a healthy individual. Despite the absence of disease, biologically meaningful variants—including rare, moderate, and high-impact changes were identified and interpreted in proper population and clinical context.



**The results highlight:**



* The prevalence of functional variation in healthy genomes



* The importance of population frequency and transcript-aware prioritization



* Why annotation alone is insufficient without downstream interpretation



* This pipeline forms a strong foundation for disease-focused or comparative genomics studies.







## 🚀 Future Extensions



1. ACMG/AMP variant classification
   
2. Pathway Enrichment Analysis
   
3. Multi-sample or cancer cohort analysis
   
4. Integration with COSMIC or TCGA
   
5. Machine learning–based prioritization















\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

## &nbsp;                       

## &nbsp;			*🚀 How to Run This Project*



##### 🖥️ System Requirements



Operating system:



Linux (Ubuntu ≥ 20.04) or WSL2 on Windows (used in this project)



CPU: ≥ 4 cores recommended



RAM: ≥ 16 GB (8 GB minimum for chr22-only analysis)



Disk space: ≥ 150 GB (reference + cache + intermediate files)





##### 📦 Software, Tools \& Dependencies



1️⃣ Core Bioinformatics Tools (Command Line)



These tools are required before running the pipeline:





| Tool       | Purpose                                  | Version (recommended) |

| ---------- | ---------------------------------------- | --------------------- |

| `bwa`      | Read alignment                           | ≥ 0.7.17              |

| `samtools` | BAM processing                           | ≥ 1.17                |

| `gatk`     | Duplicate marking, BQSR, variant calling | ≥ 4.5                 |

| `bcftools` | VCF manipulation \& filtering             | ≥ 1.18                |

| `fastqc`   | Raw read QC                              | ≥ 0.12                |

| `picard`   | (Used via GATK)                          | bundled               |

| `java`     | Required for GATK                        | OpenJDK 17            |





2️⃣ Reference \& Resource Files (GRCh38) 



Download and prepare the following (GATK Resource Bundle): 



| File                   | Purpose               |

| ---------------------- | --------------------- |

| `GRCh38 FASTA`         | Reference genome      |

| `FAI + DICT`           | Required for GATK     |

| `dbSNP VCF`            | Known variants (BQSR) |

| `Mills \& 1000G indels` | Known indels (BQSR)   |



These are placed in `/references`





3️⃣ Conda Environment (Used Only Where Needed)



Conda is not used for the entire pipeline, by design.



Create Conda Environment for Annotation \& Visualization and install dependencies: ensembl-vep, pandas, numpy, matplotlib, seaborn, and jupyterlab.





4️⃣ Ensembl VEP Cache (Offline Annotation)



Required for fast, reproducible annotation.







##### 🧾 Folder Structure in your System



wes\_chr22\_project/

├── fastq/           # Raw FASTQs

├── references/      # Genome \& known sites

├── logs/            # Tool logs

├── scripts/         # Bash scripts \& notebook analysis

├── results/         # VCFs, TSVs, plots

├── tmp/             # Temporary files







##### 📌 Reproducibility Note



System tools → installed via apt



Annotation \& visualization → isolated via conda



Offline VEP cache ensures deterministic results



Each step produces interpretable intermediate outputs.





\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*







***📌 Notes***



***This project is intended for educational and portfolio purposes and reflects realistic bioinformatics analysis practices used in research and industry settings.***







###### **👩‍🔬 Author: Sanjana Ghosh**

###### **📍 Field: Bioinformatics | Genomics | Variant Analysis \& Interpretation**

