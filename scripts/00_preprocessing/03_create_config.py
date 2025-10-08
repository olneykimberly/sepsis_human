#!/usr/bin/python3

# create a new output file
outfile = open('config.json', 'w')

# get all sample names
allSamples = list()
maleSamples = list()  
femaleSamples = list()
numSamples = 0

with open('../../metadata/metadata_merged.tsv', 'r') as infile:
    for line in infile:
        numSamples += 1

        line = line.replace("\t", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # Brain_sepsis_M_1_B1L1_Fryer-1_384_193_S46_L001_R1_001.fastq.gz_@A00124:303:HC2WMDSX2:1:1101:1271:1000 1:N:0:GGTACCGACC+CAGATACCAC
        # create a shorter sample name
        stemName = sampleAttributes[0] + '_' + sampleAttributes[1] + '_' + sampleAttributes[2] + '_' + sampleAttributes[3] # Brain_sepsis_M_1
        allSamples.append(stemName)
        
        # Determine sex and append to the corresponding list
        sex = sampleAttributes[2]
        if sex == 'M':
            maleSamples.append(stemName)
        elif sex == 'F':
            femaleSamples.append(stemName)

# create header and write to outfile
header = '''{{
    "Commment_Input_Output_Directories": "This section specifies the input and output directories for scripts",
    "fastq" : "/tgen_labs/jfryer/kolney/sepsis_human/fastq/",
    "rawReads" : "/tgen_labs/jfryer/kolney/sepsis_human/fastq/",
    "rawQC" : "../../rawQC/",
    "trimmedReads" : "../../trimmedReads/",
    "trimmedQC" : "../../trimmedQC/",
    "starAligned" : "../../starAligned/",
    "starAligned_SCC" : "../../startAligned_SCC/",
    "bamstats" : "../../bamstats/",

    "Comment_Reference" : "This section specifies the location of the human, Genocode reference genome",
    "ref_fa" : "/tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38.primary_assembly.genome.fa",
    "ref_gtf" : "/tgen_labs/jfryer/projects/references/human/GRCh38/gencode.v38.annotation.sorted.gtf",
    
    "GRCh38.fa" : "/tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38.primary_assembly.genome.fa",
    "GRCh38.star" : "/tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_def",
    "GRCh38.gtf" : "/tgen_labs/jfryer/projects/references/human/GRCh38/gencode.v38.annotation",

    "GRCh38.Ymasked.fa" : "/tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_Ymasked_XX.fa",
    "GRCh38.YPARs_masked.fa" : "/tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_YPARsmasked_XY.fa",

    "GRCh38.Ymasked.star" : "/tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_Ymasked_XX",
    "GRCh38.YPARs_masked.star" : "/tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_YPARsmasked_XY",

    "Comment_Sample_Info": "The following section lists the samples that are to be analyzed",
    "sample_names": {0},
    "male_names": {1},
    "female_names": {2},
'''
outfile.write(header.format(allSamples, maleSamples, femaleSamples))
# config formatting
counter = 0
with open('../../metadata/metadata_merged.tsv', 'r') as infile:
    for line in infile:
        counter += 1

        # make naming consistent, we will rename using only underscores (no hyphens)
        line = line.replace("\t", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # project_uniqueNum_1_tissue_group_XX_XX_sequencer_adapter_lane_read_001 
                          # Brain_sepsis_M_1_B1L1_Fryer-1_384_193_S46_L001_R1_001.fastq.gz_@A00124:303:HC2WMDSX2:1:1101:1271:1000 1:N:0:GGTACCGACC+CAGATACCAC

        base = '/tgen_labs/jfryer/kolney/sepsis_human/fastq/' + sampleAttributes[5] + '_' + sampleAttributes[6] + '_' + sampleAttributes[7] +  '_' + sampleAttributes[8] + '_' + sampleAttributes[9] + '_' + sampleAttributes[10] +'_' + sampleAttributes[11] 
        sampleName1 = base
        sampleName2 = sampleName1.replace("_R1_", "_R2_")
        sampleInfo = split[0]

        # create a shorter sample name
        stemName = sampleAttributes[0] + '_' + sampleAttributes[1] + '_' + sampleAttributes[2] + '_' + sampleAttributes[3] # Brain_sepsis_M_1

        # break down fastq file info
        # @A00127:312:HVNLJDSXY:2:1101:2211:1000
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
        sampleInfo = sampleInfo.split(':')
        instrument = sampleInfo[0]
        runNumber = sampleInfo[1]
        flowcellID = sampleInfo[2]

        lane = sampleAttributes[9]
        ID = stemName  # ID tag identifies which read group each read belongs to, so each read group's ID must be unique
        SM = sampleAttributes[3] # Brain_sepsis_M_1
        PU = lane  # Platform Unit
        LB = stemName
        SEX = sampleAttributes[2]
        TYPE = sampleAttributes[1]
        TISSUE = sampleAttributes[0]

        out = '''
    "{0}":{{
        "fq1": "{1}",
        "fq2": "{2}",
        "ID": "{3}",
        "SM": "{4}",
        "PU": "{5}",
        "LB": "{6}",
        "PL": "illumina",
        "sex": "{7}", 
        "group" : "{8}", 
        "tissue" : "{9}"
        '''
        outfile.write(out.format(stemName, sampleName1, sampleName2, ID, SM, PU, LB, SEX, TYPE, TISSUE))
        if (counter == numSamples):
            outfile.write("}\n}")
        else:
            outfile.write("},\n")
outfile.close()
