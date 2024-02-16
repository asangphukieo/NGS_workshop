
# Install pycirclize in python
python3 -m pip install pycirclize

wget https://raw.githubusercontent.com/asangphukieo/NGS_workshop_2024/main/Day5/Input_Circos/01_gatk_somatic_pass.vcf
wget https://raw.githubusercontent.com/asangphukieo/NGS_workshop_2024/main/Day5/Input_Circos/04_gatk_somatic_pass.vcf

# Import required packages
from pycirclize import Circos
from pycirclize.utils import ColorCycler
from pycirclize.utils import load_eukaryote_example_dataset
import numpy as np
import pandas as pd
import gzip
np.random.seed(0)

def get_vcf_names(vcf_path):
    with open(vcf_path, "rt") as ifile:
          for line in ifile:
            if line.startswith("#CHROM"):
                  vcf_names = [x.replace("\n","").replace("#","") for x in line.split('\t')]
                  break
    ifile.close()
    return vcf_names


# 3.2 import SNV file and create dictionary
names = get_vcf_names('01_gatk_somatic_pass.vcf')
vcf = pd.read_csv('01_gatk_somatic_pass.vcf',  comment='#', delim_whitespace=True, header=None, names=names)

names2 = get_vcf_names('04_gatk_somatic_pass.vcf')
vcf2 = pd.read_csv('04_gatk_somatic_pass.vcf',  comment='#', delim_whitespace=True, header=None, names=names2)


snv_dict={}
for i in range(1,23):
  #print(i)
  if "chr"+str(i) in set(vcf.CHROM): ## edit
  	snv_dict["chr"+str(i)] = list(vcf[vcf.CHROM == "chr"+str(i)].POS)

snv_dict["chrX"] = list(vcf[vcf.CHROM == "chrX"].POS)

snv_dict2={}
for i in range(1,23):
  #print(i)
  if "chr"+str(i) in set(vcf2.CHROM): ## edit
    snv_dict2["chr"+str(i)] = list(vcf2[vcf2.CHROM == "chr"+str(i)].POS)
snv_dict2["chrX"] = list(vcf2[vcf2.CHROM == "chrX"].POS)


#####
chr_bed_file, cytoband_file, chr_links = load_eukaryote_example_dataset("hg38") ## edit

# Initialize Circos from BED chromosomes
circos = Circos.initialize_from_bed(chr_bed_file, space=3)
circos.text("Homo sapiens\n(hg38)", deg=315, r=150, size=12)

# Add cytoband tracks from cytoband file
circos.add_cytoband_tracks((95, 100), cytoband_file)

# Create chromosome color dict
ColorCycler.set_cmap("hsv")
chr_names = [s.name for s in circos.sectors]
colors = ColorCycler.get_color_list(len(chr_names))
chr_name2color = {name: color for name, color in zip(chr_names, colors)}


###############################
for sector in circos.sectors:
    sector.text(sector.name, r=120, size=10, color=chr_name2color[sector.name])
    sector.get_track("cytoband").xticks_by_interval(40000000,
        label_size=8,
        label_orientation="vertical",
        label_formatter=lambda v: f"{v / 1000000:.0f} Mb")
    if sector.name in snv_dict:
      x = snv_dict[sector.name]
      y = np.random.randint(0, 100, size=len(x))
      # Scatter track
      track1 = sector.add_track((85, 95), r_pad_ratio=0.1)
      track1.axis()
      track1.scatter(x, y, color="red")
    if sector.name in snv_dict2:
      x = snv_dict2[sector.name]
      y = np.random.randint(0, 100, size=len(x))
      # Scatter track
      track2 = sector.add_track((70, 80), r_pad_ratio=0.1)
      track2.axis()
      track2.scatter(x, y, color="green")

fig = circos.plotfig()
fig.savefig("circos_ouput.png", dpi=100)
