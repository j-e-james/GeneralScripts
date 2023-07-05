import sys, argparse, gzip, re 
from collections import defaultdict

"""Define a few things: Reference genome, Results file, GFF file"""
Ref = '/Users/jennyjames/Desktop/AssignAncestral/TAIR10_chr_all.fas.gz'
Res = '/Users/jennyjames/Desktop/AssignAncestral/TAIR10_GFF3_genes_degeneracy_annotation.txt'
GFF = '/Users/jennyjames/Desktop/AssignAncestral/TAIR10_GFF3_genes.gff.txt.gz'


codon_table = '''
#These are DNA codons
#CODON BASE_1 BASE_2 BASE_3
#Ala
GCT 0fold 0fold 4fold
GCC 0fold 0fold 4fold
GCA 0fold 0fold 4fold
GCG 0fold 0fold 4fold
#Arg
#the first base here is ambiguous, it can be 2 fold (NGA or NGG) or it can be 0fold (CGT and CGC)
CGT 0fold 0fold 4fold
CGC 0fold 0fold 4fold
CGA 2fold 0fold 4fold
CGG 2fold 0fold 4fold
AGA 2fold 0fold 2fold
AGG 2fold 0fold 2fold
#Asn
AAT 0fold 0fold 2fold
AAC 0fold 0fold 2fold
#Asp
GAT 0fold 0fold 2fold
GAC 0fold 0fold 2fold
#Cys
TGT 0fold 0fold 2fold
TGC 0fold 0fold 2fold
#Gin
CAA 0fold 0fold 2fold
CAG 0fold 0fold 2fold
#Glu
GAA 0fold 0fold 2fold
GAG 0fold 0fold 2fold
#Gly
GGT 0fold 0fold 4fold
GGC 0fold 0fold 4fold
GGA 0fold 0fold 4fold
GGG 0fold 0fold 4fold
#His
CAT 0fold 0fold 2fold
CAC 0fold 0fold 2fold
#Ile
ATT 0fold 0fold 3fold
ATC 0fold 0fold 3fold
ATA 0fold 0fold 3fold
#Met - Start
ATG 0fold 0fold 0fold
#Leu
#this first one is ambiguous 2fold (for NTA and NTG) and 0fold (for CTT and CTC)
TTA 2fold 0fold 2fold
TTG 2fold 0fold 2fold
CTT 0fold 0fold 4fold
CTC 0fold 0fold 4fold
CTA 2fold 0fold 4fold
CTG 2fold 0fold 4fold
#Lys
AAA 0fold 0fold 2fold
AAG 0fold 0fold 2fold
#Phe
TTT 0fold 0fold 2fold
TTC 0fold 0fold 2fold
#Pro
CCT 0fold 0fold 4fold
CCC 0fold 0fold 4fold
CCA 0fold 0fold 4fold
CCG 0fold 0fold 4fold
#Ser
#These are all 0fold because the two codon types (TCN and AGN) do not share bases 
TCT 0fold 0fold 4fold 
TCC 0fold 0fold 4fold 
TCA 0fold 0fold 4fold 
TCG 0fold 0fold 4fold 
AGT 0fold 0fold 2fold 
AGC 0fold 0fold 2fold 
#Thr
ACT 0fold 0fold 4fold
ACC 0fold 0fold 4fold
ACA 0fold 0fold 4fold
ACG 0fold 0fold 4fold
#Trp
TGG 0fold 0fold 0fold
#Tyr
TAT 0fold 0fold 2fold
TAC 0fold 0fold 2fold
#Val
GTT 0fold 0fold 4fold
GTC 0fold 0fold 4fold
GTA 0fold 0fold 4fold
GTG 0fold 0fold 4fold
#STOP
#I will leave these as stop annotations for now, although technically they do have a fold, 
#because their occurence helps delineate sites
TAA stop stop stop
TGA stop stop stop
TAG stop stop stop 
'''

def degeneracy(codon):
	normal_codon = False
	for f in codon_table.split("\n"):
		if '#' not in f:
			if codon in f:
				normal_codon = True
				return(f[4:].split())
	if not normal_codon:
		return(['unknown', 'unknown', 'unknown'])	
		
		

def CDS_Ref_Regions(Ref_scaf_name, CDS_list, Ref_scaf_seq):
	"""identifies regions of the scaffold that correspond to the CDS_list"""
	CDS = [x for x in CDS_list if Ref_scaf_name in x]
	CDS.sort(key=lambda x: x[2])	
	CDS_dir = defaultdict(list)
	for x in CDS:
		### 0 indexed numbers, hence the -1 for the start positions.
		CDS_dir[x[-1]].append((x[2], x[3], x[4], Ref_scaf_seq[x[2]-1:x[3]]))
	return(CDS_dir)


def Reverse_compliment(Full_CDS):
	"""Returns reverse compliment of a sequence"""
	Rev_Full_CDS = Full_CDS[::-1]
	Complement_Full_CDS = []
	for base in Rev_Full_CDS:
		if base == 'A':
			Complement_Full_CDS.append('T')
		elif base == 'T':
			Complement_Full_CDS.append('A')
		elif base == 'C':
			Complement_Full_CDS.append('G')
		elif base == 'G':
			Complement_Full_CDS.append('C')
		else:
			Complement_Full_CDS.append('N')
	return(''.join(Complement_Full_CDS))


def degenerate_annotator(CDS_dir):
	"""Returns sites, by scaffold, position, and parent gene, annotated by their degeneracy"""
	per_scaffold_degeneracy = []
	for k, v in CDS_dir.items():
		direction = v[0][2]	
		
		"""processing the CDS_dir a little more- the list v"""
		###remove duplicates
		set_v = set(v)
		###but now turn into a list and reorder by exon starting position
		v = list(set_v)
		v.sort(key=lambda x: int(x[0]))	
		
		Full_CDS = ''.join([x[-1] for x in v])

		if direction == '-':
			Full_CDS = Reverse_compliment(Full_CDS)
			
		if len(Full_CDS) % 3 != 0:
			print('CDS incomplete! There is a problem in the GFF for this gene, it will be excluded from annotations:')
			print(k)
			print([x[0] for x in v])
			print(v)
			print(Full_CDS)		
		
		else:								
			Full_CDS_codons = re.findall('...?',Full_CDS)
			degenerate_list = [degeneracy(codon) for codon in Full_CDS_codons]
			degenerate_list = [item for sublist in degenerate_list for item in sublist]

			if direction == '-':	
				degenerate_list.reverse()

			positions = [[x for x in range(int(exon[0]),int(exon[1])+1)] for exon in v]
			positions = [item for sublist in positions for item in sublist]
		
			for codon_degeneracy, pos in zip(degenerate_list, positions):				
				per_scaffold_degeneracy.append(Ref_scaf+' '+str(pos)+' '+codon_degeneracy+' '+k)
	return(per_scaffold_degeneracy)


CDS_list = []
"""New GFF3 parser- focus entirely on coding regions, pulling out only CDS"""
with gzip.open(GFF, 'rt') as GFF3_file:
	for line in GFF3_file:
		if line.startswith("#"):
			continue
		elif 'CDS' in line:
			try:
				sline = line.split()
				scaf, type, start, end, dir, Parent_info  = sline[0], sline[2], int(sline[3]), int(sline[4]), sline[6], sline[8]
				###Check the use of this: should return the identity of the gene.
				Parent_info = Parent_info.split(";")[0].replace("ID=","").split(",")[0]
				CDS_list.append((scaf, type, start, end, dir, Parent_info))
			except (ValueError, IndexError):
				print("Line conversion failed. Skipping %s.\n" % line)
				continue  

""""Programme begins- main loop, reading in files"""			
seq = []
Ref_scaf = ''
with gzip.open(Ref, 'rt') as Ref_file, open(Res, 'a') as resfile:
	for line in Ref_file:	
	
		if line.startswith(">"):
		
			### then we have hit a new chromosome.
			### assumes the header is of the form: #>CHROM_NAME other stuff by white space
			Ref_scaf_name = line[1:].split()[0]

			if len(seq) != 0:
				print(Ref_scaf)
				full_seq = ''.join(seq)
				CDS_dir = CDS_Ref_Regions(Ref_scaf, CDS_list, full_seq)	
				
				"""We will focus on the '.1' verison of the protein, ignoring splice variants for now, so that results
				only contain one degeneracy annotation per site"""
				main_splice_var = []
				for k, v in CDS_dir.items():
					if k[-2:] == '.1':
						main_splice_var.append(k)
				CDS_dir_main_splice_var = {key: CDS_dir[key] for key in main_splice_var}
							
				###process the dictionary, annotating sites by degeneracy
				per_scaffold_degeneracy = degenerate_annotator(CDS_dir_main_splice_var)					
				resfile.write('\n'.join(per_scaffold_degeneracy))

				seq = []								
				"""for testing, to get a single chromosome:"""
 				#break
				
		else:
			Ref_scaf = Ref_scaf_name
			seq.append(line.rstrip().upper())
	
	"""and pick up our trailing scaffold"""	
	print(Ref_scaf)
	full_seq = ''.join(seq)
	CDS_dir = CDS_Ref_Regions(Ref_scaf, CDS_list, full_seq)
	main_splice_var = []
	for k, v in CDS_dir.items():
		if k[-2:] == '.1':
			main_splice_var.append(k)
	CDS_dir_main_splice_var = {key: CDS_dir[key] for key in main_splice_var}
				
	###process the dictionary, annotating sites by degeneracy
	per_scaffold_degeneracy = degenerate_annotator(CDS_dir_main_splice_var)					
	resfile.write('\n'.join(per_scaffold_degeneracy))						
						
					