import gzip, os, re, dadi, numpy
from collections import defaultdict
from collections import Counter


#### This programme is agnostic to whether data is phased or not- i.e. whether the symbol / or | should be expected for 'counts': 
#### check what is appropriate before starting


"""
SAM format- order of columns
Col	Field	Description
1	QNAME	Query (pair) NAME
2	FLAG	bitwise FLAG
3	RNAME	Reference sequence NAME
4	POS	1-based leftmost POSition/coordinate of clipped sequence
5	MAPQ	MAPping Quality (Phred-scaled)
6	CIAGR	extended CIGAR string
7	MRNM	Mate Reference sequence NaMe (‘=’ if same as RNAME)
8	MPOS	1-based Mate POSistion
9	ISIZE	Inferred insert SIZE
10	SEQ	query SEQuence on the same strand as the reference
11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
12	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE
"""

### setup to generate the line for est-sfs- order A, C, G, T
nt_order = ['A', 'C', 'G', 'T']
OGBase = '\t0,0,0,0\n'


### organise outgroup sam format files, to extract the vcf positions. 
### Return a dictionary of SAM files- keys are Reference name (Chrom/scaffold)	start:end,
### values are OG string 
def SAMfile_dict(SAMfile):
	OG1_dict = {}
	with (gzip.open if 'gz' in SAMfile else open)(SAMfile, 'rt') as OG1:
# 		OG1 = [x for x in OG1 if '@' not in x]
		count = 0
		for line in OG1:
			if '@' not in line:
				count = count + 1
				line = line.split('\t')
# 				print(line)
				try:
					Ref_name, Ref_start, Ref_end, Seq = line[2], line[3], str(int(line[3])+len(line[9]) - 1), line[9]
				except:
					break
				OG1_dict[Ref_name+'\t'+Ref_start+'\t'+Ref_end] = Seq	
	return(OG1_dict)

AssignAncestralDir = '/Users/jennyjames/Desktop/AssignAncestral/'
vcffile = '/Volumes/MyPassport/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.missense.gz'

outfile_name =  AssignAncestralDir+'Athal.missense'

###OG1 should be the closer outgroup, OG2 the further away
outgroup1file = AssignAncestralDir+'Alyrata_AthalRef_algn.sam.gz' 
outgroup2file = AssignAncestralDir+'Crubella_AthalRef_algn.sam.gz'

OG1_dict = SAMfile_dict(outgroup1file)
OG2_dict = SAMfile_dict(outgroup2file)
print(outgroup1file)
print(outgroup2file)
print('Loaded sam file information')

### search sam file dictionary for a match with the vcf location, and return the 
### outgroup nucleotide in Est-SFS format
def AncBase_ID(OG_dict, scaffold, pos, ref, alt):
	OG = [key for key in OG_dict.keys() if scaffold in key]	
	OG = [key for key in OG if int(key.split('\t')[1]) <= int(pos) < int(key.split('\t')[2])+1]

	AncBase = [0,0,0,0]	
	for OG_items in OG:
		scaf, start, stop = OG_items.split('\t')
# 		OG_range = [x for x in range(int(start), int(stop)+1)]
# 		if int(pos) in range(int(start), int(stop)+1):	
		OG_range = [x for x in range(int(start), int(stop)+1)]	
		vcf_position = OG_range.index(int(pos))	
		ancestral = OG_dict[OG_items][vcf_position].upper()
		### ignoring non-ACGT positions, e.g. Ns or Xs.
		if ancestral in 'ACGT':		
			ref_index = nt_order.index(ancestral)
			AncBase[ref_index] = 1			
		###if we hit the index, we don't need to keep iterating
		break
	return(AncBase)
	
	
	
Vcf_by_position = {}
with gzip.open(vcffile, 'rt') as vcffile:
# 	count = 0
	for line in vcffile:	
		if '#' not in line:
# 			count = count+1
# 			if count > 150000:

			splitline = line.split('\t')

			### scaffold names must match
			scaffold, pos, ref, alt = 'Chr'+splitline[0], splitline[1], splitline[3], splitline[4]
	
			### we need to ignore sites that are not bi-allelic
			if len(alt) == 1 and len(ref) == 1:	
	
				OG1_Anc = AncBase_ID(OG1_dict, scaffold, pos, ref, alt)
				OG2_Anc = AncBase_ID(OG2_dict, scaffold, pos, ref, alt)

# 				print(scaffold+' '+pos+' '+ref+' '+alt)
	
				PolBase = [0,0,0,0]
				PolBaseDown = [0,0,0,0]

				position = str(scaffold)+'\t'+str(pos)
				ref_index = nt_order.index(ref)
				alt_index = nt_order.index(alt)			

				### we don't care whether data is phased or not
				counts_strike = re.findall('[10]\|[10]', line)
				counts_dash = re.findall('[10]/[10]', line)
				counts = counts_strike + counts_dash
				totalcount = len(counts)*2
				altcount = ''.join(counts).count('1')
				refcount = totalcount-altcount			
				PolBase[alt_index] = altcount
				PolBase[ref_index] = refcount
	
				resultsstring = ','.join(str(x) for x in PolBase) + '\t'+ ','.join(str(x) for x in OG1_Anc) + '\t'+ ','.join(str(x) for x in OG2_Anc)
		
				if totalcount > 100:	
					### random resample down to 100- this might take mutations down to 0
					x = numpy.random.hypergeometric(altcount, refcount, 100)
					PolBaseDown[alt_index] = x
					PolBaseDown[ref_index] = 100-x
					downsamplestring = ','.join(str(x) for x in PolBaseDown) +  '\t'+ ','.join(str(x) for x in OG1_Anc) + '\t'+ ','.join(str(x) for x in OG2_Anc)
				else:
					### if we don't meet 100, the site can pass through without resampling down
					downsamplestring = resultsstring		
	
				print(resultsstring)
				print(downsamplestring)		

				with open(outfile_name+'.Est', 'a') as estFile, open(outfile_name+'.Est.Down.100', 'a') as estDownFile:	
					estFile.write(resultsstring+'\n')
					estDownFile.write(downsamplestring+'\n')		

		
		

		


