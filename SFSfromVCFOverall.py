import gzip, os, re, numpy, time, random
from collections import defaultdict

# 
# vcf_filename = '/Users/jennyjames/Dropbox/1001GenomesProjectData/1001genomes_snp-short-indel_only_ACGTN_v3.1.synonymous.vcf.snpeff.gz'
# vcf_full_filename = '/Users/jennyjames/Dropbox/1001GenomesProjectData/1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz'

def testprint(string):
	if test ==True:
	   print(string)
	else:
		pass
		

def SFS_convert(pop_mutations, n):
	"""converts mutation output into SFS"""
	SFS = {} # full SFS with empty sites
	for x in range(1, n):
		SFS[x] = []
	for mut in pop_mutations:
		SFS[int(mut)].append(mut)
	SFS_counts = []
	for k, v in SFS.items():
		SFS_counts.append(len(v))
	return SFS_counts

   
vcf_poly = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome/1001genomes_snp-short-indel_with_tair10_only_ACGTN.4fold.Poly.Biallelic.polarized.vcf.gz'
vcf_mono = '/Users/jennyjames/Dropbox/ArabidopsisCrossGenome/1001genomes_snp-short-indel_with_tair10_only_ACGTN.4fold.Mono.Biallelic.vcf.gz'
genefile = '/Users/jennyjames/Dropbox/AssignAncestral/TAIR10_GFF3_genes.gff.txt.gz'

"""This test value is given to VCFByGene- it throttles the number of lines read in from the VCF
   when splitting based on gene"""
test = False

### will we take a haploid sample
haploid = True

### we will randomly sample down to this size, to account for missingness due to lack of sequencing information 
subsample_number = 500

if haploid:
	res_filename = vcf_poly[:-6]+'OverallSFS.haploid.csv'
else:
	res_filename = vcf_poly[:-6]+'OverallSFS.csv'



mono_counts = 0
with gzip.open(vcf_mono, 'rt') as monofile:
	for line in monofile:
		if '#' not in line:
			if 'PASS' in line:
				mono_counts = mono_counts+1
			if test == True:
				if mono_counts > 1000:
					break				
				

testprint(mono_counts)

polymorphisms = []
count = 0
with gzip.open(vcf_poly, 'rt') as vcffile:
	for line in vcffile:
	
		if '#' not in line:
			count = count+1
			
			Ref, Alt = line.split('\t')[3], line.split('\t')[4]

			AncState = line.split('AA=')[1].split(';')[0]								
			testprint(Ref+' '+Alt+' '+AncState)
	
			### we don't care whether data is phased or not
			polymorphic = re.findall('1[/|\|]0|0[/|\|]1|1[/|\|]1', line)
			not_polymorphic = re.findall('0[/|\|]0', line)
		
			if haploid:
			
				"""" choose one chromosome randomly per individual. This is relevant here due to very high levels of inbreeding, all are homozygous, i.e. allele_count % 2 == 0 """
			
				### we don't care whether data is phased or not		
				poly_haploid_sample = [re.split('/|\|', x)[random.getrandbits(1)] for x in polymorphic]
	
				allele_count = ''.join(poly_haploid_sample).count('1')
				polymorphic = allele_count
				not_polymorphic = len(not_polymorphic)

# 			if ploidy:
	
			else:
				### don't sample per individual
				allele_count = ''.join(polymorphic).count('1')	
				polymorphic = allele_count
				not_polymorphic = len(not_polymorphic)*2				
	
			"""sample down to subsample_number, and account for missingness. Any site not sequenced in at least subsample_number individuals is removed.
			   It is possible for this subsampling to reduce the number of polymorphic sites observed to 0, so they will not contribute to the SFS."""
			if not_polymorphic+polymorphic > subsample_number:
				x = numpy.random.hypergeometric(polymorphic, not_polymorphic, subsample_number)
				polymorphic, not_polymorphic = x, subsample_number-x
	
				###if resampling hasn't reduced the number of polymorphic sites to 0
				if polymorphic != 0 and not_polymorphic != 0:
					if AncState == Alt:
						derived_state_counts, anc_state_counts = not_polymorphic, polymorphic
					elif AncState == Ref:
						derived_state_counts, anc_state_counts = polymorphic, not_polymorphic
					else:
						print("Unidentified Anc state, (Ref, Alt, AA): "+Ref, Alt, AncState)
						break
						
					testprint((anc_state_counts, derived_state_counts))
					polymorphisms.append(derived_state_counts)
				else:
					testprint("Resampling has gotten rid of the site")
					testprint((polymorphic, not_polymorphic))
			else:
				testprint("Sampled in less than half of the individuals")
				testprint(not_polymorphic+polymorphic)


		if test == True:
			if count > 1000:
				break
				

SFS = SFS_convert(polymorphisms, subsample_number)
testprint(mono_counts)
testprint(len(polymorphisms))
L = mono_counts+len(polymorphisms)
with open(res_filename, 'a') as resfile:
	resfile.write(','.join(str(x) for x in SFS)+','+str(L)+'\n')



	
