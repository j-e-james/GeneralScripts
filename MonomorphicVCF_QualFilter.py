import os, re, time, gzip, numpy

###Now filter monomorphc sites to only high quality
with open('/proj/snic2020-6-184/Allotetraploid/AthalianaReads/1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz', 'rt') as VCFFile, open('/proj/snic2020-6-184/Allotetraploid/AthalianaReads/1001genomes_snp-short-indel_with_tair10_only_ACGTN.biallellic.filtered', 'a') as FilterFile:
	for line in VCFFile:
		if '#' not in line:
			if 'PASS' in line:
			
				ref, alt = line.split('\t')[3],line.split('\t')[4] 			
				### we need to ignore sites that are not bi-allelic
				if len(alt) == 1 and len(ref) == 1:	

					Qual = float(line.split('\t')[5])
					if Qual < 30:
						print(line)				
						
					else:
						Format = (re.findall('[0|1]\|[0|1]:[0-9]+:[0-9]+', line))	
						### remove sites that are not monomorphic
						if any('1|' in x or '|1' in x for x in Format):
							print('Polymorphic, skipping...')
							
						else:
							print('Monomorphic')
							Gq = [int(x.split(':')[1]) for x in Format]
							if any(x<20 for x in Gq):
								print(line)
							else:
								Dp = line.split('\t')[7]
								Dp = int(Dp.split('DP=')[1])
								if Dp != 0:
									Qual = float(line.split('\t')[5])
									Info = line.split('\t')[7]
						#			Dp = float(Info.split('DP=')[1].split(';')[0])
									QD = Qual/Dp
									if QD > 2:
										FilterFile.write(line)
		else:
			FilterFile.write(line)
							