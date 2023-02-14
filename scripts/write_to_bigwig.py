import numpy as np
import pyBigWig as bw
import pandas as pd


def write_bw(df,outfile):
	# Write dataframe to bigwig
	# Dataframe shape: chr start end value
	if not np.all(df.columns[:3] == ['chr','start','end']):
		print('Error: df columns must be [chr start end [values]]')

	# Create bw file
	bw_fid = bw.open(outfile,'w')

	# get header
	CHR = np.unique(df.chr)
	my_head = [(c, df[df['chr']==c].iloc[-1,2]) for c in CHR]
	bw_fid.addHeader(my_head)
	for c in CHR:
		my_chr = df[df.chr==c].chr.values.astype(str)
		start = df[df.chr==c].start.values
		end = df[df.chr==c].end.values-1
		vals = df[df.chr==c].iloc[:,-1].values
		bw_fid.addEntries(my_chr, start, ends=end, values=vals)
	bw_fid.close()

def write_bed(df,outfile):

	fid_neg = open(f'{outfile}_negative_regions.bed','w')
	fid_0 = open(f'{outfile}_neutral_regions.bed','w')
	fid_pos = open(f'{outfile}_positive_regions.bed','w')
	for l in df.values:
		if l[3] == -1:
			fid_neg.write('%s\t%u\t%u\n'%(l[0],l[1],l[2]))
		elif l[3] == 0:
			fid_0.write('%s\t%u\t%u\n'%(l[0],l[1],l[2]))
		elif l[3] == 1:
			fid_pos.write('%s\t%u\t%u\n'%(l[0],l[1],l[2]))
	fid_neg.close()
	fid_0.close()
	fid_pos.close()
