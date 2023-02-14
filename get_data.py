import numpy as np
import pandas as pd
import pyBigWig
import sys
sys.path.insert(1, 'scripts')
from write_to_bigwig import *

if __name__ == '__main__':

	# bw track files to load
	track_folder = 'data/tracks'
	track_file= {
			'Lmnb2_wt':'RPE_SCII_Lmnb2_Combined_track_NOSmooth_20kb.bw',
			'Lmnb2_ko':'RPE_Top2B_Lmnb2_Combined_track_NOSmooth_20kb.bw',
			#'Delta_Lmnb2':'RPE_Top2B_Lmnb2_Differential_trackSmooth_20kb.bw',
			'Top2':'Differential_CCseq_VP16vsWT_Nosmooth.bw',
			'H3K9me2_wt':'pADamID-RPE_SCII_H3K9me2-20kb-combined.bw',
			'H3K9me3_wt':'pADamID-RPE_SCII_H3K9me3-20kb-combined.bw',
			'H3K9me2_ko':'pADamID-RPE_Top2B_H3K9me2-20kb-combined.bw',
			'H3K9me3_ko':'pADamID-RPE_Top2B_H3K9me3-20kb-combined.bw',
			'Gilbert_supercoiling':'supercoil_mean_hg38.bw',
			'Ste_supercoiling_wt':'bTMP_AVD_SC_800-20kb.bw',
			'Ste_supercoiling_ko':'bTMP_AVD_Beta_800-20kb.bw'
	}

	# get chr start end for dataframe
	CHR = [f'chr{i+1}' for i in range(22)] + ['chrX','chrY']
	bin_size = 20000
	chr=[]
	start=[]
	end=[]
	for c in CHR:
		# get longest chr across tracks
		l=0
		for f in track_file:
			infile = f'{track_folder}/{track_file[f]}'
			bw = pyBigWig.open(infile)
			l = max(l,bw.chroms()[c])
			bw.close()

		n = l//bin_size + 1
		chr.extend([c]*n)
		start.extend(np.arange(0,n*bin_size,bin_size))
		end.extend(np.arange(bin_size,n*bin_size,bin_size))
		end.append(l)
	N = len(chr)

	# initialize data frame
	data = {'chr':chr,'start':start,'end':end}
	for f in track_file:
		data[f] = np.zeros(N)
	Tracks = pd.DataFrame(data=data)

	# now fill in values
	for f in track_file:
		print(f)
		infile = f'{track_folder}/{track_file[f]}'
		bw = pyBigWig.open(infile)
		for i in Tracks.index:
			Tracks.loc[i,f] = np.nanmean( bw.values(Tracks.loc[i,'chr'],Tracks.loc[i,'start'],Tracks.loc[i,'start']+1) )
		bw.close()

	# averaging and difference H3K9me2/3
	Tracks['H3K9me2'] = Tracks[['H3K9me2_wt','H3K9me2_ko']].mean(axis=1)
	Tracks['H3K9me3'] = Tracks[['H3K9me3_wt','H3K9me3_ko']].mean(axis=1)
	Tracks['Delta_H3K9me2'] = Tracks['H3K9me2_ko'] - Tracks['H3K9me2_wt']
	Tracks['Delta_H3K9me3'] = Tracks['H3K9me3_ko'] - Tracks['H3K9me3_wt']
	Tracks['Delta_Lmnb2'] = Tracks['Lmnb2_ko'] - Tracks['Lmnb2_wt']
	Tracks['Delta_Ste_supercoiling'] = Tracks['Ste_supercoiling_ko']- Tracks['Ste_supercoiling_wt']

	# fill in RNA-seq values
	track_file={
		'RNAseq1':'RPE170_10_Scramble_rep1_stranded_PE_Signal.Unique.str1.out.bw',
		'RNAseq2':'RPE170_10_Scramble_rep2_stranded_PE_Signal.Unique.str1.out.bw',
	}
	for f in track_file:
		print(f)
		infile = f'{track_folder}/{track_file[f]}'
		bw = pyBigWig.open(infile)
		i=0
		for i in Tracks.index:
			Tracks.loc[i,f] = np.nansum(bw.values(Tracks.loc[i,'chr'],Tracks.loc[i,'start'],Tracks.loc[i,'end']))
	Tracks['RNAseq'] = np.log( Tracks['RNAseq1'] + Tracks['RNAseq2'] + 1 )
	#Tracks['RNAseq'][np.isinf(Tracks['RNAseq'])] = np.nan

	bin_size_kb = int(bin_size/1000)

	outfile = f'results/Tracks_{bin_size_kb}kb.txt'
	Tracks.to_csv(outfile,sep='\t',index=False)

	# save RNAseq and Delta Lmnb2 as bigwig
	print('write to bigwig')
	for my_track in ['RNAseq','Delta_Lmnb2']:
		outfile = f'results/bigwig/{my_track}_{bin_size_kb}kb.bw'
		write_bw(Tracks.loc[:,['chr','start','end',my_track]],outfile)

	print('done')
