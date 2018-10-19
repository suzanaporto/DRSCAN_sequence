#generate sequences
import main3 as m3
import pandas as pd
import argparse

def gen_sequence(snp_name, info_list,genome_version):

	#parser = argparse.ArgumentParser(description='')
	#parser.add_argument('--snp_name', default="rs268", type=str)
	#parser.add_argument('--snp_location', default=231322, type=int)
	#parser.add_argument('--snp_chrom', default=11, type=int)
	#parser.add_argument('--snp_charac', default='SNV', type=str)
	#parser.add_argument('--snp_al_wt', default='A', type=str)
	#parser.add_argument('--snp_al_v', default='C', type=str)
	#parser.add_argument('--verify', default='n', type=str)
	#parser.add_argument('--delete_snp', default='n', type=str)
	#parser.add_argument('--return_list', default=False, type=bool)
	#parser.add_argument('--genome_version', default=GRCh37, type=str)
	#args = vars(parser.parse_args())  
	
	res = []
	for i in info_list:
		
		snp_location = i['location']
		snp_chrom = i['chrom']
		snp_al = i['allele_wt']
		snp_al_v = i['allele_v']

		args = {"snp_name":snp_name,
				"snp_location":snp_location,
				"snp_chrom":snp_chrom,
				"snp_al_wt":snp_al,
				"snp_al_v":snp_al_v,
				"snp_charac":"SNV",
				'return_list':True,
				"verify":'n',
				"genome_version":genome_version}
		
		res.extend(m3.main(args))
	
	print(res)
		
	#TODO res to dataframe
	return res
	
