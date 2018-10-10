#generate sequences
import main3 as m3
import pandas as pd


def gen_sequence(snp_name, info_list):
	res = []
	for i in info_list:
		
		snp_location = i['location']
		snp_chrom = i['chrom']
		snp_al = i['allele_wt']
		snp_al_v = i['allele_v']

		args = {"snp_name":snp_name,
			"snp_location":snp_location,
			"snp_chrom":snp_chrom,
			"snp_charac":"SNV",
			"snp_al_wt":snp_al,
			"snp_al_v":snp_al_v,
			"verify":"n",
			"delete_snp":"n",
			'return_list':True}
		
		print(args['snp_name'])
		
		res.extend(m3.main(args))
	
	print(res)
		
	#TODO res to dataframe
	return res
	
