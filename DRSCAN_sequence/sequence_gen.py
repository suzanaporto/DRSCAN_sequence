#generate sequences
import main3 as m3
import pandas as pd


def gen_sequence(snp_name, info_list):
	res = []
	for i in info_list:
		
		snp_location = i['location']
		snp_chrom = i['chrom']
		snp_al = i['allele_wt']
		snp_v = i['allele_v']
		verf = "n"
		del_s = "n"
		type_s = "SNV"
		
		print("location ",snp_location)

		args = {"snp_name":snp_name,
			"snp_location":snp_location,
			"snp_chrom":snp_chrom,
			"snp_al_wt":snp_al,
			"snp_al_v":snp_v,
			"snp_charac":type_s,
			"verify":verf,
			"delete_snp":del_s,
			"return_list":True}
		
		print(args['snp_name'])
		
		res.append(m3.main(args))
	
	print(res)
		
	#TODO res to dataframe
	return res
	
