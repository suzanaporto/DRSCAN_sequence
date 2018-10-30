#generate sequences
import main3 as m3
from snp import Snp
from allele import Allele
import pandas as pd
import argparse

def gen_sequence(snp_names, info_list,genome_version):
	#parser = argparse.ArgumentParser(description='')
	#parser.add_argument('--snp_name', default="rs268", type=list)
	#parser.add_argument('--verify', default='n', type=str)
	#parser.add_argument('--delete_snp', default='n', type=str)
	#parser.add_argument('--return_list', default=False, type=bool)
	#parser.add_argument('--genome_version', default=GRCh37, type=str)
	#args = vars(parser.parse_args())
	snp_list_object = []
	minor_allele_list = []
	res = []
  for i,snp_name in zip(info_list,snp_names):
    snp_location = i['location']
    snp_chrom = i['chrom']
    snp_al = i['allele_wt']
    snp_al_v = i['allele_v']
    
    allele_comum_insert = Allele(nome=snp_al,local=snp_location,cromossomo=snp_chrom,is_comum=True,snp_pos=0)
    
    for i in snp_al_v:
      if i == 'A' or i == 'C' or i == 'G' or i == 'T':
        minor_allele_insert = Allele(nome=i,local=snp_location,cromossomo=snp_chrom,is_comum=False,snp_pos=0)
        minor_allele_list.append(minor_allele_insert)
        
        #cria um objeto do tipo snp e insere as informações da snp nele
        snp_insert = Snp(name=snp_name, location=snp_location, chrom=snp_chrom,
                         charact="SNV",ancestral_al= allele_comum_insert,
                         minor_al=minor_allele_list)
        #coloca a snp na lista de snps
        snp_list_object.append(snp_insert)
  args = {"snp_list":snp_list_object,
          "return_list":True,
          "verify":'n',
          "delete_snp":'n',
          "genome_version":genome_version}
  
  res.append(m3.main(args))
  #print(res)
  return res
	
