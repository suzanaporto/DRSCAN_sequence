#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
#main class
########Make snp sequences
#
#Description: program receives as input snp name and size, process the information 
#and return as output sequences
# with snps in middle position
#com as combinações das snps se necessário. Está modificado para python3
#snps assembly: 
#
########Algoritmo:
#1)Definir as variaveis
#2)Pegar informações da snp atraves do servidor e colocar a snp na lista
#3)Organizar a lista com base no local
#4)Fazer um loop na lista para ver se existe uma distancia maenor de 50 nts entre uma snp e outra
#5)Se não tiver, fazer as combinacões com uma snp e se tiver juntas snps e fazer as sequencias com as snps e as suas
#combinações
#
########O que ainda falta ser feito:
#1) Fazer tratamento de erros diversos
#
#
#
from snp import Snp
from allele import Allele
from graphDAO import GraphDAO
from snpDAO import SnpDAO
from sequencia import Sequence
from epigenome import Epigenome 
import pandas as pd
import numpy as np
import requests, sys
import argparse

#----------------------------------------------------Métodos-----------------------------------------------------///

#----------------------------------------------Programa começa aqui----------------------------------------------///

def main(args):

	#cria uma lista vazia de snps
	#lista_snp = []

	#declara o servidor
	#server = "http://rest.ensembl.org"

	#usuario define o tamanho da lista
	#lista_tamanho = int(input("Defina o tamanho da sua lista(25 ou 50): "))

	#Usuário insere quantas snps ele que fazer as sequencia
	#numero_de_snps = int(input("coloque com quantas snps você deseja fazer as sequencias: "))

	#lista de minor allele_list
	#minor_allele_list = []

	#----------------------Isso aqui dentro de um while para poder pegar as várias snps do usuario------------------///

	#for i in range(numero_de_snps):
	genome_version = args['genome_version']
	
	lista_snp = args['snp_list']

	'''
	allele1 = Allele(nome='C',local=6704583,cromossomo=1,is_comum=True,snp_pos=0)
	allele2 = Allele(nome='A',local=6704583,cromossomo=1,is_comum=False,snp_pos=0)
	allele3 = Allele(nome='T',local=6704603,cromossomo=1,is_comum=True,snp_pos=0)
	allele4 = Allele(nome='C',local=6704603,cromossomo=1,is_comum=False,snp_pos=0)
	allele5 = Allele(nome='T',local=135536917,cromossomo=5,is_comum=True,snp_pos=0)
	allele6 = Allele(nome='C',local=135536917,cromossomo=5,is_comum=False,snp_pos=0)
	allele7 = Allele(nome='T',local=6704583,cromossomo=1,is_comum=False,snp_pos=0)

	snp1 = Snp(name='rs278020',location=6704583,chrom=1,charact='snv',ancestral_al=allele1,
				minor_al=[allele2,allele7])
	snp2 = Snp(name='rs278021',location=6704603,chrom=1,charact='snv',ancestral_al=allele3,
				minor_al=[allele4])
	snp3 = Snp(name='rs123121',location=135536917,chrom=5,charact='snv',ancestral_al=allele5,
				minor_al=[allele6])
	
	lista_snp = [snp1,snp2]

	for i in lista_snp:
		print("nome: ",i.name)
		print("location: ",i.location)
		print("chromossome: ",i.chrom)
		print("wild_type",i.ancestral_al)
		for m_a in i.minor_al:
			print("variation: ",m_a.nome)
	'''
	# Perguntar se o usuário gostaria de saber se as snps estão em um elemento regulatório 
	# baseado no roadmap epigenomics

	lista_snp.sort(key=lambda x: (x.chrom, x.location), reverse=False)
	'''
	input_usuario_1 = args['verify']

	epigenome = Epigenome()
	if input_usuario_1 == 'y':
		numero_de_tecidos = int(input("Insira a quantidade tecidos para a análise: "))
		for tecidos in range(numero_de_tecidos):
			tecido_nome = input("Insira o id do tecido (ex. E071): ")
			epigenome.verificar_snps(tecido_nome,lista_snp)

	# Perguntar se o usuário quer selecionar snps
	input_usuario_2 = args['delete_snp']
	if input_usuario_2 == 'y':
		snp_number = int(input("Insira quantas snps: "))
		for indice in range(snp_number):
			nome_da_snp = input("Insira o ID da snp: ")
			for snp in lista_snp:
				contador_de_snp = 0
				if snp.name == nome_da_snp:
					lista_snp.remove(snp)
    '''

	#---------------------sort na lista/ por cromossomo e por local-------------------------///

	#sort pelo cromossomo
	#lista_snp.sort(key=lambda x: x.chrom, reverse=False)

	#sort pelo local
	#lista_snp.sort(key=lambda x: x.location, reverse=False)

	# os dois ao mesmo tempo/aparentemente o algoritmo usado para sort do python é chamado timsort
	lista_snp.sort(key=lambda x: (x.chrom, x.location), reverse=False)

	#-----------------------------------------------------------------verificar o tamanho----------------------------------------------------///


	#varias listas
	#uma lista para cada cromossomo
	lista_por_chrom = []
	#lista das snps das combinacoes
	lista_comb = []
	#nova lista que vai colocar os objetos sequencia
	#lista_das_sequencias_completas = []
	#o grafo para fazer as combinacoes
	#grafo={}
	#inicio = 0
	#fim = 0
	#ultima_pos_da_snp=0
	snp_stuff = SnpDAO()
	#contador_de_comb = 0
	# lista provisória. Acho que vou apagar depois
	#lista_de_sequencias = []

	#cria um contador
	#contador = 0
	#cont = 1
	#tamanho da lista de snps
	tamanho_lista_snp = len(lista_snp) -1
	#print ("Tamanho: ",str(tamanho_lista_snp))
	#for i in lista_snp:
		#print("NOME NO LISTA_SNP: ",i.name)
	first_write = True
	columns_df = ['Name',
				'Chromossome',
				'Location',
				'Characteristic',
				'Allele Wild Type',
				'Allele Variation']
	snp_df = pd.DataFrame(([i.name,
							i.chrom,
							i.location,
							i.charact,i.ancestral_al,i.minor_al] for i in lista_snp), columns=columns_df)
	#print(snp_df)
	#snp_df_2 = snp_df.groupby('Chromossome')['Name'].apply(list).reset_index()
	snp_df_2 = snp_df.groupby('Chromossome')
	groups = dict(list(snp_df_2))
	#iterate chromossomes groups
	for chrom_id in groups:
		#get dataframe for whole group
		chrom_df = snp_df_2.get_group(chrom_id)
		#calculate difference between locations
		chrom_df['Difference'] = chrom_df['Location'].diff()
		#null values is zero
		chrom_df.fillna(0,inplace = True)
		print(chrom_df)
		#get size until the end of iterations
		size = len(chrom_df)-1
		#check is there is only one in group
		print("SIZE: ",size)
		#make sequence if there is only one element
		if (size == 0):
			print("AQUI")
			name = chrom_df.iloc[0]['Name']
			location = chrom_df.iloc[0]['Location']
			chrom = chrom_df.iloc[0]['Chromossome']
			charact = chrom_df.iloc[0]['Characteristic']
			al_wt = chrom_df.iloc[0]['Allele Wild Type']
			al_v = chrom_df.iloc[0]['Allele Variation']
			snp = Snp(name=name,location=location,chrom=chrom,charact=charact,ancestral_al=al_wt,
						minor_al=al_v)
			first_write = snp_stuff.request_sequence(snp,genome_version,first_write)
		elif (size > 0):
			#iterate chrmossome list
			for i in range(len(chrom_df)-1):
				next_smaller = chrom_df.iloc[i + 1]['Difference'] <50
				bigger = chrom_df.iloc[i]['Difference'] >=50
				current_0 = chrom_df.iloc[i+1]['Difference'] == 0
				smaller = chrom_df.iloc[i]['Difference'] <50
				#check if it needs to be in list of combinations
				if ( (bigger and next_smaller) or (smaller and not(current_0)) or (current_0 and next_smaller) ):
					name = chrom_df.iloc[i]['Name']
					location = chrom_df.iloc[i]['Location']
					chrom = chrom_df.iloc[i]['Chromossome']
					charact = chrom_df.iloc[i]['Characteristic']
					al_wt = chrom_df.iloc[i]['Allele Wild Type']
					al_v = chrom_df.iloc[i]['Allele Variation']
					snp = Snp(name=name,location=location,chrom=chrom,charact=charact,ancestral_al=al_wt,
							minor_al=al_v)
					lista_comb.append(snp)
					#if its last iteration
					if i == size-1:
						if (next_smaller):
							name = chrom_df.iloc[i+1]['Name']
							location = chrom_df.iloc[i+1]['Location']
							chrom = chrom_df.iloc[i+1]['Chromossome']
							charact = chrom_df.iloc[i+1]['Characteristic']
							al_wt = chrom_df.iloc[i+1]['Allele Wild Type']
							al_v = chrom_df.iloc[i+1]['Allele Variation']
							snp = Snp(name=name,location=location,chrom=chrom,charact=charact,ancestral_al=al_wt,
									minor_al=al_v)
							lista_comb.append(snp)
							first_write = snp_stuff.request_sequence_combinations(lista_comb,genome_version,first_write)
							#delete combinations list
							del lista_comb[:]
				#if its not in combinations list make sequence
				else:
					#if there is more than one in combinations list make sequence
					if (len(lista_comb) > 1 ):
						first_write = snp_stuff.request_sequence_combinations(lista_comb,genome_version,first_write)
                    	#delete combinations list
						del lista_comb[:]
					name = chrom_df.iloc[i]['Name']
					location = chrom_df.iloc[i]['Location']
					chrom = chrom_df.iloc[i]['Chromossome']
					charact = chrom_df.iloc[i]['Characteristic']
					al_wt = chrom_df.iloc[i]['Allele Wild Type']
					al_v = chrom_df.iloc[i]['Allele Variation']
					snp = Snp(name=name,location=location,chrom=chrom,charact=charact,ancestral_al=al_wt,
							minor_al=al_v)
					first_write = snp_stuff.request_sequence(snp,genome_version,first_write)
	lista_comb = []
	#--------------------------------------------------------fazendo sequencias direto na verificacao---------------------------------------///
	'''
	#faz um loop na lista toda
	if tamanho_lista_snp >= 1:
                for cont_i in range(tamanho_lista_snp + 1):
                        #coloca na lista de cromossomos cada iteração de i
                        lista_por_chrom.append(lista_snp[cont_i])
                        #print("Nome: " + lista_snp[cont_i].name)
                        #print("TAmanho" + str(len(lista_snp)) )
                        #print("CONT_I: ",str(cont_i))
						#verificar se só tem um na lista e se ele é o último
						#e se o cromossomo é diferente do anterior
                        varfy_list_chrom_len = len(lista_por_chrom)==1
                        varfy_cont_i_last = cont_i==len(lista_snp)-1
                        #print('Chrom Atual' + str(lista_snp[cont_i].chrom))
                        #print('Chrom Atual' + str(lista_snp[cont_i-1].chrom))
                        varfy_before = lista_snp[cont_i].chrom != lista_snp[cont_i-1].chrom
                        #print("VERIFICAR:",varfy_before)
                        if (  (varfy_list_chrom_len) and  (varfy_cont_i_last) and (varfy_before)):
                            first_write = snp_stuff.request_sequence(lista_snp[cont_i],genome_version,first_write)
                        else:
			            		#verifica se o cromossomo atual é diferente do próximo, ou 
								#se a iteração atual é a última(penúltima snp) e 
								# o cromossomo atual é igual ao próximo ()
                                crom_dif = lista_snp[cont_i].chrom != lista_snp[cont_i+1].chrom
                                last_iteration = cont_i+1==tamanho_lista_snp
                                crom_eq = lista_snp[cont_i].chrom == lista_snp[cont_i+1].chrom
                                if (crom_dif) or (last_iteration == crom_eq) :
				                        #essa foi a mudança que eu coloquei pra poder colocar o último da lista
                                        if cont_i+1==tamanho_lista_snp and lista_snp[cont_i].chrom == lista_snp[cont_i+1].chrom:
                                                lista_por_chrom.append(lista_snp[cont_i+1])
                                        #print (str(len(lista_por_chrom)))
				                        #verifica se a lista de cromossomo for maior do que 1
                                        if len(lista_por_chrom) > 1:
					                            #faz uma iteração na lista de cromossomos até o penúltimo
                                                for i in range(len(lista_por_chrom)-1):
                                                        #se for o primerio a diferença entre primeiro e segundo for maior do que 50 ou 
														#se for diferente do primeiro a diferença entre atual e anterior é maior do que 50
														#e a diferença entre próximo e atual é maior de 50
                                                        dif_f_s = lista_por_chrom[i+1].location - lista_por_chrom[i].location >=50
                                                        dif_a_before = lista_por_chrom[i].location - lista_por_chrom[i-1].location >=50
                                                        dif_a_after = lista_por_chrom[i+1].location - lista_por_chrom[i].location >=50
                                                        if (i==0 and dif_f_s) or (i!=0 and dif_a_before and dif_f_s):
                        
														        #faz a sequência com uma snp só
                                                                print("Debug6")
                                                                first_write = snp_stuff.request_sequence(lista_por_chrom[i],genome_version,first_write)
							                                    #se a próxima iteração for a última(penultima snp) e o próximo menos o atual for maior do que 50
                                                                if i+1 == len(lista_por_chrom)-1 and dif_f_s:
                           
                                                                        #faz a sequência da última snp da lista de cromossomos
                                                                        print("Debug5")
                                                                        first_write = snp_stuff.request_sequence(lista_por_chrom[i+1],genome_version,first_write)
                    
                                                        #se tiver as condições do if
                                                        else:
                                                                #se for a primeira iteração da lista e a diferença for de menos de 50 ou 
														        #for diferente da primeira iteração e a diferenca com o proximo for menor do que 50
                                                                dif_small = lista_por_chrom[i+1].location - lista_por_chrom[i].location < 50
                                                                dif_small_2 = lista_por_chrom[i+1].location - lista_por_chrom[i].location < 50
                                                                if (i==0 and dif_small) or (i!=0 and dif_small_2):
                                                                        #coloca a snp na lista de combinações
                                                                        lista_comb.append(lista_por_chrom[i])
                                                                        dif_small_3 = lista_por_chrom[i+1].location - lista_por_chrom[i].location < 50
                                                                        #se o proximo for o penultimo da lista e a diferença entre eles for de menos de 50 
                                                                        if (i+1==len(lista_por_chrom)-1 and dif_small_3):
                                
                                                                                #coloca o proximo na lista de combinações
                                                                                lista_comb.append(lista_por_chrom[i+1])
                                                                                #faz lista de combinações
                                                                                #print (lista_comb)
                                                                                first_write = snp_stuff.request_sequence_combinations(lista_comb,genome_version,first_write)
                                                                                #deleta a lista de combinações
                                                                                del lista_comb[:]
							                                    # se não atender ao respectivo if
                                                                else:
                                                                        # se não for o primeiro e a diferença com o proximo for maior do que 50 e a diferença com o anterior for menor do que 50
                                                                        if i!=0 and lista_por_chrom[i+1].location - lista_por_chrom[i].location >=50 and lista_por_chrom[i].location - lista_por_chrom[i-1].location < 50:
            
                                                                                #coloca a snp na lista de combinacoes
                                                                                lista_comb.append(lista_por_chrom[i])
                                                                                #faz lista comb
                                                                                #print (lista_comb)
                                                                                first_write = snp_stuff.request_sequence_combinations(lista_comb,genome_version,first_write)

                                                                                #se o proximo for a ultima iteração e a diferença entre atual e ultimo for maior do que 50
                                                                                if i+1==len(lista_por_chrom)-1 and lista_por_chrom[i+1].location - lista_por_chrom[i].location >=50:

                                                                                        #faz a sequencia normal com o ultimo
                                                                                        print("Debug1")
                                                                                        first_write = snp_stuff.request_sequence(lista_por_chrom[i+1],genome_version,first_write)
                                        #lista de cromossomos é maior do que 1
                                        else:

                                                # faz a sequencia normal com a única snp da lista
                                                print("Debug3")
                                                first_write =snp_stuff.request_sequence(lista_por_chrom[0],genome_version,first_write)
                                                #() se o proximo for o ultimo da lista de snps, faça a sequencia normal do último(no caso de tiver só dois na lista)
                                                if cont_i + 1 == tamanho_lista_snp:
                                                        print("Debug2")
                                                        first_write =snp_stuff.request_sequence(lista_snp[-1],genome_version,first_write)  
                            	        #deletar a lista das snps por cromossomo 
                                        del lista_por_chrom[:]
	else:
                print("Debug4")
                first_write =snp_stuff.request_sequence(lista_snp[0],genome_version,first_write)
	'''
	return txt_2_list(args['return_list'], filename = "sequenciasdef.fna")
		
				
def txt_2_list(return_list, filename = "sequenciasdef.fna"):
  if(return_list):
    with open(filename,"r") as f:
      txt_to_list =[]
      for line in f:
        txt_to_list.append(line)
    f.close
    return txt_to_list
  else:
    return False		
	

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='')
	parser.add_argument('--snp_list', default=[], type=list)
	parser.add_argument('--verify', default='n', type=str)
	parser.add_argument('--delete_snp', default='n', type=str)
	parser.add_argument('--return_list', default=False, type=bool)
	parser.add_argument('--genome_version', default='GRCh37', type=str)
	args = vars(parser.parse_args())    
	main(args)
