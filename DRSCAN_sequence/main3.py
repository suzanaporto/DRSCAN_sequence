#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
#
#main class
########Make snp sequences
#
#Description: program receives as input snp name and size, process the information and return as output sequences
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
import requests, sys
import argparse

#----------------------------------------------------Métodos-----------------------------------------------------///

#----------------------------------------------Programa começa aqui----------------------------------------------///

def main(args):

	#cria uma lista vazia de snps
	lista_snp = []

	#declara o servidor
	server = "http://rest.ensembl.org"

	#usuario define o tamanho da lista
	#lista_tamanho = int(input("Defina o tamanho da sua lista(25 ou 50): "))

	#Usuário insere quantas snps ele que fazer as sequencia
	#numero_de_snps = int(input("coloque com quantas snps você deseja fazer as sequencias: "))

	#lista de minor allele_list
	minor_allele_list = []

	#----------------------Isso aqui dentro de um while para poder pegar as várias snps do usuario------------------///

	#for i in range(numero_de_snps):
	#pedir para o usuario inserir uma snp(nome da snp)
	snp_nome = args['snp_name']
	snp_local = int(args['snp_location'])
	#alguns dos atributos para o alelo
	snp_chrom = int(args['snp_chrom'])
	snp_charact = args['snp_charac']
	snp_ancestral =  args['snp_al_wt']
	minor_allele =  args['snp_al_v']
	genome_version = args['genome_version']
	minor_allele.upper()

	allele_comum_insert = Allele(nome=snp_ancestral,local=snp_local,cromossomo=snp_chrom,is_comum=True,snp_pos=0)

	for i in minor_allele:
		if i == 'A' or i == 'C' or i == 'G' or i == 'T':
			minor_allele_insert = Allele(nome=i,local=snp_local,cromossomo=snp_chrom,is_comum=False,snp_pos=0)
			minor_allele_list.append(minor_allele_insert)

	#cria um objeto do tipo snp e insere as informações da snp nele
	snp_insert = Snp(name=snp_nome, location=snp_local, chrom=snp_chrom,
                charact=snp_charact,ancestral_al= allele_comum_insert, minor_al=minor_allele_list)

	#coloca a snp na lista de snps
	lista_snp.append(snp_insert)

	#deleta tudo da lista de minor alleles
	minor_allele_list = []

	'''
	allele1 = Allele(nome='T',local=17282721,cromossomo=19,is_comum=True,snp_pos=0)
	allele2 = Allele(nome='A',local=17282721,cromossomo=19,is_comum=False,snp_pos=0)
	allele3 = Allele(nome='A',local=36662856,cromossomo=14,is_comum=True,snp_pos=0)
	allele4 = Allele(nome='G',local=36662856,cromossomo=14,is_comum=False,snp_pos=0)
	allele5 = Allele(nome='C',local=204515489,cromossomo=1,is_comum=True,snp_pos=0)
	allele6 = Allele(nome='G',local=204515489,cromossomo=1,is_comum=False,snp_pos=0)
	allele7 = Allele(nome='T',local=56909835,cromossomo=5,is_comum=True,snp_pos=0)
	allele8 = Allele(nome='A',local=56909835,cromossomo=5,is_comum=False,snp_pos=0)
	allele9 = Allele(nome='G',local=56909835,cromossomo=5,is_comum=False,snp_pos=0)
	allele10 = Allele(nome='T',local=36647076,cromossomo=14,is_comum=True,snp_pos=0)

	allele11 = Allele(nome='C',local=36647076,cromossomo=14,is_comum=False,snp_pos=0)
	allele12 = Allele(nome='C',local=27970429,cromossomo=12,is_comum=True,snp_pos=0)
	allele13 = Allele(nome='G',local=27970429,cromossomo=12,is_comum=False,snp_pos=0)
	allele14 = Allele(nome='G',local=105147856,cromossomo=4,is_comum=True,snp_pos=0)
	allele15 = Allele(nome='T',local=105147856,cromossomo=4,is_comum=False,snp_pos=0)
	allele16 = Allele(nome='C',local=54969306,cromossomo=17,is_comum=True,snp_pos=0)
	allele17 = Allele(nome='A',local=54969306,cromossomo=17,is_comum=False,snp_pos=0)
	allele18 = Allele(nome='G',local=21494705,cromossomo=10,is_comum=True,snp_pos=0)
	allele19 = Allele(nome='A',local=21494705,cromossomo=10,is_comum=False,snp_pos=0)
	allele20 = Allele(nome='C',local=36647066,cromossomo=14,is_comum=True,snp_pos=0)
	allele21 = Allele(nome='A',local=36647066,cromossomo=14,is_comum=False,snp_pos=0)
	allele22 = Allele(nome='T',local=18438074,cromossomo=19,is_comum=True,snp_pos=0)
	allele23 = Allele(nome='C',local=18438074,cromossomo=19,is_comum=False,snp_pos=0)
	allele24 = Allele(nome='A',local=13712867,cromossomo=6,is_comum=True,snp_pos=0)
	allele25 = Allele(nome='G',local=13712867,cromossomo=6,is_comum=False,snp_pos=0)
	allele26 = Allele(nome='G',local=52546588,cromossomo=16,is_comum=True,snp_pos=0)
	allele27 = Allele(nome='A',local=52546588,cromossomo=16,is_comum=False,snp_pos=0)
	allele28 = Allele(nome='C',local=21520483,cromossomo=10,is_comum=True,snp_pos=0)
	allele29 = Allele(nome='T',local=21520483,cromossomo=10,is_comum=False,snp_pos=0)
	allele30 = Allele(nome='T',local=56919325,cromossomo=5,is_comum=True,snp_pos=0)
	allele31 = Allele(nome='C',local=56919325,cromossomo=5,is_comum=False,snp_pos=0)

	snp1 = Snp(name='rs8108174',location=17282721,chrom=19,charact='snv',ancestral_al=allele1,minor_al=[allele2])
	snp2 = Snp(name='rs12883049',location=36662856,chrom=14,charact='snv',ancestral_al=allele3,minor_al=[allele4])
	snp3 = Snp(name='rs4951392',location=204515489,chrom=1,charact='snv',ancestral_al=allele5,minor_al=[allele6])
	snp4 = Snp(name='rs252923',location=56909835,chrom=5,charact='snv',ancestral_al=allele7,minor_al=[allele8, allele9])

	snp5 = Snp(name='rs72679711',location=36647076,chrom=14,charact='snv',ancestral_al=allele10,minor_al=[allele11])
	snp6 = Snp(name='rs141089797',location=27970429,chrom=12,charact='snv',ancestral_al=allele12,minor_al=[allele13])
	snp7 = Snp(name='rs62331150',location=105147856,chrom=4,charact='snv',ancestral_al=allele14,minor_al=[allele15])

	snp8 = Snp(name='rs7208403',location=54969306,chrom=17,charact='snv',ancestral_al=allele16,minor_al=[allele17])

	snp9 = Snp(name='rs12770228',location=21494705,chrom=10,charact='snv',ancestral_al=allele18,minor_al=[allele19])
	snp10 = Snp(name='rs67155512',location=36647066,chrom=14,charact='snv',ancestral_al=allele20,minor_al=[allele21])
	snp11 = Snp(name='rs2303696',location=18438074,chrom=19,charact='snv',ancestral_al=allele22,minor_al=[allele23])
	snp12 = Snp(name='rs405447',location=13712867,chrom=6,charact='snv',ancestral_al=allele24,minor_al=[allele25])

	snp13 = Snp(name='rs3095602',location=52546588,chrom=16,charact='snv',ancestral_al=allele26,minor_al=[allele27])

	snp14 = Snp(name='rs34853448',location=21520483,chrom=10,charact='snv',ancestral_al=allele28,minor_al=[allele29])
	snp15 = Snp(name='rs832533',location=56919325,chrom=5,charact='snv',ancestral_al=allele30,minor_al=[allele31])


	lista_snp = [snp1,snp2,snp3,snp4,snp5,snp6,snp7,snp8,snp9,snp10,snp11,snp12,snp13,snp14,snp15]
	'''
	# Perguntar se o usuário gostaria de saber se as snps estão em um elemento regulatório baseado no roadmap epigenomics

	lista_snp.sort(key=lambda x: (x.chrom, x.location), reverse=False)

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


	#--------------------------------sort na lista/ por cromossomo e por local--------------------------------------------///

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
	lista_das_sequencias_completas = []
	#o grafo para fazer as combinacoes
	grafo={}
	inicio = 0
	fim = 0
	ultima_pos_da_snp=0
	snp_stuff = SnpDAO()
	contador_de_comb = 0
	# lista provisória. Acho que vou apagar depois
	lista_de_sequencias = []

	#cria um contador
	contador = 0
	cont = 1
	#tamanho da lista de snps
	tamanho_lista_snp = len(lista_snp) - 1

	#--------------------------------------------------------fazendo sequencias direto na verificacao---------------------------------------///

	#faz um loop na lista toda
	if tamanho_lista_snp >= 1:
		for cont_i in range(tamanho_lista_snp):
			#coloca na lista de cromossomos cada iteração de i
			lista_por_chrom.append(lista_snp[cont_i])

			#verifica se o cromossomo atual é diferente do próximo, ou se a iteração atual é a última(penúltima snp) e o cromossomo atual é igual ao próximo ()
			if (lista_snp[cont_i].chrom != lista_snp[cont_i+1].chrom) or (cont_i+1==tamanho_lista_snp and lista_snp[cont_i].chrom == lista_snp[cont_i+1].chrom) :
				#essa foi a mudança que eu coloquei pra poder colocar o último da lista
				if cont_i+1==tamanho_lista_snp and lista_snp[cont_i].chrom == lista_snp[cont_i+1].chrom:
					lista_por_chrom.append(lista_snp[cont_i+1])

				print (str(len(lista_por_chrom)))
				#verifica se a lista de cromossomo for maior do que 1
				if len(lista_por_chrom) > 1:
					print ("Entrou aqui")
					#faz uma iteração na lista de cromossomos até o penúltimo
					for i in range(len(lista_por_chrom)-1):
                    
						#se for o primerio a diferença entre primeiro e segundo for maior do que 50 ou se for diferente do primeiro a diferença entre atual e anterior é maior do que 50 e a diferença entre próximo e atual é maior de 50
						if (i==0 and lista_por_chrom[i+1].location - lista_por_chrom[i].location >=50) or (i!=0 and lista_por_chrom[i].location - lista_por_chrom[i-1].location >=50 and lista_por_chrom[i+1].location - lista_por_chrom[i].location >=50):
                        
							#faz a sequência com uma snp só
							snp_stuff.request_sequence(lista_por_chrom[i],genome_version)
                        
							#se a próxima iteração for a última(penultima snp) e o próximo menos o atual for maior do que 50
							if i+1 == len(lista_por_chrom)-1 and lista_por_chrom[i+1].location - lista_por_chrom[i].location >=50:
                            
								#faz a sequência da última snp da lista de cromossomos
								snp_stuff.request_sequence(lista_por_chrom[i+1],genome_version)
                    
						#se tiver as condições do if
						else:
							#se for a primeira iteração da lista e a diferença for de menos de 50 ou for diferente da primeira iteração e a diferenca com o proximo for menor do que 50
							if (i==0 and lista_por_chrom[i+1].location - lista_por_chrom[i].location < 50) or (i!=0 and lista_por_chrom[i+1].location - lista_por_chrom[i].location < 50):
								#coloca a snp na lista de combinações
								lista_comb.append(lista_por_chrom[i])
            
								#se o proximo for o penultimo da lista e a diferença entre eles for de menos de 50 
								if (i+1==len(lista_por_chrom)-1 and lista_por_chrom[i+1].location - lista_por_chrom[i].location < 50):
                                
									#coloca o proximo na lista de combinações
									lista_comb.append(lista_por_chrom[i+1])
									#faz lista de combinações
									#print (lista_comb)
									snp_stuff.request_sequence_combinations(lista_comb)
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
									snp_stuff.request_sequence_combinations(lista_comb)

									#se o proximo for a ultima iteração e a diferença entre atual e ultimo for maior do que 50
									if i+1==len(lista_por_chrom)-1 and lista_por_chrom[i+1].location - lista_por_chrom[i].location >=50:

										#faz a sequencia normal com o ultimo
										snp_stuff.request_sequence(lista_por_chrom[i+1],genome_version)
				#lista de cromossomos é maior do que 1
				else:

					# faz a sequencia normal com a única snp da lista
					snp_stuff.request_sequence(lista_por_chrom[0],genome_version)
					#() se o proximo for o ultimo da lista de snps, faça a sequencia normal do último(no caso de tiver só dois na lista)
					if cont_i + 1 == tamanho_lista_snp:
						snp_stuff.request_sequence(lista_snp[-1],genome_version)  
				#deletar a lista das snps por cromossomo 
				del lista_por_chrom[:]
	else:
		snp_stuff.request_sequence(lista_snp[0],genome_version)
		
	return txt_2_list(args['return_list'], filename = "sequenciasdef.txt")
		
				
def txt_2_list(return_list, filename = "sequenciasdef.txt"):

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
	parser.add_argument('--snp_name', default="rs268", type=str)
	parser.add_argument('--snp_location', default=231322, type=int)
	parser.add_argument('--snp_chrom', default=11, type=int)
	parser.add_argument('--snp_charac', default='SNV', type=str)
	parser.add_argument('--snp_al_wt', default='A', type=str)
	parser.add_argument('--snp_al_v', default='C', type=str)
	parser.add_argument('--verify', default='n', type=str)
	parser.add_argument('--delete_snp', default='n', type=str)
	parser.add_argument('--return_list', default=False, type=bool)
	parser.add_argument('--genome_version', default='GRCh37', type=str)
	args = vars(parser.parse_args())    
	main(args)
