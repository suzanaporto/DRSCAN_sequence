#!usr/local/bin/python3
# -*- coding: utf-8 -*-
#
#Class: class was made to make methods that send a request to the server and get the sequences
#Attributes in methods:
#list_size <int>
#snp_info <snp>
#middle_location <int>
#
#
import requests, sys
from graphDAO import GraphDAO

class SnpDAO (object):

    #request para sequencia com a snp np meio/retorna a string
    def request_sequence_middle(self, start_location,end_location, snp_chrom,genome_version):
        #declara o servidor
        server = "http://rest.ensembl.org"

        #declara as variáveis para montar as urls para fazer o request
        x = start_location - 50
        y = end_location + 50
        ext = "/sequence/region/human/" + str(snp_chrom) +":"+ str(x) + ".."+ str(y) + ":1?coord_system_version=" + genome_version

        r1 = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

        if not r1.ok:
          r1.raise_for_status()
          sys.exit()

        return r1.text
    #sem uso no programa atualmente
    def request_sequence_head(self, list_size, head_location,snp_chrom):
        #declara o servidor
        server = "http://rest.ensembl.org"

        #declara as variáveis para montar as urls para fazer o request
        right_left_size = int(list_size/2)
        x = head_location - list_size
        y = middle_location + list_size
        ext = "/sequence/region/human/" + str(snp_chrom) +":"+ str(x) + ".."+ str(y) + ":1?"

        r1 = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

        if not r1.ok:
          r1.raise_for_status()
          sys.exit()

        return r1.text
    #sem uso no programa
    def request_sequence_tail(self, list_size, tail_location,snp_chrom):
        #declara o servidor
        server = "http://rest.ensembl.org"

        #declara as variáveis para montar as urls para fazer o request
        right_left_size = int(list_size/2)
        x = tail_location - list_size
        y = middle_location + list_size
        ext = "/sequence/region/human/" + str(snp_chrom) +":"+ str(x) + ".."+ str(y) + ":1?coord_system_version=" + genome_version

        r1 = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})

        if not r1.ok:
          r1.raise_for_status()
          sys.exit()

        return r1.text      

    def request_sequence(self,snp,genome_version,filename="sequenciasdef.fna"):
        #declara o servidor
        server = "http://rest.ensembl.org"
        x = snp.location - 50
        y = snp.location + 50
        ext = "/sequence/region/human/" + str(snp.chrom) +":"+ str(x) + ".."+ str(y) + ":1?coord_system_version=" + genome_version
        r1 = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
        if not r1.ok:
            r1.raise_for_status()
            sys.exit()
		# Coloca as sequencias de snp no meio na variavel declarada
        starting_at_one = 1
        tamanho_seq = 50*2 + 1
        seq_meio = r1.text
        seq_meio = seq_meio[:50] + snp.ancestral_al.nome + seq_meio[51:]
        print (">sequence_wild_type|"+str(snp.name) +"|"+str(snp.chrom)+"|"+str(snp.ancestral_al)+"|"+str((snp.location-x)+starting_at_one))
        print (seq_meio)
        
        f = open(filename,"w")   #create add file in write mode
        line_seq = ">sequence_wild_type|"+str(snp.name) +"|"+str(snp.chrom)+"|"+str(snp.ancestral_al)+"|"+str((snp.location-x)+starting_at_one) + '\n'
        f.write(line_seq)
        f.write(seq_meio + '\n')  #writes o/p to add.txt file
        for j in snp.minor_al:
            seq_meio_alt = seq_meio[:50] + j.nome + seq_meio[51:]
            print (">sequence_variation|"+str(snp.name) +"|"+str(snp.chrom)+"|"+j.nome+"|"+str((snp.location-x)+starting_at_one))
            print (seq_meio_alt)
            f.write(">sequence_variation|"+str(snp.name) +"|"+str(snp.chrom)+"|"+j.nome+"|"+str((snp.location-x)+starting_at_one)+ '\n')
            f.write(seq_meio_alt + '\n')  #writes o/p to add.txt file
        f.close()

    def request_sequence_combinations(self,lista_comb,genome,first_alleles=[],last_alleles=[],lista_de_alelos=[],lista_de_comb_sets=[],cont=1):
        g = GraphDAO()
        pos_relativa_lista = []
        #ver as snps da lista de combinacoes e coloca os seus alelos na lista de alelos com numeros para representar as suas posicoes e snps
        for j in lista_comb:
            j.ancestral_al.snp_pos = cont
            lista_de_alelos.append(j.ancestral_al)
            for x in j.minor_al:
                x.snp_pos = cont
                lista_de_alelos.append(x)
            cont+=1
            
        #coloca os primeiros alelos numa lista e os ultimos alelos em outra lista
        for y in lista_de_alelos:
            if y.snp_pos == 1:
                first_alleles.append(y)
            else:
                if y.snp_pos == lista_de_alelos[-1].snp_pos:
                    ultima_pos_da_snp = lista_de_alelos[-1].snp_pos
                    last_alleles.append(y)
        
        #criar grafo com a lista de alelos devidamente preenchida/por algum motivo ele precisou do numero da ultima posicao
        grafo = g.create_graph(lista_de_alelos,ultima_pos_da_snp)
        
        #fazer as combinacoes
        for k in first_alleles:
            for f in last_alleles:
                lista_de_comb_sets.append(g.find_all_paths(grafo,start=k,end=f))
                
        #fazer as sequencias
        first_location_snp = lista_comb[0].location
        first_location_chrom = lista_comb[0].chrom
        last_location_snp = lista_comb[-1].location
        last_location_chrom = lista_comb[-1].location
        #sequencia buscada do ensembl.org
        request_text_middle = self.request_sequence_middle(first_location_snp,last_location_snp,first_location_chrom,genome)
        for combinacao in lista_de_comb_sets:
            for j in combinacao:
                #quantidade de elementos
                for l in range(len(j)):
                    if (l == 0):
                        #index do primeiro elemento
                        real_index = 50
                    else:
                        #calculo do index real
                        real_index = real_index + lista_comb[l].location - lista_comb[l-1].location
                    print("É pra ser o index real da snp2: " + str(real_index+1))
                    pos_relativa_lista.append(str(real_index+1))
                    request_text_middle = request_text_middle[:real_index] + j[l] + request_text_middle[real_index+1:]
                print("combinacoes: ",j)
                string_dos_alelos = '|'.join(j)
                string_nomes_snps = '|'.join(map(str, lista_comb))
                pos_relativa = '|'.join(pos_relativa_lista)
                #string_nomes_snps = '|'.join(str(lista_comb))
                print (">sequence_combinations|"+ string_nomes_snps + "|" + str(lista_comb[0].chrom) + "|" +string_dos_alelos+ "|" + pos_relativa)
                print(request_text_middle)
                f = open("sequenciasdef.fna","w")   #create add file in write mode
                f.write(">sequence_combinations|"+ string_nomes_snps + "|" + str(lista_comb[0].chrom) + "|" +string_dos_alelos+ "|" + pos_relativa + '\n')
                f.write(request_text_middle + '\n')  #writes o/p to add.txt file
                f.close()
                pos_relativa_lista = []
        real_index = 0
        cont = 1
        del first_alleles[:]
        del last_alleles[:]
        del lista_comb[:]
        del lista_de_alelos[:]
        del lista_de_comb_sets[:]
