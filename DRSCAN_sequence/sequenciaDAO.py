
class SequenciaDAO(object):


    #esse metodo é para substituir no principal
    def sequencias_combinacao(self,lista_combinacoes):
        first_alleles = []
        last_alleles = []
        lista_comb = []
        lista_de_alelos = []
        lista_de_comb_sets = []

        for j in lista_comb:
            print (j.name)
            j.ancestral_al.snp_pos = cont
            lista_de_alelos.append(j.ancestral_al)
            for x in j.minor_al:
                x.snp_pos = cont
                lista_de_alelos.append(x)
            cont +=1

        for y in lista_de_alelos:
            if y.snp_pos == 1:
                first_alleles.append(y)
            else:
                if y.snp_pos == lista_de_alelos[-1].snp_pos:
                    ultima_pos_da_snp = lista_de_alelos[-1].snp_pos
                    last_alleles.append(y)
    
        #criar grafo com a lista de alelos devidamente preenchida/por algum motivo ele precisou do numero da ultima posicao
        grafo = g.create_graph(lista_de_alelos,ultima_pos_da_snp)
    
        print (grafo)
    
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
        request_text_middle = snp_stuff.request_sequence_middle(first_location_snp,last_location_snp,first_location_chrom)
        for k in lista_de_comb_sets:
            for j in k:
                print (j)
                for i in range(len(j)):
                    print (i)
                    print(j[i])
                    if (i == 0):
                        real_index = 50
                    else:
                        real_index = real_index + lista_comb[i].location - lista_comb[i-1].location
                    print (real_index)
                    request_text_middle = request_text_middle[:real_index] + j[i] + request_text_middle[real_index+1:]
                print(request_text_middle)
                #colocas só as sequencias nessa lista
                lista_de_sequencias.append(request_text_middle)
                #criar objetos para colocá-los numa lista de sequencias
        #deletar as listas
        real_index = 0
        cont = 1
        del first_alleles[:]
        del last_alleles[:]
        del lista_comb[:]
        del lista_de_alelos[:]
        del lista_de_comb_sets[:]