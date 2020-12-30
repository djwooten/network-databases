import pandas as pd
import networkx as nx
from itertools import combinations

def sanitize_gene_name(name, upper=True, to_remove=["-","_","."]):
    """
    Will capitalize (if upper==True) and remove characters from the gene name
    """
    name = name.strip()
    if upper: name = name.upper()
    for rem in to_remove:
        name = name.replace(rem,"")
    return name

def sanitize_genes(genes, upper=True, to_remove=["-","_","."]):
    """
    Sanitizes a list of gene names
    """
    ret = []
    for g in genes:
        ret.append(sanitize_gene_name(g, upper=upper, to_remove=to_remove))
    return ret

def build_alias_g(entrez=True, hgnc=True, ensembl=True, upper=True, to_remove=["-","_","."]):
    G = nx.Graph()

    if entrez or isinstance(entrez, str):
        if isinstance(entrez, str):
            aliases = pd.read_csv(entrez, delimiter="\t")
        else:
            aliases = pd.read_csv("../../data/aliases/Homo_sapiens.gene_info.tsv", delimiter="\t")

        for i in aliases.index:
            gene = aliases.loc[i, "Symbol"]
            gene = sanitize_gene_name(gene, upper=upper, to_remove=to_remove)
            
            if gene=="": continue

            G.add_node(gene)

            syns = aliases.loc[i, "Synonyms"].split("|")
            for syn in syns:
                if syn=="-": continue
                syn = sanitize_gene_name(syn, upper=upper, to_remove=to_remove)
                if syn=="": continue
                
                G.add_edge(gene, syn)

    if hgnc or isinstance(hgnc, str):
        if isinstance(hgnc, str):
            aliases = pd.read_csv(entrez, delimiter="\t")
        else:
            aliases = pd.read_csv("../../data/aliases/hgnc.tsv", delimiter="\t")

        for i in aliases.index:
            gene = aliases.loc[i, "Approved symbol"]
            if not isinstance(gene, str): continue
            gene = sanitize_gene_name(gene, upper=upper, to_remove=to_remove)
            
            if gene=="": continue

            G.add_node(gene)

            syns = []
            prev_syms = aliases.loc[i, "Previous symbols"]
            if isinstance(prev_syms, str):
                syns += prev_syms.split(",")
            
            alias_syms = aliases.loc[i, "Alias symbols"]
            if isinstance(alias_syms, str):
                syns += alias_syms.split(",")
            
            syns = set(syns)
            for syn in syns:
                if syn=="-" or syn=="": continue
                syn = sanitize_gene_name(syn, upper=upper, to_remove=to_remove)
                if syn=="": continue
                
                G.add_edge(gene, syn)

    if ensembl or isinstance(ensembl, str):
        if isinstance(hgnc, str):
            aliases = pd.read_csv(ensembl, delimiter="\t")
        else:
            aliases = pd.read_csv("../../data/aliases/mart_export.tsv", delimiter="\t")

        for i in aliases.index:
            gene = aliases.loc[i, "Gene stable ID"]
            gene = sanitize_gene_name(gene, upper=upper, to_remove=to_remove)
            
            if gene=="" or not isinstance(gene, str): continue

            G.add_node(gene)

            syns = []
            gname = aliases.loc[i, "Gene name"]
            if isinstance(gname, str):
                syns.append(gname)
            gsyn = aliases.loc[i, "Gene Synonym"]
            if isinstance(gsyn, str):
                syns.append(gsyn)

            for syn in set(syns):
                if syn=="-" or syn=="": continue
                syn = sanitize_gene_name(syn, upper=upper, to_remove=to_remove)
                if syn=="": continue
                
                G.add_edge(gene, syn)
    return G

def get_alias_graph():
    """
    Reads the alias graph from ../../data/aliases/alias_graph.graphml

    You should read it once, and pass it to the necessary functions
    """
    return nx.read_graphml("../../data/aliases/alias_graph.graphml")

def get_aliases(gene, alias_graph=None):
    if alias_graph is None:
        alias_graph = get_alias_graph()

    if (not alias_graph.has_node(gene)): return set()
    
    return nx.descendants(alias_graph, gene)

def get_all_aliases(genes, alias_graph=None):
    if alias_graph is None:
        alias_graph = get_alias_graph()

    aliases = set(genes)
    for gene in genes:
        aliases = aliases.union(get_aliases(gene, alias_graph=alias_graph))
    return aliases

def process_gene_file(fname, outfname=None, alias_graph=None):
    genes = set()
    with open(fname,"r") as infile:
        for line in infile:
            genes.add(line.strip())
    
    if alias_graph is None: alias_graph = get_alias_graph()
    aliases = get_all_aliases(fname, alias_graph=alias_graph)
    
    if outfname is None:
        extension = fname.split('.')[-1]
        file_prefix = '.'.join(fname.split('.')[:-1])
        outfname = "%s_aliases.%s"%(fileprefix, extension)
    outfile = open(outfname, "w")
    outfile.write("\n".join(aliases))
    outfile.close()

def relabel_network(G, master_gene_list, alias_graph=None):
    mapping = dict()
    for gene in G.nodes():
        if gene in master_gene_list: continue
        
        if alias_graph is None: alias_graph = get_alias_graph()

        aliases = get_aliases(gene, alias_graph=alias_graph)
        for alias in aliases:
            if alias in master_gene_list:
                mapping[gene] = alias
                break
    return nx.relabel_nodes(G, mapping)

if __name__ == "__main__":
    G = build_alias_g()
    nx.write_graphml(G, "../../data/aliases/alias_graph.graphml")