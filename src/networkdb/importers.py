import networkx as nx
import json
import unicodedata
import pandas as pd

def load_cytoscape_cx(fname):
    
    infile = open(fname, "r")
    data = json.load(infile)
    infile.close()

    G = nx.MultiDiGraph()
    nodes = []
    edges = []
    node_attributes = []
    edge_attributes = []
    positions = []  
    for d in data:
        if "nodes" in d:
            nodes = d["nodes"]
        if "edges" in d:
            edges = d["edges"]
        if "nodeAttributes" in d:
            node_attributes = d["nodeAttributes"]
        if "edgeAttributes" in d:
            edge_attributes = d["edgeAttributes"]
        if "cartesianLayout" in d:
            positions = d["cartesianLayout"]

    if len(nodes)==0:
        return G

    # Nodes:
    idx_to_node = dict()
    idx_to_attributes = dict()
    idx_to_position = dict()

    for position in positions:
        try:
            idx = position['node']
            idx_to_position[idx] = (position['x'], position['y'])
        except:
            pass

    for node_attribute in node_attributes:
        try:
            idx = node_attribute['po']
            attribute = node_attribute['n']
            value = node_attribute['v']
        
            if not idx in idx_to_attributes:
                idx_to_attributes[idx] = dict()
            idx_to_attributes[idx][attribute] = value
        except:
            pass

    for node in nodes:
        try:
            name = node['n']
            if name=="":
                continue
            idx = node['@id']

            idx_to_node[idx] = name

            attributes = dict()

            # Get attributes
            if idx in idx_to_attributes:
                attributes = idx_to_attributes[idx]

            # Get position
            if idx in idx_to_position:
                attributes['pos'] = idx_to_position[idx]

            # Sometimes other attributes are hidden here
            for k in node:
                if not k in ['n', '@id']:
                    attributes[k] = node[k]
            
            
            G.add_node(name, **attributes)
        except:
            pass


    # Edges:
    idx_to_edge = dict()
    idx_to_edge_attributes = dict()

    for edge_attribute in edge_attributes:
        try:
            idx = edge_attribute['po']
            attribute = edge_attribute['n']
            value = edge_attribute['v']
            
            if not idx in idx_to_edge_attributes:
                idx_to_edge_attributes[idx] = dict()

            idx_to_edge_attributes[idx][attribute] = value
        except:
            pass

    for edge in edges:
        try:
            idx = edge['@id']
            source = idx_to_node[edge['s']]
            target = idx_to_node[edge['t']]

            attributes = dict()
            attributes['fname']=fname

            # Get attributes
            if idx in idx_to_edge_attributes:
                attributes = idx_to_edge_attributes[idx]

            # Get other attributes
            for k in edge.keys():
                if not k in ['@id','s','t']:
                    attributes[k] = edge[k]
            
            G.add_edge(source, target, **attributes)
        except:
            pass

    return G

def remove_control_characters(s):
    return "".join(ch for ch in s if unicodedata.category(ch)[0]!="C")

def load_signor(fname=None, direct_only=False, min_score=0):
    
    if fname is None:
        fname = "../../data/networks/all_data_12_11_20.tsv"

    G = nx.MultiDiGraph()
    
    id_dict = dict()
    
    df = pd.read_csv(fname, delimiter="\t")
    for _i in df.index:
        line = list(df.loc[_i])

        entityA, typeA, idA, databaseA, entityB, typeB, idB, databaseB, effect, mechanism, residue, sequence, tax_id, cell_Data, tissue_Data, modulator_complex, target_complex, modificationA, modAseq, modificationB, modBseq, pmid, direct, notes, annotator, sentence, score, signor_id = line

        if pd.isna(entityA) or pd.isna(entityB):
            continue

        if direct_only:
            if not direct.upper()=="YES": continue
        try:
            score = float(score)
        except ValueError:
            score = 0
        

        if score < min_score: continue

        id_dict[entityA.upper()] = {'id':idA, 'type':typeA, 'database':databaseA}
        id_dict[entityB.upper()] = {'id':idB, 'type':typeB, 'database':databaseB}


        if effect is None: effect=""
        if mechanism is None: mechanism=""
        if pmid is None: pmid=""
        if direct is None: direct=""
        if sentence is None: sentence=""
        if score is None: score=0
        if signor_id is None: signor_id=""

        if isinstance(sentence, str):
            sentence = remove_control_characters(sentence)
        else:
            sentence=""

        G.add_edge(entityA.upper(), entityB.upper(), effect=effect, mechanism=mechanism, pmid=pmid, direct=direct, sentence=sentence, score=score, signor_id=signor_id, fname=fname.split("/")[-1])

    nx.set_node_attributes(G, id_dict)
    
    return G

def load_trrust(fname=None):

    if fname is None:
        fname = "../../data/networks/trrust_rawdata.human.tsv"

    G = nx.MultiDiGraph()
    df = pd.read_csv(fname, header=None, delimiter="\t", dtype=str)
    for _i in df.index:
        line = list(df.loc[_i])
            
        source, target, effect, pmids = line

        G.add_edge(source.upper(), target.upper(), effect=effect, pmid=pmids, fname=fname.split("/")[-1])

    return G

def load_regnetwork(fname=None, valid_evidence=None, minimum_confidence=None):

    if fname is None:
        fname = "../../data/networks/regnetwork_nov_12_20.csv"

    if minimum_confidence is None or minimum_confidence.upper() == "LOW":
        minimum_confidence = ["LOW", "MEDIUM", "HIGH"]
    elif minimum_confidence.upper() == "MEDIUM":
        minimum_confidence = ["MEDIUM", "HIGH"]
    else:
        minimum_confidence = ["HIGH"]

    G = nx.MultiDiGraph()
    df = pd.read_csv(fname, header=None, dtype=str)
    for _i in df.index:
        line = list(df.loc[_i])
        source, sid, target, tid, databases, evidence, confidence = line
        #databases = database.split(",")
        
        if valid_evidence is not None and evidence.upper()!=valid_evidence.upper(): continue
        if not confidence.upper() in minimum_confidence: continue

        
        G.add_edge(source.upper(), target.upper(), databases=databases, evidence=evidence, confidence=confidence, fname=fname.split("/")[-1])
    
    return G

if __name__=="__main__":
    pass
    print("Signor")
    G = load_signor()
    nx.write_graphml(G, "../../data/networks/signor.graphml")
    
    print("Trrust")
    G = load_trrust()
    nx.write_graphml(G,"../../data/networks/trrust.graphml")
    
    print("regnetwork")
    G = load_regnetwork()
    nx.write_graphml(G,"../../data/networks/regnetwork.graphml")
    