import networkx as nx
import pandas as pd
import itertools

def out_neighbors(G, node, domain=None):
    """
    Parameters:
    -----------
        G : nx.DiGraph
            Network
        
        node : str
            The node to get out_neighbors of
        
        domain : array-like or None
            If not None, targets of node that are not in domain are ignored

    Returns:
    --------
    targets : list
        Targets of node that are in the given domain (or all targets)
    """
    neighbors = set()
    for edge in G.edges(node):
        neighbor = edge[1]
        
        if domain is not None and neighbor in domain:
            neighbors.add(neighbor)
        
        else: neighbors.add(neighbor)
    return neighbors

def descendents_stop_at_nodes(G, start_node, stop_nodes=[]):
    """
    Parameters:
    -----------
        G : nx.DiGraph
            Network
        
        start_node : str
            The node to start with
        
        stop_nodes : array-like
            The search for descendents of start_node will not continue past stop_nodes
        
    """
    descendents = set([start_node])
    if start_node in stop_nodes: return descendents
    children = out_neighbors(G, start_node)
    
    while len(children) > 0:
        child = children.pop()
        descendents.add(child)
        
        if not child in stop_nodes:
            new_children = out_neighbors(G, child).difference(descendents)
            for x in new_children: children.add(x)
        
    return descendents

def get_regulators(G):
    regulators = set()
    for node in G.nodes():
        if G.out_degree(node) > 0:
            regulators.add(node)
    return regulators

def get_source_nodes(G):
    sources = set()
    for node in G.nodes():
        if G.in_degree(node) == 0:
            sources.add(node)
    return sources

def get_sink_nodes(G):
    sinks = set()
    for node in G.nodes():
        if G.out_degree(node) == 0:
            sinks.add(node)
    return sinks

def remove_sink_nodes(G, new_G=None):
    if new_G is None:
        new_G = G.copy()
    sink_nodes = get_sink_nodes(new_G)
    if len(sink_nodes) == 0:
        return new_G
    else:
        for node in sink_nodes:
            new_G.remove_node(node)
        return remove_sink_nodes(new_G, new_G=new_G)

def get_massive_component(G):
    sccs = dict()
    biggest = None
    for scc in nx.strongly_connected_components(G):
        sccG = G.subgraph(scc)
        sccs[sccG] = sccG.number_of_nodes()
        if biggest is None: biggest=sccG
        elif sccs[biggest] < sccs[sccG]: biggest=sccG
    return sccs, biggest


def merge_networks(*args):
    if len(args)==0: return nx.DiGraph()
    if len(args)==1: return args[0]

    G = nx.compose(args[0], args[1])
    for other_g in args[2:]:
        G = nx.compose(G, other_g)
    
    return G

def attempt_include_nodes_in_G(bigG, littleG, nodes, max_path_length=2):
    
    new_G = littleG.copy()

    current_nodes = list(littleG.nodes())

    try_node_results = dict()
    for try_node in nodes:
        if try_node in current_nodes:
            continue
        # shortest paths from littleG to node
        to_paths = []
        to_node_l = 1000000000

        # shortest path from node to littleG
        from_paths = []
        from_node_l = 1000000000

        for cur_node in current_nodes:
            # Find a path from littleG to the node
            try:
                to_path = nx.shortest_path(bigG, cur_node, try_node)
                l = len(to_path)-1
                if l < to_node_l:
                    to_paths = [to_path,]
                    to_node_l = l
                elif l == to_node_l:
                    to_paths.append(to_path)
            except:
                pass

            # Find a path from the node to littleG
            try:
                from_path = nx.shortest_path(bigG, try_node, cur_node)
                l = len(from_path)-1
                if l < from_node_l:
                    from_paths = [from_path,]
                    from_node_l = l
                elif l == from_node_l:
                    from_paths.append(from_path)
            except:
                pass
                
        if from_node_l + to_node_l <= max_path_length:
            print(from_node_l, to_node_l)
            for path in from_paths + to_paths:
                add_path(new_G, path)
    
    return new_G

def add_path(G, path):
    if len(path) < 2: return
    G.add_edge(path[0], path[1])
    for _i in range(1, len(path)-1):
        G.add_edge(path[_i], path[_i+1])

def get_sccs(G):
    if G.is_directed():
        return list(nx.components.strongly_connected_component_subgraphs(G))
    else:
        return list(nx.components.strongly_connected_component_subgraphs(G))


def remove_node(G, node):
    """
    Removes a node, but connects all of its in_neighbors to its out_neighbors
    """
    e_in = list(G.in_edges(node))
    e_out = list(G.out_edges(node))
    G.remove_node(node)
    for e in itertools.product(e_in, e_out):
        G.add_edge(e[0], e[1], indirect="Through %s"%node)
    return G
