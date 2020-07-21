# plotting.py - plotting utilities
# ASS, 21 Apr 2020
#
# This file contains some utility functions for plotting CRNs
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import random
import networkx as nx
import statistics
from bokeh.models import (BoxSelectTool, Circle,Square, EdgesAndLinkedNodes, HoverTool,
                          MultiLine, NodesAndLinkedEdges, Plot, Range1d, TapTool,PanTool,WheelZoomTool)
from bokeh.palettes import Spectral4
from bokeh.models.graphs import from_networkx
from fa2 import ForceAtlas2
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt

def updateLimits(limits,xvalues):
    for value in xvalues:
        if(value < limits[0]):
            limits[0] = value
        if(value > limits[1]):
            limits[1] = value
    return limits

def makeArrows2(graph_renderer,graph,positions,headsize=3,headangle=np.pi/6):
    """this function draws an arrow shape at the end of graph lines"""
    xs,ys = [],[]
    xbounds = [0,0]
    ybounds = [0,0]
    for edge in graph.edges:
        #iterate through all the edges
        from_node = edge[0]
        to_node = edge[1]
        from_x = positions[from_node][0]
        from_y = positions[from_node][1]
        to_x = positions[to_node][0]
        to_y = positions[to_node][1]
        updateLimits(xbounds,[from_x,to_x])
        updateLimits(ybounds,[from_y,to_y])
        #above, we get all the variables
        #below, we are assuming the "to" position is the middle of the
        #coordinate space
        ydif = from_y-to_y
        xdif = from_x-to_x
        #next we calculate the angle from the destination node 
        #to the source node
        angl = np.arctan2(ydif,xdif)
        #the arrow consists of three added points, one on either side 
        #of the line and one in the middle
        p1x = to_x+headsize*np.cos(angl+headangle) #left side of the arrow
        p1y = to_y+headsize*np.sin(angl+headangle) #left side of the arrow
        p2x = to_x+headsize*np.cos(angl-headangle) #right side of the arrow
        p2y = to_y+headsize*np.sin(angl-headangle) #right side of the arrow
        p3x = to_x+headsize*.7*np.cos(angl) #middle of the arrow
        p3y = to_y+headsize*.7*np.sin(angl) #middle of the arrow
        xs.append([from_x,p3x,p1x,to_x,p2x,p3x]) #'xs' is a list of lists which represent each line from node to node
        ys.append([from_y,p3y,p1y,to_y,p2y,p3y]) #'ys' is the same thing except the y positions
    graph_renderer.edge_renderer.data_source.data['xs'] = xs #this part replaces the lines with the ones made by this function
    graph_renderer.edge_renderer.data_source.data['ys'] = ys
    return xbounds,ybounds

def graphPlot(DG,DGspecies,DGreactions,plot,layout="force",positions=None,posscale = 1.0,layoutfunc=None):
    """given a directed graph, plot it!
    Inputs:
    DG: a directed graph of type DiGraph
    DGspecies: a directed graph which only contains the species nodes
    DGreactions: a directed graph which only contains the reaction nodes
    plot: a bokeh plot object
    layout: graph layout function. 
                'force' uses fa2 to push nodes apart
                'circle' plots the nodes and reactions in two overlapping circles, with the reactions on the inside of the circle
                'custom' allows user input "layoutfunc". Internally, layoutfunc is passed the three inputs (DG, DGspecies, DGreactions)
                                                        and should output a position dictionary with node {<node number>:(x,y)}

    positions: a dictionary of node names and x,y positions. this gets passed into the layout function
    posscale: multiply the scaling of the plot. This only affects the arrows because the arrows are a hack :("""
    if(layout=="force"):
        #below are parameters for the force directed graph visualization
        forceatlas2 = ForceAtlas2(
                            # Behavior alternatives
                            outboundAttractionDistribution=True,  # Dissuade hubs
                            linLogMode=False,  # NOT IMPLEMENTED
                            adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                            edgeWeightInfluence=1.0,

                            # Performance
                            jitterTolerance=1.0,  # Tolerance
                            barnesHutOptimize=True,
                            barnesHutTheta=1.2,
                            multiThreaded=False,  # NOT IMPLEMENTED

                            # Tuning
                            scalingRatio=2.4*posscale,
                            strongGravityMode=False,
                            gravity=1.0,

                            # Log
                            verbose=False)

        positions = forceatlas2.forceatlas2_networkx_layout(DG, pos=positions, iterations=2000) 
    elif(layout == "circle"):
        positions = nx.circular_layout(DGspecies,scale=50*posscale)
        positions.update(nx.circular_layout(DGreactions,scale=35*posscale))
    elif(layout == "custom"):
        positions = layoutfunc(DG,DGspecies,DGreactions)
    reaction_renderer = from_networkx(DGreactions, positions, center=(0, 0))
    species_renderer = from_networkx(DGspecies, positions, center=(0, 0))
    edges_renderer = from_networkx(DG, positions, center=(0, 0))

    #edges
    edges_renderer.node_renderer.glyph = Circle(size=12,line_alpha=0,fill_alpha=0, fill_color="color")
    edges_renderer.edge_renderer.glyph = MultiLine( line_alpha=0.2, line_width=4)
    edges_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
    edges_renderer.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)
    xbounds,ybounds = makeArrows2(edges_renderer,DG,positions,headsize=5) #make the arrows!
    
    #we want to find the middle of the graph and plot a square that is 1:1 aspect ratio
    
    #find the midpoint of the graph
    xmid = statistics.mean(xbounds)
    ymid = statistics.mean(ybounds)
    #now, subtract the middle from the edges
    xmiddlized = [a-xmid for a in xbounds]
    ymiddlized = [a-ymid for a in ybounds]
    #now, find the biggest dimension
    absdim = max([abs(a) for a in xmiddlized+ymiddlized])
    xlim = [xmid-absdim*1.05,xmid + absdim*1.05]
    ylim = [ymid-absdim*1.05,ymid + absdim*1.05]
    #now set it on the plot!
    plot.x_range = Range1d(xlim[0],xlim[1])
    plot.y_range = Range1d(ylim[0],ylim[1])

    #reactions
    reaction_renderer.node_renderer.glyph = Square(size=8, fill_color=Spectral4[0])
    reaction_renderer.node_renderer.selection_glyph = Square(size=8, fill_color=Spectral4[2])
    reaction_renderer.node_renderer.hover_glyph = Square(size=8, fill_color=Spectral4[1])

    #nodes
    species_renderer.node_renderer.glyph = Circle(size=12, fill_color="color")
    species_renderer.node_renderer.selection_glyph = Circle(size=15, fill_color=Spectral4[2])
    species_renderer.node_renderer.hover_glyph = Circle(size=15, fill_color=Spectral4[1])
    
    #this part adds the interactive elements that make it so that the lines are highlighted 
    #when you mouse over and click
    edge_hover_tool = HoverTool(tooltips= None,renderers=[edges_renderer])
    species_hover_tool = HoverTool(tooltips=[("name", "@species"), ("type", "@type")],\
                                        renderers=[species_renderer],attachment="right")
    rxn_hover_tool = HoverTool(tooltips=[("reaction", "@species"), ("type", "@type"),("k_f","@k"),("k_r","@k_r")],\
                                        renderers=[reaction_renderer],attachment="right")
    
    plot.add_tools(edge_hover_tool,species_hover_tool,rxn_hover_tool, TapTool(), BoxSelectTool(),PanTool(),WheelZoomTool())

    edges_renderer.selection_policy = NodesAndLinkedEdges()
    edges_renderer.inspection_policy = EdgesAndLinkedNodes()

    plot.renderers.append(edges_renderer)
    plot.renderers.append(reaction_renderer)
    plot.renderers.append(species_renderer)

def generate_networkx_graph(CRN,useweights=False,use_pretty_print=False,pp_show_material=True,
                                                    pp_show_rates=True,pp_show_attributes=True,
                                                colordict=None):
    """generates a networkx DiGraph object that represents the CRN.
    input:
    ==========================
    CRN: a CRN from mixture.get_model() for example
    useweights: this will attempt to represent the reaction rates by the length of edges.
                short edges are fast rates. It doesn't look very good usually
    use_pretty_print: this uses the "pretty print" function to represent reactions and nodes a bit cleaner
    the next three parameters are pretty_print parameters
        pp_show_material: default false because this is listed in "type"
        pp_show_rates: default true because this is useful information
        pp_show_attributes
    colordict: a dictionary containing which node types are what color based upon the following keywords:
        Keywords are chosen to match species.material_type
            {"complex": "cyan",
            "protein": "green",
            "dna": "grey",
            "rna": "orange",
            "ligand": "pink",
            "phosphate": "yellow",
            "nothing":"purple"}

    When using a custom colordict, the following attributes will be checked to find colors with the first keys taking precedence:
        repr(species): "color"
        species.name: "color"
        (species.material_type, tuple(species.attributes)): "color"
        species.material_type: "color"
        tuple(species.attributes): "color"

    output:
    ==================
    CRNgraph: the DiGraph object containing all nodes and edges
    CRNspeciesonly: a DiGraph object with only species
    CRNreactionsonly: a DiGraph object with only reactions
    """
    if not colordict:
        colordict = {"complex":"cyan","protein":"green",
                     "dna":"grey","rna":"orange",
                     "ligand":"pink","phosphate":"yellow","nothing":"purple"}
    CRNgraph = nx.DiGraph()
    allnodenum = 1 #every node has an index
    #this starts at 1 because "nothing" is node 0
    nodedict = {} #this is so that we can write out the reactions in
                #the reaction "species" field
                #it has {species:index}
    rxnlist = [] #list of numbers corresponding to only reaction nodes
    #sometimes reactions have no products. I want this to be represented in the graph with this
    #"nothing" node. However, usually we are making degradation reactions which yield the
    #degradation enzyme, so then it doesn't go to nothing. This means actually this node
    #isn't use for anything. But i think it's good to have just in case.
    defaultcolor = "grey"
    nodedict["nothing"]=0
    CRNgraph.add_node(0)
    CRNgraph.nodes[0]["type"]="nothing"
    CRNgraph.nodes[0]["species"]="nothing"
    if("nothing" in colordict):
        CRNgraph.nodes[0]["color"]=colordict["nothing"]
    for species in CRN.species:
        #add all species first

        if repr(species) in colordict:
            mycol = colordict[repr(species)]
        elif species.name in colordict:
            mycol = colordict[species.name]
        elif (species.material_type, tuple(species.attributes)) in colordict:
            mycol = colordict[(species.material_type, tuple(species.attributes))]
        elif(species.material_type in colordict):
            mycol = colordict[species.material_type]
        elif tuple(species.attributes) in colordict:
            mycol = colordict[tuple(species.attributes)]
        else:
            mycol = defaultcolor

        nodedict[species]=allnodenum
        CRNgraph.add_node(allnodenum)
        CRNgraph.nodes[allnodenum]["type"]=str(species.material_type)
        if(not use_pretty_print):
            CRNgraph.nodes[allnodenum]["species"]=str(species)
        else:
            spectxt = species.pretty_print(pp_show_material)
            CRNgraph.nodes[allnodenum]["species"]=spectxt
        CRNgraph.nodes[allnodenum]["color"]=mycol
        allnodenum +=1
    #reactions follow, allnodenum is not reset between these two loops
    for rxn in CRN.reactions:
        CRNgraph.add_node(allnodenum)
        CRNgraph.nodes[allnodenum]["type"]= str(rxn.propensity_type)
        CRNgraph.nodes[allnodenum]["k"] = str(rxn.k_forward)
        CRNgraph.nodes[allnodenum]["k_r"] = str(rxn.k_reverse)
        default_color = "blue"
        #CRNgraph.nodes[allnodenum]
        kval = rxn.k_forward
        if(not useweights):
            kval = 1
        krev_val = rxn.k_reverse
        if((krev_val is not None) and (not useweights)):
            krev_val = 1
        for reactant in rxn.inputs:
            CRNgraph.add_edge(nodedict[reactant.species],allnodenum,weight=kval)
            if(krev_val is not None):
                #if the k is 0 then the node does not exist, right?
                CRNgraph.add_edge(allnodenum,nodedict[reactant.species],weight=krev_val)
        for product in rxn.outputs:
            CRNgraph.add_edge(allnodenum,nodedict[product.species],weight=kval)
            if(krev_val is not None):
                CRNgraph.add_edge(nodedict[product.species],allnodenum,weight=krev_val)
        if(len(rxn.outputs)==0):
            #this adds an edge to the "nothing" node we made in the beginning
            CRNgraph.add_edge(allnodenum,0,weight=kval)
            if(krev_val is not None):
                CRNgraph.add_edge(0,allnodenum,weight=krev_val)
        elif(len(rxn.inputs)==0):
            #this adds an edge from the "nothing" node we made in the beginning
            CRNgraph.add_edge(0,allnodenum,weight=kval)
            if(krev_val is not None):
                CRNgraph.add_edge(allnodenum,0,weight=krev_val)
        CRNgraph.nodes[allnodenum]["color"]=default_color
        if(not use_pretty_print):
            CRNgraph.nodes[allnodenum]["species"]=str(rxn)
        else:
            rxntxt = rxn.pretty_print(show_material=pp_show_material, show_rates=pp_show_rates, show_attributes=pp_show_attributes)
            CRNgraph.nodes[allnodenum]["species"]=rxntxt #this will show up as "reaction" in the tool tip
        #the name of the reaction is the string representation
        rxnlist += [allnodenum]
        allnodenum +=1
    CRNspeciesonly = CRNgraph.copy()
    CRNspeciesonly.remove_nodes_from(rxnlist)
    CRNreactionsonly = CRNgraph.copy()
    CRNreactionsonly.remove_nodes_from(range(rxnlist[0]))
    return CRNgraph,CRNspeciesonly,CRNreactionsonly