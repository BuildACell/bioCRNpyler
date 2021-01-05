# plotting.py - plotting utilities
# ASS, 21 Apr 2020
#
# This file contains some utility functions for plotting CRNs
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import math
import random
import statistics
from warnings import warn

from .components_basic import Protein
from .dna_part_cds import CDS
from .dna_part_misc import AttachmentSite
from .dna_part_promoter import Promoter
from .dna_part_rbs import RBS
from .dna_part_terminator import Terminator
from .propensities import MassAction

HAVE_MATPLOTLIB = False
try:
    import matplotlib.pyplot as plt
    from matplotlib import cm

    HAVE_MATPLOTLIB = True
except ModuleNotFoundError:
    pass
PLOT_DNA = False
try:
    import dnaplotlib as dpl
    PLOT_DNA = True
except ModuleNotFoundError:
    pass

if(PLOT_DNA and not HAVE_MATPLOTLIB):
    PLOT_DNA = False

PLOT_NETWORK = False
try:
    import networkx as nx
    from bokeh.models import (BoxSelectTool, Circle, EdgesAndLinkedNodes,
                              HoverTool, MultiLine, NodesAndLinkedEdges,
                              PanTool, Plot, Range1d, Square, TapTool,
                              WheelZoomTool)
    from bokeh.plotting import from_networkx
    from bokeh.palettes import Spectral4
    from fa2 import ForceAtlas2
    PLOT_NETWORK = True
except ModuleNotFoundError:
    pass


def updateLimits(limits, xvalues):
    for value in xvalues:
        if(value < limits[0]):
            limits[0] = value
        if(value > limits[1]):
            limits[1] = value
    return limits


def makeArrows2(graph_renderer, graph, positions, headsize=3, headangle=math.pi/6):
    """this function draws an arrow shape at the end of graph lines"""
    xs, ys = [], []
    xbounds = [0, 0]
    ybounds = [0, 0]
    for edge in graph.edges:
        # iterate through all the edges
        from_node = edge[0]
        to_node = edge[1]
        from_x = positions[from_node][0]
        from_y = positions[from_node][1]
        to_x = positions[to_node][0]
        to_y = positions[to_node][1]
        updateLimits(xbounds, [from_x, to_x])
        updateLimits(ybounds, [from_y, to_y])
        # above, we get all the variables
        # below, we are assuming the "to" position is the middle of the
        # coordinate space
        ydif = from_y-to_y
        xdif = from_x-to_x
        # next we calculate the angle from the destination node
        # to the source node
        angl = math.atan2(ydif, xdif)
        # the arrow consists of three added points, one on either side
        # of the line and one in the middle
        p1x = to_x+headsize*math.cos(angl+headangle)  # left side of the arrow
        p1y = to_y+headsize*math.sin(angl+headangle)  # left side of the arrow
        p2x = to_x+headsize*math.cos(angl-headangle)  # right side of the arrow
        p2y = to_y+headsize*math.sin(angl-headangle)  # right side of the arrow
        p3x = to_x+headsize*.7*math.cos(angl)  # middle of the arrow
        p3y = to_y+headsize*.7*math.sin(angl)  # middle of the arrow
        # 'xs' is a list of lists which represent each line from node to node
        xs.append([from_x, p3x, p1x, to_x, p2x, p3x])
        # 'ys' is the same thing except the y positions
        ys.append([from_y, p3y, p1y, to_y, p2y, p3y])
    # this part replaces the lines with the ones made by this function
    graph_renderer.edge_renderer.data_source.data['xs'] = xs
    graph_renderer.edge_renderer.data_source.data['ys'] = ys
    return xbounds, ybounds


def graphPlot(DG, DGspecies, DGreactions, plot, layout="force", positions=None, posscale=1.0, layoutfunc=None, iterations=2000, rseed=30):
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
    random.seed(rseed)
    if(not PLOT_NETWORK):
        warn("network plotting disabled because some libraries are not found")
        return
    if(layout == "force"):
        # below are parameters for the force directed graph visualization
        forceatlas2 = ForceAtlas2(
            # Behavior alternatives
            outboundAttractionDistribution=True,  # Dissuade hubs
            linLogMode=False,  # NOT IMPLEMENTED
            # Prevent overlap (NOT IMPLEMENTED)
            adjustSizes=False,
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

        positions = forceatlas2.forceatlas2_networkx_layout(
            DG, pos=positions, iterations=iterations)
    elif(layout == "circle"):
        positions = nx.circular_layout(DGspecies, scale=50*posscale)
        positions.update(nx.circular_layout(DGreactions, scale=35*posscale))
    elif(layout == "custom"):
        positions = layoutfunc(DG, DGspecies, DGreactions)
    reaction_renderer = from_networkx(DGreactions, positions, center=(0, 0))
    species_renderer = from_networkx(DGspecies, positions, center=(0, 0))
    edges_renderer = from_networkx(DG, positions, center=(0, 0))

    # edges
    edges_renderer.node_renderer.glyph = Circle(
        size=12, line_alpha=0, fill_alpha=0, fill_color="color")
    edges_renderer.edge_renderer.glyph = MultiLine(
        line_alpha=0.2, line_width=4, line_join="round")
    edges_renderer.edge_renderer.selection_glyph = MultiLine(
        line_color=Spectral4[2], line_width=5, line_join="round")
    edges_renderer.edge_renderer.hover_glyph = MultiLine(
        line_color=Spectral4[1], line_width=5, line_join="round")
    xbounds, ybounds = makeArrows2(
        edges_renderer, DG, positions, headsize=5)  # make the arrows!

    # we want to find the middle of the graph and plot a square that is 1:1 aspect ratio

    # find the midpoint of the graph
    xmid = statistics.mean(xbounds)
    ymid = statistics.mean(ybounds)
    # now, subtract the middle from the edges
    xmiddlized = [a-xmid for a in xbounds]
    ymiddlized = [a-ymid for a in ybounds]
    # now, find the biggest dimension
    absdim = max([abs(a) for a in xmiddlized+ymiddlized])
    xlim = [xmid-absdim*1.05, xmid + absdim*1.05]
    ylim = [ymid-absdim*1.05, ymid + absdim*1.05]
    # now set it on the plot!
    plot.x_range = Range1d(xlim[0], xlim[1])
    plot.y_range = Range1d(ylim[0], ylim[1])

    # reactions
    reaction_renderer.node_renderer.glyph = Square(
        size=8, fill_color=Spectral4[0])
    reaction_renderer.node_renderer.selection_glyph = Square(
        size=8, fill_color=Spectral4[2])
    reaction_renderer.node_renderer.hover_glyph = Square(
        size=8, fill_color=Spectral4[1])

    # nodes
    species_renderer.node_renderer.glyph = Circle(size=12, fill_color="color")
    species_renderer.node_renderer.selection_glyph = Circle(
        size=15, fill_color=Spectral4[2])
    species_renderer.node_renderer.hover_glyph = Circle(
        size=15, fill_color=Spectral4[1])

    # this part adds the interactive elements that make it so that the lines are highlighted
    # when you mouse over and click
    edge_hover_tool = HoverTool(tooltips=None, renderers=[edges_renderer])
    species_hover_tool = HoverTool(tooltips=[("name", "@species"), ("type", "@type")],
                                   renderers=[species_renderer], attachment="right")
    rxn_hover_tool = HoverTool(tooltips=[("reaction", "@species"), ("type", "@type"), ("k_f", "@k"), ("k_r", "@k_r")],
                               renderers=[reaction_renderer], attachment="right")

    plot.add_tools(edge_hover_tool, species_hover_tool, rxn_hover_tool,
                   TapTool(), BoxSelectTool(), PanTool(), WheelZoomTool())

    edges_renderer.selection_policy = NodesAndLinkedEdges()
    edges_renderer.inspection_policy = EdgesAndLinkedNodes()

    plot.renderers.append(edges_renderer)
    plot.renderers.append(reaction_renderer)
    plot.renderers.append(species_renderer)


def generate_networkx_graph(CRN, useweights=False, use_pretty_print=False, pp_show_material=True,
                            pp_show_rates=True, pp_show_attributes=True,
                            pp_show_compartments=True,
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
    if(not PLOT_NETWORK):
        warn("network plotting disabled because some libraries are not found")
        return None, None, None
    if not colordict:
        colordict = {"complex": "cyan", "protein": "green",
                     "dna": "grey", "rna": "orange",
                     "ligand": "pink", "phosphate": "yellow", "nothing": "purple"}
    CRNgraph = nx.DiGraph()
    allnodenum = 1  # every node has an index
    # this starts at 1 because "nothing" is node 0
    nodedict = {}  # this is so that we can write out the reactions in
    # the reaction "species" field
    # it has {species:index}
    rxnlist = []  # list of numbers corresponding to only reaction nodes
    # sometimes reactions have no products. I want this to be represented in the graph with this
    # "nothing" node. However, usually we are making degradation reactions which yield the
    # degradation enzyme, so then it doesn't go to nothing. This means actually this node
    # isn't use for anything. But i think it's good to have just in case.
    defaultcolor = "grey"
    nodedict["nothing"] = 0
    CRNgraph.add_node(0)
    CRNgraph.nodes[0]["type"] = "nothing"
    CRNgraph.nodes[0]["species"] = "nothing"
    if("nothing" in colordict):
        CRNgraph.nodes[0]["color"] = colordict["nothing"]
    for species in CRN.species:
        # add all species first

        if repr(species) in colordict:
            mycol = colordict[repr(species)]
        elif species.name in colordict:
            mycol = colordict[species.name]
        elif (species.material_type, tuple(species.attributes)) in colordict:
            mycol = colordict[(species.material_type,
                               tuple(species.attributes))]
        elif(species.material_type in colordict):
            mycol = colordict[species.material_type]
        elif tuple(species.attributes) in colordict:
            mycol = colordict[tuple(species.attributes)]
        else:
            mycol = defaultcolor

        nodedict[species] = allnodenum
        CRNgraph.add_node(allnodenum)
        CRNgraph.nodes[allnodenum]["type"] = str(species.material_type)
        if(not use_pretty_print):
            CRNgraph.nodes[allnodenum]["species"] = str(species)
        else:
            spectxt = species.pretty_print(
                show_material=pp_show_material, show_compartment=pp_show_compartments)
            CRNgraph.nodes[allnodenum]["species"] = spectxt
        CRNgraph.nodes[allnodenum]["color"] = mycol
        allnodenum += 1
    # reactions follow, allnodenum is not reset between these two loops
    for rxn in CRN.reactions:
        CRNgraph.add_node(allnodenum)
        CRNgraph.nodes[allnodenum]["type"] = str(rxn.propensity_type)
        if isinstance(rxn.propensity_type, MassAction):
            CRNgraph.nodes[allnodenum]["k"] = str(
                rxn.propensity_type.k_forward)
            CRNgraph.nodes[allnodenum]["k_r"] = str(
                rxn.propensity_type.k_reverse)
        else:
            CRNgraph.nodes[allnodenum]["k"] = str(rxn.propensity_type.k)
            CRNgraph.nodes[allnodenum]["k_r"] = ''

        default_color = "blue"
        # CRNgraph.nodes[allnodenum]
        if isinstance(rxn.propensity_type, MassAction):
            kval = rxn.propensity_type.k_forward
            CRNgraph.nodes[allnodenum]["k"] = str(kval)
        else:
            kval = rxn.k
            CRNgraph.nodes[allnodenum]["k"] = str(rxn.propensity_type.k)

        if(not useweights):
            kval = 1
        if isinstance(rxn.propensity_type, MassAction):
            krev_val = rxn.propensity_type.k_reverse
        else:
            krev_val = None
        if((krev_val is not None) and (not useweights)):
            krev_val = 1
        for reactant in rxn.inputs:
            
            CRNgraph.add_edge(nodedict[reactant.species],allnodenum,weight=kval)
            if(krev_val is not None):
                # if the k is 0 then the node does not exist, right?
                CRNgraph.add_edge(
                    allnodenum, nodedict[reactant.species], weight=krev_val)
        for product in rxn.outputs:
            #TODO species cannot find another species in the nodedict????
            CRNgraph.add_edge(allnodenum,nodedict[product.species],weight=kval)
            if(krev_val is not None):
                CRNgraph.add_edge(
                    nodedict[product.species], allnodenum, weight=krev_val)
        if(len(rxn.outputs) == 0):
            # this adds an edge to the "nothing" node we made in the beginning
            CRNgraph.add_edge(allnodenum, 0, weight=kval)
            if(krev_val is not None):
                CRNgraph.add_edge(0, allnodenum, weight=krev_val)
        elif(len(rxn.inputs) == 0):
            # this adds an edge from the "nothing" node we made in the beginning
            CRNgraph.add_edge(0, allnodenum, weight=kval)
            if(krev_val is not None):
                CRNgraph.add_edge(allnodenum, 0, weight=krev_val)
        CRNgraph.nodes[allnodenum]["color"] = default_color
        if(not use_pretty_print):
            CRNgraph.nodes[allnodenum]["species"] = str(rxn)
        else:
            rxntxt = rxn.pretty_print(
                show_material=pp_show_material, show_rates=pp_show_rates, show_attributes=pp_show_attributes)
            # this will show up as "reaction" in the tool tip
            CRNgraph.nodes[allnodenum]["species"] = rxntxt
        # the name of the reaction is the string representation
        rxnlist += [allnodenum]
        allnodenum += 1
    CRNspeciesonly = CRNgraph.copy()
    CRNspeciesonly.remove_nodes_from(rxnlist)
    CRNreactionsonly = CRNgraph.copy()
    CRNreactionsonly.remove_nodes_from(range(rxnlist[0]))
    return CRNgraph, CRNspeciesonly, CRNreactionsonly


def make_dpl_from_construct(construct, showlabels=None):
    """ This function creats a dictionary suitable for
    input into dnaplotlib for plotting constructs.
    Inputs:
    construct: a DNA_construct object
    showlabels: list of part types to show labels for. For example, [AttachmentSite,Terminator]"""
    # TODO make showlabels more general
    if(showlabels is None):
        showlabels = []
    outdesign = []
    if(HAVE_MATPLOTLIB):
        cmap = cm.Set1(range(len(construct.parts_list)*2))
    pind = 0
    for part in construct.parts_list:
        pcolor = part.color
        pcolor2 = part.color2
        if(HAVE_MATPLOTLIB):
            if(type(pcolor) == int):
                c1 = cmap[pcolor][:-1]
            else:
                c1 = cmap[pind][:-1]
            if(type(pcolor2) == int):
                c2 = cmap[pcolor2][:-1]
            else:
                c2 = cmap[random.choice(
                    list(range(len(construct.parts_list))))][:-1]
        showlabel = False
        if(type(part) in showlabels):
            showlabel = True
        outdesign += make_dpl_from_part(part, direction=part.direction == "forward",
                                        color=c1, color2=c2, showlabel=showlabel)
        pind += 1
    return outdesign


def make_dpl_from_part(part, direction=None, color=None, color2=None, showlabel=False):
    """ This function creats a dictionary suitable for
    input into dnaplotlib for plotting constructs.
    Inputs:
    part: a DNA_part object
    direction: True for forward, False for reverse. If you leave it as None, it will take from the DNA_part object
    color: this is the color of the part. Tuple with relative rgb values. if the DNA_part has a defined color it will take that first before
                looking at this variable
    color2: this is the secondary color of the part. Only relevant for RecombinaseSite2 components. Basically the
            same idea as color, above
    showlabel: if True, the label of this part will be shown."""
    regs = []
    if(direction is None and part.direction is not None):
        direction = part.direction == "forward"
    elif(direction is None):
        direction = True
    if(type(part.color) is not int):
        color = part.color
    elif(color is not None):
        part.color = color
    if(type(part.color2) is not int):
        color2 = part.color2
    elif(color2 is not None):
        part.color2 = color2
    dpl_type = "UserDefined"  # this is the default part type
    if(hasattr(part, "dpl_type")):
        part_dpl = part.dpl_type
    else:
        part_dpl = None

    if(isinstance(part, Promoter)):
        dpl_type = "Promoter"
        if(hasattr(part, "regulators")):
            regs = part.regulators
    elif(isinstance(part, RBS)):
        dpl_type = "RBS"
    elif(isinstance(part, CDS)):
        dpl_type = "CDS"
    elif(isinstance(part, Protein)):
        dpl_type = "CDS"
    elif(isinstance(part, Terminator)):
        dpl_type = "Terminator"
    elif(isinstance(part, AttachmentSite)):
        if(part.site_type == "attP" or part.site_type == "attB"):
            dpl_type = "RecombinaseSite"
        elif(part.site_type == "attL" or part.site_type == "attR"):
            dpl_type = "RecombinaseSite2"
    if(part_dpl is not None):
        # parts can have their own pre-set dnaplotlib types
        dpl_type = part_dpl
    outdesign = [{'type': dpl_type, "name": part.name,
                  "fwd": direction, 'opts': {'color': color, 'color2': color2}}]
    for reg in regs:
        # promoters with regulators have a number of "operator" symbols on them
        outdesign += [{"type": "Operator", "name": str(reg), "fwd": direction, 'opts': {
            'color': color, 'color2': color2}}]
    if(showlabel):
        outdesign[0]["opts"].update(
            {'label': str(part.name), 'label_size': 13, 'label_y_offset': -8, })
    if(not direction):
        outdesign = outdesign[::-1]
    return outdesign

def plotDesign(design,renderer = None,part_renderers=None,\
                circular=False,title=None,outfig=None):
    """helper function for doing dnaplotlib plots. You need to set the size and min max of the
    plot, and that's what this function does"""
    if(PLOT_DNA):
        if(renderer is None):
            renderer = dpl.DNARenderer(scale=5, linewidth=3)
        if(part_renderers is None):
            part_renderers = renderer.SBOL_part_renderers()
        fig = plt.figure(figsize=(len(design)*.75, 1.1))
        ax = fig.add_axes([0, 0, 1, 1])
        try:
            start, end = renderer.renderDNA(
                ax, design, part_renderers, circular=circular)
        except TypeError:
            start, end = renderer.renderDNA(ax, design, part_renderers)
        ax.axis('off')
        if title is not None:
            ax.set_title(title)
        addedsize=1
        ax.set_xlim([start-addedsize,end+addedsize])
        ax.set_ylim([-15,15])
        if(outfig is not None):
            fig.savefig(outfig,format='png')
        else:
            plt.show()
    else:
        warn("plotting DNA has been disabled because you don't have DNAplotlib")

def plotConstruct(DNA_construct_obj,dna_renderer=None,\
                                    rna_renderer=None,\
                                    plot_rnas=False,debug=False,showlabels = None,plot_dna_test=True,outfig=None):
    """helper function for making dnaplotlib plots of a DNA_construct object. Plots the
    DNAs and the RNAs that come from that DNA, using DNA_construct.explore_txtl"""
    # TODO: make the label showing more general
    if(showlabels is None):
        showlabels = []
    if(PLOT_DNA and plot_dna_test):
        if(dna_renderer is None):
            dna_renderer = dpl.DNARenderer(scale=5, linewidth=3)
        if(rna_renderer is None):
            rna_renderer=dpl.DNARenderer(scale = 5,linewidth=3,linecolor=(1,0,0))
    
    design = make_dpl_from_construct(DNA_construct_obj,showlabels=showlabels)
    circular=DNA_construct_obj.circular
    if(PLOT_DNA and plot_dna_test):
        plotDesign(design,circular=circular,title=DNA_construct_obj.get_species(),outfig=outfig)
        if(plot_rnas):
            rnas_and_proteins = DNA_construct_obj.enumerate_constructs()
            for component in rnas_and_proteins:
                if(component.get_species().material_type=="rna"):
                    rnadesign = make_dpl_from_construct(component,showlabels=showlabels)
                else:
                    continue
                rnacolor = rna_renderer.linecolor
                for part in rnadesign:
                    if("edgecolor" not in part['opts']):
                        part['opts'].update({'edgecolor':rnacolor})
                plotDesign(rnadesign,renderer=rna_renderer,title=component.get_species())
    else:
        print(DNA_construct_obj)
