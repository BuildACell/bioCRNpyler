# pathutil.py - path utilities
# ASS, 21 Apr 2020
#
# This file contains some utility functions for plotting CRNs
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import networkx as nx
import statistics
from bokeh.models import (BoxSelectTool, Circle,Square, EdgesAndLinkedNodes, HoverTool,
                          MultiLine, NodesAndLinkedEdges, Plot, Range1d, TapTool,PanTool,WheelZoomTool)
from bokeh.palettes import Spectral4
from bokeh.models.graphs import from_networkx
from fa2 import ForceAtlas2
import numpy as np

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
        from_x = positions[edge[0]][0]
        from_y = positions[edge[0]][1]
        to_x = positions[edge[1]][0]
        to_y = positions[edge[1]][1]
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

def graphPlot(DG,DGspecies,DGreactions,plot,layout="force"):
    """given a directed graph, plot it!"""
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
                            scalingRatio=2.0,
                            strongGravityMode=False,
                            gravity=1.0,

                            # Log
                            verbose=False)

        positions = forceatlas2.forceatlas2_networkx_layout(DG, pos=None, iterations=2000) 
    elif(layout == "circle"):
        #here we would arrange the nodes in a circle
        pass
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
    node_hover_tool = HoverTool(tooltips=[("name", "@species"), ("type", "@type")],\
                                        renderers=[reaction_renderer,species_renderer])
    plot.add_tools(edge_hover_tool,node_hover_tool, TapTool(), BoxSelectTool(),PanTool(),WheelZoomTool())

    edges_renderer.selection_policy = NodesAndLinkedEdges()
    edges_renderer.inspection_policy = EdgesAndLinkedNodes()

    plot.renderers.append(edges_renderer)
    plot.renderers.append(reaction_renderer)
    plot.renderers.append(species_renderer)