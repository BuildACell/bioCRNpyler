from biocrnpyler import *
# import bioscrape
import numpy as np
import bokeh.io
import bokeh.plotting
bokeh.io.output_notebook()
from bokeh.themes import Theme
from bokeh.layouts import row

# Modules needed from Bokeh.
from bokeh.io import output_file, show
from bokeh.plotting import gridplot,figure
from bokeh.models import LinearAxis, Range1d
color=bokeh.palettes.Accent[6]
    
from bokeh.io import export_png

# Function with some standard Bokeh plot settings
def create_custom_plot(title_text, x_max=8,y_max=2, yname=None):
    custom_plot = figure(
        toolbar_location='right',
        outline_line_color=None,
        min_border_right=10,
        height=400,
        width=400,
    )

    custom_plot.title.text = title_text
    custom_plot.xaxis.axis_label = 'Time (hours)'
    custom_plot.yaxis.axis_label = yname
    custom_plot.y_range = Range1d(0, y_max)
    custom_plot.x_range = Range1d(0, x_max)
    custom_plot.outline_line_color = None

    # custom_plot.yaxis
    custom_plot.ygrid.visible = False
    custom_plot.yaxis.axis_label_text_font_size = '15pt'
    custom_plot.yaxis.major_label_text_font_size = '15pt'
    custom_plot.yaxis.major_label_text_font = 'Work Sans'
    custom_plot.yaxis.axis_label_standoff = 15
    custom_plot.yaxis.axis_label_text_font_style = 'normal'

    # custom_plot.xaxis
    custom_plot.xgrid.visible = False
    custom_plot.xaxis.axis_label_text_font_size = '15pt'
    custom_plot.xaxis.major_label_text_font_size = '15pt'
    custom_plot.xaxis.major_label_text_font = 'Work Sans'
    custom_plot.xaxis.axis_label_standoff = 15
    custom_plot.xaxis.axis_label_text_font_style = 'normal'

    # custom_plot.title
    custom_plot.title.text_font_size = '18pt'
    custom_plot.title.align = 'left'
    custom_plot.title.offset = -50.0

    return custom_plot

alphaHL_monomer = IntegralMembraneProtein('alphaHL_monomer', product='alphaHL', size=7)

#Mechanisms
mech_int = Membrane_Protein_Integration() 

alphaHL_channel=MembraneChannel(alphaHL_monomer.product, substrate='ATP')

#Mechanisms
mech_tra = Simple_Transport()
compartment_internal = alphaHL_monomer.get_species().compartment
#ActivatedPromoter
activator = Species("T7RNAP",
                    material_type = "small_molecule",
                    compartment=compartment_internal)

#Create a custom set of parameters
hill_parameters = {"k":1.0, "n":4, "K":20, "kleak":0.0001}
P_activatable = ActivatablePromoter("P_activtable", activator = activator, leak = False,
                                    parameters = hill_parameters)

#Create a DNA assembly "reporter" with P_activatable for its promoter
activatable_assembly = DNAassembly("activatable_assembly", promoter = P_activatable, rbs = "rbs", initial_concentration = 1*10**-3, protein= alphaHL_monomer.membrane_protein)
all_mechanisms = {mech_tra.mechanism_type:mech_tra, mech_int.mechanism_type:mech_int}
E = EnergyTxTlExtract(components=[activatable_assembly, alphaHL_monomer, alphaHL_channel],
                      mechanisms = all_mechanisms,
                      parameter_file = "all_parameters.txt")
CRN = E.compile_crn(compartment=compartment_internal)
# print(CRN.pretty_print())
CRN.write_sbml_file('test2.xml')

# Simulate CRN
try:
    import bioscrape
    import bokeh
except ModuleNotFoundError:
    print('please install the plotting libraries: pip install biocrnpyler[all]')
else:
    from biocrnpyler import *
    import numpy as np
    import pandas as pd
    maxtime = 300000
    timepoints = np.arange(0, maxtime, 100)
    
    #Inital conditions
    # x0_dict= {'small_molecule_T7RNAP':120, alphaHL_channel.product:5,}
    x0_dict= {'small_molecule_T7RNAP_Internal':1.20, alphaHL_channel.product:4.5,
              alphaHL_channel.substrate:0.5}

    #Run Simulation
    R = CRN.simulate_with_bioscrape_via_sbml(timepoints = timepoints,
                                             initial_condition_dict = x0_dict)

    #Plot alpha-hemolysin protien and the multiple configurations
    p2 = create_custom_plot("alpha-hemolysin Expression", x_max=8,y_max= .014, yname='monomer (uM)')
    p2.line(timepoints/3600,  R['protein_alphaHL_monomer_Internal'], line_width = 3, line_alpha=.5, color=color[1], legend_label='monomer')
    p2.line(timepoints/3600,  R['complex_protein_alphaHL_monomer_7x_Internal_'], line_width = 3, line_alpha=.5, color=color[4], legend_label='polymer')
    p2.line(timepoints/3600,  R['protein_alphaHL_Passive_Internal'], line_width = 3, line_alpha=.5, color=color[5], legend_label='channel')
    
    #Plot the transport of substrate ATP
    p3 = create_custom_plot("Passive Transport of ATP", x_max=18,y_max= 5, yname='ATP (uM)')
    p3.line(timepoints/3600,  R['ATP_Internal'], line_width = 3, legend_label = "Internal", color=color[0])
    p3.line(timepoints/3600,  R['ATP_External'], line_width = 3, legend_label = "External", color=color[2], line_dash = 'dashed')
    
    # Arrange the plots in a row layout
    layout = row(p2, p3)
    # save as png
    bokeh.io.export_png(layout, filename="transport_models.png")