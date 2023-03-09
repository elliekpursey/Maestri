from ete3 import NCBITaxa, Tree, NodeStyle, AttrFace, TreeStyle
import pandas as pd
from reportlab.lib.units import cm


# create tree of topology with system types for each genus

ncbi = NCBITaxa()

species_taxids = pd.read_csv("../results/species_accessions_refseq_bacteria.csv", header=0)
taxid_list = list(species_taxids.species_taxid)
tree = ncbi.get_topology(taxid_list)
ncbi.annotate_tree(tree, taxid_attr="name")

def my_layout(node):
    if getattr(node, "rank", None):
        if node.is_leaf():
            sciname_face = AttrFace("sci_name", fsize=7, fgcolor="black")
            node.add_face(sciname_face, column=0, position="branch-bottom")
            
nstyle = NodeStyle()
nstyle["shape"] = "square"
nstyle["size"] = 3
nstyle["fgcolor"] = "black"

for n in tree.traverse():
   n.set_style(nstyle)

ts = TreeStyle()
ts.layout_fn = my_layout
ts.mode = "c"
ts.show_leaf_name = False
ts.scale =  120 # 120 pixels per branch length unit
ts.branch_vertical_margin = 25
ts.show_scale = False

tree.render("../results/Maestri_all_species_unannotated.svg", w=1000, tree_style=ts)

