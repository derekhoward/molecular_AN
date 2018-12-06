#if you get rpy2 errors in pycharm, try running it from the command line
import data_processing as data
import svg_utils
from pathlib import Path


# define the SVGs files you want to modify
adult_svgs = ['1532_112360908.svg', '1823_112363270.svg', '1297_112282739.svg']
fetal_svgs = ['0893_101892619.svg', '1097_101892615.svg', '1352_101892610.svg']
fetal_brainstem_svgs = ['0391_102182817.svg', '0639_102182810.svg']
human_diagram = 'human_diagram.svg'

# input directories
svg_dir = Path('./data/svg')
adult_dir = svg_dir / 'slices' / 'adult'
fetal_dir = svg_dir / 'slices' / 'fetal21'
fetal_brainstem_dir = svg_dir / 'slices' / 'fetal21_brainstem'

# define output directory
figures_dir = Path('./figures')
figures_dir.mkdir(exist_ok=True)

# get data
adult_exp = data.get_dataset('adult', 'reannotator')
fetal_exp = data.get_dataset('fetal', 'reannotator')
negraes = data.get_genelist('negraes')

# create tables to match AUC values to structures
adult_lookup = svg_utils.create_auc_lookup(exp_df=adult_exp, gene_list=negraes, ontology='adult')
fetal_lookup = svg_utils.create_auc_lookup(exp_df=fetal_exp, gene_list=negraes, ontology='fetal')

adult_lookup = adult_lookup.rename(index=str, columns={"AUROC": "AUC"})
fetal_lookup = fetal_lookup.rename(index=str, columns={"AUROC": "AUC"})


svg_utils.modify_svg(svg_dir / human_diagram, figures_dir / human_diagram, graph_id='adult', lookup_table=adult_lookup)

for adult_svg in adult_svgs:
    #save a cleaned up version of original svg in figures folder
    svg_utils.clean_SVG(adult_dir/adult_svg, figures_dir/adult_svg)
    output = f'colourized_{adult_svg}'
    svg_utils.modify_svg(figures_dir/adult_svg, figures_dir/output, graph_id='adult', lookup_table=adult_lookup)

for fetal_svg in fetal_svgs:
    svg_utils.clean_SVG(fetal_dir/fetal_svg, figures_dir/fetal_svg)
    output = f'colourized_{fetal_svg}'
    svg_utils.modify_svg(figures_dir / fetal_svg, figures_dir / output, graph_id='fetal', lookup_table=fetal_lookup)

for brainstem_svg in fetal_brainstem_svgs:
    print(f'modifying brainstem svg: {brainstem_svg}')
    svg_utils.clean_SVG(fetal_brainstem_dir/brainstem_svg, figures_dir/brainstem_svg)
    output = f'colourized_{brainstem_svg}'
    svg_utils.modify_svg(figures_dir/brainstem_svg, figures_dir/output, graph_id='fetal', lookup_table=fetal_lookup)
