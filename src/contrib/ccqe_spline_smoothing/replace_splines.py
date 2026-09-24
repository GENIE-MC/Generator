import xml.etree.ElementTree as ET
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description="Replaces keys in one GENIE spline file with those from another GENIE spline file. Only checks spline names.")
parser.add_argument('-b', '--base', type=str, required=True, help="Base file to make changes to.")
parser.add_argument('-r', '--revision', type=str, required=True, help="Revision file to take changes from.")
parser.add_argument('-o', '--outpath', type=str, default="./output.xml", help="Path to output file.")
args = parser.parse_args()

base_tree = ET.parse(Path(args.base))
rev_tree  = ET.parse(Path(args.revision))

root_base = base_tree.getroot()
root_rev  = rev_tree.getroot()

# assert base and revision have same tune
base_tune = root_base.findall("./genie_tune")[0].get('name')
rev_tune = root_rev.findall("./genie_tune")[0].get('name')

assert( base_tune == rev_tune )

splines_base = {}
splines_rev  = {}

# A dictionary of models for CCQE.
# Perhaps there is a simpler way of extracting these from GENIE...
models = {
    "AR23_20i_00_000": "genie::NievesQELCCPXSec/ZExp", # same for G18_10(i, j) and G24_20* tunes
    "AR23_20i_01_000": "genie::NievesQELCCPXSec/ZExp_minerva",
    "AR23_20i_01_001": "genie::NievesQELCCPXSec/ZExp_minervaNature",
    "AR23_20i_02_000": "genie::NievesQELCCPXSec/ZExp_lqcd",
    "AR25_20i_00_000": "genie::NievesQELCCPXSec/ZExpNoRPA",
    "AR25_20i_01_000": "genie::NievesQELCCPXSec/ZExpNoRPA_minerva",
    "AR25_20i_01_001": "genie::NievesQELCCPXSec/ZExpNoRPA_minervaNature",
    "AR25_20i_02_000": "genie::NievesQELCCPXSec/ZExpNoRPA_lqcd",
    "G18_01a_00_000" : "genie::LwlynSmithQELCCPXSec/Dipole", # same for all G18_(01, 02)* and G18_10(a, b, c, d) tunes
    "G21_11a_00_000" : "genie::HybridXSecAlgorithm/SuSAv2-QEL" # same for all G21_11* tunes
}

'''
# pick out the appropriate configuration
desired_conf = None
if "G18_01" in base_tune or "G18_02" in base_tune or \
   ("G18_10" in base_tune and base_tune[6] in ['a', 'b', 'c', 'd']):
    desired_conf = models["G18_01a_00_000"]
elif "G18_10" in base_tune or "G24_20" in base_tune or base_tune == "AR23_20i_00_000":
    desired_conf = models["AR23_20i_00_000"]
elif "G21_11" in base_tune:
    desired_conf = models["G21_11a_00_000"]
else:
    desired_conf = models[base_tune]
'''

# Base splines
for spline in root_base.findall(".//spline"):
    name = spline.get('name')
    '''
    if('Weak[CC],QES' in name and 'charm' not in name and 'strange' not in name):
        # Split configuration out
        namebits = name.split('/')
        name = f"{namebits[0]}/{desired_conf}"
        for nb in namebits[2:]:
            name = f"{name}/{nb}"
        if(desired_conf == "CRPA-QEL"):
            name = "genie::HybridXSecAlgorithm/CRPA-QEL"
            for nb in namebits[2:]:
                name = f"{name}/{nb}"
    '''
    splines_base[name] = spline

# Revision splines
for spline in root_rev.findall(".//spline"):
    name = spline.get('name')
    '''
    if('Weak[CC],QES' in name and 'charm' not in name and 'strange' not in name):
        name = desired_conf
    '''
    splines_rev[name] = spline

print("Found", len(splines_base), "splines in", args.base)
print("Found", len(splines_rev), "splines in", args.revision)

for name, rsp in splines_rev.items():
    try:
        bsp = splines_base[name]
        # also replace the spline name
        old_name = bsp.get("name")
        bsp.attrib["name"] = name
        old_knots = list(bsp.findall("knot"))
        for knot in old_knots:
            bsp.remove(knot)
        # Copy knots
        for knot in rsp:
            bsp.append(knot)
    except:
        continue

base_tree.write(Path(args.outpath), encoding="ISO-8859-1", xml_declaration=True)
print("Saved output XML to", args.outpath)
