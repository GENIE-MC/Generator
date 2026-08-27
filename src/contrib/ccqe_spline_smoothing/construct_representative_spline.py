import argparse
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel

from pathlib import Path

# For conversion from number --> 1e-38 cm2
# * CONV_CONST: number --> 1e-38cm2, / CONV_CONST: 1e-38cm2 --> number
hbarc = 1.973269804e-16
meter = 1.0 / hbarc
m2    = meter*meter
cm2   = 1.0e-4 * m2
CONV_CONST = 1.0e+38 / cm2

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

# Custom formatter for XML tags to be "genie-like"
def format_text(elem):
    if elem.tag in ["E", "xsec"]:
        raw = (elem.text or "").strip()
        val  = float(raw)
        #width = 9 if elem.tag == "E" else 22
        #elem.text = f" {raw:<{width}} "
        if elem.tag == "E":
            
            elem.text = f" {val:10.5f} "
        elif elem.tag == "xsec":
            s    = f"{val:.10e}"
            m, e = s.split("e")
            ei   = int(e)
            s    = f"{m}e{ei:+03d}"
            elem.text = f" {s:<16} "
    for child in elem:
        format_text(child)
        
def genie_indent(elem, level=0):
    indent = "  "
    i = "\n" + level * indent
    child = "\n" + (level+1) * indent
    children = list(elem)
    if not children:
        return
        
    if elem.tag in ["knot"]:
        elem.text = " "
        for kid in children:
            genie_indent(kid, level+1)
            kid.tail = " "
    else:
        elem.text = child
        for kid in children:
            genie_indent(kid, level+1)
        for kid in children[:-1]:
            kid.tail = child
        children[-1].tail = i

def main(args):

    df = pd.read_parquet(Path(args.input).resolve())
    idx = df.index.to_numpy()

    mean, std, err = np.zeros_like( idx ), np.zeros_like( idx ), np.zeros_like( idx )
    for i, pt in enumerate(idx):
        mean[i] = np.mean(df.loc[pt].to_numpy())
        std [i] = np.std (df.loc[pt].to_numpy())
    err = std / np.sqrt( float(df.iloc[0].to_numpy().shape[0]) )

    zerodims = mean[mean == 0.0].shape[0]
    idx  = idx [zerodims:]
    mean = mean[zerodims:]
    std  = std [zerodims:]
    err  = err [zerodims:]

    # The points gmkspl gives you are equidistant in log space, we'll take advantage of that
    # Support is only in idx[:]
    lidx  = np.log(idx)
    ckernel = ConstantKernel(1.0) * RBF(length_scale = 0.5)
    cgp = GaussianProcessRegressor(kernel=ckernel, alpha=err**2, normalize_y=True)
    cgp.fit(lidx.reshape(-1, 1), mean)

    resp, cov = cgp.predict(lidx.reshape(-1, 1), return_cov=True)
    eresp = np.sqrt(np.diag(cov))

    out_mean, out_err, out_resp, out_eresp = (np.zeros_like(df.index.to_numpy()) for _ in range(4))
    out_mean [zerodims:] = np.clip(mean,  0.0, None)
    out_err  [zerodims:] = np.clip(err,   0.0, None)
    out_resp [zerodims:] = np.clip(resp,  0.0, None)
    out_eresp[zerodims:] = np.clip(eresp, 0.0, None)

    # Write the output to a parquet and an xml
    df_out = pd.DataFrame(index = df.index)
    df_out["Smoothed xsec [1e-38 cm2/A]"] = out_resp
    df_out["Error on smoothed xsec"] = out_eresp
    df_out[f"Mean xsec from {len(df.columns)} nucleons [1e-38 cm2/A]"] = out_mean
    df_out["Error on mean xsec"] = out_err
    df_out.to_parquet(Path(f"{args.output}.parquet").resolve())

    # careful which tune was given...
    model_choice = None
    if "G18_01" in args.tune or "G18_02" in args.tune or \
       ("G18_10" in args.tune and args.tune[6] in ['a', 'b', 'c', 'd']):
        model_choice = models["G18_01a_00_000"]
    elif "G18_10" in args.tune or "G24_20" in args.tune or args.tune == "AR23_20i_00_000":
        model_choice = models["AR23_20i_00_000"]
    elif "G21_11" in args.tune:
        model_choice = models["G21_11a_00_000"]
    else:
        model_choice = models[args.tune]

    with open(Path(f"{args.output}.xml").resolve(), 'wb') as fxout:
        root = ET.Element("genie_xsec_spline_list")
        tree = ET.ElementTree(root)
        root.set("version", "3.00")
        root.set("uselog", "1")
        tune = ET.SubElement(root, "genie_tune")
        tune.set("name", f"{args.tune}")
        spl  = ET.SubElement(tune, "spline")
        nucpdg = 2112 if int(args.nupdg) > 0 else 2212
        spl.set("name", f"{model_choice}/nu:{args.nupdg};tgt:{args.tgtpdg};N:{nucpdg};proc:Weak[CC],QES;")
        lknots = len(df.index)
        spl.set("nknots", str(lknots))
        for Enu, xs in zip(df.index, out_resp):
            knot = ET.SubElement(spl, "knot")
            Eknot = ET.SubElement(knot, "E")
            Eknot.text = f"{Enu}"
            Xknot = ET.SubElement(knot, "xsec")
            Xknot.text = f"{xs / CONV_CONST}"
        #ET.indent(tree)
        format_text(root)
        genie_indent(root)
        tree.write(fxout, encoding="utf-8", xml_declaration=True)

    print(f"Saved outputs: {args.output}.parquet, {args.output}.xml")
    
if __name__ == "__main__":    
    parser = argparse.ArgumentParser(
        description = "From a single parquet file produced from the concatenation of many, many nucleon throws for CCQE splines, I will construct a single representative for that parquet and save it as a new parquet. I will also save a single XML file with that spline. Pass the input parquet and the tune used to run it, as well as the neutrino and target used."
    )
    parser.add_argument('-i', "--input", type=str, required=True, help="Input parquet.")
    parser.add_argument('-o', "--output", type=str, required=True, help="Output basename (no extension).")
    parser.add_argument('-t', "--tune", type=str, required=True, help="GENIE tune used.")
    parser.add_argument('-p', "--nupdg", type=int, required=True, help="Neutrino PDG code.")
    parser.add_argument('-n', "--tgtpdg", type=int, required=True, help="Target PDG code.")
    args = parser.parse_args()

    main(args)
