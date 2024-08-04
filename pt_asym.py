"""
Plot 2d histogram of pt_asym vs pt_avg for ABCD regions 
of double J/Psi events.

Variables needed:
pt1: pT of the leading J/Psi particle
pt2: pT of the second J/Psi particle

Variables plotted:
pt_asym = (pt1 - pt2) / (pt1 + pt2)
pt_avg = (pt1 + pt2) / 2


For each ABCD region:
1. Read json file of the specified mass region and ABCD region.
    json files can be made using bkgest_abcd_sort.py.
    
    Examples of mass regions: 
    "JPsi": events with dimuon masses of 3.05-3.15 GeV.
    "control1": events with dimuon masses of 2.5-3.05 GeV.
    "control2": events with dimuon masses of 3.15-4.0 GeV.

2. Calculate pt_asym and pt_avg. Plot 2d histogram.

3. Draw and save to pdf file.

"""

"""
############################### Functions ###############################
1. percent_bins(a,b,p)
    Set up histogram bin edges such that each bin width is p percent
    of the central bin value.

    a: float
    first bin edge (min)

    b: float
    last bin edge (max)

    p: float
    percent for bin width, decimal in range [0, 1]

    Return: duple
    (number of bins, array of bin edges)


2. th2f_asym_create_draw(mass_region, region_abcd, 
                        bins_pt_avg, bins_pt_asym,
                        mass_region_formatted)
    
    For a specified mass region and ABCD region, create and draw 
    2d histogram of pt_asym vs pt_avg.
    Print and save to pdf file.
    Serves as the target for each concurrent.futures process.

    mass_region: str
    Name to specify the range of dimuon masses considered
    (e.g. "JPsi", "control1").

    region_abcd: str ("A", "B", "C", "D", or "ABCD")
    Name to specify the ABCD region.

    bins_pt_avg: list [number of bins, array of bin edges]
    Bins for the pt_avg axis.

    bins_pt_asym: list [number of bins, array of bin edges]
    Bins for the pt_asym axis.

    mass_region_formatted: str
    Formatted name (including capitalization and LaTeX) of the
    mass region to be displayed in histogram title.

    Return: None


############################ Main function ############################
Specify the range of dimuon masses (mass_region).
List all ABCD regions and formatted names to be used in histograms.
Set up bins for pt_asym and pt_avg.

For each ABCD region, call th2f_asym_create_draw(...) to create a parallel 
process to plot pt_asym vs pt_avg.

"""

from ROOT import *
from numpy import array, arange
import json
import concurrent.futures


def percent_bins(a,b,p):
        

    mult = 1 + (p / 2.)
    lbins = []

    i = 0
    x = a
    while x < b:
        if i % 2 == 0:      # Even number iterations.
            lbins.append(x)
        x = x * mult
        i += 1

    return (len(lbins) -1, array(lbins, float))


def th2f_asym_create_draw(mass_region, region_abcd, 
                        bins_pt_avg, bins_pt_asym,
                        mass_region_formatted):
    
    
    # Set up 2d histogram.
    nbins_pt_avg = bins_pt_avg[0]
    abins_pt_avg = bins_pt_avg[1]
    nbins_pt_asym = bins_pt_asym[0]
    abins_pt_asym = bins_pt_asym[1]
    histname = "th2f_{}_pt_asym_vs_pt_avg_{}".\
        format(mass_region, region_abcd)
    H2d = TH2F(histname, "", 
               nbins_pt_avg, abins_pt_avg, nbins_pt_asym, abins_pt_asym)


    # Fill in histogram.
    file_name = "dvar_doubleJPsi_new_{}_{}.json".\
        format(mass_region, region_abcd)

    with open(file_name, "r") as file:
        dvar = json.load(file)

    lpt1 = dvar["l{}_pt_dimlead_{}".format(mass_region, region_abcd)]
    lpt2 = dvar["l{}_pt_dimlast_{}".format(mass_region, region_abcd)]

    for event in range(len(lpt1)):
        pt1 = lpt1[event]
        pt2 = lpt2[event]
        pt_avg = float((pt1 + pt2) / 2)
        pt_asym = float((pt1 - pt2) / (pt1 + pt2))
        H2d.Fill(pt_avg, pt_asym)
    

    # Draw and print histogram.
    H2d_clone = H2d.Clone(histname+"_clone")
    H2d_clone.Scale(1.0, "width")

    C = TCanvas("c_"+histname)
    P_hist = TPad("pad_"+histname, "", 0, 0, 1, 1)
    P_hist.SetTopMargin(0.13)
    P_hist.SetBottomMargin(0.22)
    P_hist.SetLeftMargin(0.17)
    P_hist.SetRightMargin(0.15)
    P_hist.Draw()

    C.cd()
    if region_abcd == "ABCD":
         region_abcd_formatted = "all regions ABCD"
    else:
         region_abcd_formatted = "region "+region_abcd
    
    """ title_line1 = "{} asymmetry of J/#psi pairs for {} {}".\
                    format("p_{T}", region_abcd_formatted, 
                           mass_region_formatted) """
    title_line1 = "{} asymmetry of double J/#psi events".\
                    format("p_{T}")
    
    T_line1 = TPaveText(0.4, 0.92, 0.7, 0.95, "NDC")
    T_line1.AddText(title_line1)
    T_line1.SetTextSize(0.05)
    T_line1.SetFillColor(0)
    T_line1.SetBorderSize(0)
    T_line1.Draw()

    T_line2 = TPaveText(0.85, 0.12, 0.9, 0.22, "NDC")
    T_line2.AddText("Events / bin area")
    T_line2.SetTextFont(42)
    T_line2.SetTextSize(0.04)
    T_line2.SetFillColor(0)
    T_line2.SetBorderSize(0)
    T_line2.Draw()

    P_hist.cd()
    P_hist.SetLogx()
    gStyle.SetPalette(kLake)

    H2d_clone.SetStats(0)
    H2d_clone.GetXaxis().SetTitle("#hat{p}_{T,J/#psi}")
    H2d_clone.GetXaxis().SetTitleSize(0.06)
    H2d_clone.GetXaxis().SetTitleOffset(1.5)
    H2d_clone.GetXaxis().SetLabelSize(0.05)
    H2d_clone.GetYaxis().SetTitle("A_{2J/#psi}^{p_{T}}")
    H2d_clone.GetYaxis().SetTitleSize(0.06)
    H2d_clone.GetYaxis().SetTitleOffset(1.2)
    H2d_clone.GetYaxis().SetLabelSize(0.05)

    H2d_clone.Draw("colz")

    C.Print("fig_pt_asym_doubleJPsi_new_{}_{}.pdf".\
            format(mass_region, region_abcd))
    P_hist.Clear()
    P_hist.Update()
    C.Clear()
    C.Update()
    C.Close()
         


if __name__ == "__main__":

    print("################ Main function started ################")

    #mass_region = "all"
    #mass_region = "control1"
    mass_region = "JPsi"
    #mass_region = "control2"

    lregions_abcd = [
        "A",
        #"B",
        #"C",
        #"D",
        #"ABCD"
    ]

    print("mass region: {}".format(mass_region))
    print("ABCD regions: {}".format(lregions_abcd))

    bins_pt_avg = list(percent_bins(1., 1e4, 0.3))
    #bins_pt_asym = list(percent_bins(1., 1e5, 0.3))

    abins_pt_asym = arange(0., 1.05, 0.05)
    nbins_pt_asym = len(abins_pt_asym) - 1
    bins_pt_asym = [nbins_pt_asym, abins_pt_asym]

    dmass_regions_formatted = {
        "all": "",
        "JPsi": "(J/#psi)",
        "control1": "(control 1)",
        "control2": "(control 2)",
    }

    # For each ABCD region, create parallel process to make plots.
    n = len(lregions_abcd)
    gROOT.SetBatch(True)
    lfutures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers = n) \
    as executor:
        for region_abcd in lregions_abcd:
            mass_region_formatted = dmass_regions_formatted[mass_region]

            future = executor.submit(
                        th2f_asym_create_draw, mass_region, region_abcd, 
                        bins_pt_avg, bins_pt_asym, mass_region_formatted)
            
            lfutures.append(future)


    print("################ Main function ended ################")

