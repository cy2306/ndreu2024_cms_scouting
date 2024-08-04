"""
For each range of average dimuon masses (mass_region), plot 2d histogram 
of four-muon charge (qfour) vs mass asymmetry (masym) to show ABCD regions.

Variables needed:
m_asymm = abs((m1-m2) / (m1+m2)),
q_four = q_four = abs(q1) + abs(q2),

where the indices 1,2 indicate dimuons in a pair 
(not in any particular order).

Examples of mass regions: 
    "JPsi": events with dimuon masses of 3.05-3.15 GeV.
    "control1": events with dimuon masses of 2.5-3.05 GeV.
    "control2": events with dimuon masses of 3.15-4.0 GeV.

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


2. abcd_regions_plot(mass_region)
    For a specified mass region, plot 2d histogram of masym vs qfour, with
    ABCD regions outlined and labeled.
    Print and save to pdf file.
    Serves as the target for each concurrent.futures process.

    mass_region: str
    Name to specify the range of dimuon masses considered
    (e.g. "JPsi", "control1").

    Return: None


############################ Main function ############################
List all ranges of dimuon masses (mass_region) to plot.
For each mass region, call abcd_regions_plot(mass_region) to create a 
parallel process to plot the ABCD regions in qfour vs masym.

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


def abcd_regions_plot(mass_region):

    # Set up 2d histogram for qfour vs masym.
    nbins_masym_temp, abins_masym_temp = percent_bins(1e-2, 0.5, 0.1)
    lbins_masym = list(abins_masym_temp)
    lbins_masym.append(0.05)
    lbins_masym.append(0.075)
    lbins_masym.append(0.5)
    lbins_masym.sort()
    abins_masym = array(lbins_masym)
    nbins_masym = len(abins_masym) - 1

    abins_qfour = array([-0.0001, 2, 4, 6])
    nbins_qfour= len(abins_qfour) - 1

    H2d = TH2F("th2f_"+mass_region, "",
               nbins_masym, abins_masym, nbins_qfour, abins_qfour)
    
    # Fill in histogram.
    for region_abcd in ["A", "B", "C", "D"]:
        file_name = "dmasym_qfour_new_{}_{}.json".\
            format(mass_region, region_abcd)
        
        with open(file_name, "r") as file:
            dvar = json.load(file)

        list_masym = "l{}_masym_{}".format(mass_region, region_abcd)
        list_qfour = "l{}_qfour_{}".format(mass_region, region_abcd)
        
        for i in range(len(dvar[list_masym])):
            masym = dvar[list_masym][i]
            qfour = dvar[list_qfour][i]

            H2d.Fill(masym, qfour)
    

    # Draw and print histogram.
    H2d_clone = H2d.Clone("th2f_"+mass_region+"_clone")
    H2d_clone.Scale(1.0, "width")

    C = TCanvas("c_"+mass_region, "", 800, 650)
    P_hist = TPad("pad_"+mass_region, "", 0, 0, 1, 1)
    P_hist.SetTopMargin(0.05)
    P_hist.SetBottomMargin(0.2)
    P_hist.SetLeftMargin(0.15)
    P_hist.SetRightMargin(0.15)
    P_hist.Draw()

    P_hist.cd()
    P_hist.SetLogx()
    gStyle.SetPalette(kBird)

    H2d_clone.SetStats(0)
    H2d_clone.GetXaxis().SetTitle("m_{asym} (GeV/c^{2})")
    H2d_clone.GetXaxis().SetTitleSize(0.06)
    H2d_clone.GetXaxis().SetTitleOffset(1.2)
    H2d_clone.GetXaxis().SetLabelSize(0.05)
    H2d_clone.GetYaxis().SetTitle("q_{4#mu}")
    H2d_clone.GetYaxis().SetTitleSize(0.07)
    H2d_clone.GetYaxis().SetTitleOffset(0.8)
    H2d_clone.GetYaxis().SetNdivisions(306)
    H2d_clone.GetYaxis().SetLabelSize(0.05)

    H2d_clone.Draw("colz")

    # Draw lines to divide ABCD regions.
    xmin = H2d_clone.GetXaxis().GetXmin()
    xmax = H2d_clone.GetXaxis().GetXmax()

    line_m1 = TLine(0.05, 0., 0.05, 6.)     # masym = 0.05
    line_m1.SetLineWidth(2)
    line_m1.Draw()
    
    line_m2 = TLine(0.075, 0., 0.075, 6.)   # masym = 0.075
    line_m2.SetLineWidth(2)
    line_m2.Draw()

    line_m3 = TLine(0.5, 0., 0.5, 6.)       # masym = 0.5
    line_m3.SetLineWidth(2)
    line_m3.Draw()

    line_q = TLine(xmin, 2., xmax, 2.)      # qfour = 2
    line_q.SetLineWidth(2)
    line_q.Draw()

    # Draw text labels.
    T_legend = TPaveText(0.85, 0.12, 0.9, 0.22, "NDC")
    T_legend.AddText("Events / bin area")
    T_legend.SetTextFont(42)
    T_legend.SetTextSize(0.04)
    T_legend.SetFillColor(0)
    T_legend.SetFillStyle(0)
    T_legend.SetBorderSize(0)
    T_legend.Draw()

    T_regionA = TPaveText(0.25, 0.3, 0.35, 0.4, "NDC")
    T_regionA.AddText("A")
    T_regionA.SetTextSize(0.08)
    T_regionA.SetFillColor(0)
    T_regionA.SetFillStyle(0)
    T_regionA.SetBorderSize(0)
    T_regionA.Draw()
    
    T_regionB = TPaveText(0.25, 0.65, 0.35, 0.75, "NDC")
    T_regionB.AddText("B")
    T_regionB.SetTextSize(0.08)
    T_regionB.SetFillColor(0)
    T_regionB.SetFillStyle(0)
    T_regionB.SetBorderSize(0)
    T_regionB.Draw()

    T_regionC = TPaveText(0.65, 0.3, 0.75, 0.4, "NDC")
    T_regionC.AddText("C")
    T_regionC.SetTextSize(0.08)
    T_regionC.SetFillColor(0)
    T_regionC.SetFillStyle(0)
    T_regionC.SetBorderSize(0)
    T_regionC.Draw()
    
    T_regionD = TPaveText(0.65, 0.65, 0.75, 0.75, "NDC")
    T_regionD.AddText("D")
    T_regionD.SetTextSize(0.08)
    T_regionD.SetFillColor(0)
    T_regionD.SetFillStyle(0)
    T_regionD.SetBorderSize(0)
    T_regionD.Draw()


    C.Print("fig_masym_qfour_new_{}.pdf".\
            format(mass_region))
    P_hist.Clear()
    P_hist.Update()
    C.Clear()
    C.Update()
    C.Close()



if __name__ == "__main__":

    print("################ Main function started ################")

    lmass_regions = [
        "all",
        "JPsi"
    ]

    # For each mass region, create parallel process to plot ABCD regions.
    n = len(lmass_regions)
    gROOT.SetBatch(True)
    lfutures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers = n) \
    as executor:
        for mass_region in lmass_regions:
            future = executor.submit(
                        abcd_regions_plot, mass_region)
            
            lfutures.append(future)
    

    print("################ Main function ended ################")