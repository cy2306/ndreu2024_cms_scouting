"""
Plot 1d histogram of eta_asym vs abs(eta_avg) for ABCD regions 
of double J/Psi (dimuon) events.

Variables needed:
1. deta = eta1 - eta2,
2. eta_avg = (eta1 + eta2) / 2,
    where eta1 corresponds to the leading-pT J/Psi.

(Optional for additional plots)
3. eta1
4. eta2


Calculations:
1. N_tot = N(deta > 0) + N(deta < 0),
    where N is the number of events in a bin.

2. eta_asym = (N(deta > 0) - N(deta < 0)) / N_tot


For each ABCD region:
1. Read json file of the specified mass region and ABCD region.
    json files can be made using dimuon_pairs_abcd_sort.py.
    json files can be made using dimuon_pairs_abcd_sort.py.
    
    Examples of mass regions: 
    "JPsi": events with dimuon masses of 3.05-3.15 GeV.
    "control1": events with dimuon masses of 2.5-3.05 GeV.
    "control2": events with dimuon masses of 3.15-4.0 GeV.

2. (Optional) Plot 1d histograms of eta1, eta2, eta_avg, deta.

3. Plot N(deta > 0) - N(deta < 0) as a function of abs(eta_avg).

4. Plot N_tot = N(deta > 0) + N(deta < 0) as a function of abs(eta_avg).

5. Plot eta_asym as a function of abs(eta_avg) by setting each
    bin content to N(deta > 0) - N(deta < 0) / N_tot.

6. Draw and save plots to pdf files.

"""

"""
############################### Functions ###############################
1. th1f_create(mass_region, region_abcd, var_name, abins)
    Create 1d histogram for a specified mass region, ABCD region, 
    and variable.

    mass_region: str
    Name to specify the range of dimuon masses considered 
    (e.g. "JPsi", "control1").

    region_abcd: str ("A", "B", "C", "D", or "ABCD")
    Name to specify the ABCD region.

    var_name: str ("eta1", "eta2", "eta_avg", or "deta")
    Name to specify the variable.

    abins: array
    Array of bin edges for the histogram.
    (Use bins edges for eta_avg.)

    Return: TH1F
    Filled 1d histogram (not drawn).


2. th1f_draw(mass_region, region_abcd,
              H1d_eta1, H1d_eta2, H1d_eta_avg, H1d_deta,
                mass_region_formatted, dnames_formatted)
    
    Draw 1d histograms made with th1f_create(...).
    4 histograms on the same graph, plotted as colored lines "hist".
    Print and save to pdf file.

    mass_region, region_abcd: str
    Names to specify the mass region and ABCD region.

    H1d_eta1, H1d_eta2, H1d_eta_avg, H1d_deta: TH1F
    Outputs of th1f_create(...).

    mass_region_formatted: str
    Formatted name (including capitalization and LaTeX) of the
    mass region to be displayed in graph title.

    dnames_formatted: dict
    Contains formatted names (str) of variable names to be displayed
    in the legend.
    {
    "eta1": "Formatted string for eta1",
    "eta2": "Formatted string for eta2",
    ...
    }

    Return: None


3. th1f_Ndeta_diff_create(mass_region, region_abcd, abins_eta_avg)
    Create 1d histogram of N(deta > 0) - N(deta < 0) as a function
    of abs(eta_avg).

    mass_region, region_abcd: str
    Names to specify the mass region and ABCD region.

    abins_eta_avg: array
    Array of bin edges for abs(eta_avg).

    Return: TH1F
    Filled 1d histogram (not drawn).


4. th1f_Ndeta_tot_create(mass_region, region_abcd, abins_eta_avg)
    Create 1d histogram for N_tot = N(deta > 0) + N(deta < 0) as
    a function of abs(eta_avg).

    mass_region, region_abcd: str

    abins_eta_avg: array

    Return TH1F


5. th1f_asym_create_draw(mass_region, region_abcd, abins_eta_avg,
                        H1d_Ndeta_diff, H1d_Ndeta_tot,
                            mass_region_formatted)
    
    Create and draw 1d histogram of 
    eta_asym = (N(deta > 0) - N(deta < 0)) / N_tot
    as a function of abs(eta_avg).
    Plot as point with error bars "p" and "e0".
    Print and save to pdf file.

    mass_region, region_abcd: str

    abins_eta_avg: array

    H1d_Ndeta_diff: TH1F
    Output of th1f_Ndeta_diff_create(...).

    H1d_Ndeta_tot: TH1F
    Output of th1f_Ndeta_tot_create(...).

    mass_region_formatted: str
    Formatted name to be displayed in histogram title.

    Return: None


6. eta_asym(mass_region, region_abcd, lnames, dbins,
                mass_region_formatted, dnames_formatted)
    
    Call the above functions in order to plot eta_asym vs abs(eta_avg).
    Serves as the target for each concurrent.futures process.

    Return: None


############################ Main function ############################
Specify the range of dimuon masses (mass_region).
List all ABCD regions, variables ("eta1", "eta2", etc.), bin edges,
and formatted names to be used in histograms.

For each ABCD region, call eta_asym(...) to create a parallel process 
to plot eta_asym vs abs(eta_avg).

"""

from ROOT import *
from numpy import array, argmax, arange
import json
import concurrent.futures


def th1f_create(mass_region, region_abcd, var_name, abins):


    # Set up 1d histogram.
    nbins = len(abins) - 1
    histname = "th1f_{}_{}_{}".format(mass_region, var_name, region_abcd)
    H1d = TH1F(histname, "", nbins, abins)


    # Fill in histogram.
    file_name = "dvar_doubleJPsi_new_{}_{}.json".\
        format(mass_region, region_abcd)

    with open(file_name, "r") as file:
            dvar = json.load(file)

    list_name = "l{}_{}_{}".format(mass_region, var_name, region_abcd)
    for event in dvar[list_name]:
        H1d.Fill(event)
    
    return H1d


def th1f_draw(mass_region, region_abcd,
              H1d_eta1, H1d_eta2, H1d_eta_avg, H1d_deta,
                mass_region_formatted, dnames_formatted):
     

    # Draw and print 1d histogram.
    H1d_eta1_clone = H1d_eta1.Clone(H1d_eta1.GetName()+"_clone")
    H1d_eta2_clone = H1d_eta2.Clone(H1d_eta2.GetName()+"_clone")
    H1d_eta_avg_clone = H1d_eta_avg.Clone(H1d_eta_avg.GetName()+"_clone")
    H1d_deta_clone = H1d_deta.Clone(H1d_deta.GetName()+"_clone")
    
    H1d_eta1_clone.Scale(1.0, "width")
    H1d_eta2_clone.Scale(1.0, "width")
    H1d_eta_avg_clone.Scale(1.0, "width")
    H1d_deta_clone.Scale(1.0, "width")
    
    C = TCanvas("c_"+region_abcd, "", 800, 650)
    P_hist = TPad("pad_eta_"+region_abcd, "", 0, 0, 1, 1)
    P_hist.SetTopMargin(0.17)
    P_hist.SetBottomMargin(0.2)
    P_hist.SetLeftMargin(0.2)
    P_hist.SetRightMargin(0.05)
    P_hist.Draw()

    C.cd()
    if region_abcd == "ABCD":
         region_abcd_formatted = "all regions ABCD"
    else:
         region_abcd_formatted = "region "+region_abcd
    
    title_line1 = "#eta of J/#psi pairs for {} {}".\
        format(region_abcd_formatted, mass_region_formatted)
    
    T_line1 = TPaveText(0.4, 0.92, 0.7, 0.95, "NDC")
    T_line1.AddText(title_line1)
    T_line1.SetTextSize(0.06)
    T_line1.SetFillColor(0)
    T_line1.SetBorderSize(0)
    T_line1.Draw()

    P_hist.cd()

    H1d_eta1_clone.SetLineColor(kPink)
    H1d_eta2_clone.SetLineColor(kGreen+2)
    H1d_eta_avg_clone.SetLineColor(kBlack)
    H1d_deta_clone.SetLineColor(kMagenta)

    H1d_eta1_clone.SetLineStyle(1)         # solid
    H1d_eta2_clone.SetLineStyle(1)         # solid
    H1d_eta_avg_clone.SetLineStyle(3)      # dotted
    H1d_deta_clone.SetLineStyle(2)         # dashed

    H1d_eta1_clone.SetLineWidth(2)
    H1d_eta2_clone.SetLineWidth(2)
    H1d_eta_avg_clone.SetLineWidth(3)
    H1d_deta_clone.SetLineWidth(2)

    aymax_abcd = array([H1d_eta1_clone.GetMaximum(), 
                        H1d_eta2_clone.GetMaximum(),
                        H1d_eta_avg_clone.GetMaximum(), 
                        H1d_deta_clone.GetMaximum()])
    ymax_abcd = aymax_abcd[argmax(aymax_abcd)]

    H1d_eta1_clone.SetStats(0)
    H1d_eta1_clone.GetXaxis().SetTitle("#eta")
    H1d_eta1_clone.GetXaxis().SetTitleSize(0.06)
    H1d_eta1_clone.GetXaxis().SetTitleOffset(1.5)
    H1d_eta1_clone.GetXaxis().SetLabelSize(0.05)
    H1d_eta1_clone.GetYaxis().SetTitle("Events / bin width")
    H1d_eta1_clone.GetYaxis().SetTitleSize(0.06)
    H1d_eta1_clone.GetYaxis().SetTitleOffset(1.5)
    H1d_eta1_clone.GetYaxis().SetLabelSize(0.05)
    H1d_eta1_clone.SetMaximum(1.2 * ymax_abcd)

    H1d_eta1_clone.Draw("hist")
    H1d_eta2_clone.Draw("hist same")
    H1d_eta_avg_clone.Draw("hist same")
    H1d_deta_clone.Draw("hist same")

    L = TLegend(0.25, 0.57, 0.45, 0.8)
    L.SetLineWidth(0)
    L.AddEntry(H1d_eta1_clone, dnames_formatted["eta1"], "l")
    L.AddEntry(H1d_eta2_clone, dnames_formatted["eta2"], "l")
    L.AddEntry(H1d_eta_avg_clone, dnames_formatted["eta_avg"], "l")
    L.AddEntry(H1d_deta_clone, dnames_formatted["deta"], "l")
    L.Draw()

    C.Print("eta_asym_doubleJPsi_new_{}_{}_eta.pdf".\
            format(mass_region, region_abcd))
    P_hist.Clear()
    P_hist.Update()
    C.Clear()
    C.Update()
    C.Close()


def th1f_Ndeta_diff_create(mass_region, region_abcd, abins_eta_avg):
     

    # Set up 1d histogram.
    nbins_eta_avg = len(abins_eta_avg) - 1
    histname = "th1f_{}_Ndeta_diff_{}".format(mass_region, region_abcd)
    H1d = TH1F(histname, "", nbins_eta_avg, abins_eta_avg)

    # Fill in histogram.
    file_name = "dvar_doubleJPsi_new_{}_{}.json".\
        format(mass_region, region_abcd)

    with open(file_name, "r") as file:
        dvar = json.load(file)

    ldeta = dvar["l{}_deta_{}".format(mass_region, region_abcd)]
    leta_avg = dvar["l{}_eta_avg_{}".format(mass_region, region_abcd)]

    for event in range(len(ldeta)):
        eta_avg = abs(leta_avg[event])
        if ldeta[event] > 0:
            H1d.Fill(eta_avg, 1.)
        elif ldeta[event] < 0:
            H1d.Fill(eta_avg, -1.)
    
    print("Done {}".format(histname))
    
    return H1d
     

def th1f_Ndeta_tot_create(mass_region, region_abcd, abins_eta_avg):
     

    # Set up 1d histogram.
    nbins_eta_avg = len(abins_eta_avg) - 1
    histname = "th1f_{}_Ndeta_tot_{}".format(mass_region, region_abcd)
    H1d = TH1F(histname, "", nbins_eta_avg, abins_eta_avg)


    # Fill in histogram.
    file_name = "dvar_doubleJPsi_new_{}_{}.json".\
        format(mass_region, region_abcd)

    with open(file_name, "r") as file:
        dvar = json.load(file)

    leta_avg = dvar["l{}_eta_avg_{}".format(mass_region, region_abcd)]

    for event in range(len(leta_avg)):
        eta_avg = abs(leta_avg[event])
        H1d.Fill(eta_avg)
    
    print("Done {}".format(histname))
    
    return H1d


def th1f_asym_create_draw(mass_region, region_abcd, abins_eta_avg,
                        H1d_Ndeta_diff, H1d_Ndeta_tot,
                            mass_region_formatted):
    
    
    # Set up 1d histogram. 
    nbins_eta_avg = len(abins_eta_avg) - 1
    histname = "th1f_{}_asym_vs_eta_avg_{}".\
        format(mass_region, region_abcd)
    H1d_asym = TH1F(histname, "", nbins_eta_avg, abins_eta_avg)
    
    # Fill in histogram.
    for bin in range(1, nbins_eta_avg + 1):
        Ndeta_diff = float(H1d_Ndeta_diff.GetBinContent(bin))
        Ndeta_tot = float(H1d_Ndeta_tot.GetBinContent(bin))
        if Ndeta_tot == 0: continue
        asym = Ndeta_diff / Ndeta_tot
        H1d_asym.SetBinContent(bin, asym)
    

    # Draw and print histogram.
    H1d_asym_clone = H1d_asym.Clone(histname+"_clone")

    C = TCanvas("c_"+histname, "", 800, 650)
    P_hist = TPad("pad_asym_"+region_abcd, "", 0, 0, 1, 1)
    P_hist.SetTopMargin(0.05)
    P_hist.SetBottomMargin(0.2)
    P_hist.SetLeftMargin(0.17)
    P_hist.SetRightMargin(0.05)
    P_hist.Draw()

    C.cd()
    """ if region_abcd == "ABCD":
         region_abcd_formatted = "all regions ABCD"
    else:
         region_abcd_formatted = "region "+region_abcd
    
    title_line1 = "#eta asymmetry of J/#psi pairs for {} {}".\
                    format(region_abcd_formatted, mass_region_formatted) """
    """ title_line1 = "#eta asymmetry of double J/#psi events"

    T_line1 = TPaveText(0.4, 0.92, 0.7, 0.95, "NDC")
    T_line1.AddText(title_line1)
    T_line1.SetTextSize(0.05)
    T_line1.SetFillColor(0)
    T_line1.SetBorderSize(0)
    T_line1.Draw() """


    P_hist.cd()

    H1d_asym_clone.SetMaximum(1.5 * H1d_asym_clone.GetMaximum())
    H1d_asym_clone.SetMinimum(1.5 * H1d_asym_clone.GetMinimum())

    H1d_asym_clone.SetStats(0)
    H1d_asym_clone.GetXaxis().SetTitle("|#hat{#eta}_{J/#psi}|")
    H1d_asym_clone.GetXaxis().SetTitleSize(0.07)
    H1d_asym_clone.GetXaxis().SetTitleOffset(1.)
    H1d_asym_clone.GetXaxis().SetNdivisions(5, 0, 5)
    H1d_asym_clone.GetXaxis().SetLabelSize(0.05)
    H1d_asym_clone.GetYaxis().SetTitle("A_{2J/#psi}^{#eta}")
    H1d_asym_clone.GetYaxis().SetTitleSize(0.07)
    H1d_asym_clone.GetYaxis().SetTitleOffset(1.1)
    H1d_asym_clone.GetYaxis().SetLabelSize(0.05)

    """ H1d_asym_clone.SetMarkerStyle(20)
    H1d_asym_clone.SetMarkerSize(0.9)
    H1d_asym_clone.SetMarkerColor(kBlack)
    H1d_asym_clone.SetLineColor(kBlack)
    H1d_asym_clone.Draw("e0")
    H1d_asym_clone.Draw("p same") """

    H1d_asym_clone.SetLineStyle(1)                       # solid line
    H1d_asym_clone.SetFillStyle(1001)                    # solid fill
    H1d_asym_clone.SetLineColor(0)

    fill_color = TColor.GetColor(153, 196, 210)          # light blue
    H1d_asym_clone.SetFillColor(fill_color)
    H1d_asym_clone.Draw("hist")

    # Draw line at y=0.
    xmin = H1d_asym_clone.GetXaxis().GetXmin()
    xmax = H1d_asym_clone.GetXaxis().GetXmax()
    line = TLine(xmin, 0., xmax, 0.)
    line.Draw()

    T = TLatex()
    T.SetNDC()
    T.SetTextSize(0.05)
    T.SetTextFont(42)      # non-bold Helvetica
    T.SetTextAlign(13)     # align top left
    T.DrawLatex(0.37, 0.92, "3.05 < #hat{m}_{#mu#mu} < 3.15 GeV/c^{2}")


    C.Print("fig_eta_asym_doubleJPsi_new_{}_{}.pdf".\
            format(mass_region, region_abcd))
    P_hist.Clear()
    P_hist.Update()
    C.Clear()
    C.Update()
    C.Close()
         


def eta_asym(mass_region, region_abcd, lnames, dbins,
                mass_region_formatted, dnames_formatted):
    

    # Make 1d hists of eta1, eta2, eta_avg.
    eta1 = lnames[0]
    eta2 = lnames[1]
    eta_avg = lnames[2]
    deta = lnames[3]

    """ H1d_eta1 = th1f_create(
                mass_region, region_abcd, eta1, dbins[deta]
                )
    H1d_eta2 = th1f_create(
                mass_region, region_abcd, eta2, dbins[deta]
                )
    H1d_eta_avg = th1f_create(
                mass_region, region_abcd, eta_avg, dbins[deta]
                )
    H1d_deta = th1f_create(
                mass_region, region_abcd, deta, dbins[deta]
                )
    
    th1f_draw(
                mass_region, region_abcd,
                H1d_eta1, H1d_eta2, H1d_eta_avg, H1d_deta,
                    mass_region_formatted, dnames_formatted
                    ) """
    
    
    # Make 1d hist of N(deta > 0) - N(deta < 0) vs eta_avg.
    H1d_Ndeta_diff = th1f_Ndeta_diff_create(
                mass_region, region_abcd, dbins[eta_avg]
                )
    H1d_Ndeta_tot = th1f_Ndeta_tot_create(
                mass_region, region_abcd, dbins[eta_avg]
                )


    # Make 1d hist of N_tot = N(deta > 0) + N(deta < 0) vs eta_avg.
    th1f_asym_create_draw(mass_region, region_abcd, dbins[eta_avg],
                        H1d_Ndeta_diff, H1d_Ndeta_tot,
                            mass_region_formatted)

    print("################ Region {} completed".format(region_abcd))


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

    lnames = [
        "eta1",
        "eta2",
        "eta_avg",
        "deta"
    ]

    dbins = {
        "eta_avg": arange(0, 2.5, 0.1),
        "deta": arange(-5, 5.1, 0.1),
    }

    dmass_regions_formatted = {
        "all": "",
        "JPsi": "(J/#psi)",
        "control1": "(control 1)",
        "control2": "(control 2)",
    }

    dnames_formatted = {
        "eta1": "#eta_{J/#psi1}",
        "eta2": "#eta_{J/#psi2}",
        "eta_avg": "#eta_{avg}",
        "deta": "#Delta#eta"
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
                        eta_asym, mass_region, region_abcd, lnames, dbins,
                            mass_region_formatted, dnames_formatted
                            )
            
            lfutures.append(future)


    print("################ Main function ended ################")

