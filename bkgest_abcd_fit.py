"""
Create and plot background estimates using the ABCD method, using
2 variables x,y.
e.g. Abkg(mfour, alpha) = B(mfour) * R(alpha),
where R(alpha) is the best fit function of C(alpha) / D(alpha).

1. For a specified region of dimuon masses (mass region), 
    read json files for each ABCD region.
    json files can be made using bkgest_abcd_sort.py.
    
    Examples of mass regions: 
    "JPsi": events with dimuon masses of 3.05-3.15 GeV.
    "control1": events with dimuon masses of 2.5-3.05 GeV.
    "control2": events with dimuon masses of 3.15-4.0 GeV.


2. Use parallel processes from concurrent.futures to plot multiple 
    background estimates for either
    (1) 1 mass region and different combinations of x,y variables, or 
    (2) different mass regions and 1 combination of x,y variables.


For each process with a specified mass region and variables x,y:
3. Plot 1d histograms of A(x), B(x), C(x), D(x).
    Plot 1d histograms of A(y), B(y), C(y), D(y).

4. Plot R(y) by weighting each event in C(y) by 1 / D(y).
    Plot the degree 2 polynomial of best fit.

5. Plot Abkg = B(x) * R(y) by weighting each event in B(x)
    by R(y), using the best fit polynomial.

6. Plot pulls of A data vs background by calculating
    pull = (A(x) - Abkg(x,y)) / sqrt(A(x)) for each bin in x.

7. Draw and save plots to pdf files.

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

    
2. th1f_create(mass_region, region_abcd, var_name, bins)
    Create 1d histogram for a specified ABCD region, mass region,
    and variable.

    mass_region: str
    Name to specify the range of dimuon masses.

    region_abcd: str ("A", "B", "C", or "D")
    Name to specify the ABCD region.

    var_name: str
    Name to specify the variable (e.g. "alpha")

    bins: list [number of bins, array of bin edges]
    Specify bins for the histogram.
    Can be the output of percent_bins(a,b,p).

    Return: TH1F
    Filled histogram (not drawn).


3. th1f_draw(H1d_A, H1d_B, H1d_C, H1d_D,
              mass_region, var_name,
                name_formatted, mass_region_formatted, 
                axis_formatted)
    
    Draw 1d histograms of ABCD regions for a specified mass region
    and variable.
    All 4 histograms on the same graph, plotted as colored lines "hist".
    Events / bin width to normalize bin contents.
    Print and save to pdf file.

    H1d_A, H1d_B, H1d_C, H1d_D: TH1F
    1d histograms for ABCD regions.
    Should be outputs of th1f_create(...).

    mass_region, var_name: str
    Names to specify the range of dimuon masses and variable
    (e.g. mass_region = "JPsi", var_name = "alpha").

    name_formatted, mass_region_formatted: str
    Formatted names (including capitalization and LaTeX) to be displayed
    in histogram title.

    axis_formatted: str
    Formatted x-axis title (including capitalization and LaTeX).

    Return: None


4. th1f_ratiofit_create_draw(file_nameC, H1d_D, mass_region, 
                              var_namey, binsy,
                        namey_formatted, mass_region_formatted, 
                        yaxis_formatted)
    
    Create 1d histogram for C(y) / D(y) of a specified mass region 
    and variable y (e.g. "alpha").
    Find degree 2 polynomial of best fit.
    Plot C/D histogram as points with error bars "p" and "e0", along with
    best fit curve as red line.
    Print and save to pdf file.

    file_nameC: str
    Name of json file containing the region C data for a specified
    mass region and variable y.

    H1d_D: TH1F
    Histogram of region D for the same mass region and variable.
    Should be an output of th1f_create(...).

    mass_region, var_namey: str
    Names to specify mass region and variable y (e.g. "alpha" or "dR").
    
    binsy: list
    [number of bins, array of bin edges]

    namey_formatted, mass_region_formatted, yaxis_formatted: str

    Return: list [p0, p1, p2]
    Parameters for the best fit curve.


5. th1f_bkgfit_create(file_nameB, lpar1, mass_region, 
                    var_namex, var_namey, binsx)
    
    Create 1d histogram background estimate for a specified mass region 
    and variables x,y.
    Abkg(x,y) = B(x) * R(y), where R(y) = C(y) / D(y).
    e.g. Abkg(mfour, alpha) = B(mfour) * R(alpha).

    file_nameB: str
    Name of json file containing the region B data for a specified
    mass region and variable x.

    lpar1: list [p0, p1, p2]
    Parameters of C(y) / D(y) best fit curve for the same mass region,
    variable y. Should be the output of th1f_ratiofit_create_draw(...).

    mass_region, var_namex, var_namey: str

    binsx: list [number of bins, array of bin edges]

    Return: TH1F
    Filled histogram of Abkg(x,y).


6. th1f_pull_create(H1d_A, H1d_bkg, mass_region,
                     var_namex, var_namey, binsx)
    
    Create 1d pull plot of region A data vs background for a specified
    mass region and variables x,y.

    H1d_A: TH1F
    1d histogram of region A for a specified mass region, 
    variable x (e.g. "mfour").
    Should be an output of th1f_create(...).

    H1d_bkg: TH1F
    1d histogram of background estmiate for the same mass region,
    variable y (e.g. "alpha").
    Should be the output of th1f_bkgfit_create(...).

    mass_region, var_namex, var_namey: str

    binsx: list [number of bins, array of bin edges]

    Return: TH1F
    1d pull plot of data vs background.


7. th1f_pullhist_draw(H1d_pull, mass_region,
                       var_namex, var_namey, binsx,
                    namex_formatted, namey_formatted,
                    mass_region_formatted)
    
    Draw 1d histogram of pull values for a specified mass region and
    variable x,y.
    Events / bin width vs pull values (not pull plot).
    Plotted as points with error bars "p" and "e0".
    Print and save to pdf file.

    H1d_pull: TH1F
    1d pull plot of a specified mass region and variables x,y.
    Should be the output of th1f_pull_create(...).

    mass_region, var_namex, var_namey: str

    binsx: list [number of bins, array of bin edges]

    namex_formatted, namey_formatted, mass_region_formatted: str

    Return: None


8. th1f_bkgest_draw(H1d_A, H1d_bkg, H1d_pull, 
                     mass_region, var_namex, var_namey,
                        namex_formatted, namey_formatted,
                        mass_region_formatted, xaxis_formatted)
    
    Draw 1d histograms of region A data and background estimate on the 
    same graph, and 1d pull plot underneath.
    Region A data as points and error bars "p" and "e0".
    Background estimate and pull plot as colored lines with fill "hist".
    Print and save to pdf file.

    H1d_A: TH1F
    1d histogram of region A for a specified mass region, variable x.
    Should be an output of th1f_create(...).

    H1d_bkg: TH1F
    1d histogram of background estimate for the same mass region,
    variable y. Should be the output of th1f_bkgfit_create(...).

    H1d_pull: TH1F
    1d pull plot for the same mass region and variables x,y.
    Should be the output of th1f_pull_create(...).

    mass_region, var_namex, var_name: str

    namex_formatted, namey_formatted, mass_region_formatted,
    xaxis_formatted: str

    Return: None


9. abcd_bkg_estimate(mass_region, var_namex, var_namey, binsx, binsy,
                        namex_formatted, namey_formatted,
                        mass_region_formatted,
                        xaxis_formatted, yaxis_formatted)
    
    Call the above functions 1-8 in order to perform background
    estimate for a specified mass region and variables x,y.
    Serves as the target for each concurrent.futures process.

    Return: None


############################ Main function ############################
List all mass regions, variables x and y, bins x and y,
formatted names, and formatted axis labels.

Two options for parallel processes:
1. For one mass region, create a process for each combination of
    variables x,y.

2. For one combination of variables x,y, create a process for each
    mass region.

Call abcd_bkg_estimate(...) for each process.

"""

from ROOT import *
from numpy import array, empty, argmin, argmax, arange
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


def th1f_create(mass_region, region_abcd, var_name, bins):

    # Set up 1d histogram.
    nbins = bins[0]
    abins = bins[1]

    histname = "th1f_{}_{}_{}".format(mass_region, var_name, region_abcd)
    H1d = TH1F(histname, "", nbins, abins)

    # Fill in histogram with values from file.
    file_name = "dvar_doubleJPsi_new_{}_{}.json".\
        format(mass_region, region_abcd)

    with open(file_name, "r") as file:
            dvar = json.load(file)
    
    list_name = "l{}_{}_{}".format(mass_region, var_name, region_abcd)
    for event in dvar[list_name]:
        H1d.Fill(event)
    
    """ print("Done 1d histogram for region {} of {} ({})".\
                    format(region_abcd, var_name, mass_region)) """
    
    return H1d
    

def th1f_draw(H1d_A, H1d_B, H1d_C, H1d_D,
              mass_region, var_name,
                name_formatted, mass_region_formatted, 
                axis_formatted):


    # Draw and print 1d histogram.
    histname = "th1f_{}_{}".format(mass_region, var_name)

    H1d_A_clone = H1d_A.Clone(histname+"_A_clone")
    H1d_B_clone = H1d_B.Clone(histname+"_B_clone")
    H1d_C_clone = H1d_C.Clone(histname+"_C_clone")
    H1d_D_clone = H1d_D.Clone(histname+"_D_clone")
    
    H1d_A_clone.Scale(1.0, "width")
    H1d_B_clone.Scale(1.0, "width")
    H1d_C_clone.Scale(1.0, "width")
    H1d_D_clone.Scale(1.0, "width")
    
    C = TCanvas()
    P_hist = TPad("pad_"+histname, "", 0, 0, 1, 1)
    P_hist.SetTopMargin(0.22)
    P_hist.SetBottomMargin(0.2)
    P_hist.SetLeftMargin(0.2)
    P_hist.SetRightMargin(0.05)
    P_hist.Draw()

    C.cd()
    title_line1 = "Regions ABCD of {}".format(name_formatted)
    title_line2 = "{}".format(mass_region_formatted)
    
    T_line1 = TPaveText(0.4, 0.92, 0.7, 0.95, "NDC")
    T_line2 = TPaveText(0.4, 0.82, 0.7, 0.85, "NDC")
    T_line1.AddText(title_line1)
    T_line2.AddText(title_line2)
    T_line1.SetTextSize(0.06)
    T_line2.SetTextSize(0.06)
    T_line1.SetFillColor(0)
    T_line2.SetFillColor(0)
    T_line1.SetBorderSize(0)
    T_line2.SetBorderSize(0)
    T_line1.Draw()
    T_line2.Draw()

    P_hist.cd()
    P_hist.SetLogx()

    H1d_A_clone.SetLineColor(kBlack)
    H1d_B_clone.SetLineColor(kBlue)
    H1d_C_clone.SetLineColor(kGreen+2)
    H1d_D_clone.SetLineColor(kPink)

    H1d_A_clone.SetLineStyle(1)         # solid
    H1d_B_clone.SetLineStyle(2)         # dashed
    H1d_C_clone.SetLineStyle(3)         # dotted
    H1d_D_clone.SetLineStyle(4)         # dash-dotted

    H1d_A_clone.SetLineWidth(2)
    H1d_B_clone.SetLineWidth(2)
    H1d_C_clone.SetLineWidth(2)
    H1d_D_clone.SetLineWidth(2)

    aymax_abcd = array([H1d_A_clone.GetMaximum(), H1d_B_clone.GetMaximum(),
                        H1d_C_clone.GetMaximum(), H1d_D_clone.GetMaximum()])
    ymax_abcd = aymax_abcd[argmax(aymax_abcd)]

    H1d_A_clone.SetStats(0)
    H1d_A_clone.GetXaxis().SetTitle(axis_formatted)
    H1d_A_clone.GetXaxis().SetTitleSize(0.06)
    H1d_A_clone.GetXaxis().SetTitleOffset(1.5)
    H1d_A_clone.GetXaxis().SetLabelSize(0.05)
    H1d_A_clone.GetYaxis().SetTitle("Events / bin width")
    H1d_A_clone.GetYaxis().SetTitleSize(0.06)
    H1d_A_clone.GetYaxis().SetTitleOffset(1.5)
    H1d_A_clone.GetYaxis().SetLabelSize(0.05)
    H1d_A_clone.SetMaximum(1.2 * ymax_abcd)

    H1d_A_clone.Draw("hist")
    H1d_B_clone.Draw("hist same")
    H1d_C_clone.Draw("hist same")
    H1d_D_clone.Draw("hist same")

    L = TLegend(0.23, 0.58, 0.46, 0.75)
    L.SetLineWidth(0)
    L.AddEntry(H1d_A_clone, "Signal region A", "l")
    L.AddEntry(H1d_B_clone, "Control region B", "l")
    L.AddEntry(H1d_C_clone, "Control region C", "l")
    L.AddEntry(H1d_D_clone, "Control region D", "l")
    L.Draw()

    C.Print("bkgest_doubleJPsi_new_{}_{}_ABCD.pdf".\
            format(mass_region, var_name))
    P_hist.Clear()
    P_hist.Update()
    C.Clear()
    C.Update()
    C.Close()


def th1f_ratiofit_create_draw(file_nameC, H1d_D, mass_region, 
                              var_namey, binsy,
                        namey_formatted, mass_region_formatted, 
                        yaxis_formatted):
    
    # Set up 1d histogram.
    nbinsy = binsy[0]
    abinsy = binsy[1]

    histname = "th1f_{}_{}_ratio".\
        format(mass_region, var_namex)
    H1d_ratio = TH1F(histname, "", nbinsy, abinsy)

    
    # Fill in histogram with C/D.
    with open(file_nameC, "r") as file:
        dvar = json.load(file)
    
    list_namey = "l{}_{}_C".format(mass_region, var_namey)

    for valuey in dvar[list_namey]:
        eventsD = H1d_D.GetBinContent(H1d_ratio.FindBin(valuey))
        if eventsD == 0:
            weight = 0
        else:
            weight = float(1 / eventsD)
        H1d_ratio.Fill(valuey, weight)


    # Fit polynomial function.
    polyfit1 = TF1("polyfit1", "pol2")

    H1d_ratio_clone = H1d_ratio.Clone(histname+"_clone")
    H1d_ratio_clone.Fit(polyfit1)
    
    # Draw and print histogram.
    C = TCanvas("c_"+histname, "", 800, 650)
    P_hist = TPad("pad_"+histname, "", 0, 0, 1, 1)
    P_hist.SetTopMargin(0.05)
    P_hist.SetBottomMargin(0.15)
    P_hist.SetLeftMargin(0.17)
    P_hist.SetRightMargin(0.05)
    P_hist.Draw()

    C.cd()
    """ title_line1 = "q_{4#mu} = 0 pass-fail ratio"
    title_line2 = "of {} {}".\
                    format(namey_formatted, mass_region_formatted)
    
    T_line1 = TPaveText(0.4, 0.92, 0.7, 0.95, "NDC")
    T_line2 = TPaveText(0.4, 0.82, 0.7, 0.85, "NDC")
    T_line1.AddText(title_line1)
    T_line2.AddText(title_line2)
    T_line1.SetTextSize(0.06)
    T_line2.SetTextSize(0.06)
    T_line1.SetFillColor(0)
    T_line2.SetFillColor(0)
    T_line1.SetBorderSize(0)
    T_line2.SetBorderSize(0)
    T_line1.Draw()
    T_line2.Draw() """


    P_hist.cd()
    P_hist.SetLogx()

    """ H1d_ratio_clone.SetLineStyle(1)                       # solid line
    H1d_ratio_clone.SetFillStyle(1001)                    # solid fill
    H1d_ratio_clone.SetLineColor(0)

    fill_color = TColor.GetColor(206, 181, 206)           # light purple
    H1d_ratio_clone.SetFillColor(fill_color) """

    ytitle = "q_{4#mu} = 0 pass-fail ratio of events"

    H1d_ratio_clone.SetStats(0)
    H1d_ratio_clone.GetXaxis().SetTitle(yaxis_formatted)
    H1d_ratio_clone.GetXaxis().SetTitleSize(0.06)
    H1d_ratio_clone.GetXaxis().SetTitleOffset(1.)
    H1d_ratio_clone.GetXaxis().SetLabelSize(0.06)
    H1d_ratio_clone.GetYaxis().SetTitle(ytitle)
    H1d_ratio_clone.GetYaxis().SetTitleSize(0.06)
    H1d_ratio_clone.GetYaxis().SetTitleOffset(1.)
    H1d_ratio_clone.GetYaxis().SetLabelSize(0.05)

    H1d_ratio_clone.SetMarkerStyle(20)
    H1d_ratio_clone.SetMarkerSize(0.8)
    H1d_ratio_clone.SetMarkerColor(kBlack)
    H1d_ratio_clone.SetLineColor(kBlack)

    #H1d_ratio_clone.Draw("hist")
    H1d_ratio_clone.Draw("e0")
    H1d_ratio_clone.Draw("p same")
    polyfit1.Draw("same")

    T = TLatex()
    T.SetNDC()
    T.SetTextSize(0.04)
    T.SetTextFont(42)      # non-bold Helvetica
    T.SetTextAlign(13)     # align top left
    T.DrawLatex(0.6, 0.9, mass_region_formatted)

    L = TLegend(0.6, 0.72, 0.9, 0.82)
    L.SetLineWidth(0)
    L.SetTextSize(0.04)
    L.AddEntry(H1d_ratio_clone, "Data", "lep")
    L.AddEntry(polyfit1, "Polynomial fit R_{C/D}", "l")
    L.Draw()

    C.Print("fig_bkgest_doubleJPsi_new_{}_fit_{}_ratio.pdf".\
            format(mass_region, var_namey))
    P_hist.Clear()
    P_hist.Update()
    C.Clear()
    C.Update()
    C.Close()

    # Return parameters of fit equations.
    lpar1 = []
    for i in range(polyfit1.GetNpar()):
        lpar1.append(polyfit1.GetParameter(i))

    return lpar1


def th1f_bkgfit_create(file_nameB, lpar1, mass_region, 
                    var_namex, var_namey, binsx):

    
    def polyfit(x, lpar1):
        if len(lpar1) == 3:
            lpar1.append(0)
            lpar1.append(0)
        return lpar1[0] + lpar1[1]*x + lpar1[2]*x**2 \
            + lpar1[3]*x**3 + lpar1[4]*x**4
    
    
    # Set up 1d histogram.
    nbinsx = binsx[0]
    abinsx = binsx[1]

    histname = "th1f_{}_{}_{}_Abkg".\
        format(mass_region, var_namex, var_namey)
    H1d_bkg = TH1F(histname, "", nbinsx, abinsx)

    
    # Fill in histogram with B(x)*R(y).
    with open(file_nameB, "r") as file:
        dvar = json.load(file)
    
    list_namex = "l{}_{}_B".format(mass_region, var_namex)
    list_namey = "l{}_{}_B".format(mass_region, var_namey)
    events = len(dvar[list_namex])

    for i in range(0, events):
        valuex = dvar[list_namex][i]
        valuey = dvar[list_namey][i]
        weight = polyfit(valuey, lpar1)
        H1d_bkg.Fill(valuex, weight)

    print("Done 1d histogram for bkg of {} as function of {} ({})".\
          format(var_namex, var_namey, mass_region))

    return H1d_bkg


def th1f_pull_create(H1d_A, H1d_bkg, mass_region,
                     var_namex, var_namey, binsx):
    

    # Set up 1d histogram.
    nbinsx = binsx[0]
    abinsx = binsx[1]

    pullname = "th1f_{}_{}_{}_Apull".\
                format(mass_region, var_namex, var_namey)
    H1d_pull = TH1F(pullname, "", nbinsx, abinsx)

    # Calculate pulls and fill in histogram.
    for binx in range(1, nbinsx + 1):
        
        """ # If plotting m_avg, skip J/Psi bin.
        if var_namex == "mavg" \
        and mass_region == "all" \
        and H1d_pull.GetBinCenter(binx) > 3.0 \
        and H1d_pull.GetBinCenter(binx) < 3.2:
            continue """

        data = H1d_A.GetBinContent(binx)
        if data == 0: continue 
        bkg = H1d_bkg.GetBinContent(binx)
        #if bkg == 0: continue
        uncert = TMath.Sqrt(data)
        pull = (data - bkg) / uncert
        
        H1d_pull.SetBinContent(binx, pull)

    print("Done 1d histogram for pulls of {} as function of {} ({})".\
          format(var_namex, var_namey, mass_region))

    return H1d_pull
    

def th1f_pullhist_draw(H1d_pull, mass_region,
                       var_namex, var_namey, binsx,
                    namex_formatted, namey_formatted,
                    mass_region_formatted):
        

    # Set up and fill 1d histogram.
    nbinsx = binsx[0]
    apulls = empty(nbinsx, float)
    for bin in range(1, nbinsx +1):
        apulls[bin -1] = H1d_pull.GetBinContent(bin)
    
    min_pull = int(apulls[argmin(apulls)]-1)
    max_pull = int(apulls[argmax(apulls)]+1)
    nbins_pull = max_pull - min_pull

    histname = "th1f_{}_{}_pullhist".format(var_namex, var_namey)
    H1d_pullhist = TH1F(histname, "", 
                        nbins_pull, float(min_pull), float(max_pull))

    for pull in apulls:
        H1d_pullhist.Fill(pull)

    # Draw and print histogram.
    H1d_pullhist_clone = H1d_pullhist.Clone(histname+"_clone")

    C = TCanvas()
    P_hist = TPad("pad_"+histname, "", 0, 0, 1, 1)
    P_hist.SetTopMargin(0.22)
    P_hist.SetBottomMargin(0.2)
    P_hist.SetLeftMargin(0.15)
    P_hist.SetRightMargin(0.05)
    P_hist.Draw()

    C.cd()
    title_line1 = "Pulls of data vs. background of {}".\
        format(namex_formatted)
    title_line2 = "as function of {} {}".\
        format(namey_formatted, mass_region_formatted)
    
    T_line1 = TPaveText(0.4, 0.92, 0.7, 0.95, "NDC")
    T_line2 = TPaveText(0.4, 0.82, 0.7, 0.85, "NDC")
    T_line1.AddText(title_line1)
    T_line2.AddText(title_line2)
    T_line1.SetTextSize(0.06)
    T_line2.SetTextSize(0.06)
    T_line1.SetFillColor(0)
    T_line2.SetFillColor(0)
    T_line1.SetBorderSize(0)
    T_line2.SetBorderSize(0)
    T_line1.Draw()
    T_line2.Draw()

    P_hist.cd()

    H1d_pullhist_clone.SetStats(0)
    H1d_pullhist_clone.GetXaxis().SetTitle("Pull")
    H1d_pullhist_clone.GetXaxis().SetTitleSize(0.06)
    H1d_pullhist_clone.GetXaxis().SetTitleOffset(1.5)
    H1d_pullhist_clone.GetXaxis().SetLabelSize(0.05)
    H1d_pullhist_clone.GetYaxis().SetTitle("Bins")
    H1d_pullhist_clone.GetYaxis().SetTitleSize(0.06)
    H1d_pullhist_clone.GetYaxis().SetTitleOffset(1.2)
    H1d_pullhist_clone.GetYaxis().SetLabelSize(0.05)

    H1d_pullhist_clone.SetMarkerStyle(20)
    H1d_pullhist_clone.SetMarkerSize(0.7)
    H1d_pullhist_clone.SetMarkerColor(kBlue)
    H1d_pullhist_clone.SetLineColor(kBlue)

    H1d_pullhist_clone.Draw("e0")
    H1d_pullhist_clone.Draw("p same")

    # Display mean and standard deviation.
    C.cd()
    mean = H1d_pullhist_clone.GetMean()
    std_dev = H1d_pullhist_clone.GetStdDev()
    S = TPaveText(0.75, 0.75, 0.85, 0.65)
    #S = TPaveText(0.25, 0.75, 0.35, 0.65)
    S.AddText("Mean: {:.3f}".format(mean))
    S.AddText("Std dev: {:.3f}".format(std_dev))
    S.SetTextFont(42)       # unbold
    S.SetTextSize(0.04)
    S.SetTextSize(0.04)
    S.SetFillColor(0)
    S.SetFillColor(0)
    S.SetBorderSize(0)
    S.SetBorderSize(0)
    S.Draw()

    C.Print("bkgest_doubleJPsi_new_{}_{}_{}_pullhist.pdf".\
            format(mass_region, var_namex, var_namey))
    P_hist.Clear()
    P_hist.Update()
    C.Clear()
    C.Update()
    C.Close()


def th1f_bkgest_draw(H1d_A, H1d_bkg, H1d_pull, 
                     mass_region, var_namex, var_namey,
                        namex_formatted, namey_formatted,
                        mass_region_formatted, xaxis_formatted):
    
    
    # Draw and print 1d histograms of A, Abkg, Apull.
    H1d_A_clone = H1d_A.Clone(H1d_A.GetName()+"_clone")
    H1d_bkg_clone = H1d_bkg.Clone(H1d_bkg.GetName()+"_clone")
    H1d_pull_clone = H1d_pull.Clone(H1d_pull.GetName()+"_clone")
    
    H1d_bkg_clone.Scale(1.0, "width")
    H1d_A_clone.Scale(1.0, "width")


    C = TCanvas("canvas_{}_{}".format(var_namex, var_namey), "",
                800, 650)
    P_hist = TPad("th1f_x_data_bkg_{}_{}".format(var_namex, var_namey), "",
                    0, 0.35, 1, 1)     # top 65% of canvas
    P_pull = TPad("th1f_x_pull_{}_{}".format(var_namex, var_namey), "",
                    0, 0, 1, 0.34)     # bottom 35% of canvas

    P_hist.SetTopMargin(0.05)
    P_hist.SetBottomMargin(0.02)
    P_hist.SetLeftMargin(0.15)
    P_hist.SetRightMargin(0.05)
    P_pull.SetTopMargin(0.02)
    P_pull.SetBottomMargin(0.35)
    P_pull.SetLeftMargin(0.15)
    P_pull.SetRightMargin(0.05)

    P_hist.Draw()
    P_pull.Draw()

    C.cd()
    
    """ hist_title_line1 = "Background estimate of {}".\
        format(namex_formatted)
    hist_title_line2 = "as function of {} {}".\
        format(namey_formatted, mass_region_formatted)
    
    hist_title_line1 = "Background estimate of {} as function of #alpha {}".\
                        format("m_{4#mu}", mass_region_formatted)
    
    T_line1 = TPaveText(0.3, 0.93, 0.7, 0.96, "NDC")
    T_line1.AddText(hist_title_line1)
    T_line1.SetTextSize(0.045)
    T_line1.SetFillColor(0)
    T_line1.SetBorderSize(0)
    T_line1.Draw()

    T_line2 = TPaveText(0.3, 0.86, 0.7, 0.89, "NDC")
    T_line2.AddText(hist_title_line2)
    T_line2.SetTextSize(0.045)
    T_line2.SetFillColor(0)
    T_line2.SetBorderSize(0)
    T_line2.Draw() """

    # Region A background and data.
    P_hist.cd()
    P_hist.SetLogx()

    amax = array([H1d_A_clone.GetMaximum(), H1d_bkg_clone.GetMaximum()])
    ymax = amax[argmax(amax)]

    H1d_bkg_clone.SetLineStyle(1)              # solid line
    H1d_bkg_clone.SetFillStyle(1001)           # solid fill
    H1d_bkg_clone.SetLineColor(0)
    H1d_bkg_clone.SetFillColor(kOrange)
    H1d_A_clone.SetMarkerStyle(20)
    H1d_A_clone.SetMarkerSize(0.9)
    H1d_A_clone.SetLineColor(kBlack)
    H1d_A_clone.SetMarkerColor(kBlack)

    H1d_bkg_clone.SetTitle("")
    H1d_bkg_clone.GetXaxis().SetTitle("")
    H1d_bkg_clone.GetXaxis().SetLabelOffset(100)
    H1d_bkg_clone.GetYaxis().SetTitle("Events / bin width")
    H1d_bkg_clone.GetYaxis().SetTitleSize(0.09)
    H1d_bkg_clone.GetYaxis().SetLabelSize(0.05)
    H1d_bkg_clone.GetYaxis().SetTitleOffset(0.7)
    H1d_bkg_clone.SetMaximum(1.2 * ymax)
    H1d_bkg_clone.SetStats(0)

    H1d_bkg_clone.Draw("hist")
    H1d_A_clone.Draw("p same")
    H1d_A_clone.Draw("e0 same")

    T = TLatex()
    T.SetNDC()
    T.SetTextSize(0.07)
    T.SetTextFont(42)      # non-bold Helvetica
    T.SetTextAlign(13)     # align top left
    T.DrawLatex(0.55, 0.9, mass_region_formatted)
    
    L = TLegend(0.55, 0.6, 0.9, 0.82)
    L.SetLineWidth(0)
    L.SetTextSize(0.07)
    L.AddEntry(H1d_A_clone, "Data", "lep")   #3.15 < #hat{m}_{#mu#mu} < 4.0 GeV
    L.AddEntry(H1d_bkg_clone, "Background", "f")
    L.Draw()
    
    
    # Pull plot.
    P_pull.cd()
    P_pull.SetLogx()

    H1d_pull_clone.SetLineStyle(1)              # solid line
    H1d_pull_clone.SetFillStyle(1001)           # solid fill
    H1d_pull_clone.SetLineColor(0)

    fill_color = TColor.GetColor(153, 196, 210)  # light blue
    H1d_pull_clone.SetFillColor(fill_color)

    H1d_pull_clone.SetStats(0)
    #H1d_pull_clone.GetXaxis().SetTitle(xaxis_formatted)
    H1d_pull_clone.GetXaxis().SetTitle("m_{4#mu} (GeV/c^{2})")
    H1d_pull_clone.GetXaxis().SetTitleSize(0.13)
    H1d_pull_clone.GetXaxis().SetTitleOffset(1)
    H1d_pull_clone.GetXaxis().SetLabelSize(0.1)
    H1d_pull_clone.GetYaxis().SetTitle("Pull")
    H1d_pull_clone.GetYaxis().SetTitleSize(0.15)
    H1d_pull_clone.GetYaxis().SetTitleOffset(0.4)
    H1d_pull_clone.GetYaxis().SetNdivisions(504)
    H1d_pull_clone.GetYaxis().SetLabelSize(0.11)

    H1d_pull_clone.Draw("hist")

    # Draw line at y=0.
    xmin = H1d_pull_clone.GetXaxis().GetXmin()
    xmax = H1d_pull_clone.GetXaxis().GetXmax()
    line = TLine(xmin, 0., xmax, 0.)
    line.Draw()
    

    C.Print("fig_bkgest_doubleJPsi_new_{}_fit_{}_{}.pdf".\
            format(mass_region, var_namex, var_namey))
    L.Clear()
    P_hist.Clear()
    P_hist.Update()
    P_pull.Clear()
    P_pull.Update()
    C.Clear()
    C.Update()
    C.Close()



def abcd_bkg_estimate(mass_region, var_namex, var_namey, binsx, binsy,
                        namex_formatted, namey_formatted,
                        mass_region_formatted,
                        xaxis_formatted, yaxis_formatted):
    
    print("##### Started {} as function of {} ({}) #####").\
            format(var_namex, var_namey, mass_region)
    
    # Create and fill 1d hists of A,B,C,D regions.
    H1d_A = th1f_create(mass_region, "A", var_namex, binsx)
    H1d_B = th1f_create(mass_region, "B", var_namex, binsx)
    #H1d_B = th1f_create(mass_region, "B", var_namey, binsy)
    H1d_C = th1f_create(mass_region, "C", var_namey, binsy)
    H1d_D = th1f_create(mass_region, "D", var_namey, binsy)

    """ H1d_A_y = th1f_create(mass_region, "A", var_namey, binsy)
    H1d_B_y = th1f_create(mass_region, "B", var_namey, binsy)
    H1d_C_x = th1f_create(mass_region, "C", var_namex, binsx)
    H1d_D_x = th1f_create(mass_region, "D", var_namex, binsx) """

    print("Done 1d histograms for regions ABCD of {} and {} {}").\
            format(var_namex, var_namey, mass_region)


    # Draw and print 1d hists.
    """ th1f_draw(H1d_A, H1d_B, H1d_C_x, H1d_D_x,
              mass_region, var_namex,
                namex_formatted, mass_region_formatted, 
                xaxis_formatted) """
    
    """ th1f_draw(H1d_A_y, H1d_B_y, H1d_C, H1d_D,
              mass_region, var_namey,
                namey_formatted, mass_region_formatted, 
                yaxis_formatted) """

    
    # Draw and print 1d hist of the ratio R = C/D.
    lpar1 = th1f_ratiofit_create_draw(
                    "dvar_doubleJPsi_new_{}_C.json".format(mass_region), H1d_D, 
                    mass_region, var_namey, binsy,
                        namey_formatted, mass_region_formatted, 
                        yaxis_formatted
                        )

    # Create and fill 1d hists of Abkg = B(x)*R(y) and A pull plot.
    H1d_bkg = th1f_bkgfit_create(
                    "dvar_doubleJPsi_new_{}_B.json".format(mass_region),
                    lpar1, mass_region, var_namex, var_namey, binsx
                    )
    
    H1d_pull = th1f_pull_create(
                    H1d_A, H1d_bkg, mass_region,
                     var_namex, var_namey, binsx
                     )
    
    
    # Draw and print 1d hists of A, Abkg, Apull.
    """ th1f_pullhist_draw(
                    H1d_pull, mass_region,
                    var_namex, var_namey, binsx,
                        namex_formatted, namey_formatted,
                        mass_region_formatted
                        ) """
    
    th1f_bkgest_draw(H1d_A, H1d_bkg, H1d_pull, 
                     mass_region, var_namex, var_namey,
                        namex_formatted, namey_formatted,
                        mass_region_formatted, xaxis_formatted
                        )
    

    print("##### Completed {} as function of {} ({}) #####").\
            format(var_namex, var_namey, mass_region)


if __name__ == "__main__":

    print("################ Main function started ################")

    #mass_region = "all" 
    #mass_region = "JPsi" 
    #mass_region = "other"
    
    lmass_regions = [#"control1", 
                     "JPsi", 
                     "control2",
                     #"2gev", 
                     #"4gev", 
                     #"5gev",
                     #"8gev",
                     #"10gev",
                     #"20gev"
                     ]

    lnames_x = [#"mavg", 
                "mfour",
                #"alpha", 
                #"dR", #"deta",
                #"pt_lead", #"pt_last",
                #"pt_dimlead", "pt_dimlast"
                ]
    
    lnames_y = [#"mavg", "mfour",
                "alpha", 
                #"dR", #"deta",
                #"pt_lead", #"pt_last",
                #"pt_dimlead", #"pt_dimlast"
                ]

    print("x variables: {}".format(lnames_x))
    print("y variables: {}".format(lnames_y))


    # Set bins for mavg to match J/Psi window.
    nbins_mavg_temp, abins_mavg_temp = percent_bins(1e-1, 1e2, 0.2)
    lbins_mavg = [bin for bin in abins_mavg_temp \
                  if bin < 3.0 or bin > 3.2]
    lbins_mavg.append(3.0)
    lbins_mavg.append(3.2)
    lbins_mavg.sort()
    bins_mavg = [len(lbins_mavg) -1, array(lbins_mavg)]

    """ abins_mavg = arange(1., 5.1, 0.1)
    bins_mavg = [len(abins_mavg) -1, abins_mavg] """


    dbins_x = {
        #"mavg": list(percent_bins(1e-1, 1e2, 0.1)),
        "mavg": bins_mavg,
        #"mfour": list(percent_bins(1e-2, 1e3, 0.2)),
        "mfour": list(percent_bins(1., 5e2, 0.2)),
        "alpha": list(percent_bins(1e-3, 1., 0.2)),
        "dR": list(percent_bins(1e-5, 50., 0.3)),
        #"dR": list(percent_bins(1e-4, 1e-3, 0.3)),
        "deta": list(percent_bins(1e-8, 50., 0.1)),
        "pt_lead": list(percent_bins(1e-1, 1e5, 0.3)),
        "pt_last": list(percent_bins(1., 1e2, 0.1)),
        "pt_dimlead": list(percent_bins(5e-2, 1e7, 0.3)),
        "pt_dimlast": list(percent_bins(1e-2, 1e2, 0.1))
    }

    dbins_y = dbins_x

    dmass_regions_formatted = {
        "JPsi": "3.05 < #hat{m}_{#mu#mu} < 3.15 GeV/c^{2}",
        "control1": "2.5 < #hat{m}_{#mu#mu} < 3.05 GeV/c^{2}",
        "control2": "3.15 < #hat{m}_{#mu#mu} < 4.0 GeV/c^{2}",
        "2gev": "2.0 < #hat{m}_{#mu#mu} < 3.0 GeV/c^{2}",
        "4gev": "4.0 < #hat{m}_{#mu#mu} < 5.0 GeV/c^{2}",
        "5gev": "5.0 < #hat{m}_{#mu#mu} < 6.0 GeV/c^{2}",
        "8gev": "8.0 < #hat{m}_{#mu#mu} < 9.0 GeV/c^{2}",
        "10gev": "10.0 < #hat{m}_{#mu#mu} < 11.0 GeV/c^{2}",
        "20gev": "20.0 < #hat{m}_{#mu#mu} < 21.0 GeV/c^{2}",
    }

    dnames_formatted = {
        "mavg": "dimuon m_{avg}",
        "mfour": "m_{4#mu}",
        "alpha": "#alpha",
        "dR": "#DeltaR",
        "deta": "|#Delta#eta|",
        "pt_lead": "leading #mu p_{T}",
        "pt_last": "slowest #mu p_{T}",
        "pt_dimlead": "leading di-#mu p_{T}",
        "pt_dimlast": "slowest di-#mu p_{T}"
    }

    daxis_titles = {
        "mavg": "Dimuon m_{avg} (GeV)",
        "mfour": "m_{4#mu} (GeV)",
        "alpha": "#alpha",
        "dR": "#DeltaR",
        "deta": "|#Delta#eta|",
        "pt_lead": "Leading #mu p_{T} (GeV)",
        "pt_last": "Slowest #mu p_{T} (GeV)",
        "pt_dimlead": "Leading di-#mu p_{T} (GeV)",
        "pt_dimlast": "Slowest di-#mu p_{T} (GeV)"
    }

    # List all possible pairs of variables x,y.
    l2dpairs = []
    for y in range(0, len(lnames_y)):
        for x in range(0, len(lnames_x)):
            namey = lnames_y[y]
            namex = lnames_x[x]
            if namex != namey:
                l2dpairs.append([namex, namey])

    
    """ # Plot background estimates for 1 mass region and multiple xy pairs.
    # For each xy pair, create parallel process to make background estimate.
    n = len(l2dpairs)
    gROOT.SetBatch(True)
    lfutures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers = n) \
    as executor:
        for lpair in l2dpairs:

            var_namex = lpair[0]
            var_namey = lpair[1]
            binsx = dbins_x[var_namex]
            binsy = dbins_y[var_namey]

            namex_formatted = dnames_formatted[var_namex]
            namey_formatted = dnames_formatted[var_namey]
            mass_region_formatted = dmass_regions_formatted[mass_region]
            xaxis_formatted = daxis_titles[var_namex]
            yaxis_formatted = daxis_titles[var_namey]

            future = executor.submit(
                        abcd_bkg_estimate, 
                        mass_region, var_namex, var_namey, binsx, binsy,
                        namex_formatted, namey_formatted,
                        mass_region_formatted,
                        xaxis_formatted, yaxis_formatted)
            
            lfutures.append(future)
        
        for future in concurrent.futures.as_completed(lfutures):
            try:
                None
            except Exception as exc:
                print("Exception in {} as function of {} ({}): {} ").\
                        format(var_namex, var_namey, mass_region, exc) """
            

    # Plot background estimates for multiple mass regions and 1 xy pair.
    # For each mass region, create parallel process to make background estimate.
    n = len(lmass_regions)
    gROOT.SetBatch(True)
    lfutures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers = n) \
    as executor:

        lpair = l2dpairs[0]
        var_namex = lpair[0]
        var_namey = lpair[1]
        binsx = dbins_x[var_namex]
        binsy = dbins_y[var_namey]

        namex_formatted = dnames_formatted[var_namex]
        namey_formatted = dnames_formatted[var_namey]
        xaxis_formatted = daxis_titles[var_namex]
        yaxis_formatted = daxis_titles[var_namey]

        for mass_region in lmass_regions:
            mass_region_formatted = dmass_regions_formatted[mass_region]

            future = executor.submit(
                        abcd_bkg_estimate, 
                        mass_region, var_namex, var_namey, binsx, binsy,
                        namex_formatted, namey_formatted,
                        mass_region_formatted,
                        xaxis_formatted, yaxis_formatted)
            
            lfutures.append(future)


    print("################ Main function ended ################")

