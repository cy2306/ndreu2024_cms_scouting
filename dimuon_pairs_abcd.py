# ABCD method on dimuon pairs.
import ROOT
from ROOT import *
from numpy import *

#########################################################################
#########  Functions for creating and drawing histograms.  ##############

# percent_bins creates an array of bin edges, where the width of each bin
# is p% of the central value.

# Arguments: 
# a: lower limit of histogram range (float)
# b: upper limit of histogram range (float)
# p: percent of central value, determines width of each bin (float, 0<p<1)

# Return:
# number of bins (int), array of bin edges (array of floats)

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

############################################################################

# th1f_abcd_make creates histograms of ABCD regions for a given variable.

# Arguments:
# name: name of variable to be used in histogram names (string)
# nbins: number of bins (int)
# abins: array of bin edges (array of floats)

# Returns:
# list of empty histograms for ABCD regions

def th1f_abcd_make(name, nbins, abins):
    H_A = TH1F("h_"+name+"_A", "", \
              nbins, abins)
    H_B = TH1F("h_"+name+"_B", "", \
                nbins, abins)
    H_C = TH1F("h_"+name+"_C", "", \
                nbins, abins)
    H_D = TH1F("h_"+name+"_D", "", \
                nbins, abins)

    return [H_A, H_B, H_C, H_D]

############################################################################

# th1f_apull_make creates histogram and pull plot of region A data vs bkg.

# Arguments:
# name: name of variable to be used in histogram names (string)
# lhist: list of *filled* histograms of ABCD regions.
# nbins: number of bins (int)
# abins: array of bin edges (array of floats)

# Returns:
# list of filled histograms for region A data, background, and pull

def th1f_apull_make(name, lhist, nbins, abins):
    H_A = lhist[0]
    H_B = lhist[1]
    H_C = lhist[2]
    H_D = lhist[3]

    # Creating histograms for region A data and background.
    H_Adata = H_A.Clone("h_"+name+"_Adata")      # Clone A.
    H_Abkg = H_C.Clone("h_"+name+"_Abkg")        # Clone C.
    H_Abkg.Divide(H_D)                          # Divide by D.
    H_Abkg.Multiply(H_B)                        # Multiply by B.

    H_pull = TH1F("pull_"+name, "", \
                   nbins, abins)

    # Creating pull plot for region A data and background.
    for bin in range(1, nbins + 1):
        data = H_Adata.GetBinContent(bin)
        if data == 0: continue 
        bkg = H_Abkg.GetBinContent(bin)
        uncert = TMath.Sqrt(data)
        pull = (data - bkg) / uncert
        H_pull.SetBinContent(bin, pull)
    
    return [H_Adata, H_Abkg, H_pull]

############################################################################

# th1f_abcd_draw draws ABCD regions on a single histogram.

# Arguments:
# name: name of variable to be used in histogram names (string)
# title: name of variable to be displayed in histogram title (string)
# xlabel: x-axis label
# lhist: list of filled ABCD histograms

# Draws one plot with all four histograms, and prints in pdf file.

def th1f_abcd_draw(name, title, xlabel, lhist):
    H_A = lhist[0]
    H_B = lhist[1]
    H_C = lhist[2]
    H_D = lhist[3]

    C = TCanvas()

    C.cd()
    C.SetLogx()
    H_A.Scale(1.0, "width")
    H_B.Scale(1.0, "width")
    H_C.Scale(1.0, "width")
    H_D.Scale(1.0, "width")

    H_A.SetLineColor(kBlack)
    H_B.SetLineColor(kBlue)
    H_C.SetLineColor(kGreen+2)
    H_D.SetLineColor(kPink)

    H_A.SetLineStyle(1)         # solid
    H_B.SetLineStyle(2)         # dashed
    H_C.SetLineStyle(3)         # dotted
    H_D.SetLineStyle(4)         # dash-dotted

    # For setting histogram height.
    ymax_A = H_A.GetMaximum()
    aymax_abcd = array([H_A.GetMaximum(), H_B.GetMaximum(), \
                        H_C.GetMaximum(), H_D.GetMaximum()])
    ymax_abcd = aymax_abcd[argmax(aymax_abcd)]
    """ print()
    print("ymax of A for "+name+" = ", H_A.GetMaximum()) """

    H_A.SetTitle("ABCD regions of "+title+" (projection)")      # +" for J/Psi"
    H_A.GetXaxis().SetTitle(xlabel)
    H_A.GetYaxis().SetTitle("Events / bin width")
    H_A.GetXaxis().SetTitleSize(0.05)
    H_A.GetYaxis().SetTitleSize(0.05)
    H_A.GetXaxis().SetLabelSize(0.06)
    H_A.GetYaxis().SetLabelSize(0.06)
    H_A.SetMaximum(1.2 * ymax_abcd)
    H_A.SetStats(0)

    H_A.Draw("hist")
    H_B.Draw("hist same")
    H_C.Draw("hist same")
    H_D.Draw("hist same")

    L = TLegend(0.7, 0.75, 0.9, 0.9)
    L.AddEntry(H_A, "Signal region A", "l")
    L.AddEntry(H_B, "Control region B", "l")
    L.AddEntry(H_C, "Control region C", "l")
    L.AddEntry(H_D, "Control region D", "l")
    L.Draw()

    C.Print("h_proj_abcd_0.05_alpha0.1_"+name+"_2.pdf")
    H_A.SetMaximum(ymax_A)                  # Resets maximum of H_A.
    L.Clear()
    C.Clear()
    C.Update()

    # Number of events in ABCD regions.
    eventsA = H_A.Integral(1, H_A.GetNbinsX())
    eventsB = H_B.Integral(1, H_B.GetNbinsX())
    eventsC = H_C.Integral(1, H_C.GetNbinsX())
    eventsD = H_D.Integral(1, H_D.GetNbinsX())
    eventstot = eventsA + eventsB + eventsC + eventsD

    print("Done m_asymm 0.05 ABCD histogram of "+name)
    print("events shown in A of "+name+" = ", eventsA)
    print("events shown in B of "+name+" = ", eventsB)
    print("events shown in C of "+name+" = ", eventsC)
    print("events shown in D of "+name+" = ", eventsD)
    print("total events shown in ABCD of "+name+" = ", eventstot)
 
############################################################################

# th1f_apull_draw draws histogram of region A data vs background,
# along with pull plot.

# Arguments:
# name: name of variable to be used in histogram names (string)
# title: name of variable to be displayed in histogram title (string)
# xlabel: x-axis label of histogram and pull plot
# lplots: list of *filled* histograms for region A data, bkg, and pull

# Draws histogram with pull plot, and prints in pdf file.

def th1f_apull_draw(name, title, xlabel, lplots):
    H_Adata = lplots[0]
    H_Abkg = lplots[1]
    H_Apull = lplots[2]

    C = TCanvas()
    P_hist = TPad("h_"+name+"_data_bkg", "", \
                    0, 0.4, 1, 1)     # top 60% of canvas
    P_pull = TPad("p_"+name, "", \
                    0, 0, 1, 0.4)     # bottom 40% of canvas

    P_hist.SetBottomMargin(0)
    P_hist.SetLeftMargin(0.15)
    P_pull.SetTopMargin(0)
    P_pull.SetBottomMargin(0.2)
    P_pull.SetLeftMargin(0.15)

    P_hist.Draw()
    P_pull.Draw()

    # Region A background and data.
    P_hist.cd()
    P_hist.SetLogx()
    H_Abkg.Scale(1.0, "width")
    H_Adata.Scale(1.0, "width")
    
    H_Abkg.SetLineStyle(1)              # solid line
    H_Abkg.SetFillStyle(1001)           # solid fill
    H_Abkg.SetLineColor(0)
    H_Abkg.SetFillColor(kOrange)
    H_Adata.SetLineColor(kBlack)
    H_Adata.SetMarkerColor(kBlack)

    # For setting histogram height.
    amax_apull = array([H_Adata.GetMaximum(), H_Abkg.GetMaximum()])
    ymax_apull = amax_apull[argmax(amax_apull)]
    """ print()
    print("ymax of A data and bkg for "+name+" = ", ymax_apull)
    print("ymax of A data for "+name+" = ", H_Adata.GetMaximum())
    print("ymax of A bkg for "+name+" = ", H_Abkg.GetMaximum()) """

    H_Abkg.SetTitle(title+" (projection)")                # +" for J/Psi"
    H_Abkg.GetYaxis().SetTitle("Events / bin width")
    H_Abkg.GetYaxis().SetTitleSize(0.05)
    H_Abkg.GetYaxis().SetLabelSize(0.04)
    H_Abkg.SetMaximum(1.2 * ymax_apull)
    H_Abkg.SetStats(0)

    H_Abkg.Draw("hist")
    H_Adata.Draw("e0 same")
  
    L = TLegend(0.8, 0.7, 0.9, 0.9)
    L.AddEntry(H_Adata, "Data", "lep")
    L.AddEntry(H_Abkg, "Background", "f")
    L.Draw()


    # Pull plot.
    P_pull.cd()
    P_pull.SetLogx()

    H_Apull.SetStats(0)
    H_Apull.GetXaxis().SetTitle(xlabel)
    H_Apull.GetYaxis().SetTitle("Pull")
    H_Apull.GetXaxis().SetTitleSize(0.07)
    H_Apull.GetYaxis().SetTitleSize(0.07)
    H_Apull.GetXaxis().SetLabelSize(0.06)
    H_Apull.GetYaxis().SetLabelSize(0.06)
    H_Apull.GetXaxis().SetTitleOffset(1.2)
    H_Apull.GetYaxis().SetTitleOffset(0.7)

    H_Apull.Draw("e0")

    # Drawing line at y=0.
    xmin = H_Apull.GetXaxis().GetXmin()
    xmax = H_Apull.GetXaxis().GetXmax()
    line = TLine(xmin, 0., xmax, 0.)
    line.Draw()

    C.Print("h_proj_apull_0.05_alpha0.1_"+name+"_2.pdf")
    L.Clear()
    P_hist.Clear()
    P_hist.Update()
    P_pull.Clear()
    P_pull.Update()
    C.Clear()
    C.Update()

    # Number of events in A data and background.
    eventsA_data = H_Adata.Integral(1, H_Adata.GetNbinsX())
    eventsA_bkg = H_Abkg.Integral(1, H_Abkg.GetNbinsX())

    print("Done m_asymm 0.05 A pull plot of "+name)
    print("events in A data of "+name+" = ", eventsA_data)
    print("events in A bkg of "+name+" = ", eventsA_bkg)

############################################################################

# th2f_abcd_make creates 2d histograms of ABCD regions for given variables.

# Arguments:
# namex: name of x variable to be used in histogram names (string)
# namey: name of y variable to be used in histogram names (string)
# nbinsx: number of bins along x (int)
# abinsx: array of bin edges along x (array of floats)
# nbinsy: number of bins along y (int)
# abinsy: array of bin edges along y (array of floats)

# Returns:
# list of empty 2d histograms for ABCD regions

def th2f_abcd_make(namex, namey, nbinsx, abinsx, nbinsy, abinsy):
    H_A = TH2F("h_"+namey+"_vs_"+namex+"_A", "", \
              nbinsx, abinsx, nbinsy, abinsy)
    H_B = TH2F("h_"+namey+"_vs_"+namex+"_B", "", \
              nbinsx, abinsx, nbinsy, abinsy)
    H_C = TH2F("h_"+namey+"_vs_"+namex+"_C", "", \
              nbinsx, abinsx, nbinsy, abinsy)
    H_D = TH2F("h_"+namey+"_vs_"+namex+"_D", "", \
              nbinsx, abinsx, nbinsy, abinsy)

    return [H_A, H_B, H_C, H_D]

############################################################################

# th2f_abcd_projx_make projects 2d histograms of ABCD regions 
# onto the x axis.

# Arguments:
# namex: name of x variable to be used in histogram names (string)
# lhist2d: list of *filled* 2d histograms of ABCD regions.

# Returns:
# list of filled 1d histograms of the x-projections of ABCD regions.

def th2f_abcd_projx_make(namex, lhist2d):
    H_A = lhist2d[0]
    H_B = lhist2d[1]
    H_C = lhist2d[2]
    H_D = lhist2d[3]

    H_Aprojx = H_A.ProjectionX("h_proj_"+namex+"_A")
    H_Bprojx = H_B.ProjectionX("h_proj_"+namex+"_B")
    H_Cprojx = H_C.ProjectionX("h_proj_"+namex+"_C")
    H_Dprojx = H_D.ProjectionX("h_proj_"+namex+"_D")

    return [H_Aprojx, H_Bprojx, H_Cprojx, H_Dprojx]

############################################################################

# th2f_apull_projx_make creates 2d histogram and pull plot of 
# region A data vs bkg, projected onto the x-axis.

# Arguments:
# namex: name of x variable to be used in histogram names (string)
# namey: name of y variable to be used in histogram names (string)
# lhist2d: list of *filled* 2d histograms of ABCD regions.
# nbinsx: number of bins along x (int)
# abinsx: array of bin edges along x (array of floats)

# Returns:
# list of filled 2d histograms of the x-projections of 
# region A data, background, and pull.

def th2f_apull_projx_make(namex, namey, lhist2d, nbinsx, abinsx):
    H_A = lhist2d[0]
    H_B = lhist2d[1]
    H_C = lhist2d[2]
    H_D = lhist2d[3]

    # Creating 2d histograms for region A data and background.
    H_Adata = H_A.Clone("h_"+namey+"_vs_"+namex+"_Adata")    # Clone A.
    H_Abkg = H_C.Clone("h_"+namey+"_vs_"+namex+"_Abkg")      # Clone C.
    H_Abkg.Divide(H_D)                                      # Divide by D.
    H_Abkg.Multiply(H_B)                                    # Multiply by B.

    # Projecting region A data and bkg onto x-axis.
    H_Adata_projx = H_Adata.ProjectionX("h_proj_"+namex+"_Adata")
    H_Abkg_projx = H_Abkg.ProjectionX("h_proj_"+namex+"_Abkg")

    # Creating x-projection pull plot for region A data and background.
    H_pull_projx = TH1F("pull_"+namex, "", \
                   nbinsx, abinsx)

    for binx in range(1, nbinsx + 1):
        data = H_Adata_projx.GetBinContent(binx)
        if data == 0: continue 
        bkg = H_Abkg_projx.GetBinContent(binx)
        uncert = TMath.Sqrt(data)
        pull = (data - bkg) / uncert
        H_pull_projx.SetBinContent(binx, pull)
    
    return [H_Adata_projx, H_Abkg_projx, H_pull_projx]


#########################################################################
######  Algorithm for filtering events and filling in histograms.  ######

# Large dataset (ScoutingMultiEta_NoBias.root)
F = TFile("/afs/crc.nd.edu/user/c/cyang22/CMSSW_10_6_19_patch2/src/dimuon_pairs/ScoutingMultiEta_NoBias.root")
T = F.Get("Events")

m = 0.1057					# Mass of muon (GeV/c^2)


# Creating 1d histograms.
#nbins_mavg, abins_mavg = percent_bins(2.5, 3.5, 0.01)     # J/Psi (3.1 GeV).
""" nbins_mavg, abins_mavg = percent_bins(1e-1, 5e2, 0.05)
labcd_mavg = th1f_abcd_make("mavg", nbins_mavg, abins_mavg)
print("Done labcd_mavg, 0.05 binning")

nbins_mfour, abins_mfour = percent_bins(5e-2, 1e3, 0.05)
labcd_mfour = th1f_abcd_make("mfour", nbins_mfour, abins_mfour)
print("Done labcd_mfour, 0.05 binning")

nbins_alpha, abins_alpha = percent_bins(1e-3, 5., 0.05)
labcd_alpha = th1f_abcd_make("alpha", nbins_alpha, abins_alpha)
print("Done labcd_alpha, 0.05 binning")

nbins_ptlead, abins_ptlead = percent_bins(1., 1e3, 0.05)
labcd_ptlead = th1f_abcd_make("ptlead", nbins_ptlead, abins_ptlead)
print("Done labcd_ptlead, 0.05 binning")

nbins_ptlast, abins_ptlast = percent_bins(1., 1e2, 0.05)
labcd_ptlast = th1f_abcd_make("ptlast", nbins_ptlast, abins_ptlast)
print("Done labcd_ptlast, 0.05 binning")

nbins_pt_dimlead, abins_pt_dimlead = percent_bins(5e-2, 5e2, 0.05)
labcd_pt_dimlead = th1f_abcd_make("pt_dimlead", \
                              nbins_pt_dimlead, abins_pt_dimlead)
print("Done labcd_pt_dimlead, 0.05 binning")

nbins_pt_dimlast, abins_pt_dimlast = percent_bins(1e-2, 1e2, 0.05)
labcd_pt_dimlast = th1f_abcd_make("pt_dimlast", \
                              nbins_pt_dimlast, abins_pt_dimlast)
print("Done labcd_pt_dimlast, 0.05 binning")

nbins_dR, abins_dR = percent_bins(1e-5, 50., 0.1)
labcd_dR = th1f_abcd_make("dR", nbins_dR, abins_dR)
print("Done labcd_dR, 0.1 binning")

nbins_deta, abins_deta = percent_bins(1e-8, 50., 0.1)
labcd_deta = th1f_abcd_make("deta", nbins_deta, abins_deta)
print("Done labcd_deta, 0.1 binning") """

# Creating 2d histograms.
abins_alpha = arange(1e-3, 1.1, 0.1)
nbins_alpha = len(abins_alpha) -1

#nbins_mavg, abins_mavg = percent_bins(2.5, 3.5, 0.01)     # J/Psi (3.1 GeV).
nbins_mavg, abins_mavg = percent_bins(1e-1, 5e2, 0.05)
labcd_alpha_vs_mavg = th2f_abcd_make("mavg", "alpha", nbins_mavg, \
                            abins_mavg, nbins_alpha, abins_alpha)
print("Done labcd_alpha_vs_mavg, 0.05 binning")

nbins_mfour, abins_mfour = percent_bins(5e-2, 1e3, 0.05)
labcd_alpha_vs_mfour = th2f_abcd_make("mfour", "alpha", nbins_mfour, \
                            abins_mfour, nbins_alpha, abins_alpha)
print("Done labcd_alpha_vs_mfour, 0.05 binning")

nbins_ptlead, abins_ptlead = percent_bins(1., 1e3, 0.05)
labcd_alpha_vs_ptlead = th2f_abcd_make("ptlead", "alpha", nbins_ptlead, \
                            abins_ptlead, nbins_alpha, abins_alpha)
print("Done labcd_alpha_vs_ptlead, 0.05 binning")

nbins_ptlast, abins_ptlast = percent_bins(1., 1e2, 0.05)
labcd_alpha_vs_ptlast = th2f_abcd_make("ptlast", "alpha", nbins_ptlast, \
                            abins_ptlast, nbins_alpha, abins_alpha)
print("Done labcd_alpha_vs_ptlast, 0.05 binning")

nbins_pt_dimlead, abins_pt_dimlead = percent_bins(5e-2, 5e2, 0.05)
labcd_alpha_vs_pt_dimlead = th2f_abcd_make("pt_dimlead", "alpha", \
                            nbins_pt_dimlead, abins_pt_dimlead, \
                                nbins_alpha, abins_alpha)
print("Done labcd_alpha_vs_pt_dimlead, 0.05 binning")

nbins_pt_dimlast, abins_pt_dimlast = percent_bins(1e-2, 1e2, 0.05)
labcd_alpha_vs_pt_dimlast = th2f_abcd_make("pt_dimlast", "alpha", \
                            nbins_pt_dimlast, abins_pt_dimlast, \
                                nbins_alpha, abins_alpha)
print("Done labcd_alpha_vs_pt_dimlast, 0.05 binning")

nbins_dR, abins_dR = percent_bins(1e-5, 50., 0.1)
labcd_alpha_vs_dR = th2f_abcd_make("dR", "alpha", \
                            nbins_dR, abins_dR, nbins_alpha, abins_alpha)
print("Done labcd_alpha_vs_dR, 0.1 binning")

nbins_deta, abins_deta = percent_bins(1e-8, 50., 0.1)
labcd_alpha_vs_deta = th2f_abcd_make("deta", "alpha", \
                            nbins_deta, abins_deta, \
                                nbins_alpha, abins_alpha)
print("Done labcd_alpha_vs_deta, 0.1 binning")

#########################################################################

k = -1
for e in T:
    k += 1
    """ if k == 2000: break    #20000
    if k % 5000 == 0:
        print("k = ", k) """

    n_muons = len(T.muon_pt)
    if n_muons < 4: continue
    #if n_muons < 6: continue

    # Listing all muon four-vectors and charges in the event.
    afourv = empty(n_muons, object)
    aq = empty(n_muons, int)  
    apt = empty(n_muons, float)      
    for n in range(n_muons):
        n_pt = T.muon_pt[n]
        n_eta = T.muon_eta[n]
        n_phi = T.muon_phi[n]
        n_q = T.muon_q[n]
        
        V = TLorentzVector()
        V.SetPtEtaPhiM(n_pt, n_eta, n_phi, m)
        afourv[n] = V
        aq[n] = n_q
        apt[n] = n_pt
        # Indices of all arrays above should match.

    """ print("n_muons = ", n_muons)
    print("afourv")
    print(afourv)
    print("aq")
    print(aq)
    print("apt")
    print(apt) """

    # Listing indices for all possible dimuons. ("pairs")
    lpairs = []
    for i in range(0, len(afourv)):
        for j in range(i+1, len(afourv)):
            lpairs.append([i,j])        

    """ print("lpairs")
    print(lpairs) """

    # Listing indices for all possible pairs of dimuons. ("pairs of pairs")
    lpairs_pairs = []
    for p in range(0, len(lpairs)):
        for q in range(p, len(lpairs)):
            # If no two elements in the pairs p,q are the same.
            if lpairs[q][0] != lpairs[p][0] \
                and lpairs[q][0] != lpairs[p][1] \
                and lpairs[q][1] != lpairs[p][0] \
                    and lpairs[q][1] != lpairs[p][1]:
                 lpairs_pairs.append([p,q]) 

    """ print("lpairs_pairs")
    print(lpairs_pairs) """

    # Listing the four-vectors of dimuon pairs according to their indices.
    avectors = empty((len(lpairs_pairs),2), object)
    i = 0
    for pp in lpairs_pairs:
        p11 = lpairs[pp[0]][0]          # index of pair 1 muon 1
        p12 = lpairs[pp[0]][1]          # index of pair 1 muon 2

        p21 = lpairs[pp[1]][0]
        p22 = lpairs[pp[1]][1]

        """ print()
        print("p11 vector:", afourv[p11])
        print("p12 vector:", afourv[p12])
        print("p21 vector:", afourv[p21])
        print("p22 vector:", afourv[p22]) """
        
        v1 = afourv[p11] + afourv[p12]  # dimuon four-vector 1
        v2 = afourv[p21] + afourv[p22]  # dimuon four-vector 2
        avectors[i, 0] = v1
        avectors[i, 1] = v2

        i += 1

    """ print("avectors")
    print(avectors) """

    # Finding dimuon pair with the lowest mass asymmetry.
    # Indices of avectors and am_asymm should match.
    am_asymm = empty(len(avectors), float)
    i = 0
    for v_dimuon in avectors:
        m1 = v_dimuon[0].M()
        m2 = v_dimuon[1].M()
        am_asymm[i] = abs((m1-m2)/(m1+m2))
        
        i += 1

    """ print("am_asymm")
    print(am_asymm) """

    i_min = argmin(am_asymm)   # index of lowest m_asymm.

    """ print("i_min = ", i_min) """

    # Four-vectors of dimuon pair.
    vdimuon1 = avectors[i_min, 0]
    vdimuon2 = avectors[i_min, 1]

    # Mass asymmetry of dimuon pair
    m_asymm = am_asymm[i_min]

    # Total charge of the four muons.
    q1 = aq[p11] + aq[p12]              # Charge of dimuon 1.
    q2 = aq[p21] + aq[p22]              # Charge of dimuon 2.
    q_four = abs(q1) + abs(q2)
    # Absolute values make q_four >= 0.
    """ print()
    print("q1 = ", q1)
    print(aq[p11], aq[p12])
    print("q2 = ", q2)
    print(aq[p21], aq[p22])
    print("q_four = ", q_four) """

    # Individual dimuon mass.
    m_dimuon1 = vdimuon1.M()
    m_dimuon2 = vdimuon2.M()

    # Average mass of dimuon pair.
    m_avg = (m_dimuon1 + m_dimuon2)/2

    # Four-muon mass.
    V_fourmuon = vdimuon1 + vdimuon2
    m_four = V_fourmuon.M()

    # Alpha (ratio of dimuon mass over four-muon mass).
    alpha = float(m_avg) / float(m_four)

    # pT of leading and last muons.
    apt_four = [apt[p11], apt[p12], apt[p21], apt[p22]]
    pt_lead = max(apt_four)
    pt_last = min(apt_four)

    # pT of leading and last dimuons.
    pt_dimuon = [vdimuon1.Pt(), vdimuon2.Pt()]
    pt_dimlead = max(pt_dimuon)
    pt_dimlast = min(pt_dimuon)

    # deltaR of dimuon pair.
    deltaR = vdimuon1.DeltaR(vdimuon2)

    # delta_eta of dimuon pair.
    eta1 = vdimuon1.Eta()
    eta2 = vdimuon2.Eta()
    delta_eta = abs(eta1 - eta2)

    """ print()
    print("m_asymm = ", m_asymm)
    print("m_avg = ", m_avg)
    print("m_four = ", m_four)
    print("q1 = ", q1)
    print("q2 = ", q2)
    print("q_four = ", q_four)
    print()
    print("apt_four")
    print(apt_four)
    print("pt_lead = ", pt_lead) """

    # Sorting event into A,B,C,D regions and saving to arrays.
    #if m_avg < 3.0 or m_avg > 3.2: continue     # J/Psi (3.1 GeV)
    
    if m_asymm <= 0.05:    
        if q_four < 1:
            # Filling 1d histograms.
            """ labcd_mavg[0].Fill(m_avg)
            labcd_mfour[0].Fill(m_four)
            labcd_alpha[0].Fill(alpha)
            labcd_ptlead[0].Fill(pt_lead)
            labcd_ptlast[0].Fill(pt_last)
            labcd_pt_dimlead[0].Fill(pt_dimlead)
            labcd_pt_dimlast[0].Fill(pt_dimlast)
            labcd_dR[0].Fill(deltaR)
            labcd_deta[0].Fill(delta_eta) """

            # Filling 2d histograms.
            labcd_alpha_vs_mavg[0].Fill(m_avg, alpha)
            labcd_alpha_vs_mfour[0].Fill(m_four, alpha)
            labcd_alpha_vs_ptlead[0].Fill(pt_lead, alpha)
            labcd_alpha_vs_ptlast[0].Fill(pt_last, alpha)
            labcd_alpha_vs_pt_dimlead[0].Fill(pt_dimlead, alpha)
            labcd_alpha_vs_pt_dimlast[0].Fill(pt_dimlast, alpha)
            labcd_alpha_vs_dR[0].Fill(deltaR, alpha)
            labcd_alpha_vs_deta[0].Fill(delta_eta, alpha)

            """ print("Region A")
            print() """

        elif q_four >= 1:
            # Filling 1d histograms.
            """ labcd_mavg[1].Fill(m_avg)
            labcd_mfour[1].Fill(m_four)
            labcd_alpha[1].Fill(alpha)
            labcd_ptlead[1].Fill(pt_lead)
            labcd_ptlast[1].Fill(pt_last)
            labcd_pt_dimlead[1].Fill(pt_dimlead)
            labcd_pt_dimlast[1].Fill(pt_dimlast)
            labcd_dR[1].Fill(deltaR)
            labcd_deta[1].Fill(delta_eta) """

            # Filling 2d histograms.
            labcd_alpha_vs_mavg[1].Fill(m_avg, alpha)
            labcd_alpha_vs_mfour[1].Fill(m_four, alpha)
            labcd_alpha_vs_ptlead[1].Fill(pt_lead, alpha)
            labcd_alpha_vs_ptlast[1].Fill(pt_last, alpha)
            labcd_alpha_vs_pt_dimlead[1].Fill(pt_dimlead, alpha)
            labcd_alpha_vs_pt_dimlast[1].Fill(pt_dimlast, alpha)
            labcd_alpha_vs_dR[1].Fill(deltaR, alpha)
            labcd_alpha_vs_deta[1].Fill(delta_eta, alpha)

            """ print("Region B")
            print() """
    
    elif m_asymm > 0.05:
        if q_four < 1:
            # Filling 1d histograms.
            """ labcd_mavg[2].Fill(m_avg)
            labcd_mfour[2].Fill(m_four)
            labcd_alpha[2].Fill(alpha)
            labcd_ptlead[2].Fill(pt_lead)
            labcd_ptlast[2].Fill(pt_last)
            labcd_pt_dimlead[2].Fill(pt_dimlead)
            labcd_pt_dimlast[2].Fill(pt_dimlast)
            labcd_dR[2].Fill(deltaR)
            labcd_deta[2].Fill(delta_eta) """

            # Filling 2d histograms.
            labcd_alpha_vs_mavg[2].Fill(m_avg, alpha)
            labcd_alpha_vs_mfour[2].Fill(m_four, alpha)
            labcd_alpha_vs_ptlead[2].Fill(pt_lead, alpha)
            labcd_alpha_vs_ptlast[2].Fill(pt_last, alpha)
            labcd_alpha_vs_pt_dimlead[2].Fill(pt_dimlead, alpha)
            labcd_alpha_vs_pt_dimlast[2].Fill(pt_dimlast, alpha)
            labcd_alpha_vs_dR[2].Fill(deltaR, alpha)
            labcd_alpha_vs_deta[2].Fill(delta_eta, alpha)

            """ print("Region C")
            print() """
        
        elif q_four >= 1:
            # Filling 1d histograms.
            """ labcd_mavg[3].Fill(m_avg)
            labcd_mfour[3].Fill(m_four)
            labcd_alpha[3].Fill(alpha)
            labcd_ptlead[3].Fill(pt_lead)
            labcd_ptlast[3].Fill(pt_last)
            labcd_pt_dimlead[3].Fill(pt_dimlead)
            labcd_pt_dimlast[3].Fill(pt_dimlast)
            labcd_dR[3].Fill(deltaR)
            labcd_deta[3].Fill(delta_eta) """

            # Filling 2d histograms.
            labcd_alpha_vs_mavg[3].Fill(m_avg, alpha)
            labcd_alpha_vs_mfour[3].Fill(m_four, alpha)
            labcd_alpha_vs_ptlead[3].Fill(pt_lead, alpha)
            labcd_alpha_vs_ptlast[3].Fill(pt_last, alpha)
            labcd_alpha_vs_pt_dimlead[3].Fill(pt_dimlead, alpha)
            labcd_alpha_vs_pt_dimlast[3].Fill(pt_dimlast, alpha)
            labcd_alpha_vs_dR[3].Fill(deltaR, alpha)
            labcd_alpha_vs_deta[3].Fill(delta_eta, alpha)

            """ print("Region D")
            print() """
    

#########################################################################

# Arrays of filled 1d histograms.
""" labcd_mavg_filled = [labcd_mavg[0], labcd_mavg[1], labcd_mavg[2], \
                     labcd_mavg[3]]
print("Done labcd_mavg_filled")

labcd_mfour_filled = [labcd_mfour[0], labcd_mfour[1], labcd_mfour[2], \
                      labcd_mfour[3]]
print("Done labcd_mfour_filled")

labcd_alpha_filled = [labcd_alpha[0], labcd_alpha[1], labcd_alpha[2], \
                      labcd_alpha[3]]
print("Done labcd_alpha_filled")

labcd_ptlead_filled = [labcd_ptlead[0], labcd_ptlead[1], labcd_ptlead[2], \
                      labcd_ptlead[3]]
print("Done labcd_ptlead_filled")

labcd_ptlast_filled = [labcd_ptlast[0], labcd_ptlast[1], labcd_ptlast[2], \
                      labcd_ptlast[3]]
print("Done labcd_ptlast_filled")

labcd_pt_dimlead_filled = [labcd_pt_dimlead[0], labcd_pt_dimlead[1], \
                       labcd_pt_dimlead[2], labcd_pt_dimlead[3]]
print("Done labcd_pt_dimlead_filled")

labcd_pt_dimlast_filled = [labcd_pt_dimlast[0], labcd_pt_dimlast[1], \
                       labcd_pt_dimlast[2], labcd_pt_dimlast[3]]
print("Done labcd_pt_dimlast_filled")

labcd_dR_filled = [labcd_dR[0], labcd_dR[1], labcd_dR[2], labcd_dR[3]]
print("Done labcd_dR_filled")

labcd_deta_filled = [labcd_deta[0], labcd_deta[1], labcd_deta[2], \
                     labcd_deta[3]]
print("Done labcd_deta_filled") """

# Arrays of filled 2d histograms.
labcd_alpha_vs_mavg_filled = [labcd_alpha_vs_mavg[0], \
                              labcd_alpha_vs_mavg[1], \
                                labcd_alpha_vs_mavg[2], \
                                    labcd_alpha_vs_mavg[3]]
print("Done labcd_alpha_vs_mavg_filled")

labcd_alpha_vs_mfour_filled = [labcd_alpha_vs_mfour[0], \
                              labcd_alpha_vs_mfour[1], \
                                labcd_alpha_vs_mfour[2], \
                                    labcd_alpha_vs_mfour[3]]
print("Done labcd_alpha_vs_mfour_filled")

labcd_alpha_vs_ptlead_filled = [labcd_alpha_vs_ptlead[0], \
                              labcd_alpha_vs_ptlead[1], \
                                labcd_alpha_vs_ptlead[2], \
                                    labcd_alpha_vs_ptlead[3]]
print("Done labcd_alpha_vs_ptlead_filled")

labcd_alpha_vs_ptlast_filled = [labcd_alpha_vs_ptlast[0], \
                              labcd_alpha_vs_ptlast[1], \
                                labcd_alpha_vs_ptlast[2], \
                                    labcd_alpha_vs_ptlast[3]]
print("Done labcd_alpha_vs_ptlast_filled")

labcd_alpha_vs_pt_dimlead_filled = [labcd_alpha_vs_pt_dimlead[0], \
                              labcd_alpha_vs_pt_dimlead[1], \
                                labcd_alpha_vs_pt_dimlead[2], \
                                    labcd_alpha_vs_pt_dimlead[3]]
print("Done labcd_alpha_vs_pt_dimlead_filled")

labcd_alpha_vs_pt_dimlast_filled = [labcd_alpha_vs_pt_dimlast[0], \
                              labcd_alpha_vs_pt_dimlast[1], \
                                labcd_alpha_vs_pt_dimlast[2], \
                                    labcd_alpha_vs_pt_dimlast[3]]
print("Done labcd_alpha_vs_pt_dimlast_filled")

labcd_alpha_vs_dR_filled = [labcd_alpha_vs_dR[0], \
                              labcd_alpha_vs_dR[1], \
                                labcd_alpha_vs_dR[2], \
                                    labcd_alpha_vs_dR[3]]
print("Done labcd_alpha_vs_dR_filled")


labcd_alpha_vs_deta_filled = [labcd_alpha_vs_deta[0], \
                              labcd_alpha_vs_deta[1], \
                                labcd_alpha_vs_deta[2], \
                                    labcd_alpha_vs_deta[3]]
print("Done labcd_alpha_vs_deta_filled")


# Projecting 2d histograms onto x-axis.
""" labcd_proj_mavg = th2f_abcd_projx_make("proj_mavg", \
                                labcd_alpha_vs_mavg_filled) """

############################################################################
#########################  Drawing histograms.  ############################

# 1d histograms of ABCD regions.
""" th1f_abcd_draw("mavg", "average dimuon mass", \
                       "Average dimuon mass (GeV)", \
                        labcd_mavg_filled)

th1f_abcd_draw("mfour", "four-muon mass", \
                       "Four-muon mass (GeV)", \
                        labcd_mfour_filled)

th1f_abcd_draw("alpha", "#alpha", \
                       "#alpha (ratio of average dimuon mass to four-muon mass)", \
                        labcd_alpha_filled)

th1f_abcd_draw("ptlead", "leading muon pT", \
                       "pT of leading muon (GeV)", \
                        labcd_ptlead_filled)

th1f_abcd_draw("ptlast", "slowest muon pT", \
                       "pT of slowest muon (GeV)", \
                        labcd_ptlast_filled)

th1f_abcd_draw("pt_dimlead", "leading dimuon pT", \
                       "pT of leading dimuon (GeV)", \
                        labcd_pt_dimlead_filled)

th1f_abcd_draw("pt_dimlast", "slowest dimuon pT", \
                       "pT of slowest dimuon (GeV)", \
                        labcd_pt_dimlast_filled)

th1f_abcd_draw("dR", "#Delta R of dimuon pair", \
                       "#Delta R of dimuon pair", \
                        labcd_dR)

th1f_abcd_draw("deta", "| #Delta #eta | of dimuon pair", \
                       "| #Delta #eta | of dimuon pair", \
                        labcd_deta) """


# 1d histograms of region A data, bkg, and pull.
""" lapull_mavg = th1f_apull_make("mavg", labcd_mavg_filled, \
                              nbins_mavg, abins_mavg)
th1f_apull_draw("mavg", "Average dimuon mass", \
                "Average dimuon mass (GeV)", \
                    lapull_mavg)


lapull_mfour = th1f_apull_make("mfour", labcd_mfour_filled, \
                               nbins_mfour, abins_mfour)
th1f_apull_draw("mfour", "Four-muon mass", \
                "Four-muon mass (GeV)", \
                    lapull_mfour)


lapull_alpha= th1f_apull_make("alpha", labcd_alpha_filled, \
                               nbins_alpha, abins_alpha)
th1f_apull_draw("alpha", "#alpha", \
                "#alpha (ratio of average dimuon mass to four-muon mass)", \
                    lapull_alpha)


lapull_ptlead= th1f_apull_make("ptlead", labcd_ptlead_filled, \
                               nbins_ptlead, abins_ptlead)
th1f_apull_draw("ptlead", "Leading muon pT", \
                "pT of leading muon (GeV)", \
                    lapull_ptlead)


lapull_ptlast= th1f_apull_make("ptlast", labcd_ptlast_filled, \
                               nbins_ptlast, abins_ptlast)
th1f_apull_draw("ptlast", "Slowest muon pT", \
                "pT of slowest muon (GeV)", \
                    lapull_ptlast)


lapull_pt_dimlead= th1f_apull_make("pt_dimlead", labcd_pt_dimlead_filled, \
                               nbins_pt_dimlead, abins_pt_dimlead)
th1f_apull_draw("pt_dimlead", "Leading dimuon pT", \
                "pT of leading dimuon (GeV)", \
                    lapull_pt_dimlead)


lapull_pt_dimlast= th1f_apull_make("pt_dimlast", labcd_pt_dimlast_filled, \
                               nbins_pt_dimlast, abins_pt_dimlast)
th1f_apull_draw("pt_dimlast", "Slowest dimuon pT", \
                "pT of slowest dimuon (GeV)", \
                    lapull_pt_dimlast)


lapull_dR = th1f_apull_make("dR", labcd_dR_filled, nbins_dR, abins_dR)
th1f_apull_draw("dR", "#Delta R of dimuon pair", \
                "#Delta R of dimuon pair", \
                    lapull_dR)


lapull_deta = th1f_apull_make("deta", labcd_deta_filled, \
                               nbins_deta, abins_deta)
th1f_apull_draw("deta", "| #Delta #eta | of dimuon pair", \
                "| #Delta #eta | of dimuon pair", \
                    lapull_deta) """



# Histograms of x-projections of ABCD regions.
""" th1f_abcd_draw("mavg", "average dimuon mass (projection)", \
                       "Average dimuon mass (GeV)", \
                        labcd_proj_mavg) """


# Histograms of x-projections of region A data, bkg, and pull.
lapull_proj_mavg = th2f_apull_projx_make("mavg", "alpha", \
                            labcd_alpha_vs_mavg_filled, \
                              nbins_mavg, abins_mavg)
th1f_apull_draw("mavg", "Average dimuon mass", \
                "Average dimuon mass (GeV)", \
                    lapull_proj_mavg)


lapull_proj_mfour = th2f_apull_projx_make("mfour", "alpha", \
                            labcd_alpha_vs_mfour_filled, \
                              nbins_mfour, abins_mfour)
th1f_apull_draw("mfour", "Four-muon mass", \
                "Four-muon mass (GeV)", \
                    lapull_proj_mfour)


lapull_proj_ptlead= th2f_apull_projx_make("ptlead", "alpha", \
                            labcd_alpha_vs_ptlead_filled, \
                               nbins_ptlead, abins_ptlead)
th1f_apull_draw("ptlead", "Leading muon pT", \
                "pT of leading muon (GeV)", \
                    lapull_proj_ptlead)


lapull_proj_ptlast= th2f_apull_projx_make("ptlast", "alpha", \
                            labcd_alpha_vs_ptlast_filled, \
                               nbins_ptlast, abins_ptlast)
th1f_apull_draw("ptlast", "Slowest muon pT", \
                "pT of slowest muon (GeV)", \
                    lapull_proj_ptlast)


lapull_proj_pt_dimlead= th2f_apull_projx_make("pt_dimlead", "alpha", \
                            labcd_alpha_vs_pt_dimlead_filled, \
                               nbins_pt_dimlead, abins_pt_dimlead)
th1f_apull_draw("pt_dimlead", "Leading dimuon pT", \
                "pT of leading dimuon (GeV)", \
                    lapull_proj_pt_dimlead)


lapull_proj_pt_dimlast= th2f_apull_projx_make("pt_dimlast", "alpha", \
                            labcd_alpha_vs_pt_dimlast_filled, \
                               nbins_pt_dimlast, abins_pt_dimlast)
th1f_apull_draw("pt_dimlast", "Slowest dimuon pT", \
                "pT of slowest dimuon (GeV)", \
                    lapull_proj_pt_dimlast)


lapull_proj_dR = th2f_apull_projx_make("dR", "alpha", \
                            labcd_alpha_vs_dR_filled, \
                                nbins_dR, abins_dR)
th1f_apull_draw("dR", "#Delta R of dimuon pair", \
                "#Delta R of dimuon pair", \
                    lapull_proj_dR)


lapull_proj_deta = th2f_apull_projx_make("deta", "alpha", \
                            labcd_alpha_vs_deta_filled, \
                               nbins_deta, abins_deta)
th1f_apull_draw("deta", "| #Delta #eta | of dimuon pair", \
                "| #Delta #eta | of dimuon pair", \
                    lapull_proj_deta)

