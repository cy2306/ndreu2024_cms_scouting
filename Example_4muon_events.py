import ROOT
from ROOT import *

F = TFile("/afs/crc.nd.edu/user/c/cyang22/CMSSW_10_6_19_patch2/src/ScoutingMultiEta_NoBias.root")
T = F.Get("Events") # loads the TTree we want to use

m = 0.1057					# mass of muon (GeV/c^2)

# Creating 2d histograms.
""" # Dimuon mass assymetry vs average dimuon mass of all possible pairs.
H1_avg_asymm = TH2F("dimuon_pairs1", "Dimuon pairs 1", \
                    500, 0., 50., 500, 0., 1.)
H2_avg_asymm = TH2F("dimuon_pairs2", "Dimuon pairs 2", \
                    500, 0., 50., 500, 0., 1.)
H3_avg_asymm = TH2F("dimuon_pairs3", "Dimuon pairs 3", \
                    500, 0., 50., 500, 0., 1.) """

# Histograms with lowest dimuon mass assymmetry.
""" H_asymm_avg = TH2F("lowest_masymm_vs_mavg", \
                   "Lowest dimuon m_asymm vs m_avg", \
                   100, 0., 50., 100, 0., 1.)
H_asymm_fourm = TH2F("lowest_masymm_vs_fourm", \
                     "Lowest dimuon m_asymm vs four-muon mass", \
                   100, 0., 100., 100, 0., 1.)
H_avg_fourm = TH2F("closest_mavg_vs_fourm", \
                   "Closest dimuon m_avg vs four-muon mass", \
                    100, 0., 100., 100, 0., 50.) """

# Higgs mass vs dimuon mass assymetry or average dimuon mass.
""" H_asymm_mH = TH2F("mH_vs_asymm", "mH vs dimuon mass assym", \
                  500, 0., 1., 500, 0., 100.)
H_avg_mH = TH2F("mH_vs_avg", "mH vs dimuon avg mass", \
                  500, 0., 1., 500, 0., 100.) """

# 1d histograms for pairs with lowest mass assymetry.
H_m = TH1F("individual_dimuon_mass", \
           "Individual dimuon mass in J/Psi region", \
           50, 2., 4.)
H_mavg = TH1F("average_dimuon_mass", \
              "Average mass of dimuon pair in J/Psi region", \
              50, 2., 4.)
H_fourm = TH1F("four_muon_mass", \
               "Four muon mass in J/Psi region", \
               50, 5., 100.)
H_deltaR = TH1F("double_J_Psi_angle,", "Angle between double J/Psi", \
                50, 0., 5.)


j = 0
for e in T:
    j += 1
    if j == 10000: break

    n_muons = T.muon_num 
    if n_muons != 4: continue   # Only looking at events with 4 muons.

    muon_fourv = []        # List for all 4 muon four-vectors in an event.
    for n in range(n_muons):
        n_pt = T.muon_pt[n]
        n_eta = T.muon_eta[n]
        n_phi = T.muon_phi[n]
        
        V = TLorentzVector()
        V.SetPtEtaPhiM(n_pt, n_eta, n_phi, m)
        muon_fourv.append(V)
    
    # Total mass of all four muons (mass of exotic Higgs).
    V_four = muon_fourv[0] + muon_fourv[1] + muon_fourv[2] + muon_fourv[3]
    m_four = V_four.M()

    # Mother particles "a" of all possible dimuon pairs.
    a01 = muon_fourv[0] + muon_fourv[1]
    a23 = muon_fourv[2] + muon_fourv[3]

    a02 = muon_fourv[0] + muon_fourv[2]
    a13 = muon_fourv[1] + muon_fourv[3]

    a03 = muon_fourv[0] + muon_fourv[3]
    a12 = muon_fourv[1] + muon_fourv[2]

    # Mass asymmetry and average mass of all possible dimuon pairs.
    m_asymm1 = abs((a01.M()-a23.M()) / (a01.M()+a23.M()))
    m_asymm2 = abs((a02.M()-a13.M()) / (a02.M()+a13.M()))
    m_asymm3 = abs((a03.M()-a12.M()) / (a03.M()+a12.M()))
    list_m_asymm = [m_asymm1, m_asymm2, m_asymm3]

    m_avg1 = (a01.M()+a23.M()) / 2
    m_avg2 = (a02.M()+a13.M()) / 2
    m_avg3 = (a03.M()+a12.M()) / 2
    list_m_avg = [m_avg1, m_avg2, m_avg3]

    i = list_m_asymm.index(min(list_m_asymm)) # Index of lowest m_asymm.
    
    # Setting restriction on mass asymmetry.
    if list_m_asymm[i] > 0.02: continue

    # Looking only at J/Psi events (3.1 GeV).
    if list_m_avg[i] < 3.0 or list_m_avg[i] > 3.2: continue

    # Filling histograms.
    """ H_asymm_avg.Fill(list_m_avg[i], list_m_asymm[i])
    H_asymm_fourm.Fill(m_four, list_m_asymm[i])
    H_avg_fourm.Fill(m_four, list_m_avg[i]) """

    """ H1_avg_asymm.Fill(m_avg1, m_asymm1)
    H2_avg_asymm.Fill(m_avg2, m_asymm2)
    H3_avg_asymm.Fill(m_avg3, m_asymm3) """

    # Individual dimuon masses.
    if i == 0:
        H_m.Fill(a01.M())
        H_m.Fill(a23.M())
    elif i == 1:
        H_m.Fill(a02.M())
        H_m.Fill(a13.M())
    elif i == 2:
        H_m.Fill(a03.M())
        H_m.Fill(a12.M())

    # Average dimuon mass and four-muon mass.
    H_mavg.Fill(list_m_avg[i])
    H_fourm.Fill(m_four)

    # Angle between double J/Psi.
    if i == 0:
        deltaR = a01.DeltaR(a23)
    elif i == 1:
        deltaR = a02.DeltaR(a13)
    elif i == 2:
        deltaR = a03.DeltaR(a12)

    H_deltaR.Fill(deltaR)


# Drawing and saving histograms.
""" C1 = TCanvas()
C1.Divide(1,3)
C1.cd(1)
H1_avg_asymm.Draw("colz")
H1_avg_asymm.GetXaxis().SetTitle("Average dimuon mass")
H1_avg_asymm.GetYaxis().SetTitle("Dimuon mass asymmetry")


C1.cd(2)
H2_avg_asymm.Draw("colz")
H2_avg_asymm.GetXaxis().SetTitle("Average dimuon mass")
H2_avg_asymm.GetYaxis().SetTitle("Dimuon mass asymmetry")

C1.cd(3)
H3_avg_asymm.Draw("colz")
H3_avg_asymm.GetXaxis().SetTitle("Average dimuon mass")
H3_avg_asymm.GetYaxis().SetTitle("Dimuon mass asymmetry")

C1.Print("hist_dimuon_masymm_vs_mavg.pdf") """

""" C2 = TCanvas()
C2.Divide(1,3)
gStyle.SetPalette(56)       # Color palette.

C2.cd(1)
H_asymm_avg.Draw("colz")
H_asymm_avg.GetXaxis().SetTitle("Average dimuon mass")
H_asymm_avg.GetYaxis().SetTitle("Dimuon mass asymmetry")

C2.cd(2)
H_asymm_fourm.Draw("colz")
H_asymm_fourm.GetXaxis().SetTitle("Four-muon mass")
H_asymm_fourm.GetYaxis().SetTitle("Dimuon mass asymmetry")

C2.cd(3)
H_avg_fourm.Draw("colz")
H_avg_fourm.GetXaxis().SetTitle("Four-muon mass")
H_avg_fourm.GetYaxis().SetTitle("Average dimuon mass")

C2.Print("hist_lowest_masymm.pdf") """


""" C3 = TCanvas()
C3.Divide(1,2)
C3.cd(1)
H_avg_mH.Draw("hist")
H_avg_mH.GetXaxis().SetTitle("Average dimuon mass")
H_avg_mH.GetYaxis().SetTutke("Higgs mass")

C3.cd(2)
H_asymm_mH.Draw("hist")
H_asymm_mH.GetXaxis().SetTitle("Dimuon mass asymmetry")
H_asymm_mH.GetYaxis().SetTutke("Higgs mass")

C3.Print("hist_mhiggs_vs_mdimuon.pdf") """


C4 = TCanvas()
C4.Divide(2,2)

C4.cd(1)
#H_m.SetStats(0)
H_m.Draw("hist")
H_m.GetXaxis().SetTitle("Individual dimuon mass (GeV)")
H_m.GetYaxis().SetTitle("Number of dimuon pairs")

C4.cd(2)
#H_mavg.SetStats(0)
H_mavg.Draw("hist")
H_mavg.GetXaxis().SetTitle("Average mass of dimuon pair (GeV)")
H_mavg.GetYaxis().SetTitle("Number of events")

C4.cd(3)
#H_fourm.SetStats(0)
H_fourm.Draw("hist")
H_fourm.GetXaxis().SetTitle("Mass of all four muons per event (GeV)")
H_fourm.GetYaxis().SetTitle("Number of events")

C4.cd(4)
#H_deltaR.SetStats(0)
H_deltaR.Draw("hist")
H_deltaR.GetXaxis().SetTitle("Angle (deltaR) between double J/Psi")
H_deltaR.GetYaxis().SetTitle("Number of events")

C4.Print("hist_double_J_psi_fulldata_test.pdf")