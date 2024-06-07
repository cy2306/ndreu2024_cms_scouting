# Looking at pairs of dimuons (from events with n >= 4 muons) 
#with the lowest mass asymmetry.
import ROOT
from ROOT import *

# Large dataset (ScoutingMultiEta_NoBias.root)
F = TFile("/afs/crc.nd.edu/user/c/cyang22/CMSSW_10_6_19_patch2/src/dimuon_pairs/ScoutingMultiEta_NoBias.root")
T = F.Get("Events")

# Small dataset (scouting_ntuple_1.root)
""" F = TFile("/afs/crc.nd.edu/user/c/cyang22/CMSSW_10_6_19_patch2/src/Example/scouting_ntuple_1.root")
T = F.Get("scoutingntuplizer") """


m = 0.1057					# mass of muon (GeV/c^2)

# Creating histograms for dimuon mass, m_avg, and four-muon mass.
H_m = TH1F("individual_dimuon_mass", \
           "Individual dimuon mass in Rho region", \
           100, 0., 3.)
H_mavg = TH1F("average_dimuon_mass", \
              "Average mass of dimuon pair in Rho region", \
              100, 0., 3.)
H_fourm = TH1F("four_muon_mass", \
               "Four muon mass in Rho region", \
               100, 0., 100.)
H_deltaR = TH1F("parent_deltaR", \
                "DeltaR between dimuon parents in Rho region", \
                50, 0., 5.)


k = 0
for e in T:
    k += 1
    #if k == 20000: break    #10000

    n_muons = len(T.muon_pt)
    if n_muons < 4: continue   # Only looking at events with n>=4 muons.
    #if n_muons != 4: continue   # For checking with previous results.

    #if n_muons < 6: continue    # For checking code with a single event.

    """ print("n_muons: ", n_muons)
    print("k", k) """

    # Listing all muon four-vectors in an event.
    lfourv = []
    lq = []        
    for n in range(n_muons):
        n_pt = T.muon_pt[n]
        n_eta = T.muon_eta[n]
        n_phi = T.muon_phi[n]
        n_q = T.muon_q[n]
        
        V = TLorentzVector()
        V.SetPtEtaPhiM(n_pt, n_eta, n_phi, m)
        lfourv.append(V)
        lq.append(n_q)
        # Indices of lfourv and lq should match.

    """ print("lfourv")
    print(lfourv) """

    # Listing indices for all possible dimuons. ("pairs")
    lpairs = []
    for i in range(0, len(lfourv)):
        for j in range(i+1, len(lfourv)):
            
            """ print("q sum: ", lq[i]+lq[j]) """

            # If dimuon has total charge of 0.
            if abs((lq[i] + lq[j]) - 0) < 1e-3:
                lpairs.append([i,j])
                
    """ print()
    print("lpairs")
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
    
    if len(lpairs_pairs) == 0: continue
    """ print()
    print("lpairs_pairs")
    print(lpairs_pairs) """


    # Listing the four-vectors of dimuon pairs according to their indices.
    lvectors = []
    for pp in lpairs_pairs:
        p11 = lpairs[pp[0]][0]          # index of pair 1 muon 1
        p12 = lpairs[pp[0]][1]          # index of pair 1 muon 2

        p21 = lpairs[pp[1]][0]
        p22 = lpairs[pp[1]][1]

        """ print()
        print("p11 vector:", lfourv[p11])
        print("p12 vector:", lfourv[p12])
        print("p21 vector:", lfourv[p21])
        print("p22 vector:", lfourv[p22]) """
        
        v1 = lfourv[p11] + lfourv[p12]  # dimuon four-vector 1
        v2 = lfourv[p21] + lfourv[p22]  # dimuon four-vector 2
        lvectors.append([v1,v2])
    """ print()
    print("lvectors")
    print(lvectors) """

    # Finding dimuon pair with the lowest mass asymmetry.
    # Indices of lvectors and lm_asymm should match.
    lm_asymm = []
    for v_dimuon in lvectors:
        m1 = v_dimuon[0].M()
        m2 = v_dimuon[1].M()
        m_asymm = abs((m1-m2)/(m1+m2))
        lm_asymm.append(m_asymm)

        """ print()
        print("dimuon m1: ", m1)
        print("dimuon m2: ", m2) """

    i = lm_asymm.index(min(lm_asymm))   # index of lowest m_asymm.

    # Setting restriction on mass asymmetry.
    if lm_asymm[i] > 0.01: continue

    """ print()
    print("m_asymm:")
    print(lm_asymm)
    print("lowest: ", lm_asymm[i]) """

    # Individual dimuon mass.
    m_dimuon1 = lvectors[i][0].M()
    m_dimuon2 = lvectors[i][1].M()

    # Average mass of dimuon pair.
    m_avg = (m_dimuon1 + m_dimuon2)/2

    # Four-muon mass.
    V_fourmuon = lvectors[i][0] + lvectors[i][1]
    m_four = V_fourmuon.M()

    """ print()
    print("dimuon1: ", lvectors[i][0])
    print("dimuon2: ", lvectors[i][1])
    print("four muon:", V_fourmuon) """

    # deltaR between dimuon four-vectors.
    deltaR = lvectors[i][0].DeltaR(lvectors[i][1])

    # Restricting m_avg.
    #if m_avg < 3.5 or m_avg > 4.0: continue     # J/Psi (3.1 GeV)
    #if m_avg < 10.0 or m_avg > 10.5: continue     # Upsilon (9.46 GeV)
    if m_avg < 0.7 or m_avg > 0.9: continue

    # Filling in histograms.
    H_m.Fill(m_dimuon1)
    H_m.Fill(m_dimuon2)
    H_mavg.Fill(m_avg)
    H_fourm.Fill(m_four)
    H_deltaR.Fill(deltaR)
    
    """ print()
    print("m_dimuon1: ", m_dimuon1)
    print("m_dimuon2: ", m_dimuon2)
    print("m_avg: ", m_avg)
    print("m_fourm:", m_four) """

    """ k = 1                   # For checking code with a single event.
    if k == 1: break """        


# Drawing histograms.
C1 = TCanvas()
C1.Divide(2,2)

C1.cd(1)
H_m.Draw("hist")
H_m.GetXaxis().SetTitle("Individual dimuon mass (GeV)")
H_m.GetYaxis().SetTitle("Number of dimuons")

C1.cd(2)
H_mavg.Draw("hist")
H_mavg.GetXaxis().SetTitle("Average mass of dimuon pair (GeV)")
H_mavg.GetYaxis().SetTitle("Number of events")

C1.cd(3)
H_fourm.Draw("hist")
H_fourm.GetXaxis().SetTitle("Mass of all four muons per event (GeV)")
H_fourm.GetYaxis().SetTitle("Number of events")

C1.cd(4)
H_deltaR.Draw("hist")
H_deltaR.GetXaxis().SetTitle("DeltaR between double J/Psi")
H_deltaR.GetYaxis().SetTitle("Number of events")

C1.Print("hist_largedata_0.01_q0_rho.pdf")
    
