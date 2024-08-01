"""
Sort events in a root file into ABCD regions for background estimates 
using the ABCD method.

1. For each event in the TTree, identify the dimuon pair of lowest 
    mass asymmetry.

2. Select events based on restrictions in number of muons, 
    mass asymmetry, dimuon charge, dimuon mass, etc. 
    Events are included if all of the following are met.
    
    n_muons >= 4
    m_asymm <= 0.05 or 0.075 <= m_asymm <= 0.5
    abs(q1 + q2 + q3 + q4) < 4 
    (all four muons do not have the same charge)


3. Sort selected events into ABCD regions based on four-muon charge
    and mass asymmetry of dimuon pair.

    A: m_asymm <= 0.05 and q_four == 0
    B: m_asymm <= 0.05 and q_four != 0
    C: 0.075 <= m_asymm <= 0.5 and q_four == 0
    D: 0.075 <= m_asymm <= 0.5 and q_four != 0


4. Further sort events in mass regions based on average dimuon mass,
    and append relevant variables to the correct list to be saved.

    "control1": 2.5 <= mavg < 3.05
    "JPsi": 3.05 <= mavg < 3.15
    "control2": 3.15 <= mavg <= 4.0
    ...
    

5. Create a json file for each ABCD region and mass region. 
    Each file is a dictionary of lists containing relevant variables, 
    formatted as follows.
    {
        "variable x": [event 1 x, event 2 x, event 3 x, ...],
        "variable y": [event 1 y, event 2 y, event 3 y, ...],
        "variable z": [event 1 z, event 2 z, event 3 z, ...],
        ...
    }

"""
"""
############################### Functions ###############################

1. dvar_create(mass_region, abcd_region, lvariables)
    Creates an empty dictionary to hold variables for an ABCD region.

    mass_region: str
        Name to specify the range of dimuon masses considered.

    abcd_region: str ("A", "B", "C", "D", or "ABCD")
        Name to specify the ABCD region.

    lvariables: list of str
        Names of all the relevant variables.

    Return: dict


2. abcd_sort(m_asymm, q_four)
    Assigns an event to an ABCD region.

    m_asymm: float
        Mass asymmetry between dimuons. |m12 - m34| / (m12 + m34)

    q_four: int
        Charge of four muons. |q1 + q2| + |q3 + q4|

    Return: str ("A", "B", "C", or "D")

"""

import ROOT
from ROOT import *
from numpy import empty, argmin
import json



def dvar_create(mass_region, abcd_region):

    dvar = {
    "l"+mass_region+"_mavg_"+abcd_region: [],
    "l"+mass_region+"_mfour_"+abcd_region: [],
    "l"+mass_region+"_alpha_"+abcd_region: [],
    "l"+mass_region+"_dR_"+abcd_region: [],
    "l"+mass_region+"_deta_"+abcd_region: [],
    "l"+mass_region+"_pt_lead_"+abcd_region: [],
    "l"+mass_region+"_pt_last_"+abcd_region: [],
    "l"+mass_region+"_pt_dimlead_"+abcd_region: [],
    "l"+mass_region+"_pt_dimlast_"+abcd_region: [],
}
    return dvar


def abcd_sort(m_asymm, q_four):
    if m_asymm <= 0.05:
        if q_four < 1:
            return "A"
        
        elif q_four >= 1:
            return "B"

    elif m_asymm >= 0.075:
        if q_four < 1:
            return "C"
        
        elif q_four >= 1:
            return "D"


if __name__ == "__main__":

    # Dataset: ScoutingMultiEta_NoBias.root
    F = TFile("/afs/crc.nd.edu/user/c/cyang22/CMSSW_10_6_19_patch2/src/dimuon_pairs/ScoutingMultiEta_NoBias.root")
    T = F.Get("Events")

    m = 0.1057					# Mass of muon (GeV/c^2)


    # Create empty dictionaries for ABCD regions of all variables.
    # All in J/Psi, control 1, and control 2.
    """ dvar_all_A = dvar_create("all", "A")
    dvar_all_B = dvar_create("all", "B")
    dvar_all_C = dvar_create("all", "C")
    dvar_all_D = dvar_create("all", "D")

   # Other masses excluding J/Psi region.
    dvar_other_A = dvar_create("other", "A")
    dvar_other_B = dvar_create("other", "B")
    dvar_other_C = dvar_create("other", "C")
    dvar_other_D = dvar_create("other", "D") """

    """ # Control 1 mass region (2.5-3.0 GeV)
    dvar_control1_A = dvar_create("control1", "A")
    dvar_control1_B = dvar_create("control1", "B")
    dvar_control1_C = dvar_create("control1", "C")
    dvar_control1_D = dvar_create("control1", "D")

    # J/Psi mass region (3.0-3.2 GeV)
    dvar_JPsi_A = dvar_create("JPsi", "A")
    dvar_JPsi_B = dvar_create("JPsi", "B")
    dvar_JPsi_C = dvar_create("JPsi", "C")
    dvar_JPsi_D = dvar_create("JPsi", "D")

    # Control 2 mass region (3.2-4.0 GeV)
    dvar_control2_A = dvar_create("control2", "A")
    dvar_control2_B = dvar_create("control2", "B")
    dvar_control2_C = dvar_create("control2", "C")
    dvar_control2_D = dvar_create("control2", "D") """

    # Other mass regions.
    """ # 2.0-3.0 GeV
    dvar_2gev_A = dvar_create("2gev", "A")
    dvar_2gev_B = dvar_create("2gev", "B")
    dvar_2gev_C = dvar_create("2gev", "C")
    dvar_2gev_D = dvar_create("2gev", "D")

    # 4.0-5.0 GeV
    dvar_4gev_A = dvar_create("4gev", "A")
    dvar_4gev_B = dvar_create("4gev", "B")
    dvar_4gev_C = dvar_create("4gev", "C")
    dvar_4gev_D = dvar_create("4gev", "D")

    # 5.0-6.0 GeV
    dvar_5gev_A = dvar_create("5gev", "A")
    dvar_5gev_B = dvar_create("5gev", "B")
    dvar_5gev_C = dvar_create("5gev", "C")
    dvar_5gev_D = dvar_create("5gev", "D") """

    # 8-9 GeV
    dvar_8gev_A = dvar_create("8gev", "A")
    dvar_8gev_B = dvar_create("8gev", "B")
    dvar_8gev_C = dvar_create("8gev", "C")
    dvar_8gev_D = dvar_create("8gev", "D")

    # 10-11 GeV
    dvar_10gev_A = dvar_create("10gev", "A")
    dvar_10gev_B = dvar_create("10gev", "B")
    dvar_10gev_C = dvar_create("10gev", "C")
    dvar_10gev_D = dvar_create("10gev", "D")

    # 20-21 GeV
    dvar_20gev_A = dvar_create("20gev", "A")
    dvar_20gev_B = dvar_create("20gev", "B")
    dvar_20gev_C = dvar_create("20gev", "C")
    dvar_20gev_D = dvar_create("20gev", "D")


    ######################################################################
    # Loop through TTree and fill in dictionaries.
    
    k = -1
    for e in T:
        k += 1
        #if k == 20000: break    #20000
        if k % 100000 == 0:
            print("k = ", k)

        n_muons = len(T.muon_pt)
        if n_muons < 4: continue
        #if n_muons < 6: continue

        # List all muon four-vectors and charges in the event.
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

        # List indices for all possible dimuons. ("pairs")
        lpairs = []
        for i in range(0, len(afourv)):
            for j in range(i+1, len(afourv)):
                lpairs.append([i,j])        

        """ print("lpairs")
        print(lpairs) """

        # List indices for all possible pairs of dimuons. ("pairs of pairs")
        lpairs_pairs = []
        for p in range(0, len(lpairs)):
            for q in range(p, len(lpairs)):
                # If no two elements in the pairs p,q are the same,
                # append to list.
                if lpairs[q][0] != lpairs[p][0] \
                    and lpairs[q][0] != lpairs[p][1] \
                    and lpairs[q][1] != lpairs[p][0] \
                        and lpairs[q][1] != lpairs[p][1]:
                    lpairs_pairs.append([p,q]) 

        """ print("lpairs_pairs")
        print(lpairs_pairs) """

        # List the four-vectors of dimuon pairs according to their indices.
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

        # Find dimuon pair with the lowest mass asymmetry.
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


        # Get properties of the lowest m_asymm dimuon pair.

        vdimuon1 = avectors[i_min, 0]       # Dimuon four-vectors.
        vdimuon2 = avectors[i_min, 1]

        m_asymm = am_asymm[i_min]           # Mass asymmetry.
        if m_asymm > 0.05 \
        and m_asymm < 0.075: continue       # Exclude events in buffer
                                            # zone between regions A, C.
        if m_asymm > 0.5: continue          # Exclude events with 
                                            # large m_asymm.

        q1 = aq[p11] + aq[p12]              # Charge of dimuon 1.
        q2 = aq[p21] + aq[p22]              # Charge of dimuon 2.
        if abs(q1 + q2) == 4: continue      # Exclude events with 
                                            # 4 muons of the same charge.
        q_four = abs(q1) + abs(q2)          # Total charge of four muons.
                                            # q_four >= 0.
        
        """ print()
        print("q1 = ", q1)
        print(aq[p11], aq[p12])
        print("q2 = ", q2)
        print(aq[p21], aq[p22])
        print("q_four = ", q_four) """

        m_dimuon1 = vdimuon1.M()            # Individual dimuon masses.
        m_dimuon2 = vdimuon2.M()

        mavg = (m_dimuon1 + m_dimuon2)/2    # Average mass of dimuon pair.
        """ if mavg < 2.5 \
        or mavg > 4.0: continue             # Exclude events outside J/Psi
                                            # region and control regions. """

        V_fourmuon = vdimuon1 + vdimuon2
        mfour = V_fourmuon.M()              # Four-muon mass.

        alpha = float(mavg) / float(mfour)  # Alpha ratio
                                            # (dimuon mass : four-muon mass).

        dR = vdimuon1.DeltaR(vdimuon2)      # deltaR of dimuon pair.

        eta1 = vdimuon1.Eta()
        eta2 = vdimuon2.Eta()
        deta = abs(eta1 - eta2)             # delta_eta of dimuon pair.

        apt_four = [apt[p11], apt[p12], 
                    apt[p21], apt[p22]]
        pt_lead = max(apt_four)             # pT of leading muon.
        pt_last = min(apt_four)             # pT of last muon.

        pt_dimuon = [vdimuon1.Pt(), 
                     vdimuon2.Pt()]
        pt_dimlead = max(pt_dimuon)         # pT of leading dimuon.
        pt_dimlast = min(pt_dimuon)         # pT of last dimuon.



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

        # Sort event into ABCD regions.
        # Temporary dictionary for holding the values of each event.
        dvalues = {
            "mavg": mavg, 
            "mfour": mfour, 
            "alpha": alpha, 
            "dR": dR, 
            "deta": deta, 
            "pt_lead": pt_lead,
            "pt_last": pt_last,
            "pt_dimlead": pt_dimlead,
            "pt_dimlast": pt_dimlast}


        region = abcd_sort(m_asymm, q_four)

        """ if region == "A":
            if mavg >= 2.5 and mavg < 3.0:        # Control 1 mass region.
                for key, value in dvalues.items():
                    dvar_all_A["lall_"+key+"_A"].append(value)
                    dvar_control1_A["lcontrol1_"+key+"_A"].append(value)

            elif mavg >= 3.0 and mavg <= 3.2:     # J/Psi mass region.
                for key, value in dvalues.items():
                    dvar_all_A["lall_"+key+"_A"].append(value)
                    dvar_JPsi_A["lJPsi_"+key+"_A"].append(value)
            
            elif mavg > 3.2 and mavg <= 4.0:      # Control 2 mass region.
                for key, value in dvalues.items():
                    dvar_all_A["lall_"+key+"_A"].append(value)
                    dvar_control2_A["lcontrol2_"+key+"_A"].append(value)

    
        elif region == "B":
            if mavg >= 2.5 and mavg < 3.0:        # Control 1 mass region.
                for key, value in dvalues.items():
                    dvar_all_B["lall_"+key+"_B"].append(value)
                    dvar_control1_B["lcontrol1_"+key+"_B"].append(value)

            elif mavg >= 3.0 and mavg <= 3.2:     # J/Psi mass region.
                for key, value in dvalues.items():
                    dvar_all_B["lall_"+key+"_B"].append(value)
                    dvar_JPsi_B["lJPsi_"+key+"_B"].append(value)
            
            elif mavg > 3.2 and mavg <= 4.0:      # Control 2 mass region.
                for key, value in dvalues.items():
                    dvar_all_B["lall_"+key+"_B"].append(value)
                    dvar_control2_B["lcontrol2_"+key+"_B"].append(value)

        elif region == "C":
            if mavg >= 2.5 and mavg < 3.0:        # Control 1 mass region.
                for key, value in dvalues.items():
                    dvar_all_C["lall_"+key+"_C"].append(value)
                    dvar_control1_C["lcontrol1_"+key+"_C"].append(value)

            elif mavg >= 3.0 and mavg <= 3.2:     # J/Psi mass region.
                for key, value in dvalues.items():
                    dvar_all_C["lall_"+key+"_C"].append(value)
                    dvar_JPsi_C["lJPsi_"+key+"_C"].append(value)
            
            elif mavg > 3.2 and mavg <= 4.0:      # Control 2 mass region.
                for key, value in dvalues.items():
                    dvar_all_C["lall_"+key+"_C"].append(value)
                    dvar_control2_C["lcontrol2_"+key+"_C"].append(value)

        elif region == "D":
            if mavg >= 2.5 and mavg < 3.0:        # Control 1 mass region.
                for key, value in dvalues.items():
                    dvar_all_D["lall_"+key+"_D"].append(value)
                    dvar_control1_D["lcontrol1_"+key+"_D"].append(value)

            elif mavg >= 3.0 and mavg <= 3.2:     # J/Psi mass region.
                for key, value in dvalues.items():
                    dvar_all_D["lall_"+key+"_D"].append(value)
                    dvar_JPsi_D["lJPsi_"+key+"_D"].append(value)
            
            elif mavg > 3.2 and mavg <= 4.0:      # Control 2 mass region.
                for key, value in dvalues.items():
                    dvar_all_D["lall_"+key+"_D"].append(value)
                    dvar_control2_D["lcontrol2_"+key+"_D"].append(value) """
        

        """ if region == "A":
            if mavg >= 2.0 and mavg < 3.0:        # 2-3 GeV.
                for key, value in dvalues.items():
                    dvar_2gev_A["l2gev_"+key+"_A"].append(value)

            elif mavg >= 4.0 and mavg < 5.0:     # 4-5 GeV.
                for key, value in dvalues.items():
                    dvar_4gev_A["l4gev_"+key+"_A"].append(value)
            
            elif mavg >= 5.0 and mavg < 6.0:      # 5-6 GeV.
                for key, value in dvalues.items():
                    dvar_5gev_A["l5gev_"+key+"_A"].append(value)

    
        elif region == "B":
            if mavg >= 2.0 and mavg < 3.0:        # 2-3 GeV.
                for key, value in dvalues.items():
                    dvar_2gev_B["l2gev_"+key+"_B"].append(value)

            elif mavg >= 4.0 and mavg < 5.0:     # 4-5 GeV.
                for key, value in dvalues.items():
                    dvar_4gev_B["l4gev_"+key+"_B"].append(value)
            
            elif mavg >= 5.0 and mavg < 6.0:      # 5-6 GeV.
                for key, value in dvalues.items():
                    dvar_5gev_B["l5gev_"+key+"_B"].append(value)

        elif region == "C":
            if mavg >= 2.0 and mavg < 3.0:        # 2-3 GeV.
                for key, value in dvalues.items():
                    dvar_2gev_C["l2gev_"+key+"_C"].append(value)

            elif mavg >= 4.0 and mavg < 5.0:     # 4-5 GeV.
                for key, value in dvalues.items():
                    dvar_4gev_C["l4gev_"+key+"_C"].append(value)
            
            elif mavg >= 5.0 and mavg < 6.0:      # 5-6 GeV.
                for key, value in dvalues.items():
                    dvar_5gev_C["l5gev_"+key+"_C"].append(value)

        elif region == "D":
            if mavg >= 2.0 and mavg < 3.0:        # 2-3 GeV.
                for key, value in dvalues.items():
                    dvar_2gev_D["l2gev_"+key+"_D"].append(value)

            elif mavg >= 4.0 and mavg < 5.0:     # 4-5 GeV.
                for key, value in dvalues.items():
                    dvar_4gev_D["l4gev_"+key+"_D"].append(value)
            
            elif mavg >= 5.0 and mavg < 6.0:      # 5-6 GeV.
                for key, value in dvalues.items():
                    dvar_5gev_D["l5gev_"+key+"_D"].append(value) """

        if region == "A":
            if mavg >= 8.0 and mavg < 9.0:         # 8-9 GeV.
                for key, value in dvalues.items():
                    dvar_8gev_A["l8gev_"+key+"_A"].append(value)
            elif mavg >= 10.0 and mavg < 11.0:     # 10-11 GeV.
                for key, value in dvalues.items():
                    dvar_10gev_A["l10gev_"+key+"_A"].append(value)

            elif mavg >= 20.0 and mavg < 21.0:     # 20-21 GeV.
                for key, value in dvalues.items():
                    dvar_20gev_A["l20gev_"+key+"_A"].append(value)

    
        elif region == "B":
            if mavg >= 8.0 and mavg < 9.0:         # 8-9 GeV.
                for key, value in dvalues.items():
                    dvar_8gev_B["l8gev_"+key+"_B"].append(value)
            elif mavg >= 10.0 and mavg < 11.0:     # 10-11 GeV.
                for key, value in dvalues.items():
                    dvar_10gev_B["l10gev_"+key+"_B"].append(value)

            elif mavg >= 20.0 and mavg < 21.0:     # 20-21 GeV.
                for key, value in dvalues.items():
                    dvar_20gev_B["l20gev_"+key+"_B"].append(value)

        elif region == "C":
            if mavg >= 8.0 and mavg < 9.0:         # 8-9 GeV.
                for key, value in dvalues.items():
                    dvar_8gev_C["l8gev_"+key+"_C"].append(value)
            elif mavg >= 10.0 and mavg < 11.0:     # 10-11 GeV.
                for key, value in dvalues.items():
                    dvar_10gev_C["l10gev_"+key+"_C"].append(value)

            elif mavg >= 20.0 and mavg < 21.0:     # 20-21 GeV.
                for key, value in dvalues.items():
                    dvar_20gev_C["l20gev_"+key+"_C"].append(value)

        elif region == "D":
            if mavg >= 8.0 and mavg < 9.0:         # 8-9 GeV.
                for key, value in dvalues.items():
                    dvar_8gev_D["l8gev_"+key+"_D"].append(value)
            elif mavg >= 10.0 and mavg < 11.0:     # 10-11 GeV.
                for key, value in dvalues.items():
                    dvar_10gev_D["l10gev_"+key+"_D"].append(value)

            elif mavg >= 20.0 and mavg < 21.0:     # 20-21 GeV.
                for key, value in dvalues.items():
                    dvar_20gev_D["l20gev_"+key+"_D"].append(value)
    


    #####################################################################    
    
    # Save dictionaries to json files.
    """ ddvar = {
        "all": [dvar_all_A, dvar_all_B, dvar_all_C, dvar_all_D],
        "control1": [dvar_control1_A, dvar_control1_B, 
                     dvar_control1_C, dvar_control1_D],
        "JPsi": [dvar_JPsi_A, dvar_JPsi_B, dvar_JPsi_C, dvar_JPsi_D],
        "control2": [dvar_control2_A, dvar_control2_B, 
                     dvar_control2_C, dvar_control2_D]
    } """

    """ ddvar = {
        "2gev": [dvar_2gev_A, dvar_2gev_B, 
                     dvar_2gev_C, dvar_2gev_D],
        "4gev": [dvar_4gev_A, dvar_4gev_B, 
                     dvar_4gev_C, dvar_4gev_D],
        "5gev": [dvar_5gev_A, dvar_5gev_B, 
                     dvar_5gev_C, dvar_5gev_D],
    } """

    ddvar = {
        "8gev": [dvar_8gev_A, dvar_8gev_B, 
                     dvar_8gev_C, dvar_8gev_D],
        "10gev": [dvar_10gev_A, dvar_10gev_B, 
                     dvar_10gev_C, dvar_10gev_D],
        "20gev": [dvar_20gev_A, dvar_20gev_B, 
                     dvar_20gev_C, dvar_20gev_D],
    }


    for mass_region in ddvar.keys():
        i = 0
        for abcd_region in ["A", "B", "C", "D"]:
            filename = "dvar_new_{}_{}.json".\
                        format(mass_region, abcd_region)
            with open(filename, "w") as file:
                json.dump(ddvar[mass_region][i],
                          file, indent=4)
            
            events = len((ddvar[mass_region][i]["l{}_mavg_{}".\
                        format(mass_region, abcd_region)]))
            print("Number of events in {} = {}".\
                  format(filename, events))
            i += 1

        print("Done dvar_new_{}_A/B/C/D.json".format(mass_region))
