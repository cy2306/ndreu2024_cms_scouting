"""
Sort events in a root file into ABCD regions for background estimates 
using the ABCD method for paired dimuon events.
Save the mass asymmetry and four-muon charge of each event.

    m_asymm = abs((m1-m2)/(m1+m2)),
    q_four = q_four = abs(q1) + abs(q2),

where the indices 1,2 indicate dimuons in a pair 
(not in any particular order).

    
For each event in the TTree:
1. Identify the dimuon pair of lowest mass asymmetry.

2. Select events based on restrictions in number of muons, 
    mass asymmetry, dimuon charge, dimuon mass, etc. 
    Events are included if all of the following are met.
    
    a. n_muons >= 4
    b. m_asymm <= 0.05 or 0.075 <= m_asymm <= 0.5
    c. abs(q1 + q2 + q3 + q4) < 4 
        (all four muons do not have the same charge) 


3. Get the m_asymm and q_four.        


4. Sort selected events into ABCD regions based on four-muon charge
    and mass asymmetry of dimuon pair.
    
    A: m_asymm <= 0.05 and q_four == 0
    B: m_asymm <= 0.05 and q_four != 0
    C: 0.075 <= m_asymm <= 0.5 and q_four == 0
    D: 0.075 <= m_asymm <= 0.5 and q_four != 0


5. Further sort events in mass regions based on average dimuon mass,
    and append relevant variables to the correct list to be saved.

    "control1": 2.5 <= mavg < 3.05
    "JPsi": 3.05 <= mavg < 3.15
    "control2": 3.15 <= mavg <= 4.0
    ...

    
After looping through all events: 
6. Create a json file for each ABCD region and mass region. 
    Each file is a dictionary of lists containing relevant variables, 
    formatted as follows.
    {
        "m_asymm": [event 1 m_asymm, event 2 m_asymm, ...],
        "q_four": [event 1 q_four, event 2 q_four, ...],
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
    "l"+mass_region+"_masym_"+abcd_region: [],
    "l"+mass_region+"_qfour_"+abcd_region: [],
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
    # All events.
    dvar_all_A = dvar_create("all", "A")
    dvar_all_B = dvar_create("all", "B")
    dvar_all_C = dvar_create("all", "C")
    dvar_all_D = dvar_create("all", "D")

    # J/Psi mass region (3.0-3.2 GeV)
    dvar_JPsi_A = dvar_create("JPsi", "A")
    dvar_JPsi_B = dvar_create("JPsi", "B")
    dvar_JPsi_C = dvar_create("JPsi", "C")
    dvar_JPsi_D = dvar_create("JPsi", "D")


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

        # List indices for all possible dimuons. ("pairs")
        lpairs = []
        for i in range(0, len(afourv)):
            for j in range(i+1, len(afourv)):
                lpairs.append([i,j])  
        
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

        # List the four-vectors of dimuon pairs according to their indices.
        avectors = empty((len(lpairs_pairs),2), object)
        i = 0
        for pp in lpairs_pairs:
            p11 = lpairs[pp[0]][0]          # index of pair 1 muon 1
            p12 = lpairs[pp[0]][1]          # index of pair 1 muon 2

            p21 = lpairs[pp[1]][0]
            p22 = lpairs[pp[1]][1]

            v1 = afourv[p11] + afourv[p12]  # dimuon four-vector 1
            v2 = afourv[p21] + afourv[p22]  # dimuon four-vector 2
            avectors[i, 0] = v1
            avectors[i, 1] = v2

            i += 1
        
        # Find dimuon pair with the lowest mass asymmetry.
        # Indices of avectors and am_asymm should match.
        am_asymm = empty(len(avectors), float)
        i = 0
        for v_dimuon in avectors:
            m1 = v_dimuon[0].M()
            m2 = v_dimuon[1].M()
            am_asymm[i] = abs((m1-m2)/(m1+m2))
            
            i += 1
        
        i_min = argmin(am_asymm)   # index of lowest m_asymm.

        
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
        
        m_dimuon1 = vdimuon1.M()            # Individual dimuon masses.
        m_dimuon2 = vdimuon2.M()

        mavg = (m_dimuon1 + m_dimuon2)/2    # Average mass of dimuon pair.
        

        # Sort event into ABCD regions.
        # Temporary dictionary for holding the values of each event.
        dvalues = {
            "masym": m_asymm, 
            "qfour": q_four
        }

        region = abcd_sort(m_asymm, q_four)

        if region == "A":
            if mavg >= 3.05 and mavg <= 3.15:           # J/Psi mass region.
                for key, value in dvalues.items():
                    dvar_JPsi_A["lJPsi_"+key+"_A"].append(value)
            for key, value in dvalues.items():          # all events.
                    dvar_all_A["lall_"+key+"_A"].append(value)
        
        elif region == "B":
            if mavg >= 3.05 and mavg <= 3.15:           # J/Psi mass region.
                for key, value in dvalues.items():
                    dvar_JPsi_B["lJPsi_"+key+"_B"].append(value)
            for key, value in dvalues.items():          # all events
                    dvar_all_B["lall_"+key+"_B"].append(value)

        elif region == "C":
            if mavg >= 3.05 and mavg <= 3.15:           # J/Psi mass region.
                for key, value in dvalues.items():
                    dvar_JPsi_C["lJPsi_"+key+"_C"].append(value)
            for key, value in dvalues.items():          # all events
                    dvar_all_C["lall_"+key+"_C"].append(value)

        elif region == "D":
            if mavg >= 3.05 and mavg <= 3.15:           # J/Psi mass region.
                for key, value in dvalues.items():
                    dvar_JPsi_D["lJPsi_"+key+"_D"].append(value)
            for key, value in dvalues.items():          # all events
                    dvar_all_D["lall_"+key+"_D"].append(value)

    
    #####################################################################    
    
    # Save dictionaries to json files.
    ddvar = {
        "all": [dvar_all_A, dvar_all_B, dvar_all_C, dvar_all_D],
        "JPsi": [dvar_JPsi_A, dvar_JPsi_B, dvar_JPsi_C, dvar_JPsi_D]
    }

    for mass_region in ddvar.keys():
        i = 0
        for abcd_region in ["A", "B", "C", "D"]:
            filename = "dmasym_qfour_new_{}_{}.json".\
                        format(mass_region, abcd_region)
            with open(filename, "w") as file:
                json.dump(ddvar[mass_region][i],
                          file, indent=4)
            
            events = len((ddvar[mass_region][i]["l{}_masym_{}".\
                        format(mass_region, abcd_region)]))
            print("Number of events in {} = {}".\
                  format(filename, events))
            i += 1

        print("Done dmasym_qfour_new_{}_A/B/C/D.json".format(mass_region))