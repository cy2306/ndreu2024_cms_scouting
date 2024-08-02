"""
Sort events in a root file into ABCD regions (according to ABCD method
of background estimation for paired dimuon events) 
and save the momentum vector (p-vector) components of

    1. plead_mu1: muon 1 of the leading-pT dimuon,
    2. psec_mu1: muon 2 of the leading-pT dimuon,
    3. psec_mu1: muon 1 of the second dimuon,
    4. psec_mu2: muon 2 of the second dimuon,
    5. plead_dimu: leading-pT dimuon,
    6. psec_dimu: second dimuon,
    7. pfourmu: four-muon particle,

given in the four-muon rest frame (except pfourmu is given in lab frame).


For each event in the TTree:
1. Identify the dimuon pair of lowest mass asymmetry.

2. Select events based on restrictions in number of muons, 
    mass asymmetry, dimuon charge, dimuon mass, etc. 
    Events are included if all of the following are met.
    
    a. n_muons >= 4
    b. m_asymm <= 0.05 or 0.075 <= m_asymm <= 0.5
    c. abs(q1 + q2 + q3 + q4) < 4 
        (all four muons do not have the same charge


3. Get the p-vector components of all particles (listed above).
    a. Create four-vectors (pT, eta, phi, m) for all particles.
    b. Transform all four-vectors to the four-muon rest frame.
    c. Get the momentum components (px, py, pz) of each particle.


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
5. Create a json file for each ABCD region and mass region. 
    Each file is a dictionary of lists containing relevant variables, 
    formatted as follows.
    {
        "lplead_mu1": [ [p1x, p1y, p1z],    [p2x, p2y, p2z], ...],
        "lpsec_mu1": [ [p1x, p1y, p1z],    [p2x, p2y, p2z], ...], 
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
from numpy import array, empty, argmax, argmin
from math import cos, sin, sinh, acos, sqrt
import json


def dpvectors_create(mass_region, abcd_region, lvariables):

    dpvectors = {}
    for variable in lvariables:
        key = "l{}_{}_{}".\
            format(mass_region, variable, abcd_region)
        dpvectors[key] = []

    return dpvectors


def lpx_py_pz(fourv):

    # Get the momentum three-vector (px, py, pz).
    px = fourv.Px()
    py = fourv.Py()
    pz = fourv.Pz()
    
    return [px, py, pz]


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


    # Create empty dictionaries for ABCD regions of all four-vectors.
    # List of four-vectors to be saved.
    lvariables = [
        "plead_mu1", "psec_mu1", "plead_mu2", "psec_mu2",
        "plead_dimu", "psec_dimu", "pfourmu"
    ]
    
    # All in J/Psi, control 1, and control 2.
    dpvectors_all_ABCD = dpvectors_create("all", "ABCD", lvariables)
    dpvectors_all_A = dpvectors_create("all", "A", lvariables)
    dpvectors_all_B = dpvectors_create("all", "B", lvariables)
    dpvectors_all_C = dpvectors_create("all", "C", lvariables)
    dpvectors_all_D = dpvectors_create("all", "D", lvariables)

    # Control 1 mass region (2.5-3.05 GeV)
    dpvectors_control1_ABCD = dpvectors_create("control1", "ABCD", lvariables)
    dpvectors_control1_A = dpvectors_create("control1", "A", lvariables)
    dpvectors_control1_B = dpvectors_create("control1", "B", lvariables)
    dpvectors_control1_C = dpvectors_create("control1", "C", lvariables)
    dpvectors_control1_D = dpvectors_create("control1", "D", lvariables)

    # J/Psi mass region (3.05-3.15 GeV)
    dpvectors_JPsi_ABCD = dpvectors_create("JPsi", "ABCD", lvariables)
    dpvectors_JPsi_A = dpvectors_create("JPsi", "A", lvariables)
    dpvectors_JPsi_B = dpvectors_create("JPsi", "B", lvariables)
    dpvectors_JPsi_C = dpvectors_create("JPsi", "C", lvariables)
    dpvectors_JPsi_D = dpvectors_create("JPsi", "D", lvariables)

    # Control 2 mass region (3.15-4.0 GeV)
    dpvectors_control2_ABCD = dpvectors_create("control2", "ABCD", lvariables)
    dpvectors_control2_A = dpvectors_create("control2", "A", lvariables)
    dpvectors_control2_B = dpvectors_create("control2", "B", lvariables)
    dpvectors_control2_C = dpvectors_create("control2", "C", lvariables)
    dpvectors_control2_D = dpvectors_create("control2", "D", lvariables)
    
    
    ######################################################################
    # Loop through TTree and fill in dictionaries.
    
    k = -1
    for e in T:
        k += 1
        #if k == 20000: break    # For testing code.
        if k % 500000 == 0:
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
 

        # For each possible pair of dimuons, list the four-vectors of 
        # the muons and dimuons according to their indices.
        avectors = empty((len(lpairs_pairs),6), object)
        i = 0
        for pp in lpairs_pairs:
            p11 = lpairs[pp[0]][0]          # index of pair 1 muon 1
            p12 = lpairs[pp[0]][1]          # index of pair 1 muon 2
            p21 = lpairs[pp[1]][0]          # index of pair 2 muon 1
            p22 = lpairs[pp[1]][1]          # index of pair 2 muon 2

            vmu11 = afourv[p11]             # pair 1 muon 1 four-vector
            vmu12 = afourv[p12]             # pair 1 muon 2 four-vector
            vmu21 = afourv[p21]             # pair 2 muon 1 four-vector
            vmu22 = afourv[p22]             # pair 2 muon 2 four-vector
            avectors[i, 0] = vmu11
            avectors[i, 1] = vmu12
            avectors[i, 2] = vmu21
            avectors[i, 3] = vmu22

            vdimu1 = vmu11 + vmu12          # dimuon four-vector 1
            vdimu2 = vmu21 + vmu22          # dimuon four-vector 2
            avectors[i, 4] = vdimu1
            avectors[i, 5] = vdimu2

            i += 1
        
        # Find dimuon pair with the lowest mass asymmetry.
        # Indices of avectors and am_asymm should match.
        am_asymm = empty(len(avectors), float)
        i = 0
        for v_dimuon in avectors:
            m1 = v_dimuon[4].M()
            m2 = v_dimuon[5].M()
            am_asymm[i] = abs((m1-m2)/(m1+m2))
            
            i += 1

        i_min = argmin(am_asymm)   # index of lowest m_asymm.


        # Get properties of the lowest m_asymm dimuon pair.
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
        
        m_dimu1 = avectors[i_min, 4].M()    # Individual dimuon masses.
        m_dimu2 = avectors[i_min, 5].M()

        mavg = (m_dimu1 + m_dimu2)/2        # Average mass of dimuon pair.
        if mavg < 2.5 \
        or mavg > 4.0: continue             # Exclude events outside J/Psi
                                            # region and control regions.


        # Get four-vectors of the lowest m_asymm dimuon pair.
        vlowest_mu11 = avectors[i_min, 0]        # muon four-vectors
        vlowest_mu12 = avectors[i_min, 1]
        vlowest_mu21 = avectors[i_min, 2]
        vlowest_mu22 = avectors[i_min, 3]

        vlowest_dimu1 = avectors[i_min, 4]       # dimuon four-vectors
        vlowest_dimu2 = avectors[i_min, 5]

        alowest_fourv = array([
            vlowest_mu11, vlowest_mu12, vlowest_mu21, vlowest_mu22,
            vlowest_dimu1, vlowest_dimu2
        ])

        # Get pT of the muons and dimuons.
        alowest_pt = empty(6, float)
        i = 0
        for fourv in alowest_fourv:
            pt = fourv.Pt()
            alowest_pt[i] = pt
            i += 1
        
        """ print("pt of muon1, muon2, muon3, muon4, dimuon1, dimuon2")
        print(alowest_pt) """
        

        # Order muon and dimuon four-vectors by pT.
        muons1 = [alowest_pt[0], alowest_pt[1]]
        muons2 = [alowest_pt[2], alowest_pt[3]]
        dimuons = [alowest_pt[4], alowest_pt[5]]
        
        """ print("muon pair 1: {}".format(muons1))
        print("muon pair 2: {}".format(muons2))
        print("dimuon pair: {}".format(dimuons)) """

        lead_mu1 = argmax(muons1)
        sec_mu1 = argmin(muons1)
        lead_mu2 = argmax(muons2) + 2
        sec_mu2 = argmin(muons2) + 2
        lead_dimu = argmax(dimuons) + 4
        sec_dimu = argmin(dimuons) + 4

        """ print("lead muon1 pt: ({}) {}".\
              format(lead_mu1, alowest_pt[lead_mu1]))
        print("second muon1 pt: ({}) {}".\
              format(sec_mu1, alowest_pt[sec_mu1]))
        print("lead muon2 pt: ({}) {}".\
              format(lead_mu2, alowest_pt[lead_mu2]))
        print("second muon2 pt: ({}) {}".\
              format(sec_mu2, alowest_pt[sec_mu2]))
        print("lead dimuon pt: ({}) {}".\
              format(lead_dimu, alowest_pt[lead_dimu]))
        print("second dimuon pt: ({}) {}".\
              format(sec_dimu, alowest_pt[sec_dimu])) """
        
        # Get the four-vectors in the four-muon rest frame.
        vfourmu = alowest_fourv[lead_dimu] \
                    + alowest_fourv[sec_dimu]
        vlead_mu1 = alowest_fourv[lead_mu1] - vfourmu
        vsec_mu1 = alowest_fourv[sec_mu1] - vfourmu
        vlead_mu2 = alowest_fourv[lead_mu2] - vfourmu
        vsec_mu2 = alowest_fourv[sec_mu2] - vfourmu
        vlead_dimu = alowest_fourv[lead_dimu] - vfourmu
        vsec_dimu = alowest_fourv[sec_dimu] - vfourmu
        
        # Get the momentum three-vector components [px, py, pz].
        # Pair 1 leading and second muon.
        plead_mu1 = lpx_py_pz(vlead_mu1)
        psec_mu1 = lpx_py_pz(vsec_mu1)

        # Pair 2 leading and second muon.
        plead_mu2 = lpx_py_pz(vlead_mu2)
        psec_mu2 = lpx_py_pz(vsec_mu2)

        # Leading and second dimuon.
        plead_dimu = lpx_py_pz(vlead_dimu)
        psec_dimu = lpx_py_pz(vsec_dimu)

        # Four-muon.
        pfourmu = lpx_py_pz(vfourmu)


        # Sort event into ABCD regions.
        # Temporary dictionary for holding the values of each event.
        dvalues = {
            "plead_mu1": plead_mu1,
            "psec_mu1": psec_mu1,
            "plead_mu2": plead_mu2,
            "psec_mu2": psec_mu2,
            "plead_dimu": plead_dimu,
            "psec_dimu": psec_dimu,
            "pfourmu": pfourmu
            }
        

        # All events.
        if mavg >= 2.5 and mavg < 3.05:        # Control 1 mass region.
            for key, value in dvalues.items():
                dpvectors_all_ABCD["lall_"+key+"_ABCD"].append(value)
                dpvectors_control1_ABCD["lcontrol1_"+key+"_ABCD"].append(value)

        elif mavg >= 3.05 and mavg < 3.15:     # J/Psi mass region.
            for key, value in dvalues.items():
                dpvectors_all_ABCD["lall_"+key+"_ABCD"].append(value)
                dpvectors_JPsi_ABCD["lJPsi_"+key+"_ABCD"].append(value)
        
        elif mavg >= 3.15 and mavg <= 4.0:      # Control 2 mass region.
            for key, value in dvalues.items():
                dpvectors_all_ABCD["lall_"+key+"_ABCD"].append(value)
                dpvectors_control2_ABCD["lcontrol2_"+key+"_ABCD"].append(value)

        # Sort into ABCD regions.
        region = abcd_sort(m_asymm, q_four)

        if region == "A":
            if mavg >= 2.5 and mavg < 3.05:        # Control 1 mass region.
                for key, value in dvalues.items():
                    dpvectors_all_A["lall_"+key+"_A"].append(value)
                    dpvectors_control1_A["lcontrol1_"+key+"_A"].append(value)

            elif mavg >= 3.05 and mavg < 3.15:     # J/Psi mass region.
                for key, value in dvalues.items():
                    dpvectors_all_A["lall_"+key+"_A"].append(value)
                    dpvectors_JPsi_A["lJPsi_"+key+"_A"].append(value)
            
            elif mavg >= 3.15 and mavg <= 4.0:      # Control 2 mass region.
                for key, value in dvalues.items():
                    dpvectors_all_A["lall_"+key+"_A"].append(value)
                    dpvectors_control2_A["lcontrol2_"+key+"_A"].append(value)

    
        elif region == "B":
            if mavg >= 2.5 and mavg < 3.05:        # Control 1 mass region.
                for key, value in dvalues.items():
                    dpvectors_all_B["lall_"+key+"_B"].append(value)
                    dpvectors_control1_B["lcontrol1_"+key+"_B"].append(value)

            elif mavg >= 3.05 and mavg < 3.15:     # J/Psi mass region.
                for key, value in dvalues.items():
                    dpvectors_all_B["lall_"+key+"_B"].append(value)
                    dpvectors_JPsi_B["lJPsi_"+key+"_B"].append(value)
            
            elif mavg >= 3.15 and mavg <= 4.0:      # Control 2 mass region.
                for key, value in dvalues.items():
                    dpvectors_all_B["lall_"+key+"_B"].append(value)
                    dpvectors_control2_B["lcontrol2_"+key+"_B"].append(value)

        elif region == "C":
            if mavg >= 2.5 and mavg < 3.05:        # Control 1 mass region.
                for key, value in dvalues.items():
                    dpvectors_all_C["lall_"+key+"_C"].append(value)
                    dpvectors_control1_C["lcontrol1_"+key+"_C"].append(value)

            elif mavg >= 3.05 and mavg < 3.15:     # J/Psi mass region.
                for key, value in dvalues.items():
                    dpvectors_all_C["lall_"+key+"_C"].append(value)
                    dpvectors_JPsi_C["lJPsi_"+key+"_C"].append(value)
            
            elif mavg >= 3.15 and mavg <= 4.0:      # Control 2 mass region.
                for key, value in dvalues.items():
                    dpvectors_all_C["lall_"+key+"_C"].append(value)
                    dpvectors_control2_C["lcontrol2_"+key+"_C"].append(value)

        elif region == "D":
            if mavg >= 2.5 and mavg < 3.05:        # Control 1 mass region.
                for key, value in dvalues.items():
                    dpvectors_all_D["lall_"+key+"_D"].append(value)
                    dpvectors_control1_D["lcontrol1_"+key+"_D"].append(value)

            elif mavg >= 3.05 and mavg < 3.15:     # J/Psi mass region.
                for key, value in dvalues.items():
                    dpvectors_all_D["lall_"+key+"_D"].append(value)
                    dpvectors_JPsi_D["lJPsi_"+key+"_D"].append(value)
            
            elif mavg >= 3.15 and mavg <= 4.0:      # Control 2 mass region.
                for key, value in dvalues.items():
                    dpvectors_all_D["lall_"+key+"_D"].append(value)
                    dpvectors_control2_D["lcontrol2_"+key+"_D"].append(value)



    #####################################################################    
    
    # Save dictionaries to json files.
    ddpvectors = {
        "all": [dpvectors_all_A, dpvectors_all_B,
                dpvectors_all_C, dpvectors_all_D, 
                dpvectors_all_ABCD],
        "control1": [dpvectors_control1_A, dpvectors_control1_B,
                     dpvectors_control1_C, dpvectors_control1_D,
                     dpvectors_control1_ABCD],
        "JPsi": [dpvectors_JPsi_A, dpvectors_JPsi_B,
                     dpvectors_JPsi_C, dpvectors_JPsi_D,
                     dpvectors_JPsi_ABCD],
        "control2": [dpvectors_control2_A, dpvectors_control2_B,
                     dpvectors_control2_C, dpvectors_control2_D,
                     dpvectors_control2_ABCD]
    }

    for mass_region in ddpvectors.keys():
        i = 0
        for abcd_region in ["A", "B", "C", "D", "ABCD"]:
            filename = "dpvectors_doubleJPsi_new_{}_{}.json".\
                        format(mass_region, abcd_region)
            with open(filename, "w") as file:
                json.dump(ddpvectors[mass_region][i],
                            file, indent=4)
            
            events = len((ddpvectors[mass_region][i]["l{}_pfourmu_{}".\
                        format(mass_region, abcd_region)]))
            print("Number of events in {} = {}".\
                    format(filename, events))
            i += 1

        print("Done dpvectors_doubleJPsi_new_{}_A/B/C/D.json".format(mass_region))
