"""
Plot 1d histogram of the CPV-sensitive angle deltaPhi for ABCD regions 
of paired dimuon events.

Variables needed: 
Components of momentum vectors (px, py, pz) of

    1. plead_mu1: muon 1 of the leading-pT dimuon,
    2. psec_mu1: muon 2 of the leading-pT dimuon,
    3. psec_mu1: muon 1 of the second dimuon,
    4. psec_mu2: muon 2 of the second dimuon,
    5. plead_dimu: leading-pT dimuon,

given in the four-muon rest frame.


Calculations:
    x : cross product
    * : dot product
    abs(): norm of a vector

    1. Unit vectors normal to dimuon planes:
        n1 = (plead_mu1 x psec_mu1) / abs(plead_mu1 x psec_mu1),
        n2 = (plead_mu2 x psec_mu2) / abs(plead_mu2 x psec_mu2)

    2. a = (plead_dimu * (n1 x n2)) / abs(plead_dimu * (n1 x n2))

    3. deltaPhi = a(arccos(n1 * n2))


For each ABCD region:
1. Read json file of the specified mass region and ABCD region.
    json files can be made using pvectors_abcd_sort.py.
    Files should be a dictionaries that contain the p-vector components
    of all particles needed (listed above), in the form
    {
        "lplead_mu1": [ [p1x, p1y, p1z],    [p2x, p2y, p2z], ...],
        "lpsec_mu1": [ [p1x, p1y, p1z],    [p2x, p2y, p2z], ...], 
        ...
    }
    
    Examples of mass regions: 
    "JPsi": events with dimuon masses of 3.05-3.15 GeV.
    "control1": events with dimuon masses of 2.5-3.05 GeV.
    "control2": events with dimuon masses of 3.15-4.0 GeV.

2. Convert each list of vector components into an array of 
    TVector 3 p-vector objects.

3. Use the arrays of p-vectors to calculate deltaPhi for each event.

4. Plot 1d histogram of deltaPhi distribution. Draw and save to pdf file.

"""
"""
############################### Functions ###############################
1. array_pvectors(lvectors, n)
    Create array of the momentum vectors (p-vectors) of one type of 
    particle (plead_mu1, plead_mu2, psec_mu1, psec_mu2, or plead_dimu)
    in all events.

    lvectors: list (e.g. lplead_mu1 or lplead_mu2)
    Contains p-vector components for one type of particle in all events.
    e.g. lplead_mu1 is a list of the form
    [ [p1x, p1y, p1z],  [p2x, p2y, p2z], ... ]

    n: int
    Number of events.

    Return: array of TVector3 objects (p-vector of each event).
    

2. array_interplane_angle(aplead_mu1, apsec_mu1, 
                           aplead_mu2, apsec_mu2, 
                           aplead_dimu, n)
    
    Calculate the angle deltaPhi for each event and save in an array.
    
    aplead_mu1, apsec_mu1, aplead_mu2, apsec_mu2, aplead_dimu: arrays
    p-vectors of all the particles needed in all events.
    Should be outputs of array_pvectors(lvectors, n).

    n: int
    Number of events.

    Return: array of floats (deltaPhi of each event).

    
3. th1f_angle_create_draw(array_angles, bins, mass_region, region_abcd)
    Draw 1d histogram of deltaPhi distribution.
    Events/width vs deltaPhi, plotted as points with error bars "p", "e0". 
    Print and save to pdf file.

    array_angle: array
    deltaPhi of each event. 
    Should be the output of array_interplane_angle(...).

    bins: list
    [number of bins, array of bin edges]

    mass_region: str
    Name for the mass region

    region_abcd: str ("A", "B", "C", "D", or "ABCD")
    Name for the ABCD region

    Return: None

    
4. aplanarity(mass_region, region_abcd)
    Call the above functions in order to calculate and plot deltaPhi.
    Serves as the target for each concurrent.futures process.

    Return: None


############################ Main function ############################
Specify the range of dimuon masses (mass_region) and list ABCD regions.

For each ABCD region, call aplanarity(mass_region, region_abcd) to 
create a parallel process to plot deltaPhi.

"""

from ROOT import *
from numpy import empty, arange
from math import acos
import json
import concurrent.futures


def array_pvectors(lvectors, n):

    # Get the momentum vectors (px, py, pz).
    apvectors = empty(n, object)
    i = 0
    for pvector_list in lvectors:
        px = pvector_list[0]
        py = pvector_list[1]
        pz = pvector_list[2]
        
        pvector = TVector3(px, py, pz)
        apvectors[i] = pvector
        
        """ print("pvector: {}".format([
            pvector.X(), pvector.Y(), pvector.Z()])) """

        i += 1
        if i == n: break
    
    return apvectors


def array_interplane_angle(aplead_mu1, apsec_mu1, 
                           aplead_mu2, apsec_mu2, 
                           aplead_dimu, n):

    array_angles = empty(n, float)
    for i in range(n):

        # Take the cross product of each muon pair to get the normal
        # vector to the dimuon plane.
        crossv1 = aplead_mu1_[i].Cross(apsec_mu1[i])
        crossv2 = aplead_mu2[i].Cross(apsec_mu2[i])

        normv1 = crossv1.Unit()
        normv2 = crossv2.Unit()

        """ print("lead muon 1: {}".format([
            aplead_mu1[i].X(), aplead_mu1[i].Y(), 
            aplead_mu1[i].Z()]))
        print("second muon 1: {}".format([
            apsec_mu1[i].X(), apsec_mu1[i].Y(),
            apsec_mu1[i].Z()]))
        print("cross product 1: {}".format([
            crossv1.X(), crossv1.Y(), crossv1.Z()]))
        print("normal 1: {}".format([
            normv1.X(), normv1.Y(), normv1.Z()]))
        print("lead muon 2: {}".format([
            aplead_mu2[i].X(), aplead_mu2[i].Y(),
            aplead_mu2[i].Z()]))
        print("second muon 2: {}".format([
            apsec_mu2[i].X(), apsec_mu2[i].Y(),
            apsec_mu2[i].Z()]))
        print("cross product 2: {}".format([
            crossv2.X(), crossv2.Y(), crossv2.Z()]))
        print("normal 2: {}".format([
            normv2.X(), normv2.Y(), normv2.Z()])) """

        # Calculate the directional term a.
        adot = float(aplead_dimu[i].Dot(normv1.Cross(normv2)))
        a = adot / abs(adot)

        # Calculate the interplane angle.
        angle = a * acos(normv1.Dot(normv2))
        array_angles[i] = angle

        """ print("a: {}".format(a))
        print("angle: {}".format(angle)) """
    
    return array_angles


def th1f_angle_create_draw(array_angles, bins, mass_region, region_abcd):
    
    # Set up 1d histogram.
    nbins = bins[0]
    abins = bins[1]

    histname = "th1f_{}_angle_{}".format(mass_region, region_abcd)
    H1d = TH1F(histname, "", nbins, abins)
    
    # Fill in histogram.
    for angle in array_angles:
        H1d.Fill(angle)
    
    # Format names to be displayed on histogram.
    if mass_region == "all":
        #mass_region_formatted = "(2.5 < #hat{m}_{#mu#mu} < 4.0 GeV)"
        mass_region_formatted = ""
    elif mass_region == "control1":
        #mass_region_formatted = "(2.5 < #hat{m}_{#mu#mu} < 3.05 GeV)"
        mass_region_formatted = "(control 1)"
    elif mass_region == "JPsi":
        #mass_region_formatted = "(3.05 < #hat{m}_{#mu#mu} < 3.15 GeV)"
        mass_region_formatted = "(J/#psi)"
    elif mass_region == "control2":
        #mass_region_formatted = "(3.15 < #hat{m}_{#mu#mu} < 4.0 GeV)"
        mass_region_formatted = "(control 2)"

    if region_abcd == "ABCD":
        region_abcd_formatted = "all ABCD"
    else:
        region_abcd_formatted = "region {}".format(region_abcd)

    # Draw and print histogram.
    H1d_clone = H1d.Clone(histname+"_clone")

    C = TCanvas("c_"+histname, "", 800, 650)
    P_hist = TPad("pad_"+histname, "", 0, 0, 1, 1)
    P_hist.SetTopMargin(0.05)
    P_hist.SetBottomMargin(0.17)
    P_hist.SetLeftMargin(0.15)
    P_hist.SetRightMargin(0.05)
    P_hist.Draw()

    C.cd()
    """ title_line1 = "Angle between dimuon planes for {} {}".\
        format(region_abcd_formatted, mass_region_formatted)
    
    T_line1 = TPaveText(0.4, 0.92, 0.7, 0.95, "NDC")
    T_line1.AddText(title_line1)
    T_line1.SetTextSize(0.05)
    T_line1.SetFillColor(0)
    T_line1.SetBorderSize(0)
    T_line1.Draw() """

    P_hist.cd()

    H1d_clone.SetStats(0)
    H1d_clone.GetXaxis().SetTitle("#Delta#Phi")
    H1d_clone.GetXaxis().SetTitleSize(0.07)
    H1d_clone.GetXaxis().SetTitleOffset(1.)
    H1d_clone.GetXaxis().SetLabelSize(0.05)
    H1d_clone.GetYaxis().SetTitle("Events")
    H1d_clone.GetYaxis().SetTitleSize(0.07)
    H1d_clone.GetYaxis().SetTitleOffset(1.)
    H1d_clone.GetYaxis().SetLabelSize(0.05)

    H1d_clone.SetMarkerStyle(20)
    H1d_clone.SetMarkerSize(0.7)
    H1d_clone.SetMarkerColor(kBlack)
    H1d_clone.SetLineColor(kBlack)

    H1d_clone.Draw("e0")
    H1d_clone.Draw("p same")

    T = TLatex()
    T.SetNDC()
    T.SetTextSize(0.05)
    T.SetTextFont(42)      # non-bold Helvetica
    T.SetTextAlign(13)     # align top left
    T.DrawLatex(0.37, 0.92, "3.05 < #hat{m}_{#mu#mu} < 3.15 GeV/c^{2}")

    C.Print("fig_aplan_doubleJPsi_new_{}_{}.pdf".\
            format(mass_region, region_abcd))
    P_hist.Clear()
    P_hist.Update()
    C.Clear()
    C.Update()
    C.Close()



def aplanarity(mass_region, region_abcd):

    file_name = "dpvectors_doubleJPsi_new_{}_{}.json".\
                format(mass_region, region_abcd)
    
    with open(file_name, "r") as file:
        dpvectors = json.load(file)
    
    lplead_mu1 = dpvectors["l{}_plead_mu1_{}".\
                       format(mass_region, region_abcd)]
    lpsec_mu1 = dpvectors["l{}_psec_mu1_{}".\
                       format(mass_region, region_abcd)]
    lplead_mu2 = dpvectors["l{}_plead_mu2_{}".\
                       format(mass_region, region_abcd)]
    lpsec_mu2 = dpvectors["l{}_psec_mu2_{}".\
                       format(mass_region, region_abcd)]
    lplead_dimu = dpvectors["l{}_plead_dimu_{}".\
                       format(mass_region, region_abcd)]
    """ lpsec_dimu = dpvectors["l{}_psec_dimu_{}".\
                       format(mass_region, region_abcd)]
    lpfourmu = dpvectors["l{}_pfourmu_{}".\
                       format(mass_region, region_abcd)] """
    
    n = len(lvfourmu)
    #n = 1

    # List the momentum vectors (px, py, pz) in the four-muon rest frame.
    aplead_mu1 = array_pvectors(lplead_mu1, n)
    apsec_mu1 = array_pvectors(lpsec_mu1, n)
    aplead_mu2 = array_pvectors(lplead_mu2, n)
    apsec_mu2 = array_pvectors(lpsec_mu2, n)
    aplead_dimu = array_pvectors(lplead_dimu, n)
    """ apsec_dimu = array_pvectors(lpsec_dimu, n)
    avfourmu = array_pvectors(lpfourmu, n) """

    print("- Done lab frame p vectors for region {} ({})".\
        format(region_abcd, mass_region))

    print("-- Done four-muon rest frame p vectors for region {} ({})".\
        format(region_abcd, mass_region))

    # Calculate and plot the angle between dimuon planes.
    array_angles = array_interplane_angle(
                    aplead_mu1, apsec_mu1, 
                    aplead_mu2, apsec_mu2, 
                    aplead_dimu, n
                    )
    
    print("--- Done interplane angles for region {} ({})".\
        format(region_abcd, mass_region))
    
    abins = arange(-3.5, 3.6, 0.2)
    nbins = len(abins) - 1
    bins = [nbins, abins]

    th1f_angle_create_draw(
                    array_angles, bins, mass_region, region_abcd,
                    )



if __name__ == "__main__":

    print("################ Main function started ################")

    #mass_region = "all"
    #mass_region = "control1" 
    mass_region = "JPsi" 
    #mass_region = "control2"
    
    lregions_abcd = [
        "A",
        #"B", "C", "D",
        #"ABCD"
    ]

    # For each region ABCD, create parallel process to plot aplanarity.
    n = len(lregions_abcd)
    gROOT.SetBatch(True)
    lfutures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers = n) \
    as executor:
        
        for region_abcd in lregions_abcd:
            future = executor.submit(
                aplanarity,
                mass_region, region_abcd
            )
            lfutures.append(future)
     
    print("################# Main function ended #################")