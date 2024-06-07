import ROOT # tells python to use ROOT commands 
from ROOT import * # tells python to just remember the names of all the 
					# ROOT functions, instead of havign to call ROOT first
from pylab import *

F = TFile("/afs/crc.nd.edu/user/g/gziemyt2/Public/scouting_ntuple_1.root") 
								# loads the file containing the data
T = F.Get("scoutingntuplizer") # loads the TTree we want to use
							# (some files might have more than one TTree)

m = 0.1057					# mass of muon (GeV/c^2)

data_pt = []
data_eta = []
data_phi = [] 
#data_four_muon = []			# muon four vectors
data_angle_diff = []		# angle difference between first two muons
							# in each event
data_four_mother = []		# mother particle four vectors
data_dimuon_mass = []

# Creating histograms for pT, eta, phi.
H_pt = TH1F("muon_pT","Muon pT", 500, 0., 40.)
H_eta = TH1F("muon_eta","Muon eta", 500, -20., 20.)
H_phi = TH1F("muon_phi", "Muon phi", 500, -10., 10.)

# Histograms using data from four-vectors.
H_angle_diff = TH1F("muon_angle_diff", \
					"Angle difference between muon four-vectors", \
						100, 0., 10.)
H_dimuon_mass = TH1F("dimuon_mass", "Dimuon mass", 500, 0., 20.)

i = 0 						# integer we'll use to count events
for e in T: 				# loop over all the events in the tree
	i+=1 					# increment i by 1
	#if i > 1000: break 		# stop the loop once we've looked at 100 events
	n_muons = T.muon_num 	
	# makes an integer "n_muons" which is the number of muons in that 
	# events. The "T dot word" is where ROOT actually looks this value up 
	# in the tree. The for loop over events automatically loads the next 
	# events in to T, so that when you do T.something, it's looking at the 
	# something for that particualr event.

	""" print("This event has", n_muons, " muons") 
	# sanity check (you should see the output appear in your termianl.
	print("Their charges are:")
	for i_q in range(n_muons): 
		# we now know how many muons there are, so we can loop over them
		print(T.muon_q[i_q]) 
		# just printing the charge for each muon. 
		# T.muon_q in this case is an array! """
	
	if n_muons != 4: continue

	list_four_muon = []
	# Creates an empty list to hold muon four vectors for each new event.

	for n in range(n_muons):
		n_pt = T.muon_pt[n]
		n_eta = T.muon_eta[n]
		n_phi = T.muon_phi[n]

		# Filling histograms with data.
		H_pt.Fill(n_pt)
		H_eta.Fill(n_eta)
		H_phi.Fill(n_phi)


		# Listing p_T, eta, phi of each muon.
		data_pt.append(n_pt)
		data_eta.append(n_eta)
		data_phi.append(n_phi)

		# Listing muon four-vectors in the same event.
		V = TLorentzVector()
		V.SetPtEtaPhiM(n_pt, n_eta, n_phi, m)
		list_four_muon.append(V)
		
		""" data_four_muon.append(V)
		# List of muon four vectors of all events. """
	
	# Angle difference between first two muon four-vectors.
	angle_diff = list_four_muon[0].DeltaR(list_four_muon[1])
	H_angle_diff.Fill(angle_diff)

	# Listing the angle between the first two muon four vectors
	# in each event.
	data_angle_diff.append(list_four_muon[0].DeltaR(list_four_muon[1]))

	
	# Mother particle four vector obtained as the sum of the 
	# first two muon four vectors in each event.
	V_mother = list_four_muon[0] + list_four_muon[1]

	# Invariant mass of mother particles.
	mother_mass = V_mother.M()
	H_dimuon_mass.Fill(mother_mass)

	# Listing mass of mother particles.
	data_four_mother.append(mother_mass)

# Drawing and saving histograms.
C = TCanvas()
C.Divide(2,3)
C.cd(1)						# Changes directory to canvas section 1.
H_pt.Draw("hist")			# Draws histogram on canvas.
#C.Print("hist_pT_4muon.pdf")	# Saves pdf of canvas.

C.cd(2)
H_eta.Draw("hist")
#C.Print("hist_eta_4muon.pdf")

C.cd(3)
H_phi.Draw("hist")
#C.Print("hist_phi_4muon.pdf")

C.cd(4)
H_angle_diff.Draw("hist")
#C.Print("hist_four_muon_angle_diff_4muon.pdf")

C.cd(5)
H_dimuon_mass.Draw("hist")
#C.Print("hist_dimuon_mass_4muon.pdf")

C.Print("hist_4muon_events.pdf")

# Plotting
fig, ((ax_pt, ax_eta_phi), (ax_eta, ax_phi), \
	  (ax_angle_diff, ax_dimuon_mass)) = subplots(3,2)

# (1) histogram of p_T.
ax_pt.hist(data_pt, bins=500)
ax_pt.set_xlabel("p_T")

# (2) phi vs eta.
ax_eta_phi.plot(data_eta, data_phi, "k.")
ax_eta_phi.set_xlabel("Eta")
ax_eta_phi.set_ylabel("Phi")

# (3) histogram of eta.
ax_eta.hist(data_eta, bins=500)
ax_eta.set_xlabel("Eta")

# (4) histogram of phi.
ax_phi.hist(data_phi, bins=500)
ax_phi.set_xlabel("Phi")

# (5) histogram of angle differences between muons four vectors.
ax_angle_diff.hist(data_angle_diff, bins=500)
ax_angle_diff.set_xlabel("Angle difference between muon four vectors")

# (6) histogram of dimuon masses.
ax_dimuon_mass.hist(data_dimuon_mass, bins=500)
ax_dimuon_mass.set_xlabel("Dimuon mass")

show()
