**Project abstract**

This project presents a search for rare processes in paired dimuon events using 
CMS scouting data. We search for exotic Higgs bosons $\Phi$ associated with the 
$pp \rightarrow \Phi \rightarrow \phi\phi \rightarrow 4\mu$ process by comparing the 
four-muon mass $m_{4\mu}$ spectrum to estimated background. No new particles are 
observed. As a proof of principle, we also measure the asymmetry in pseudorapidity 
$\eta$ and the charge-parity violation (CPV) sensitive angle $\Delta\Phi$ in a known 
resonance of $J/\psi$ pair production events. We show that background events in the 
scouting dataset can be well-modeled, and our analysis is sensitive to signs of new 
physics.

![image](https://github.com/user-attachments/assets/a4982b45-2569-4580-baef-0cf587e62c8b)

----------

**Background estimation using the ABCD method**

Events with 4 or more muons and at least 1 pair of opposite charges are selected for 
analysis. The selected events are categorized into the signal region A and three 
control regions B, C, and D according to mass asymmetry $m_{asym}$ and four-muon 
charge $q_{4\mu}$. These are defined

$m_{asym}  = \frac{|m_{\mu\mu} - m_{\mu\mu}^{'}|} {|m_{\mu\mu} + m_{\mu\mu}^{'}|},
\qquad q_{4\mu} = |q_{\mu\mu}| + |q_{\mu\mu}^{'}|,$

where $m_{\mu\mu}$, $m_{\mu\mu}^{'}$ are the two dimuon masses, and $q_{\mu\mu}$,
$q_{\mu\mu}^{'}$ are the two dimuon charges. In each event, the dimuon pair is 
identified as the combination with lowest $m_{asym}$.

The distributions of $m_{asym}$ and $q_{4\mu}$ are each divided into a signal region 
and sideband. Since the expected $\Phi \rightarrow \phi\phi$ process involves the 
production of identical charge-neutral particles, we define the signal regions to be 
$m_{asym} < 0.05$ and $q_{4\mu} = 0$, with sidebands $0.075 < m_{asym} < 0.5$ and 
$q_{4\mu} \neq 0$ . Region A consists of events in the intersection of both signal 
regions, $\{m_{asym} < 0.05\ \cap\ q_{4\mu} = 0\}$. Region B consists of 
$\{m_{asym} < 0.05\ \cap\ q_{4\mu} \neq 0 \}$, while region C consists of 
$\{0.075 < m_{asym} < 0.5\ \cap\ q_{4\mu} = 0 \}$. Region D consists of events in the 
intersection of both sidebands, $\{0.075 < m_{asym} < 0.5\ \cap\ q_{4\mu} \neq 0 \}$. 
We expect the control regions to contain negligible amounts of signal and assume the 
relationship

$N_{A}^{bkg} = N_{B}\frac{N_{C}} {N_{D}}$

between the number $N_{A}^{bkg}$ of background events in region A and the total number
of events $N_{B}$, $N_{C}$, and $N_{D}$ in each of the control regions.

![image](https://github.com/user-attachments/assets/7de26254-27aa-4ee5-9684-a36eb9766719)

The event selection and division into ABCD regions is done using 
**bkgest_abcd_sort.py**.

In a given window of $\hat{m}\_{\mu\mu}$ values, the $m_{4\mu}$ distribution is 
analyzed for possible resonances corresponding to exotic Higgs bosons. In each 
$m_{4\mu}$ distribution, the number of background events is estimated on a bin-by-bin 
basis. Due to correlations between $m_{4\mu}$ and $\hat{m}\_{\mu\mu}$, the ratio 
$N_{C} / N_{D}$ is plotted with respect to $\alpha = \hat{m}\_{\mu\mu} / m_{4\mu}$. 
A second-degree polynomial function $R_{C/D}$ is fitted to the $\alpha$ distribution. 
For each bin of $m_{4\mu}$, the estimated background $N_{A}^{bkg}(m_{4\mu}, \alpha)$ 
is then plotted by weighting each event in region B by $R_{C/D}(\alpha)$.

Plots of the ratio $N_{C} / N_{D}$ and the background estimation $N_{A}^{bkg}$ are 
made using **bkgest_abcd_fit.py**.

----------

**Measurements of angular production distributions**

As a proof of principle, $\eta$ asymmetry and the CPV sensitive angle $\Delta\Phi$ are
measured for $J/\psi$ pair production events. The $\eta$ asymmetry 
$A_{2J/\psi}^{\eta}$ is a measure of directional bias in the $J/\psi$ pair production 
and is plotted as a function of the average pseudorapidity $|\hat{\eta}\_{J/\psi}|$ 
of each $J/\psi$ pair. For each bin of $|\hat{\eta}_{J/\psi}|$,

$A_{2J/\psi}^{\eta} = \frac{N(\Delta\eta > 0) - N(\Delta\eta < 0)} 
{N(\Delta\eta > 0) + N(\Delta\eta < 0)}$,

where $\Delta\eta$ is the difference in $\eta$ between the higher- $p_{T}$ and 
lower- $p_{T}$ $J/\psi$ mesons, and $N(\Delta\eta > 0)$, $N(\Delta\eta < 0)$ are the 
number of events with positive and negative $\Delta\eta$.

The plot of $A_{2J/\psi}^{\eta}$ vs. $|\hat{\eta}\_{J/\psi}|$ is made with
**eta_asym.py**, using the data sorted through **bkgest_abcd_sort.py**.


$\Delta\Phi$ is the angle between the two $J/\psi$ planes, defined in the four-muon 
rest frame as

$\Delta\Phi = a\arccos(\hat{n}_{1} \cdot \hat{n}\_{2})$,

where $\hat{n}\_{1}$, $\hat{n}\_{2}$ are the unit vectors normal to the 
higher- $p_{T}$ and lower- $p_{T}$ $J/\psi$ planes, and $a = \pm 1$ determines the 
positive/negative direction of $\Delta\Phi$. These are defined

$\hat{n} = \frac{\overrightarrow{p}\_{\mu^{-}} \times \overrightarrow{p}\_{\mu^{+}}}
{|\overrightarrow{p}\_{\mu^{-}} \times \overrightarrow{p}\_{\mu^{+}}|}, 
\qquad a = \frac{\overrightarrow{p}\_{1} \cdot (\hat{n}\_{1} \times \hat{n}\_{2})}
{|\overrightarrow{p}\_{1} \cdot (\hat{n}\_{1} \times \hat{n}_{2})|}$,

where $\overrightarrow{p}\_{\mu^{-}}$, $\overrightarrow{p}\_{\mu^{+}}$ are the 
momentum vectors of the $\mu^{-}\mu^{+}$ pair associated with a $J/\psi$ meson, 
and $\overrightarrow{p}\_{1}$ is the momentum vector of the higher- $p_{T}$ $J/\psi$.

![image](https://github.com/user-attachments/assets/7d4b9931-ce73-4495-811f-5cecfedacbe0)
(Figure from T. Agatonovic-Jovin et. al.)

Momentum vectors of the muons, dimuons, and four-muon particles of each event are
obtained using **pvectors_abcd_sort.py**. Then, plot of the $\Delta\Phi$ distribution
is made using **aplanarity.py**.

----------

**Comments on the code**

The files used in this project are

1. **bkgest_abcd_sort.py** for selecting events, obtaining relevant quantities
   (e.g. $\hat{m}\_{\mu\mu}$, $\alpha$), and sorting into ABCD regions,
2. **bkgest_abcd_fit.py** for plotting the ABCD background estimation,
3. **eta_asym.py** for plotting $A_{2J/\psi}^{\eta}$ vs. $|\hat{\eta}\_{J/\psi}|$,
4. **pvectors_abcd_sort.py** for obtaining the momentum vectors of particles and
   sorting events into ABCD regions,
5. **aplanarity.py** for plotting $\Delta\Phi$.

All plots are produced using PyROOT. Python files are designed to produce multiple
plots of different variables (e.g. $R_{C/D}(p_{T})$ instead of $R_{C/D}(\alpha)$ ), 
$\hat{m}\_{\mu\mu}$ windows, or ABCD regions in parallel. 
The variables and regions of interest can be easily specified/changed in the code.

----------

**References**

[1] A. M. Sirunyan et al., JINST 13 (P06015), DOI: 10.1088/1748-0221/13/06/P06015.
[2] S. Mukherjee, in LLC Workshop (Amsterdam Science Park, 2018).
[3] W. Buttinger, Background Estimation with the ABCD Method Featuring the TRooFit 
Toolkit (2018).
[4] T. Agatonovic-Jovin et al., in PoS(EPS-HEP2021), Vol. 398 (2022) p. 495.
[5] The CMS collaboration, J. High Energ. Phys. 94, DOI: 10.1007/JHEP09(2014)094.

