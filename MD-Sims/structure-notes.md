## PCBP1

- domain bisector residue: ALA 93
- Unstructured NTD: LEU 13
- Unstructured CTD: (less so) SER 174

### notes on MSMs (2025-12-02)

- L100Q-short-rad-0.30:
  - lag: 70
  - reason: lower Gini of eq-probs and flatlining of second and third ITS.
  - rationale: ITS plot looks like hell, sampling issues are obvious here.
- wt-short-rad-0.30:
  - lag: 150
  - rationale: lower Gini, flatlining of all three its without preposterous longest timescale behavior.
  - comment: Ditto on L100Q. There are clearly points where the eigensolver for the transition matrix NaNs out.

## UP1

- domain bisector residue: GLU 93
- Unstructured NTD: PRO 10
- Unstructured CTD: LEU 182

### notes on MSMS (2025-12-02)
- M72I-short-rad-0.30:
  - lag: 160
  - rationale: Also works for the WT, longest timescale never really turns over.
  - comment: obviously not converged (ITS drift), but counts are less sensitive to lag-length. No weirdness in Gini.
- WT-short-rad-0.30:
  - lag: 160
  - rationale: avoids jogs in Gini and longest ITS, looks pretty turned over.
  - comment: probably the best looking of the whole bunch, until the jog at lag 180.

## hnRNPA2B1

- domain bisector residue: GLY 91
- Unstructured NTD: ARG 2
- Unstructured CTD: LEU 174

### notes on MSMs (2025-12-02)

- D76H-rad-0.30:
  - lag: 100
  - rationale: shortest lagtime available with apparently level ITS before jags in Gini.
  - comment: This looks really converged until longer lags cut out some key transition. Surprising because of the state count (171).
- D76N-rad-0.30:
  - lag: 80
  - rationale: Lowest gini range in dataset, also just before large jogs emerge in ITS and Gini.
  - comment: Something important for counts matrix connectivity is getting cut out by having lags much longer than 80.
- D76V-rad-0.30:
  - lag: 50
  - rationale: less irregularity after lower-lag jog in this region, lower gini than the rest of the lag series.
  - comment: this looks hideous; more like I would expect from a 222 center model with around 8 microseconds of data.
- WT-rad-0.30:
  - lag: 150
  - rationale: Basically the only lag with a Gini below 0.95. ITS are flat from 150 to very high values.
  - comment: The gini being pinned at 1 for most lag choices suggests this molecule either doesn't move very much, or sampling is more poor than the ITS make it look.
  