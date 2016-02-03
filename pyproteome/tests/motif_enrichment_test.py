
from unittest import TestCase

import pyproteome as pyp


FOREGROUND = """
KQTDLELsPLTKEEK
SLSTKRSsPDDGNDV
DVSPYSLsPVSNKSQ
SSRDRWIsENQDSAD
ARRNDDIsELEDLSE
ISELEDLsELEDLKD
TAFGREHsPYGPSPL
EHSPYGPsPLGWPSS
EHSPYGPsPLGWPSS
SYDPYDFsDTEEEMP
SYDPYDFsDTEEEMP
EEMPQVHtPKTADSQ
SKSDLRKsPVFSDED
LRKSPVFsDEDSDLD
PVFSDEDsDLDFDIS
EDENGDItPIKAKKM
DDDLVEFsDLESEDD
VEFSDLEsEDDERPR
KLTDEDFsPFGSGGG
TAADMYLsPVRSPKK
MYLSPVRsPKKKGST
REKEAVItPVASATQ
KGSTLDLsDLEAEKL
STEDGDGtDDFLTDK
DGTDDFLtDKEDEKA
"""

BACKGROUND = """
KKMPLDLsPLATPII
SSPAQNWtPPQPRTL
PIRSSAFsPLGGCTP
FFRQRMFsPMEEKEL
DINTFVGtPVEKLDL
GKKTKFAsDDEHDEH
QEALDYFsDKESGKQ
PDKQFLIsPPASPPV
FLISPPAsPPVGWKQ
HDDFFSTtPLQHQRI
GEDEFVPsDGLDKDE
DEIFYTLsPVNGKIT
DELFYTLsPINGKIS
GKQPLLLsEDEEDTK
PPGDYSTtPGGTLFS
GGTLFSTtPGGTRII
PGGTLFStTPGGTRI
GGTLFSTtPGGTRII
LPHDYCTtPGGTLFS
GGTLFSTtPGGTRII
IDDDFFPsSGEEAEA
DDDFFPSsGEEAEAA
KQTDLELsPLTKEEK
ASELACPtPKEDGLA
TDSSSYPsPCASPSP
SYPSPCAsPSPPSSG
DVSPYSLsPVSNKSQ
SLSTKRSsPDDGNDV
DVSPYSLsPVSNKSQ
SSRDRWIsENQDSAD
NDIYNFFsPLNPVRV
ARRNDDIsELEDLSE
ARRNDDIsELEDLSE
ISELEDLsELEDLKD
ASRYYVPsYEEVMNT
GEDGLILtPLGRYQI
APEVWGLsPKNPEPD
DVEDMELsDVEDDGS
SEDGQVFsPKKGQKK
KEHRGCDsPDPDTSY
PDTSYVLtPHTEEKY
TAFGREHsPYGPSPL
EHSPYGPsPLGWPSS
EHSPYGPsPLGWPSS
GHRELVLsSPEDLTQ
HRELVLSsPEDLTQD
RSVNFSLtPNEIKVS
EQYGFLTtPTKQLGA
KSDISPLtPRESSPL
PLTPRESsPLYSPTF
RESSPLYsPTFSDST
SYDPYDFsDTEEEMP
SYDPYDFsDTEEEMP
EEMPQVHtPKTADSQ
SKSDLRKsPVFSDED
LRKSPVFsDEDSDLD
PVFSDEDsDLDFDIS
EDENGDItPIKAKKM
RLPPKVEsLESLYFT
SLESLYFtPIPARSQ
LNNSNLFsPVNRDSE
RDSENLAsPSEYPEN
KKEQMPLtPPRFDHD
DDDLVEFsDLESEDD
VEFSDLEsEDDERPR
STRRGTFsDQELDAQ
YPAVNRFsPSPRNSP
AVNRFSPsPRNSPRP
KDEILPTtPISEQKG
KLTDEDFsPFGSGGG
RPQDSEFsPVDNCLQ
DRSSGTAsSVAFTPL
TASSVAFtPLQGLEI
GEEPSEYtDEEDTKD
GEEPTVYsDEEEPKD
VLIVYELtPTAEQKA
TEDIFPVtPPELEET
TAADMYLsPVRSPKK
MYLSPVRsPKKKGST
REKEAVItPVASATQ
DPQQLQLsPLKGLSL
GGKRSRLtPVSPESS
GGKRSRLtPVSPESS
RSRLTPVsPESSSTE
LRNPYLLsEEEDDDV
VVRRRSFsISPVRLR
RRRSFSIsPVRLRRS
DRGEFSAsPMLKSGM
EHKELSNsPLRENSF
LRENSFGsPLEFRNS
KGSTLDLsDLEAEKL
DGTDDFLtDKEDEKA
STEDGDGtDDFLTDK
DGTDDFLtDKEDEKA
EEWDPEYtPKSKKYY
"""

OUTPUT = """
| .....D.x....... |     8 /    25 |    10 /    95 | 2.73E-04 |
| .....-.x....... |    12 /    25 |    20 /    95 | 3.27E-04 |
| .....-.s....... |    10 /    25 |    15 /    95 | 3.98E-04 |
| .....-.x-...... |     8 /    25 |    11 /    95 | 8.24E-04 |
| .....-.s....E.. |     5 /    25 |     5 /    95 | 9.17E-04 |
| ..-..D.x....... |     5 /    25 |     5 /    95 | 9.17E-04 |
| .....-.s-.-.... |     7 /    25 |     9 /    95 | 1.05E-03 |
| .....-.s....-.. |     6 /    25 |     7 /    95 | 1.17E-03 |
| .....D.x-...... |     6 /    25 |     7 /    95 | 1.17E-03 |
| .......x-...... |    12 /    25 |    23 /    95 | 2.09E-03 |
| .....D.s....... |     6 /    25 |     8 /    95 | 3.80E-03 |
| .....-.s-.E.... |     6 /    25 |     8 /    95 | 3.80E-03 |
| ..-..-.x....... |     6 /    25 |     8 /    95 | 3.80E-03 |
| .....-.xD...... |     6 /    25 |     8 /    95 | 3.80E-03 |
| .....-.sD.E.E.. |     4 /    25 |     4 /    95 | 3.97E-03 |
| .....-.s.L..-.. |     4 /    25 |     4 /    95 | 3.97E-03 |
| .......sD.-.-O. |     4 /    25 |     4 /    95 | 3.97E-03 |
| ..-..D.x-...... |     4 /    25 |     4 /    95 | 3.97E-03 |
| .....D.xD...... |     4 /    25 |     4 /    95 | 3.97E-03 |
| ...-...xP..S... |     4 /    25 |     4 /    95 | 3.97E-03 |
| .....D.s-.E.... |     5 /    25 |     6 /    95 | 4.48E-03 |
| .....-.sD.-.-.. |     5 /    25 |     6 /    95 | 4.48E-03 |
| .....-.s-L-.... |     5 /    25 |     6 /    95 | 4.48E-03 |
| .......xD.E.E.. |     5 /    25 |     6 /    95 | 4.48E-03 |
| .......s-...... |    10 /    25 |    19 /    95 | 5.81E-03 |
| .......xD.-.-.. |     7 /    25 |    11 /    95 | 6.47E-03 |
| .....-.s.L..... |     6 /    25 |     9 /    95 | 9.30E-03 |
| .......s-.-.... |     9 /    25 |    17 /    95 | 9.30E-03 |
| .......sD.-.-.. |     6 /    25 |     9 /    95 | 9.30E-03 |
| ....-..x....... |     6 /    25 |     9 /    95 | 9.30E-03 |
| .......x-.-.... |    10 /    25 |    20 /    95 | 9.60E-03 |
"""


class MotifEnrichmentFullTest(TestCase):
    """
    Test pyproteome.motif.motif_enrichment using a full data set, provided
    by Brian Joughin. Compares the results to those produced by his set of
    perl scripts.
    """
    def setUp(self):
        self.foreground = [
            i.strip()
            for i in FOREGROUND.split("\n")
            if i.strip()
        ]
        self.background = [
            i.strip()
            for i in BACKGROUND.split("\n")
            if i.strip()
        ]
        self.output = []

        for line in OUTPUT.split("\n"):
            line = line.strip()

            if not line:
                continue

            tokens = line.split("|")

            motif = pyp.Motif(tokens[1].strip())

            fore_hits = int(tokens[2].split("/")[0])
            fore_size = int(tokens[2].split("/")[1])

            back_hits = int(tokens[3].split("/")[0])
            back_size = int(tokens[3].split("/")[1])

            p_val = float(tokens[4])

            self.output.append(
                (
                    motif,
                    fore_hits, fore_size,
                    back_hits, back_size,
                    p_val,
                )
            )

    def test_motif_enrichment(self):
        hits = pyp.motif.motif_enrichment(
            self.foreground, self.background,
        )

        true_positives = [i[0] for i in self.output]
        true_hits = list(hits["Motif"])

        # Check for false positives
        for motif in hits["Motif"]:
            self.assertIn(motif, true_positives)

        # Check for false-negatives
        for motif in true_positives:
            self.assertIn(motif, true_hits)

        # Check every row of output is the same
        for pd_row, out_row in zip(hits.iterrows(), self.output):
            index, row = pd_row

            self.assertEqual(row["Motif"], out_row[0])
            self.assertEqual(row["Foreground Hits"], out_row[1])
            self.assertEqual(row["Foreground Size"], out_row[2])
            self.assertEqual(row["Background Hits"], out_row[3])
            self.assertEqual(row["Background Size"], out_row[4])
            self.assertLess(abs(row["p-value"] - out_row[5]), 0.001)
