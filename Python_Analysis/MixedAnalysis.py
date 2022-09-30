#!/afs/cern.ch/user/j/jowulff/miniconda3/envs/shep/bin/python
import uproot
import numpy as np
import awkward as ak

import vector

import hist
from argparse import ArgumentParser
import correctionlib


def make_parser():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", type=str, 
                        default=None, help='Input Root File')
    parser.add_argument("-o", "--output", type=str,
                        default=None, help='Output File')
    parser.add_argument("-x", "--cross_section", type=float, nargs=1,
                        default=-1, help="Cross Section of MC Dataset")
    parser.add_argument("-l", "--int_luminosity", type=float, nargs=1,
                        default=-1, help="Integrated Lumi to scale to")
    parser.add_argument("-m", "--mc", action='store_true',
                        help='MC Flag. If set, mc==True')
    parser.add_argument("--first_data", help="Is this the first time", 
                        action='store_true',)
    #parser.add_argument("-s", "--signal", nargs=1, action='store_true',
    #help='Signal Flag. If set, signal==True. Only allowed if mc == True.')
    return parser


def load_corrector(file,):
    if file.endswith(".json.gz"):
        import gzip
        with gzip.open(file, 'rt') as file:
            data = file.read().strip()
            ceval = correctionlib.CorrectionSet.from_string(data)
    else:
        ceval = correctionlib.CorrectonSet.from_file(file)
    return ceval



def spacial_invert(vec):
    return -vec.to_Vector3D().to_Vector4D()+vector.obj(px=0, py=0, pz=0, E=vec.E)


def skimming(filename, ofilename, xs=None, lumi=None, mc_flag=False, first_data=False):
    correctionfiles = {
        'muon': "/afs/cern.ch/user/j/jowulff/Condor/TTbar/corrections/muon_Z.json.gz",
        'electron': "/afs/cern.ch/user/j/jowulff/Condor/TTbar/corrections/electron.json.gz",
        'pileup': "/afs/cern.ch/user/j/jowulff/Condor/TTbar/corrections/puWeights.json.gz",
        'jets': "/afs/cern.ch/user/j/jowulff/Condor/TTbar/corrections/jet_jerc.json.gz",
        'b_tag': "/afs/cern.ch/user/j/jowulff/Condor/TTbar/corrections/btagging.json.gz"
    }

    outfile = uproot.recreate(ofilename)
    out_dict = {'Muon_pt': np.float64,
                'Muon_eta': np.float64,
                'Muon_phi': np.float64,
                'Muon_mass': np.float64,
                'Electron_pt':  np.float64,
                'Electron_eta': np.float64,
                'Electron_phi': np.float64,
                'Electron_mass': np.float64,
                'Jet_pt': np.float64,
                'Jet_eta': np.float64,
                'Jet_phi': np.float64,
                'Jet_mass': np.float64,
                'mu_e_inv_mass': np.float64,
                'leading_lepton_pt': np.float64,
                'muon_corrections': np.float64,
                'electron_corrections': np.float64,
                'pu_corrections': np.float64,
                'b_tag_corrections': np.float64,
                'genWeight': np.float64,
                'weight': np.float64,
                'N_jet_loose': np.int32,
                'N_jet_tight': np.int32,
                'N_jet_medium': np.int32, }
    if mc_flag:
        out_dict['N_gen'] = np.int32
        out_dict['Sum_w'] = np.float64
    outfile.mktree("tout", out_dict)
    ## define histograms
    h_Muon_pt = hist.Hist(hist.axis.Regular(
        bins=100, start=0, stop=200, name="Muon pt"))
    h_Muon_eta = hist.Hist(hist.axis.Regular(
        bins=100, start=-5, stop=5, name="Muon eta"))
    h_Electron_pt = hist.Hist(hist.axis.Regular(
        bins=100, start=0, stop=200, name="Electron pt",))
    h_Electron_eta = hist.Hist(hist.axis.Regular(
        bins=100, start=-5, stop=5, name="Electron eta",))

    h_Muon_pt_weighted = hist.Hist(hist.axis.Regular(
        bins=100, start=0, stop=200, name="Muon pt"))
    h_Muon_eta_weighted = hist.Hist(hist.axis.Regular(
        bins=100, start=-5, stop=5, name="Muon eta"))
    h_Electron_pt_weighted = hist.Hist(hist.axis.Regular(
        bins=100, start=0, stop=200, name="Electron pt"))
    h_Electron_eta_weighted = hist.Hist(hist.axis.Regular(
        bins=100, start=-5, stop=5, name="Electron eta"))

    h_Muon_Electron_invariant_mass = hist.Hist(hist.axis.Regular(
        bins=100, start=12, stop=412, name="Muon Electron Inv. Mass"))
    h_Muon_Electron_invariant_mass_weighted = hist.Hist(hist.axis.Regular(
        bins=100, start=12, stop=412, name="Muon Electron Inv. Mass"))

    h_leading_lepton_pt = hist.Hist(hist.axis.Regular(
        bins=45, start=20, stop=200, name="leading lepton pt"))
    h_leading_lepton_pt_weighted = hist.Hist(hist.axis.Regular(
        bins=45, start=20, stop=200, name="leading lepton pt"))

    filter_names = ['/(Electron|Muon|Jet)_(pt|eta|phi|mass)/',
                    '/n(Muon|Electron|Jet)/',
                    '/(Electron|Muon)_charge/',
                    '/Muon_(triggerIdLoose)/',
                    '/Muon_(pfRelIso04_all|tightId)/',
                    '/HLT_(IsoMu24|Ele32_WPTight_Gsf)/',
                    'Electron_mvaFall17V2Iso_WP90',
                    '/Jet_(btagDeepB|btagDeepFlavB)/',
                    'Pileup_nTrueInt',
                    'genWeight']
    mc_filter_names = ['/(Electron|Muon)_genPart(Idx|Flav)/',
                       '/GenPart_(pdgId|genPartIdxMother)/',
                       'nGenPart']
    file = uproot.open(filename,)
    tree = file['Events']
    # ^ is xor in python
    trigger_cut = "HLT_IsoMu24 | HLT_Ele32_WPTight_Gsf"
    btag_deepflav_wp = {'loose': 0.0490, 'medium': 0.2783, 'tight': 0.71}
    # get the sumW from the runs tree
    if mc_flag:
        n_gen = file['Runs']['genEventCount'].array()
        Sum_W = file['Runs']['genEventSumw'].array()
    for events in tree.iterate(
            filter_name=filter_names+mc_filter_names, cut=trigger_cut,
            entry_stop=10000):
        if not mc_flag:
            if not first_data:
                cut = (events['HLT_IsoMu24']) & (events['HLT_Ele32_WPTight_Gsf'])
                events = events[~(cut)]
        # apply muon and electron cuts
        m_pt_cut = events['Muon_pt']>27
        m_eta_cut = np.abs(events['Muon_eta'])<2.4
        m_iso_cut = events['Muon_pfRelIso04_all']<0.15
        m_tight_cut = events["Muon_tightId"]

        muon_cuts = (m_pt_cut) &(m_eta_cut) & (m_iso_cut)&m_tight_cut
        events['muon_cuts'] = muon_cuts
        events = events[ak.any(muon_cuts, axis=1)]

        e_pt_cut =  events['Electron_pt']>35
        e_eta_cut = np.abs(events['Electron_eta'])<2.4
        e_iso_cut =  events["Electron_mvaFall17V2Iso_WP90"]

        electron_cuts = (e_pt_cut)&(e_eta_cut)&(e_iso_cut)
        events['electron_cuts'] = electron_cuts
        events = events[ak.any(electron_cuts, axis=1)]

        b_cut = events['Jet_btagDeepFlavB']>btag_deepflav_wp['medium']
        b_tag_eta = events['Jet_eta']
        events = events[ak.any(b_cut, axis=1)]
        # make sure that muon and electron are of opposite charge
        # this cut is applied after the others because you can only 
        # multiply two arrays if they both have at least one entry along the first dim
        # the pt cuts make sure of this: we only have events left where we have 
        # at least one electron and at leats one muon
        charge_cut = (events['Muon_charge', events['muon_cuts']][:, 0] * events['Electron_charge', events['electron_cuts']][:, 0]) < 0
        events = events[charge_cut]
        # create four-vectors
        events['Muon_pt'] = events['Muon_pt', events['muon_cuts']][:, 0]
        events['Muon_eta'] = events['Muon_eta', events['muon_cuts']][:, 0]
        events['Muon_phi'] = events['Muon_phi', events['muon_cuts']][:, 0]
        events['Muon_mass'] = events['Muon_mass', events['muon_cuts']][:, 0]
        events['Electron_pt'] = events['Electron_pt', events['electron_cuts']][:, 0]
        events['Electron_eta'] = events['Electron_eta', events['electron_cuts']][:, 0]
        events['Electron_phi'] = events['Electron_phi', events['electron_cuts']][:, 0]
        events['Electron_mass'] = events['Electron_mass', events['electron_cuts']][:, 0]
        # objects are ordered by pt. If there was more than one mu/e that passed all
        # requirements we take index 0 for highest pt
        muon_4d = vector.array({'pt':events['Muon_pt'],
                                'eta': events['Muon_eta'], 
                                'phi':events['Muon_phi'],
                                'mass': events['Muon_mass']})
        electron_4d = vector.array({'pt':events['Electron_pt'],
                                    'eta': events['Electron_eta'],
                                    'phi':events['Electron_phi'],
                                    'mass': events['Electron_mass']})
        # calculate the deltaR
        # cut out events with deltaR < 0.4
        delta_r_cut = muon_4d.deltaR(electron_4d) > 0.4
        muon_4d = muon_4d[delta_r_cut]; electron_4d=electron_4d[delta_r_cut]
        events = events[delta_r_cut]
        # how many jets would we have if we took other wp's?
        events['N_jet_loose'] = ak.sum(events['Jet_btagDeepFlavB']>btag_deepflav_wp['loose'], axis=1)
        events['N_jet_medium'] = ak.sum(events['Jet_btagDeepFlavB']>btag_deepflav_wp['medium'], axis=1)
        events['N_jet_tight'] = ak.sum(events['Jet_btagDeepFlavB']>btag_deepflav_wp['tight'], axis=1)
        # first drop events where there is no jet that passes the wp
        b_tag_medium_cut = events['Jet_btagDeepFlavB']>btag_deepflav_wp['medium']
        b_tag_eta = np.abs(events['Jet_eta'])<2.4
        events['b_tag_cut'] = b_tag_medium_cut&b_tag_eta
        muon_4d=muon_4d[ak.any(events['b_tag_cut'], axis=1)]
        electron_4d=electron_4d[ak.any(events['b_tag_cut'], axis=1)]
        events = events[ak.any(events['b_tag_cut'], axis=1)]
        # now our event has at least one jet passing the wp
        # event could look like this: [not_passing, not_passing, passing, not_passing, passing]
        # we want the first one that passes (index = 2)

        events['Jet_pt'] = events['Jet_pt', events['b_tag_cut']][:, 0]
        events['Jet_eta'] = events['Jet_eta', events['b_tag_cut']][:, 0]
        events['Jet_phi'] = events['Jet_phi', events['b_tag_cut']][:, 0]
        events['Jet_mass'] = events['Jet_mass', events['b_tag_cut']][:, 0]
        events['Jet_btagDeepFlavB'] = events['Jet_btagDeepFlavB', events['b_tag_cut']][:, 0]
        jet_4d = vector.array({'pt':events['Jet_pt'],
                                'eta': events['Jet_eta'],
                                'phi': events['Jet_phi'],
                                'mass': events['Jet_mass'],})
        #opposite_jet = spacial_invert(jet_4d)
        dR_mu_jet = muon_4d.deltaR(jet_4d)
        dR_e_jet = electron_4d.deltaR(jet_4d)
        dr_cut = (dR_mu_jet > 0.4) & (dR_e_jet > 0.4)
        events = events[dr_cut]
        muon_4d = muon_4d[dr_cut]
        electron_4d = electron_4d[dr_cut]
        jet_4d = jet_4d[dr_cut]
        #dR_mu_e = muon_4d.deltaR(electron_4d)

        ## coreections
        if mc_flag:
            muon_eval = load_corrector(correctionfiles['muon'])
            mu_c_trigger = muon_eval['NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight'].evaluate(
                '2018_UL', np.abs(muon_4d.eta), muon_4d.pt, 'sf')
            #mu_c_reco = muon_eval['NUM_GlobalMuons_DEN_genTracks'].evaluate(
            #    '2018_UL', np.abs(muon_4d.eta), muon_4d.pt, 'sf')
            mu_c_id = muon_eval['NUM_TightID_DEN_genTracks'].evaluate(
                '2018_UL', np.abs(muon_4d.eta), muon_4d.pt, 'sf')
            mu_c_iso = muon_eval['NUM_TightRelIso_DEN_TightIDandIPCut'].evaluate(
                '2018_UL', np.abs(muon_4d.eta), muon_4d.pt, 'sf')
            muon_c = mu_c_trigger*mu_c_id*mu_c_iso
            events['muon_corrections'] = muon_c
            electron_eval = load_corrector(correctionfiles['electron'])
            ele_c = electron_eval['UL-Electron-ID-SF'].evaluate(
                '2018', 'sf', 'wp90iso', electron_4d.eta, electron_4d.pt
            )
            events['electron_corrections'] = ele_c
            pu_eval = load_corrector(correctionfiles['pileup'])
            pu_c = pu_eval['Collisions18_UltraLegacy_goldenJSON'].evaluate(
                events['Pileup_nTrueInt'], 'nominal')
            events['pu_corrections'] = pu_c
            #jet_eval = load_corrector(correctionfiles['jets'])
            #jet_c = jet_eval['deepJet_mujets'].evaluate(
            b_tag_eval = load_corrector(correctionfiles['b_tag'])
            b_tag_c_wp = b_tag_eval['deepJet_mujets'].evaluate(
                'central', "M", 5, np.abs(jet_4d.eta), jet_4d.pt
            )
            b_tag_c_shape = b_tag_eval['deepJet_shape'].evaluate(
                'central', 5, np.abs(
                    jet_4d.eta), jet_4d.pt, events['Jet_btagDeepFlavB']
            )
            b_tag_c = b_tag_c_shape*b_tag_c_wp
            events['b_tag_corrections'] = b_tag_c
            corrections = muon_c * ele_c * pu_c * b_tag_c

            events['weight'] = events['genWeight']*lumi*xs*corrections

        # fill the histograms
        h_Muon_pt.fill(muon_4d.pt)
        h_Muon_eta.fill(muon_4d.eta)
        h_Electron_pt.fill(electron_4d.pt)
        h_Electron_eta.fill(electron_4d.eta)

        events['mu_e_inv_mass'] = (muon_4d + electron_4d).M
        events['leading_lepton_pt'] = np.max(
            [muon_4d.pt, electron_4d.pt], axis=0)
        h_Muon_Electron_invariant_mass.fill(events['mu_e_inv_mass'],)
        h_leading_lepton_pt.fill(events['leading_lepton_pt'])

        tout_dict = {'Muon_pt': events["Muon_pt"],
                     'Muon_eta': events["Muon_eta"],
                     'Muon_phi': events["Muon_phi"],
                     'Muon_mass': events["Muon_mass"],
                     'Electron_pt': events["Electron_pt"],
                     'Electron_eta': events["Electron_eta"],
                     'Electron_phi': events["Electron_phi"],
                     'Electron_mass': events["Electron_mass"],
                     'N_jet_loose': events["N_jet_loose"],
                     'N_jet_medium': events["N_jet_medium"],
                     'N_jet_tight': events["N_jet_tight"],
                     'Jet_pt': events["Jet_pt"],
                     'Jet_eta': events["Jet_eta"],
                     'Jet_phi': events["Jet_phi"],
                     'Jet_mass': events["Jet_mass"],
                     'mu_e_inv_mass': events["mu_e_inv_mass"],
                     'leading_lepton_pt': events["leading_lepton_pt"]}
        if mc_flag:
            events['N_gen'] = np.repeat(n_gen, len(events['Muon_pt']))
            events['Sum_w'] = np.repeat(Sum_W, len(events['Muon_pt']))
            tout_dict['genWeight'] = events["genWeight"]
            tout_dict['muon_corrections'] = events["muon_corrections"]
            tout_dict['electron_corrections'] = events["electron_corrections"]
            tout_dict['pu_corrections'] = events["pu_corrections"]
            tout_dict['b_tag_corrections'] = events["b_tag_corrections"]
            tout_dict['weight'] = events["weight"]
            tout_dict['N_gen'] = events['N_gen']
            tout_dict['Sum_w'] = events['Sum_w']
            outfile['tout'].extend(tout_dict)
            h_Muon_Electron_invariant_mass_weighted.fill(
                events['mu_e_inv_mass'], weight=events['weight'])
            h_leading_lepton_pt_weighted.fill(
                events['leading_lepton_pt'], weight=events['weight'])

            h_Muon_eta_weighted.fill(muon_4d.eta, weight=events['weight'])
            h_Muon_pt_weighted.fill(muon_4d.pt, weight=events['weight'])
            h_Electron_pt_weighted.fill(electron_4d.pt, weight=events['weight'])
            h_Electron_eta_weighted.fill(electron_4d.eta, weight=events['weight'])
            outfile['h_Muon_pt_weighted'] = h_Muon_pt_weighted
            outfile['h_Muon_eta_weighted'] = h_Muon_eta_weighted
            outfile['h_Electron_pt_weighted'] = h_Electron_pt_weighted
            outfile['h_Muon_eta_weighted'] = h_Muon_eta_weighted
            outfile['h_Muon_Electron_invariant_mass_weighted'] = h_Muon_Electron_invariant_mass_weighted
            outfile['h_leading_lepton_pt_weighted'] = h_leading_lepton_pt_weighted
            outfile['h_Muon_pt'] = h_Muon_pt
            outfile['h_Muon_eta'] = h_Muon_eta
            outfile['h_Electron_pt'] = h_Electron_pt
            outfile['h_Electron_eta'] = h_Electron_eta
            outfile['h_Muon_Electron_invariant_mass'] = h_Muon_Electron_invariant_mass
            outfile['h_leading_lepton_pt'] = h_leading_lepton_pt

        else:
            outfile['tout'].extend(tout_dict)
            outfile['h_Muon_pt'] = h_Muon_pt
            outfile['h_Muon_eta'] = h_Muon_eta
            outfile['h_Electron_pt'] = h_Electron_pt
            outfile['h_Electron_eta'] = h_Electron_eta
            outfile['h_Muon_Electron_invariant_mass'] = h_Muon_Electron_invariant_mass
            outfile['h_leading_lepton_pt'] = h_leading_lepton_pt


def main():
    parser = make_parser()
    args = parser.parse_args()
    skimming(filename=args.input,
             ofilename=args.output,
             xs=args.cross_section,
             lumi=args.int_luminosity,
             mc_flag=args.mc, 
             first_data=args.first_data)


if __name__ == "__main__":
    main()
