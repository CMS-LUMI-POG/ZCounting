import json
import pdb
import os
import pickle

os.sys.path.append(os.path.expandvars('$CMSSW_BASE/src/ZCounting/'))
from ZUtils.python.utils import exp, pexp, plinear, plinear_step, pquad

def apply_muon_prefire(df, apply_on="recZCount", region="I"):

    resource = 'res/prefiring.json'

    # --- get prefire corrections
    with open(resource) as file_prefire:
        res_prefire = json.load(file_prefire)

    if apply_on not in df.keys():
        print("WARNING: Requested key was not found in input dataframe")
        return

    print("apply muon prefire corrections")

    for lo, hi, era in (
        (272007, 278769, "2016preVFP"),
        (278769, 280919, "2016postVFP"),
        (280919, 284045, "2016H"),
        (297020, 307083, "2017"),
        (315252, 325274, "2018"),
        ):

        if region not in res_prefire[era]["prefMuon"].keys():
            print("WARNING: Requested key was not found in resource `{0}`".format(resource))
            continue

        loc = (df['run'] >= lo) & (df['run'] < hi)
        if sum(loc) == 0:
            continue

        factor = float(res_prefire[era]["prefMuon"][region]['nominal'])

        # for debugging
        print("Apply on {0} factors of: {1}".format(era, factor))

        df.loc[loc, apply_on] *= factor

def apply_ECAL_prefire(df, apply_on="recZCount", region="I"):

    resource = 'res/prefiring.json'

    # --- get prefire corrections
    with open(resource) as file_prefire:
        res_prefire = json.load(file_prefire)

    if apply_on not in df.keys():
        print("WARNING: Requested key was not found in input dataframe")
        return

    print("apply ECAL prefire corrections")

    for lo, hi, era, func in (
        (272007, 278769, "2016preVFP", pexp),
        (278769, 284045, "2016postVFP", pexp),
        (297020, 306926, "2017", pexp),
        (306926, 307083, "2017H", pquad),
        ):

        if region not in res_prefire[era]["prefECAL"].keys():
            print("WARNING: Requested key was not found in resource `{0}`".format(resource))
            continue

        loc = (df['run'] >= lo) & (df['run'] < hi)
        if sum(loc) == 0:
            continue

        params = res_prefire[era]["prefECAL"][region]['nominal']
        factor = func(df.loc[loc,'pileUp'], *params)
        # for debugging
        print("Apply on {0} factors of (mean): {1}".format(era, sum(factor)/len(factor)))

        df.loc[loc, apply_on] *= factor