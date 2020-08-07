import pickle

# -=- OPTIM IE -=-

# Load optim_ie for all but HR 0.75"
with open('/Users/nflagey/PycharmProjects/MSE/mse_injeff/results/' + 'no_segments_injeff_optim_curve_allbuthr075.pkl', 'rb') as f:
    ie_optim_all = pickle.load(f)
# Load optim_ie for HR 0.75"
with open('/Users/nflagey/PycharmProjects/MSE/mse_injeff/results/' + 'no_segments_injeff_optim_curve_hr075.pkl', 'rb') as f:
    ie_optim_hr075 = pickle.load(f)

ie_optim_comb = {}

for iq_key in ie_optim_all:
    print(iq_key)
    ie_optim_comb[iq_key] = {}

    for zd_key in ie_optim_all[iq_key]:
        print(zd_key)
        ie_optim_comb[iq_key][zd_key] = {}

        for field_key in ie_optim_all[iq_key][zd_key]:
            print(field_key)
            ie_optim_comb[iq_key][zd_key][field_key] = {}

            for spec_key in ie_optim_all[iq_key][zd_key][field_key]:
                print(spec_key)
                ie_optim_comb[iq_key][zd_key][field_key][spec_key] = ie_optim_all[iq_key][zd_key][field_key][spec_key]

            for spec_key in ie_optim_hr075[iq_key][zd_key][field_key]:
                print(spec_key)
                ie_optim_comb[iq_key][zd_key][field_key][spec_key] = ie_optim_hr075[iq_key][zd_key][field_key][spec_key]

# Save the IE dictionary (do it at every field to have something to play with while it's running)
with open('/Users/nflagey/PycharmProjects/MSE/mse_injeff/results/' + 'no_segments_injeff_optim_curve_all.pkl', 'wb') as f:
    pickle.dump(ie_optim_comb, f, pickle.HIGHEST_PROTOCOL)




# Load simu_ie for all but HR 0.75"
with open('/Users/nflagey/PycharmProjects/MSE/mse_injeff/results/' + 'no_segments_injeff_curve_allbuthr075.pkl', 'rb') as f:
    ie_simu_all = pickle.load(f)

# Load simu_ie for HR 0.75"
with open('/Users/nflagey/PycharmProjects/MSE/mse_injeff/results/' + 'no_segments_injeff_curve_hr075.pkl', 'rb') as f:
    ie_simu_hr075 = pickle.load(f)

ie_simu_comb = {}
ie_simu_comb_sig = {}

for iq_key in ie_simu_all[0]:
    ie_simu_comb[iq_key] = {}
    ie_simu_comb_sig[iq_key] = {}

    for zd_key in ie_simu_all[0][iq_key]:
        ie_simu_comb[iq_key][zd_key] = {}
        ie_simu_comb_sig[iq_key][zd_key] = {}

        for field_key in ie_simu_all[0][iq_key][zd_key]:
            ie_simu_comb[iq_key][zd_key][field_key] = {}
            ie_simu_comb_sig[iq_key][zd_key][field_key] = {}

            for spec_key in ie_simu_all[0][iq_key][zd_key][field_key]:
                ie_simu_comb[iq_key][zd_key][field_key][spec_key] = ie_simu_all[0][iq_key][zd_key][field_key][spec_key]
                ie_simu_comb_sig[iq_key][zd_key][field_key][spec_key] = ie_simu_all[1][iq_key][zd_key][field_key][spec_key]

            for spec_key in ie_simu_hr075[0][iq_key][zd_key][field_key]:
                ie_simu_comb[iq_key][zd_key][field_key][spec_key] = ie_simu_hr075[0][iq_key][zd_key][field_key][spec_key]
                ie_simu_comb_sig[iq_key][zd_key][field_key][spec_key] = ie_simu_hr075[1][iq_key][zd_key][field_key][spec_key]

# Save the IE dictionary (do it at every field to have something to play with while it's running)
with open('/Users/nflagey/PycharmProjects/MSE/mse_injeff/results/' + 'no_segments_injeff_curve_all.pkl', 'wb') as f:
    pickle.dump([ie_simu_comb, ie_simu_comb_sig], f, pickle.HIGHEST_PROTOCOL)

