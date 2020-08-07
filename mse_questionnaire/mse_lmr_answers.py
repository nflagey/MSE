from astropy.io import ascii
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

# Read all answers from CSV file
file = 'MSECapabilitiesScienceTeamQuestionnaireResponses.csv'
ans = ascii.read(file, header_start=0, data_start=1, encoding='utf-8-sig')


# Number of valid answers
n_ans = len(ans[(ans['swg'] != 'None') & (ans['swg'] != '0')])


# How many need HR and LMR?
n_need_both = len(ans[(ans['hr_need'] == 'Yes') & (ans['lmr_need'] == 'Yes')])
print(f"Fraction of users who need HR and LMR: {n_need_both / n_ans * 100:4.1f} %\n")


# How many need HR?
n_need_hr = len(ans[ans['hr_need'] == 'Yes'])
print(f"Fraction of users who need HR: {n_need_hr / n_ans * 100:4.1f} %\n")
hr = ans[ans['hr_need'] == 'Yes']


# How many need HR resolution > 2000? 4000?
n_hr_blue_resol = len(hr[hr['hr_blue_resol'] > 0])
n_hr_green_resol = len(hr[hr['hr_green_resol'] > 0])
n_hr_red_resol = len(hr[hr['hr_red_resol'] > 0])
print(f"Number of valid answers for the HR (blue): {n_hr_blue_resol} ")
print(f"Number of valid answers for the HR (green): {n_hr_green_resol} ")
print(f"Number of valid answers for the HR (red): {n_hr_red_resol} \n")

for i in [20000, 30000, 40000, 50000, 60000]:
    n_r20000 = len(hr[(hr['hr_blue_resol'] <= i) & (hr['hr_blue_resol'] > 0)])
    print(f"Fraction of HR users who need R<={i} (blue): {n_r20000 / n_hr_blue_resol * 100:4.1f} %")
    n_r20000 = len(hr[(hr['hr_green_resol'] <= i) & (hr['hr_green_resol'] > 0)])
    print(f"Fraction of HR users who need R<={i} (green): {n_r20000 / n_hr_green_resol * 100:4.1f} %")
    n_r20000 = len(hr[(hr['hr_red_resol'] <= i) & (hr['hr_red_resol'] > 0)])
    print(f"Fraction of HR users who need R<={i} (red): {n_r20000 / n_hr_red_resol * 100:4.1f} %\n")


# Cumulative plot for resolution needs in HR
cumul_blue = cumul_green = cumul_red = []
for i in range(10100, 60100, 100):
    cumul_blue = cumul_blue + [len(hr[(hr['hr_blue_resol'] <= i) & (hr['hr_blue_resol'] > 0)]) / n_hr_blue_resol]
    cumul_green = cumul_green + [len(hr[(hr['hr_green_resol'] <= i) & (hr['hr_green_resol'] > 0)]) / n_hr_green_resol]
    cumul_red = cumul_red + [len(hr[(hr['hr_red_resol'] <= i) & (hr['hr_red_resol'] > 0)]) / n_hr_red_resol]
plt.figure()
plt.plot(range(10100, 60100, 100), cumul_blue, c='b')
plt.plot(range(10100, 60100, 100), cumul_green, c='g')
plt.plot(range(10100, 60100, 100), cumul_red, c='xkcd:orange red')
plt.xlabel('Spectral resolution')
plt.ylabel('Fraction of valid participants')
plt.savefig('hr_resolution_cumul.png')


# Plots for HR
#   Resolution and wmin/wmax
plt.figure()
for i in range(len(hr)):
    plt.plot([hr['hr_blue_wmin'][i],hr['hr_blue_wmax'][i]],
             [hr['hr_blue_resol'][i], hr['hr_blue_resol'][i]], alpha=.25, c='b')
    plt.plot([hr['hr_green_wmin'][i],hr['hr_green_wmax'][i]],
             [hr['hr_green_resol'][i], hr['hr_green_resol'][i]], alpha=.25, c='g')
    plt.plot([hr['hr_red_wmin'][i],hr['hr_red_wmax'][i]],
             [hr['hr_red_resol'][i], hr['hr_red_resol'][i]], alpha=.25, c='xkcd:orange red')
plt.xlabel('Wavelength')
plt.ylabel('Spectral resolution')
plt.xlim([300, 1000])
plt.savefig('hr_resolution.png')

#   Resolution vs wmin/wmax
plt.figure()
sel = (hr['hr_blue_resol'] != 0) & (hr['hr_blue_wmax'] - hr['hr_blue_wmin'] !=0)
plt.scatter(hr[sel]['hr_blue_wmax'] - hr[sel]['hr_blue_wmin'], hr[sel]['hr_blue_resol'], alpha=.25, c='b')
sel = (hr['hr_green_resol'] != 0) & (hr['hr_green_wmax'] - hr['hr_green_wmin'] != 0)
plt.scatter(hr[sel]['hr_green_wmax'] - hr[sel]['hr_green_wmin'], hr[sel]['hr_green_resol'], alpha=.25, c='g')
sel = (hr['hr_red_resol'] != 0) & (hr['hr_red_wmax'] - hr['hr_red_wmin'] != 0)
plt.scatter(hr[sel]['hr_red_wmax'] - hr[sel]['hr_red_wmin'], hr[sel]['hr_red_resol'], alpha=.25, c='xkcd:orange red')
plt.xlabel('Wavelength')
plt.ylabel('Spectral resolution')
plt.savefig('hr_resolution_vs_width.png')

#  Cumulative distribution of spectral coverage
sel = (hr['hr_blue_resol'] != 0) & (hr['hr_blue_wmax'] - hr['hr_blue_wmin'] !=0)
res_blue = stats.cumfreq(hr[sel]['hr_blue_wmax'] - hr[sel]['hr_blue_wmin'], numbins=10)
x_blue = res_blue.lowerlimit + np.linspace(0, res_blue.binsize*res_blue.cumcount.size, res_blue.cumcount.size)
sel = (hr['hr_green_resol'] != 0) & (hr['hr_green_wmax'] - hr['hr_green_wmin'] != 0)
res_green = stats.cumfreq(hr[sel]['hr_green_wmax'] - hr[sel]['hr_green_wmin'], numbins=18)
x_green = res_green.lowerlimit + np.linspace(0, res_green.binsize*res_green.cumcount.size, res_green.cumcount.size)
sel = (hr['hr_red_resol'] != 0) & (hr['hr_red_wmax'] - hr['hr_red_wmin'] != 0)
res_red = stats.cumfreq(hr[sel]['hr_red_wmax'] - hr[sel]['hr_red_wmin'], numbins=30)
x_red = res_red.lowerlimit + np.linspace(0, res_red.binsize*res_red.cumcount.size, res_red.cumcount.size)
plt.figure()
plt.step(x_blue, res_blue.cumcount / res_blue.cumcount[-1], c='b')
plt.step(x_green, res_green.cumcount / res_green.cumcount[-1], c='g')
plt.step(x_red, res_red.cumcount / res_red.cumcount[-1], c='xkcd:orange red')
plt.ylim([0, 1])
plt.ylabel('Fraction of respondents')
plt.xlabel('Wavelength range (nm)')
plt.savefig('hr_width_cumul.png')

#   HR arms
hr_wav = np.arange(360, 901)
hr_arm = np.zeros_like(hr_wav)
hr_arm_n = np.zeros_like(hr_wav)
for i in range(len(hr)):
    hr_arm[(hr_wav >= hr['hr_blue_wmin'][i]) & (hr_wav <= hr['hr_blue_wmax'][i])] += 1
    hr_arm[(hr_wav >= hr['hr_green_wmin'][i]) & (hr_wav <= hr['hr_green_wmax'][i])] += 1
    hr_arm[(hr_wav >= hr['hr_red_wmin'][i]) & (hr_wav <= hr['hr_red_wmax'][i])] += 1
    if (hr['hr_blue_wmin'][i] != 360) | (hr['hr_blue_wmax'][i] != 460):
        hr_arm_n[(hr_wav >= hr['hr_blue_wmin'][i]) & (hr_wav <= hr['hr_blue_wmax'][i])] += 1
    if (hr['hr_green_wmin'][i] != 440) | (hr['hr_green_wmax'][i] != 620):
        hr_arm_n[(hr_wav >= hr['hr_green_wmin'][i]) & (hr_wav <= hr['hr_green_wmax'][i])] += 1
    if (hr['hr_red_wmin'][i] != 600) | (hr['hr_red_wmax'][i] != 900):
        hr_arm_n[(hr_wav >= hr['hr_red_wmin'][i]) & (hr_wav <= hr['hr_red_wmax'][i])] += 1
plt.figure()
plt.step(hr_wav, hr_arm, where='mid', label='All answers')
plt.step(hr_wav, hr_arm_n, where='mid', label='Only narrow')
plt.ylabel('Number of respondents')
plt.ylim([0, 19])
plt.xlabel('Wavelength (nm)')
plt.legend()
plt.savefig('hr_arms.png')
# More than 12 respondents
print(hr_wav[hr_arm > 12])
print(hr_wav[hr_arm_n > 7])

# HR multiplexing
res_blue = stats.cumfreq(hr[sel]['hr_blue_wmax'] - hr[sel]['hr_blue_wmin'], numbins=10)


# How many need LMR?
n_need_lmr = len(ans[ans['lmr_need'] == 'Yes'])
print(f"Fraction of users who need LMR: {n_need_lmr / n_ans * 100:4.1f} %\n")
lmr = ans[ans['lmr_need'] == 'Yes']


# Plot needs in LMR mode
#    Simple histogram
plt.figure()
n_blue = len(lmr[lmr['lmr_blue_need'] == 'Yes'])
n_green = len(lmr[lmr['lmr_green_need'] == 'Yes'])
n_red = len(lmr[lmr['lmr_red_need'] == 'Yes'])
n_J = len(lmr[lmr['lmr_j_need'] == 'Yes'])
n_H = len(lmr[lmr['lmr_h_need'] == 'Yes'])
print(n_blue, n_green, n_red, n_J, n_H)
plt.step(['', 'blue', 'green', 'red', 'J', 'H', ' '], [0, 0, 0, 0, 0, 0, 0], where='mid', c='k')
plt.step(['', 'blue', 'green'], [0, n_blue, 0], where='mid', c='b', linestyle=':')
plt.step(['blue', 'green', 'red'], [0, n_green, 0], where='mid', c='g', linestyle=':')
plt.step(['green', 'red', 'J'], [0, n_red, 0], where='mid', c='xkcd:orange red', linestyle=':')
plt.step(['red', 'J', 'H'], [0, n_J, 0], where='mid', c='xkcd:red', linestyle=':')
plt.step(['J', 'H', ' '], [0, n_H, 0], where='mid', c='xkcd:burgundy', linestyle=':')
plt.title('LMR arms need')
plt.savefig('lmr_needs_hist.png')
#   Connecting needs per user (sorting them nicely)
dic_lmr_score = {'lmr_blue_need': 455., 'lmr_green_need': 650., 'lmr_red_need': 875., 'lmr_j_need': 1150., 'lmr_h_need': 1550.}
lmr['lmr_score'] = 0.
for i in range(len(lmr)):
    count = 0
    for key in dic_lmr_score:
        if lmr[key][i] == 'Yes':
            lmr['lmr_score'][i] += dic_lmr_score[key]
            count += 1.
    lmr['lmr_score'][i] /= count
lmr.sort('lmr_score')
for text in ['with', 'without']:
    plt.figure()
    for i in range(len(lmr)):
        if lmr['lmr_blue_need'][i] == 'Yes':
            if lmr['lmr_blue_wmin'][i] == 0:
                plt.plot([360, 550], [i, i], c='b', ls=':')
            else:
                plt.plot([lmr['lmr_blue_wmin'][i], lmr['lmr_blue_wmax'][i]], [i, i], c='b')
            if text == 'with':
                plt.text(455, i+.1, lmr['lmr_blue_resol'][i], size=6, horizontalalignment='center')
        if lmr['lmr_green_need'][i] == 'Yes':
            if lmr['lmr_green_wmin'][i] == 0:
                plt.plot([550, 750], [i, i], c='g', ls=':')
            else:
                plt.plot([lmr['lmr_green_wmin'][i], lmr['lmr_green_wmax'][i]], [i, i], c='g')
            if text == 'with':
                plt.text(650, i+.1, lmr['lmr_green_resol'][i], size=6, horizontalalignment='center')
        if lmr['lmr_red_need'][i] == 'Yes':
            if lmr['lmr_red_wmin'][i] == 0:
                plt.plot([750, 1000], [i, i], c='xkcd:orange red', ls=':')
            else:
                plt.plot([lmr['lmr_red_wmin'][i], lmr['lmr_red_wmax'][i]], [i, i], c='xkcd:orange red')
            if text == 'with':
                plt.text(875, i+.1, lmr['lmr_red_resol'][i], size=6, horizontalalignment='center')
        if lmr['lmr_j_need'][i] == 'Yes':
            if lmr['lmr_j_wmin'][i] == 0:
                plt.plot([1000, 1300], [i, i], c='xkcd:red', ls=':')
            else:
                plt.plot([lmr['lmr_j_wmin'][i], lmr['lmr_j_wmax'][i]], [i, i], c='xkcd:red')
            if text == 'with':
                plt.text(1150, i+.1, lmr['lmr_j_resol'][i], size=6, horizontalalignment='center')
        if lmr['lmr_h_need'][i] == 'Yes':
            if lmr['lmr_h_wmin'][i] == 0:
                plt.plot([1450, 1800], [i, i], c='xkcd:burgundy', ls=':')
            else:
                plt.plot([lmr['lmr_h_wmin'][i], lmr['lmr_h_wmax'][i]], [i, i], c='xkcd:burgundy')
            if text == 'with':
                plt.text(1550, i+.1, lmr['lmr_h_resol'][i], size=6, horizontalalignment='center')
    plt.title('LMR arms need')
    plt.xlabel('Wavelength')
    plt.xlim([350, 1800])
    plt.savefig('lmr_needs_'+text+'_resolution.png')
#   Resolution and wmin/wmax
plt.figure()
for i in range(len(lmr)):
    plt.plot([lmr['lmr_blue_wmin'][i], lmr['lmr_blue_wmax'][i]],
             [lmr['lmr_blue_resol'][i], lmr['lmr_blue_resol'][i]], alpha=.25, c='b')
    plt.plot([lmr['lmr_green_wmin'][i], lmr['lmr_green_wmax'][i]],
             [lmr['lmr_green_resol'][i], lmr['lmr_green_resol'][i]], alpha=.25, c='g')
    plt.plot([lmr['lmr_red_wmin'][i], lmr['lmr_red_wmax'][i]],
             [lmr['lmr_red_resol'][i], lmr['lmr_red_resol'][i]], alpha=.25, c='xkcd:orange red')
    plt.plot([lmr['lmr_j_wmin'][i], lmr['lmr_j_wmax'][i]],
             [lmr['lmr_j_resol'][i], lmr['lmr_j_resol'][i]], alpha=.25, c='xkcd:red')
    plt.plot([lmr['lmr_h_wmin'][i], lmr['lmr_h_wmax'][i]],
             [lmr['lmr_h_resol'][i], lmr['lmr_h_resol'][i]], alpha=.25, c='xkcd:burgundy')
plt.xlabel('Wavelength')
plt.ylabel('Spectral resolution')
plt.xlim([300, 2000])
plt.savefig('lmr_resolution.png')
#   Resolution vs wrange
plt.figure()
for i in range(len(lmr)):
    if lmr['lmr_blue_wmax'][i] - lmr['lmr_blue_wmin'][i] != 0 and lmr['lmr_blue_resol'][i] != 0:
        plt.scatter((lmr['lmr_blue_wmax'][i] - lmr['lmr_blue_wmin'][i])/190, lmr['lmr_blue_resol'][i], alpha=1, c='b')
    if lmr['lmr_green_wmax'][i] - lmr['lmr_green_wmin'][i] != 0 and lmr['lmr_green_resol'][i] != 0:
        plt.scatter((lmr['lmr_green_wmax'][i] - lmr['lmr_green_wmin'][i])/200, lmr['lmr_green_resol'][i], alpha=1, c='g')
    if lmr['lmr_red_wmax'][i] - lmr['lmr_red_wmin'][i] != 0 and lmr['lmr_red_resol'][i] != 0:
        plt.scatter((lmr['lmr_red_wmax'][i] - lmr['lmr_red_wmin'][i])/250, lmr['lmr_red_resol'][i], alpha=1, c='xkcd:orange red')
    if lmr['lmr_j_wmax'][i] - lmr['lmr_j_wmin'][i] != 0 and lmr['lmr_j_resol'][i] != 0:
        plt.scatter((lmr['lmr_j_wmax'][i] - lmr['lmr_j_wmin'][i])/300, lmr['lmr_j_resol'][i], alpha=1, c='xkcd:red')
    if lmr['lmr_h_wmax'][i] - lmr['lmr_h_wmin'][i] != 0 and lmr['lmr_h_resol'][i] != 0:
        plt.scatter((lmr['lmr_h_wmax'][i] - lmr['lmr_h_wmin'][i])/350, lmr['lmr_h_resol'][i], alpha=1, c='xkcd:burgundy')
plt.xlabel('Wavelength range')
plt.ylabel('Spectral resolution')
plt.savefig('lmr_resolution_vs_width.png')


# How many need R>3000 and >50% of spectral coverage?
sel_blue = np.where(((lmr['lmr_blue_wmax'] - lmr['lmr_blue_wmin'])/190 > .50) & (lmr['lmr_blue_resol'] > 3000))
sel_green = np.where(((lmr['lmr_green_wmax'] - lmr['lmr_green_wmin'])/200 > .50) & (lmr['lmr_green_resol'] > 3000))
sel_red = np.where(((lmr['lmr_red_wmax'] - lmr['lmr_red_wmin'])/250 > .50) & (lmr['lmr_red_resol'] > 3000))
sel_j = np.where(((lmr['lmr_j_wmax'] - lmr['lmr_j_wmin'])/300 > .50) & (lmr['lmr_j_resol'] > 3000))
sel_h = np.where(((lmr['lmr_h_wmax'] - lmr['lmr_h_wmin'])/350 > .50) & (lmr['lmr_h_resol'] > 3000))
print(f"Number of MR users with wide spectral coverage needs (blue): {len(sel_blue[0])}")
print(f"Number of MR users with wide spectral coverage needs (green): {len(sel_green[0])}")
print(f"Number of MR users with wide spectral coverage needs (red): {len(sel_red[0])}")
print(f"Number of MR users with wide spectral coverage needs (j): {len(sel_j[0])}")
print(f"Number of MR users with wide spectral coverage needs (h): {len(sel_h[0])} \n")
sel = np.concatenate((sel_blue[0], sel_green[0], sel_red[0], sel_j[0], sel_h[0]))
sel = np.unique(sel)
print(f"Number of MR users with wide spectral coverage needs: {len(sel)} \n")
# print(np.sort(lmr["id"][sel_blue[0]]))
# print(np.sort(lmr["id"][sel_green[0]]))
# print(np.sort(lmr["id"][sel_red[0]]))
# print(np.sort(lmr["id"][sel_j[0]]))
# print(np.sort(lmr["id"][sel_h[0]]))


# Cumulative plot for spectral coverage needs in LMR
cumul_blue = cumul_green = cumul_red = cumul_j = cumul_h = []
for i in range(101):
    cumul_blue = cumul_blue + [len(lmr[((lmr['lmr_blue_wmax'] - lmr['lmr_blue_wmin']) / 190 * 100 <= i)
                                       & (lmr['lmr_blue_wmax'] - lmr['lmr_blue_wmin'] != 0)])]
    cumul_green = cumul_green + [len(lmr[((lmr['lmr_green_wmax'] - lmr['lmr_green_wmin']) / 200 * 100 <= i)
                                         & (lmr['lmr_green_wmax'] - lmr['lmr_green_wmin'] != 0)])]
    cumul_red = cumul_red + [len(lmr[((lmr['lmr_red_wmax'] - lmr['lmr_red_wmin']) / 250 * 100 <= i)
                                     & (lmr['lmr_red_wmax'] - lmr['lmr_red_wmin'] != 0)])]
    cumul_j = cumul_j + [len(lmr[((lmr['lmr_j_wmax'] - lmr['lmr_j_wmin']) / 300 * 100 <= i)
                                 & (lmr['lmr_j_wmax'] - lmr['lmr_j_wmin'] != 0)])]
    cumul_h = cumul_h + [len(lmr[((lmr['lmr_h_wmax'] - lmr['lmr_h_wmin']) / 350 * 100 <= i)
                                 & (lmr['lmr_h_wmax'] - lmr['lmr_h_wmin'] != 0)])]
plt.figure()
plt.plot(range(101), np.array(cumul_blue) / len(lmr[(lmr['lmr_blue_wmax'] - lmr['lmr_blue_wmin'] != 0)]), c='b')
plt.plot(range(101), np.array(cumul_green) / len(lmr[(lmr['lmr_green_wmax'] - lmr['lmr_green_wmin'] != 0)]), c='g')
plt.plot(range(101), np.array(cumul_red) / len(lmr[(lmr['lmr_red_wmax'] - lmr['lmr_red_wmin'] != 0)]), c='xkcd:orange red')
plt.plot(range(101), np.array(cumul_j) / len(lmr[(lmr['lmr_j_wmax'] - lmr['lmr_j_wmin'] != 0)]), c='xkcd:red')
plt.plot(range(101), np.array(cumul_h) / len(lmr[(lmr['lmr_h_wmax'] - lmr['lmr_h_wmin'] != 0)]), c='xkcd:burgundy')
plt.xlabel('Spectral coverage')
plt.ylim([0, 1])
plt.ylabel('Number of valid participants')
plt.savefig('lmr_coverage_cumul.png')


# How many need LMR resolution >= 2000? 4000?
n_blue_resol = len(lmr[lmr['lmr_blue_resol'] > 0])
n_green_resol = len(lmr[lmr['lmr_green_resol'] > 0])
n_red_resol = len(lmr[lmr['lmr_red_resol'] > 0])
n_j_resol = len(lmr[lmr['lmr_j_resol'] > 0])
n_h_resol = len(lmr[lmr['lmr_h_resol'] > 0])

print(f"Number of valid answers for the LMR (blue): {n_blue_resol} ")
print(f"Number of valid answers for the LMR (green): {n_green_resol} ")
print(f"Number of valid answers for the LMR (red): {n_red_resol} ")
print(f"Number of valid answers for the LMR (J): {n_j_resol} ")
print(f"Number of valid answers for the LMR (H): {n_h_resol} \n")

n_r2000 = len(lmr[(lmr['lmr_blue_resol'] <= 2000) & (lmr['lmr_blue_resol'] > 0)])
print(f"Fraction of LMR users who need R<=2000 (blue): {n_r2000 / n_blue_resol * 100:4.1f} %")
n_r2000 = len(lmr[(lmr['lmr_green_resol'] <= 2000) & (lmr['lmr_green_resol'] > 0)])
print(f"Fraction of LMR users who need R<=2000 (green): {n_r2000 / n_green_resol * 100:4.1f} %")
n_r2000 = len(lmr[(lmr['lmr_red_resol'] <= 2000) & (lmr['lmr_red_resol'] > 0)])
print(f"Fraction of LMR users who need R<=2000 (red): {n_r2000 / n_red_resol * 100:4.1f} %")
n_r2000 = len(lmr[(lmr['lmr_j_resol'] <= 2000) & (lmr['lmr_j_resol'] > 0)])
print(f"Fraction of LMR users who need R<=2000 (J): {n_r2000 / n_j_resol * 100:4.1f} %")
n_r2000 = len(lmr[(lmr['lmr_h_resol'] <= 2000) & (lmr['lmr_h_resol'] > 0)])
print(f"Fraction of LMR users who need R<=2000 (H): {n_r2000 / n_h_resol * 100:4.1f} %\n")

for i in [2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000]:
    n_r2000 = len(lmr[lmr['lmr_blue_resol'] > i])
    print(f"Fraction of LMR users who need R>{i} (blue): {n_r2000 / n_blue_resol * 100:4.1f} %")
    n_r2000 = len(lmr[lmr['lmr_green_resol'] > i])
    print(f"Fraction of LMR users who need R>{i} (green): {n_r2000 / n_green_resol * 100:4.1f} %")
    n_r2000 = len(lmr[lmr['lmr_red_resol'] > i])
    print(f"Fraction of LMR users who need R>{i} (red): {n_r2000 / n_red_resol * 100:4.1f} %")
    n_r2000 = len(lmr[lmr['lmr_j_resol'] > i])
    print(f"Fraction of LMR users who need R>{i} (J): {n_r2000 / n_j_resol * 100:4.1f} %")
    n_r2000 = len(lmr[lmr['lmr_h_resol'] > i])
    print(f"Fraction of LMR users who need R>{i} (H): {n_r2000 / n_h_resol * 100:4.1f} %\n")


# Plot for resolution needs in LMR
cumul_blue = cumul_green = cumul_red = cumul_j = cumul_h = []
for i in range(100, 10100, 100):
    cumul_blue = cumul_blue + [len(lmr[(lmr['lmr_blue_resol'] <= i) & (lmr['lmr_blue_resol'] > 0)]) / n_blue_resol]
    cumul_green = cumul_green + [len(lmr[(lmr['lmr_green_resol'] <= i) & (lmr['lmr_green_resol'] > 0)]) / n_green_resol]
    cumul_red = cumul_red + [len(lmr[(lmr['lmr_red_resol'] <= i) & (lmr['lmr_red_resol'] > 0)]) / n_red_resol]
    cumul_j = cumul_j + [len(lmr[(lmr['lmr_j_resol'] <= i) & (lmr['lmr_j_resol'] > 0)]) / n_j_resol]
    cumul_h = cumul_h + [len(lmr[(lmr['lmr_h_resol'] <= i) & (lmr['lmr_h_resol'] > 0)]) / n_h_resol]
plt.figure()
plt.plot(range(100, 10100, 100), cumul_blue, c='b', label='Blue')
plt.plot(range(100, 10100, 100), cumul_green, c='g', label='Green')
plt.plot(range(100, 10100, 100), cumul_red, c='xkcd:orange red', label='Red')
plt.plot(range(100, 10100, 100), cumul_j, c='xkcd:red', label='J')
plt.plot(range(100, 10100, 100), cumul_h, c='xkcd:burgundy', label='H')
plt.legend()
plt.xlabel('Spectral resolution')
plt.ylabel('Fraction of valid participants')
plt.savefig('lmr_resolution_cumul.png')


# How many need H-band?
n_need_h = len(ans[ans['lmr_h_need'] == 'Yes'])
print(f"Fraction of users who need H band: {n_need_h / n_ans * 100:4.1f} %")
n_need_h = len(lmr[lmr['lmr_h_need'] == 'Yes'])
print(f"Fraction of LMR users who need H band: {n_need_h / n_need_lmr * 100:4.1f} %\n")


# How many need J-band?
n_need_j = len(ans[ans['lmr_j_need'] == 'Yes'])
print(f"Fraction of users who need J band: {n_need_j / n_ans * 100:4.1f} %")
n_need_j = len(lmr[lmr['lmr_j_need'] == 'Yes'])
print(f"Fraction of LMR users who need J band: {n_need_j / n_need_lmr * 100:4.1f} %\n")


# How many need what band?
n_noneed_ir = len(lmr[(lmr['lmr_h_need'] == 'No') & (lmr['lmr_j_need'] == 'No')])
n_noneed_vis = len(lmr[(lmr['lmr_blue_need'] == 'No') & (lmr['lmr_green_need'] == 'No') & (lmr['lmr_red_need'] == 'No')])
print(f"Fraction of LMR users who do *not* need IR bands: {n_noneed_ir / n_need_lmr * 100:4.1f} %")
print(f"Fraction of LMR users who do *not* need VIS bands: {n_noneed_vis / n_need_lmr * 100:4.1f} %\n")


# How many would be impacted without H-band?
n_h_no_imp = len(lmr[lmr['h_impact'] == 0])
n_h_insig_imp = len(lmr[lmr['h_impact'] == 1])
n_h_sign_imp = len(lmr[lmr['h_impact'] == 2])
n_h_comp_imp = len(lmr[lmr['h_impact'] == 3])
print(f"Fraction of LMR users who are significantly impacted without H band: {(n_h_sign_imp + n_h_comp_imp) / n_need_lmr * 100:4.1f} %")
print(f"Fraction of LMR users who are impacted without H band: {(n_h_insig_imp + n_h_sign_imp + n_h_comp_imp) / n_need_lmr * 100:4.1f} %")
