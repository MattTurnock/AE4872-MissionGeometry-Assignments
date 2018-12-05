from astropy import units as u
from Ass2.repeating_orbit_utils import do_a_iterations
from json_to_dict import constants
from Ass1.kep_orbit_utils import get_orbit_period
import numpy as np

inc = 28
input_data = [[14, inc],
            [43/3, inc],
            [29/2, inc],
            [59/4, inc],
            [74/5, inc],
            [15, inc]]

k = 1
calc=True
prnt=False
approach='2'
extra = '_2_veri'
if calc:
    #input_data = list_perms(js,i_all)
    if prnt: print('data in: \n',input_data, '\n')
    out_data = []
    for combo in input_data:

        jk = combo[0]/k
        i = combo[1]*u.deg
        a = do_a_iterations(jk, i, approach=approach)
        h = a - constants["RE"]
        T = get_orbit_period(a, constants["muEarth"], astropy_units=u.min)
        combo=[h.value, i.value, T.value, jk, 1.0]
        out_data.append(combo)
    out_data = np.array(out_data)
    np.save('out_data%s' %extra, out_data)

out_data = np.load('out_data%s.npy' %extra)
if prnt: print('data out: \n',out_data,'\n')

tabulated = out_data
#print and save data to put in report
print('tabulated data: \n', tabulated)
np.savetxt('out_data%s.txt' %extra, tabulated, delimiter=' & ', newline=' \\\\\n\hline\n', fmt='%1.2f & %i & %1.2f & %1.2f & %i')

inc = 28*u.deg
input_data = [[14, inc],
            [43/3, inc],
            [29/2, inc],
            [59/4, inc],
            [74/5, inc],
            [15, inc]]
for pair in input_data:
    j_k, i = pair[0:2]
    a, iters = do_a_iterations(j_k, i, da_precision=100, return_its=True, approach='2')
    T = round(get_orbit_period(a, constants["muEarth"], astropy_units=u.min).value, 2)
    #print(T)
    print(round((a - constants["RE"]).value, 2), iters)

original_heights = np.array([817.14, 701.34, 645.06, 562.55, 546.31, 482.25])
diff_heights = tabulated[:,0] - original_heights
print('differences [km]: \n', abs(np.around(diff_heights, 2)))




