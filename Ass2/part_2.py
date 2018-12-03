from json_to_dict import constants
from Ass1.misc_utils import list_perms
from Ass2.repeating_orbit_utils import do_a_iterations, get_H0
import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt
from astropy.table import Table
pi = np.pi

##########################################################################################
#Approach 1 (Omega effect only)




##########################################################################################
#Approach 2 (both)

k=3

js = np.arange(39,49)
jks = js/k
i_all = np.linspace(0,180,180)

# jks = [14, 43/3]
# i_all = [28, 28]
calc=False
if calc:
    input_data = list_perms(js,i_all)
    print(input_data)
    out_data = []
    for combo in input_data:

        jk = combo[0]/k
        i = combo[1]*u.deg
        a = do_a_iterations(jk, i)
        h = a - constants["RE"]
        if 200*u.km < h < 1200*u.km:
            combo.append(h.value)
            out_data.append(combo)
    out_data = np.array(out_data)
    np.save('out_data', out_data)

out_data = np.load('out_data.npy')
print(out_data)

split_data = np.split(out_data, np.where(np.diff(out_data[:,0]))[0]+1)

plt.figure()
for data in split_data:
    xs = data[:,2]
    ys = data[:,1]
    plt.scatter(xs,ys, s=1.0)
legend=[]
for j in js:
    string = "(%s, 3)" %j
    legend.append(string)
plt.legend(legend)
plt.grid()
plt.ylabel('Inclination, i [deg]')
plt.xlabel('Orbital altitude, h [km]')
plt.savefig('repeaters.pdf')
plt.show()


