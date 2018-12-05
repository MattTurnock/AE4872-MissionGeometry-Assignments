from json_to_dict import constants
from Ass1.misc_utils import list_perms
import sympy as sp
from sympy import latex
import pandas as pd
pd.set_option('precision', 5)

from Ass2.repeating_orbit_utils import do_a_iterations
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
i_all = np.linspace(0,180,181)
#i_all = np.linspace(0,180,10)

# jks = [14, 43/3]
# i_all = [28, 28]
#def docalc():

#Function definition to do all calculations
def docalc_plot(i_all, js, k, approach='1', prnt=True, plot=True, plotshow=True, calc=True, extra='', interval=30, linewidth=1.5):
    if calc:
        input_data = list_perms(js,i_all)
        if prnt: print('data in: \n',input_data, '\n')
        out_data = []
        for combo in input_data:

            jk = combo[0]/k
            i = combo[1]*u.deg
            a = do_a_iterations(jk, i, approach=approach)
            h = a - constants["RE"]
            if 200*u.km < h < 1200*u.km:
                combo.append(h.value)
                out_data.append(combo)
        out_data = np.array(out_data)
        np.save('out_data%s' %extra, out_data)

    out_data = np.load('out_data%s.npy' %extra)
    if prnt: print('data out: \n',out_data,'\n')

    if plot:
        split_data = np.split(out_data, np.where(np.diff(out_data[:,0]))[0]+1)
        plt.figure()
        for data in split_data:
            xs = data[:,2]
            ys = data[:,1]
            plt.plot(xs, ys, linewidth=linewidth)
        legend=[]
        for j in js:
            string = "(%s, 3)" %j
            legend.append(string)
        plt.legend(legend, markerscale=5., bbox_to_anchor=(1.05, 1))
        plt.grid()
        plt.ylabel('Inclination, i [deg]')
        plt.xlabel('Orbital altitude, h [km]')
        plt.savefig('repeaters%s.pdf' %extra, bbox_inches="tight")
        if plotshow: plt.show()

    tabulated = out_data[~(out_data[:,1]%interval!=0.0), :]
    tabulated = np.around(tabulated, decimals=2)
    # tabulated = tabulated.astype('int')
    if prnt: print('tabulated data: \n', tabulated)
    np.savetxt('out_data%s.csv' %extra, tabulated)
    #tabulated = np.around(tabulated, decimals=2)
    # tabulated = pd.DataFrame(data=tabulated, columns=['j', 'i', 'a'])
    # tabulated = tabulated.astype({'j': 'int', 'i':'int', 'a':'f'})
    # tabulated = tabulated.to_latex(index=False)
    print(tabulated)

    np.savetxt("mydata.txt", tabulated, delimiter=' & ', fmt='%i %i %1.3f')#, newline=' \\\\\n')

#Do for approach 1 and approach 2
docalc_plot(i_all, js, k, approach='2', prnt=True, plot=True, plotshow=False, calc=False, extra='_2')
docalc_plot(i_all, js, k, approach='1', prnt=True, plot=True, plotshow=False, calc=False, extra='_1')


#Make a combined plot
out_data_1 = np.load('out_data_1.npy')
out_data_2 = np.load('out_data_2.npy')
split_data_1 = np.split(out_data_1, np.where(np.diff(out_data_1[:, 0]))[0] + 1)
split_data_2 = np.split(out_data_2, np.where(np.diff(out_data_2[:, 0]))[0] + 1)

linewidth=1.0
plot=True
plotshow=True
if plot:
    plt.figure()
    for data in split_data_1:
        xs = data[:, 2]
        ys = data[:, 1]
        plt.plot(xs, ys, linewidth=linewidth )#, s=1.0)
    for data in split_data_2:
        xs = data[:, 2]
        ys = data[:, 1]
        plt.plot(xs, ys, linestyle='dashed', linewidth=linewidth)#, s=1.0)
    legend = []
    for j in js:
        string = "(%s, 3)" % j
        legend.append(string)
    plt.legend(legend, markerscale=5.,bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    #plt.legend(loc="lower left", markerscale=2., scatterpoints=1, fontsize=10)
    plt.grid()
    plt.ylabel('Inclination, i [deg]')
    plt.xlabel('Orbital altitude, h [km]')
    plt.savefig('repeaters_both.pdf',bbox_inches="tight")
    if plotshow: plt.show()
