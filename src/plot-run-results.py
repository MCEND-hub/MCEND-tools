import os, sys, re
import pandas as pd
import matplotlib.pyplot as plt
from numpy import *
from matplotlib.gridspec import GridSpec
from scipy.special import comb, factorial
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as manimation
from matplotlib.animation import FuncAnimation

# plt.style.use(r'~/Dropbox/Grad-School-Things/solardarkmpl2.py')
fig = plt.figure(figsize=(6*2,4*2.5))

def plot_mcend_results(n=0):
    plot_acf = True
    lb_list = array(
        ['x', 'y', 'z', 'R', 'T_{{e}}', 'T_{{n}}',
         'V_{{ee}}', 'V_{{nn}}', 'V_{{en}}',
         'H_{{T}}', 'V_{\mathrm{BO}}', '\sigma'
        ])

    cmapcycler = ['viridis', 'plasma', 'cividis', 'Purples', 'Greys', 'Reds']

    gs = GridSpec(8, 6)
    ax1 =  plt.subplot(gs[0:2, 0:2])
    ax2 =  plt.subplot(gs[2:4, 0:2])
    ax3 =  plt.subplot(gs[4:6, 0:2])
    ax4 =  plt.subplot(gs[6:,  0:2])
    ax5 =  plt.subplot(gs[0:2, 2:4])
    ax6 =  plt.subplot(gs[0:2, 4:6])
    ax7 =  plt.subplot(gs[2:4, 2:4])
    ax8 =  plt.subplot(gs[2:4, 4:6])
    ax9 =  plt.subplot(gs[4:6, 2:4])
    ax10 = plt.subplot(gs[4:6, 4:6])
    ax11 = plt.subplot(gs[6:,  2:4])
    ax12 = plt.subplot(gs[6:,  4: ])

    ev = ['x', 'y', 'z',
          'R', 'Te', 'Tn',
          'Vee', 'Vnn', 'Ven',
          'Htot', 'Vbo']

    expec = pd.read_csv('expec.t', delimiter='\s+')
    acf_re = expec.loc[:, 'Re(Acf)'].values
    acf_im = expec.loc[:, 'Im(Acf)'].values
#     acf_abs = abs(complex(acf_re,acf_im))
    acf_abs = abs(acf_re + 1j*acf_im)

    for i, ax in enumerate(fig.axes[:-1]):
        ax.clear()
        if ev[i] == 'z':
            ax.plot(expec.time[:], expec.loc[:,ev[i]], ls='-', label=r'$\langle {:}\rangle$'.format(lb_list[i]))
            ax.set_xlim(-0.05)

        elif ev[i] == 'Vbo':
#             ax.plot(expec.R[:], expec.Ven[:], label=ev[i], alpha=0.1)
            V_BO = expec.Te[:] + expec.Ven[:] + expec.Vee[:] + expec.Vnn[:]
            ax.plot(expec.time[:], expec.R[:], label=r'$\langle {:}\rangle $'.format(lb_list[i]), alpha=0.01)

            ax.scatter(expec.time[:], expec.R[:],
                     c=V_BO.values, cmap=cmapcycler[0], label=None,
                     s=1, zorder=3, alpha=0.8)
#             ax.legend(#title=r'$\langle {:}\rangle $'.format(lb_list[i]),
#                   loc='best',
#                   handletextpad=0.2, handlelength=1,
#                   frameon=False)

            ax.set_xlabel('$t$ (fs)')
            ax.set_ylabel('$R$ ($a_0$)')

        else:
#             ax.plot(expec.time[n:], expec.loc[n:,ev[i]], ls='-')
            ax.plot(expec.time[:], expec.loc[:,ev[i]], ls='-',label=r'$\langle {:}\rangle $'.format(lb_list[i]),)
            ax.set_xlabel('$t$ (fs)')

        ax.tick_params(labelsize=9)
#         ax.legend(title=r'$\langle {:}\rangle $'.format(lb_list[i]), fontsize=8, handletextpad=0.2,
#                             handlelength=1, ncol=1, frameon=False)
        ax.legend(fontsize=10, handletextpad=0.2, handlelength=1, ncol=1, frameon=False)

        if not ev[i] == 'Vbo':
            ax.set_xlim(-0.05)

    if not plot_acf:
        ax12.plot(expec.R[:], expec.Ven[:], alpha=0.5)
        ax12.set_xlabel('$R$ ($a_0$)')

    if plot_acf:
        ax12.plot(expec.time[:],  acf_re[:], label=r'$\Re$',    alpha=0.5)
        ax12.plot(expec.time[:],  acf_im[:], label=r'$\Im$',    alpha=0.5)
        ax12.plot(expec.time[:], acf_abs[:], label=r'$\sigma$', alpha=0.8)

        ax12.set_xlabel('$t$ (fs)')
        ax12.tick_params(labelsize=9)
        ax12.set_xlim(-0.05)
        ax12.legend(fontsize=10, handletextpad=0.2, handlelength=1, ncol=3, frameon=False)

    plt.tight_layout()

ani = manimation.FuncAnimation(fig, plot_mcend_results, interval=1000)
# ani = manimation.FuncAnimation(fig, plot_mcend_results, interval=100)
plt.show()
