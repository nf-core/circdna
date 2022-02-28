import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import pi
# Say, "the default sans-serif font is COMIC SANS"
plt.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
plt.rcParams['font.family'] = "sans-serif"


#DERIVED FROM
#https://python-graph-gallery.com/391-radar-chart-with-several-individuals/
def make_classification_radar(categories,dvalueList,fname,snames):
    colors = ['cornflowerblue','orange','green','purple']
    n = len(categories)
    angles =  [i / float(n) * 2 * pi for i in range(n)]
    angles += angles[:1]

    ax = plt.subplot(111, polar=True)
    ax.set_theta_direction(-1)
    ax.set_theta_offset(120.1/360.0*2*pi)
    plt.xticks(angles[:-1], categories, size = 10)
    ax.tick_params(axis='x', pad=18)

    for ind,dvalues in enumerate(dvalueList):
        dvalues += dvalues[:1]
        ax.plot(angles,dvalues,linewidth=2,label=snames[ind],alpha=1)
        ax.fill(angles,dvalues,colors[ind],alpha=0.15)

    ax.set_rlabel_position(10)
    plt.yticks([0,0.25,0.5,0.75,1],["","0.25","0.50","0.75",""],color="grey",size=7)
    plt.ylim(-0.025,1)

    plt.legend(loc='upper right', bbox_to_anchor = (0.09,0.01),fontsize=7)
    plt.savefig(fname + "_radar.png",dpi=300)
    plt.savefig(fname + "_radar.pdf",format='pdf')
    plt.close()

