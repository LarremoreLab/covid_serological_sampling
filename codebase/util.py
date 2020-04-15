#!/usr/bin/env python

__author__ = "Sam Way"
__copyright__ = "Copyright 2017, The Clauset Lab"
__license__ = "BSD"
__maintainer__ = "Sam Way"
__email__ = "samfway@gmail.com"
__status__ = "Development"


import numpy as np
from matplotlib import rcParams


# Constants
SINGLE_FIG_SIZE = (6,4)
BAR_WIDTH = 0.6
TICK_SIZE = 15
XLABEL_PAD = 10
LABEL_SIZE = 14
TITLE_SIZE = 16
LEGEND_SIZE = 12
LINE_WIDTH = 2
LIGHT_COLOR = '0.8'
LIGHT_COLOR_V = np.array([float(LIGHT_COLOR) for i in range(3)])
DARK_COLOR = '0.4'
DARK_COLOR_V = np.array([float(DARK_COLOR) for i in range(3)])
ALMOST_BLACK = '0.125'
ALMOST_BLACK_V = np.array([float(ALMOST_BLACK) for i in range(3)])
ACCENT_COLOR_1 = np.array([255., 145., 48.]) / 255.

# Configuration
#rcParams['text.usetex'] = True #Let TeX do the typsetting
rcParams['pdf.use14corefonts'] = True
rcParams['ps.useafm'] = True
#rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}', r'\sansmath'] #Force sans-serif math mode (for axes labels)
rcParams['font.family'] = 'sans-serif' # ... for regular text
rcParams['font.sans-serif'] = ['Helvetica','Helvetica Neue', 'HelveticaNeue'] #, Avant Garde, Computer Modern Sans serif' # Choose a nice font here
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['text.color'] = ALMOST_BLACK
rcParams['axes.unicode_minus'] = False

rcParams['xtick.major.pad'] = '8'
rcParams['axes.edgecolor']  = ALMOST_BLACK
rcParams['axes.labelcolor'] = ALMOST_BLACK
rcParams['lines.color']     = ALMOST_BLACK
rcParams['xtick.color']     = ALMOST_BLACK
rcParams['ytick.color']     = ALMOST_BLACK
rcParams['text.color']      = ALMOST_BLACK
rcParams['lines.solid_capstyle'] = 'butt'

# Imports
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D


def single_fig(figsize=SINGLE_FIG_SIZE):
    return plt.subplots(1,1,figsize=figsize)


def color_bp(bp, color):
    """ Helper function for making prettier boxplots """
    c = np.array(color) # * 0.5
    c = tuple(c)

    for x in bp['boxes']:
        plt.setp(x, color=c)
        x.set_facecolor(color)
    for x in bp['medians']:
        plt.setp(x, color='w')
    for x in bp['whiskers']:
        plt.setp(x, color=c)
    for x in bp['fliers']:
        plt.setp(x, color=c)
    for x in bp['caps']:
        plt.setp(x, color=c)


def adjust_spines(ax, spines):
    """ From http://matplotlib.org/examples/pylab_examples/spine_placement_demo.html """
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])


def hide_right_top_axis(ax):
    """ Remove the top and right axis """
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


def finalize(ax, fontsize=LABEL_SIZE, labelpad=7,ignoreLegend=False):
    """ Apply final adjustments """ 
    ax.tick_params(direction='out')
    hide_right_top_axis(ax)
    ax.yaxis.label.set_size(fontsize)
    ax.xaxis.label.set_size(fontsize)
    if ignoreLegend==False:
        ax.legend(frameon=False)
    ax.tick_params(axis='both', which='major', labelsize=fontsize, pad=labelpad)


def lineswap_axis(fig, ax, zorder=-1000, lw=1, alpha=0.2, skip_zero=False):
    """ Replace y-axis ticks with horizontal lines running through the background.
        Sometimes this looks really cool. Worth having in the bag 'o tricks.
    """
    fig.canvas.draw()  # Populate the tick vals/labels. Required for get_[xy]ticklabel calls.
    
    ylabels = [str(t.get_text()) for t in ax.get_yticklabels()]
    yticks = [t for t in ax.get_yticks()]
    xlabels = [str(t.get_text()) for t in ax.get_xticklabels()]
    xticks = [t for t in ax.get_xticks()]
        
    x_draw = [tick for label, tick in zip(ylabels, yticks) if label != '']  # Which ones are real, though?
    y_draw = [tick for label, tick in zip(ylabels, yticks) if label != '']
    
    xmin = x_draw[0]
    xmax = x_draw[-1]
    
    # Draw all the lines
    for val in y_draw:
        if val == 0 and skip_zero:
            continue  # Don't draw over the bottom axis
        ax.plot([xmin, xmax], [val, val], color=ALMOST_BLACK, zorder=zorder, lw=lw, alpha=alpha)
    
    ax.spines["left"].set_visible(False)  # Remove the spine
    ax.tick_params(axis=u'y', which=u'both',length=0)  # Erase ticks by setting length=0
    ax.set_xlim(xmin, xmax)  # Retain original xlims

