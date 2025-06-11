# paper plotting formats

# test

import matplotlib.pylab as plt
import numpy as np

################################################################
# settings for axes/ticks widths, etc. -> Called for any plotting
################################################################
# axis border widths (I tend to like bolder than the default)
plt.rc('axes', linewidth=0.75)
# tick widths (I like them the same width as the border)
plt.rc('xtick.major', width=0.75)
plt.rc('ytick.major', width=0.75)
plt.rc('xtick.minor', width=0.75)
plt.rc('ytick.minor', width=0.75)
# size of markers,no outline
plt.rc('lines', markersize=4, markeredgewidth=0.0)

# If you want a range of colors, this is a useful function to generate an equally
# spaced range
# Note, this no longer works for matplotlib past version 2.2
#import matplotlib.cm as cm
#cmap = cm.get_cmap("nipy_spectral")
#colors = cmap(a / b)
#num_colors = 9
# colors = np.zeros([num_colors,4])#the four is constant
# for i in np.arange(num_colors):
#    #c = cm.spectral(i/float(num_colors),1)
#    c=cmap(i/float(num_colors))
#    colors[i,:]=c
# then, call color=colors[x,:] in the plot routine

# These are possible marker styles
points = ['o', 'v', 's', 'p', '*', 'h', '^', 'D', '+', '>', 'H', 'd', 'x', '<']

# Legends
#leg=plt.legend(loc='upper left',fontsize=5,borderaxespad=1,ncol=5,frameon=False,fancybox=False,labelspacing=.5,columnspacing=.5)
# leg.set_title('Threshold',prop={'size':5})
#fr = leg.get_frame()
# fr.set_lw(0.2) #set frame width
# other mods: numpoints=3,markerscale=0.75,handlelength=5 (number of markers in legend, relative size of markers to plot, width of legend lines)

###########################################################
# Two subplotting options: call subplot or manually set axes
###########################################################

####Using Subplot####
# default settings for margin (can be tweaked accordingly
left = 0.2  # the left side of the subplots of the figure
right = 0.85    # the right side of the subplots of the figure
bottom = 0.15   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.2   # the amount of width reserved for blank space between subplots
hspace = 0.1   # the amount of height reserved for white space between subplots

####################
###Single Subplot###
####################
# call figure #3.5x3.5 is a good number for 300dpi
fig = plt.figure(num=1, figsize=(3.5, 3.5), dpi=300,
                 facecolor='w', edgecolor='k')
# apply settings for margin listed above to the figure
plt.subplots_adjust(left=left, bottom=bottom, right=right,
                    top=top, wspace=wspace, hspace=hspace)
# call subplot object
ax = plt.subplot(1, 1, 1)  # (num rows, num columns, subplot position)

# plot curves
x = np.arange(100)/5.0
y = np.sin(x)
z = np.cos(x)
# normal line plot
plt.plot(x, y, color='black', linewidth=1.5, label='Sin')
# dashed line plot with circle markers
plt.plot(x, z, marker='o', color='blue',
         linestyle='dashed', linewidth=1.5, label='Cos')

# set labels, labels sizes, ticks, ticks sizes
plt.xlabel('Time [s]', fontsize=9)
plt.ylabel('Power [arb]', fontsize=9)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)

# display different scale for the same curve
ax2 = ax.twinx()  # make second y axis (while duplicating/twin x)
y2 = 20*np.sin(x)
# don't actually display the plot, use linestyle='None'
plt.plot(x, y2, linestyle='None')
plt.ylabel('20*Power [arb]', fontsize=7)
plt.yticks(fontsize=7)


# saving the plot
# for the paper draft, its best to use png. When we actually submit a paper
# we'll need to save the plot as a .eps file instead.
savefile = 'Figure1_1subplot.png'
plt.savefig(savefile, dpi=300, facecolor='w', edgecolor='k')

##################################
###Two Subplots (same x-axis)#####
##################################
# call figure #3.5x6 is a good number for 300dpi
fig = plt.figure(num=2, figsize=(3.5, 6), dpi=300,
                 facecolor='w', edgecolor='k')
# apply settings for margin listed above to the figure
plt.subplots_adjust(left=left, bottom=bottom, right=right,
                    top=top, wspace=wspace, hspace=hspace)

# call first subplot object
ax = plt.subplot(2, 1, 1)  # (num rows, num columns, subplot position)
# apply letter label: coordinates in subplot object space (x,y) where (0,0) = bottom left
# also make sure transform uses correct subplot object)
plt.text(0.07, 0.92, '(a)', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes)
# plot curve
plt.plot(x, y, color='black', linewidth=1.5, label='Sin')
# set labels, labels sizes, ticks, ticks sizes (Note: only y label, but both x,y ticks)
plt.ylabel('Power [arb]', fontsize=9)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
plt.xlim(0, 20)
plt.ylim(-1.5, 1.5)
# suppress xtick values (modify subplot object ax)
ax.set_xticklabels([])
# Add Zero line (y position, xmin, xmax)
plt.hlines(0, 0, 20, color='red', linestyle='dotted', linewidth=0.5)
# (x,y) is in data coordinates when transform is not specified
plt.text(5, 0.05, 'y=0', fontsize=4, color='red', horizontalalignment='center')

# cal second subplot object
ax2 = plt.subplot(2, 1, 2)  # (num rows, num columns, subplot position)
# apply letter label: coordinates in subplot object space (x,y) where (0,0) = bottom left
# also make sure transform uses correct subplot object)
plt.text(0.07, 0.92, '(b)', horizontalalignment='center',
         verticalalignment='center', transform=ax2.transAxes)
# plot curve
plt.plot(x, z, color='black', linewidth=1.5, label='Cos')
# set labels, labels sizes, ticks, ticks sizes
plt.xlabel('Time [s]', fontsize=9)
plt.ylabel('Power [arb]', fontsize=9)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)
plt.xlim(0, 20)
plt.ylim(-1.5, 1.5)
# Add Horizontal lines (x position, ymin,ymax)
plt.vlines(5, -1.5, 1.5, color='gray', linestyle='dashed', linewidth=2.0)
plt.vlines(15, -0.75, 0.75, color='gray', linestyle='dashed', linewidth=2.0)

# saving the plot
# for the paper draft, its best to use png. When we actually submit a paper
# we'll need to save the plot as a .eps file instead.
savefile = 'Figure2_2subplots.png'
plt.savefig(savefile, dpi=300, facecolor='w', edgecolor='k')

############################################
###Three or more Subplots (same x-axis)#####
############################################
# I've found using a figsize of (3.5,8) for higher numbers of subplots
# is usually sufficient.

######################################################
###Two Subplots (Unequal sizes, differe t x-axes)#####
######################################################
# call figure #3.5x6 is a good number for 300dpi
fig = plt.figure(num=3, figsize=(3.5, 6), dpi=300,
                 facecolor='w', edgecolor='k')

# First subplot
# settings for margin for first subplot
left = 0.20  # the left side of the subplots of the figure
right = 0.95    # the right side of the subplots of the figure
bottom = 0.70  # the bottom of the subplots of the figure
top = 0.98      # the top of the subplots of the figure

# make axes object (rather than subplot object)
ax = plt.axes([left, bottom, right-left, top-bottom])
# apply letter label: coordinates in subplot object space (x,y) where (0,0) = bottom left
# also make sure transform uses correct subplot object)
plt.text(0.07, 0.92, '(a)', horizontalalignment='center',
         verticalalignment='center', transform=ax.transAxes)
# plot dotted curve with triangle markers
plt.plot(x, y, marker='^', color='red',
         linestyle='dotted', linewidth=1.5, label='Sin')
# set labels, labels sizes, ticks, ticks sizes (Note: manually set x ticks)
plt.xlabel('Time [s]', fontsize=9)
plt.ylabel('Power [arb]', fontsize=9)
# second argumennt needed if plotting loglog or semilog
plt.xticks(np.array([3, 5, 6, 10, 15, 16, 17, 18]), [
           3, 5, 6, 10, 15, 16, 17, 18], fontsize=9)
plt.yticks(fontsize=9)
plt.xlim(0, 20)
plt.ylim(-1.5, 1.5)

# Second subplot
# settings for margin for first subplot
left = 0.20  # the left side of the subplots of the figure
right = 0.95    # the right side of the subplots of the figure
bottom = 0.10  # the bottom of the subplots of the figure
top = 0.60      # the top of the subplots of the figure
# make axes object (rather than subplot object)
ax2 = plt.axes([left, bottom, right-left, top-bottom])
# apply letter label: coordinates in subplot object space (x,y) where (0,0) = bottom left
# also make sure transform uses correct subplot object)
plt.text(0.07, 0.92, '(b)', horizontalalignment='center',
         verticalalignment='center', transform=ax2.transAxes)
# plot curve with errorbars, dashed with circle markers)
x2 = np.arange(10)/10.
y2 = np.exp(x2)
yerr = 0.1*y2
plt.errorbar(x2, y2, yerr=yerr, marker='o', linestyle='dashed',
             capsize=2, capthick=0.5, linewidth=1.5, label='Exp')
# set labels, labels sizes, ticks, ticks sizes (Note: manually set x ticks)
plt.xlabel('Distance [m]', fontsize=9)
# multiline with Latex
plt.ylabel(
    r'Speed [m/s]' '\n' r'speed $= e^{d}$', fontsize=9, multialignment='center')
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)

# saving the plot
# for the paper draft, its best to use png. When we actually submit a paper
# we'll need to save the plot as a .eps file instead.
savefile = 'Figure3_2subplots_diffsize.png'
plt.savefig(savefile, dpi=300, facecolor='w', edgecolor='k')
