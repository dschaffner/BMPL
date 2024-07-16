# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 18:21:25 2022

@author: dschaffner
"""

from matplotlib import pyplot as plt

plt.rcParams['backend'] = 'TkAgg'
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True

# Function to print mouse click event coordinates
def onclick(event):
    if event.key=='shift':
        print([event.xdata, event.ydata])

# Create a figure and a set of subplots
fig, ax = plt.subplots()

# Plot a line in the range of 10
ax.plot(range(10))

# Bind the button_press_event with the onclick() method
fig.canvas.mpl_connect('button_press_event', onclick)

# Display the plot
plt.show()