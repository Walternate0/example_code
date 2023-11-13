# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 17:15:02 2023

@author: spyridp
"""

import os
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dtime
from typing import Iterable
from swmm_api import read_inp_file, read_out_file, read_rpt_file, SwmmInput, swmm5_run, input_file

#%% func

def read_inps(directory: str) -> list[str]:
    return [name for name in os.listdir(directory)
            if (os.path.isfile(os.path.join(directory, name)) and name.rsplit('.')[1] == 'inp')]

def clear_section(inp, name: str, section: str) -> None:
    try:
        del inp[section]
    except:
        print(f"{name} has no {section} section")


#%% get .inp file names

folders = read_inps(os.getcwd()) # get list of model names

#%% rewrite .inp files to initialize with design storm 

TS_name = "TS_1"

for stud in folders:
    inp = read_inp_file(stud)
    clear_section(inp, stud, "TIMESERIES")
    clear_section(inp, stud, "RAINGAGES")
    #point to the storm .dat file
    inp.add_obj(input_file.sections.TimeseriesFile(TS_name, 'put path here'))
    #pay attention to timestep
    inp.add_obj(input_file.sections.RainGage(inp["SUBCATCHMENTS"].frame.iloc[0, 0], 'VOLUME', '0:05', 1, "TIMESERIES", timeseries = TS_name))
    inp["OPTIONS"].set_start(dtime.datetime(2000, 1, 1))
    inp["OPTIONS"].set_report_start(dtime.datetime(2000, 1, 1))
    inp["OPTIONS"].set_simulation_duration(dtime.timedelta(minutes = 240))
    inp["OPTIONS"].set_report_step(dtime.datetime(2000, 1, 1, 0, 5).time())
    inp.write_file(f"{stud}") # write modified .inp file
    
#%% run swmm for design storm

os.chdir('put path here')

start_t = time.time()
for stud in folders:
    swmm5_run(f"{stud}", swmm_lib_path=r'C:/Program Files/EPA SWMM 5.2.2 (64-bit)/runswmm.exe')
end_t = time.time()
print((end_t - start_t))

#%% simple map plot

poly = inp.POLYGONS.get_dataframe()
coords = inp.COORDINATES.get_dataframe()
conds = inp.CONDUITS.get_dataframe()
diams = inp.XSECTIONS.get_dataframe().loc[:, "height"]
# you can also use this for a color scheme on the subcatchments
# see how to use the same colorbar for 2 scales in the next cell
imper = inp.SUBCATCHMENTS.get_dataframe().loc[:, "imperviousness"] 

norm_min = 0.5
norm_max = 1

fig, ax = plt.subplots()

c_norm = mpl.colors.Normalize(vmin = norm_min, vmax = norm_max)
c_map = mpl.cm.RdBu

s_map = mpl.cm.ScalarMappable(cmap = c_map, norm = c_norm)
s_map.set_array([])

for j, sub in enumerate(poly.index):
    x, y = zip(*poly.loc[sub].iloc[0])
    ax.fill(x, y,
            color = "slategrey",
            alpha = 0.15,
            zorder = 1)
    ax.plot(x, y, 
            color = "black",
            lw = 0.15,
            alpha = 0.15,
            zorder = 1)

for cond in conds.index:
    from_node = conds.loc[cond, "from_node"]
    to_node = conds.loc[cond, "to_node"]
    ax.plot([coords.loc[from_node, "x"], coords.loc[to_node, "x"]], [coords.loc[from_node, "y"], coords.loc[to_node, "y"]],
            c = s_map.to_rgba(diams[cond]),
            lw = 3,
            zorder = 2)

fig.colorbar(s_map, ax = ax)

#%% function for more advanced plot

# this function needs to be passed a class instance with a specific structure in order to work, so you can't really use it
# but you can go through the code if you want to see one way to set up a more complex plot
# it's not the most elegant code ever, so if it confuses more than it helps, just delete it :)

def map_figs_sub(self, designs: tuple[tuple, tuple, tuple, tuple], norm_min = [-450, -4500], norm_max = [450, 4500]):
    
    poly = self.inp.POLYGONS.get_dataframe()
    coords = self.inp.COORDINATES.get_dataframe()
    conds = self.inp.CONDUITS.get_dataframe()
    
    
    # Create a figure and a grid of subplots
    fig_width = 11
    gap1 = 0.025*fig_width
    gap2 = 0.039*fig_width
    cbar_height = 0.02*fig_width
    fig_height = fig_width + 2*gap2 - gap1 + cbar_height
    
    fig, axs = plt.subplots(3, 2, figsize=(fig_width, fig_height))

    c_norm = mpl.colors.Normalize(vmin = norm_min[0], vmax = norm_max[0])
    c_map = mpl.cm.RdBu

    s_map = mpl.cm.ScalarMappable(cmap = c_map, norm = c_norm)
    s_map.set_array([])


    # ticks
    tick_cleanup = {"left" : False, "bottom" : False, "labelleft" : False, "labelbottom" : False}
    axs[0, 0].tick_params(**tick_cleanup)
    axs[0, 1].tick_params(**tick_cleanup)
    axs[1, 0].tick_params(**tick_cleanup)
    axs[1, 1].tick_params(**tick_cleanup)

    # Customize each subplot as needed
    axs[0, 0].set_title("Subplot 1")
    axs[0, 1].set_title("Subplot 2")
    axs[1, 0].set_title("Subplot 3")
    axs[1, 1].set_title("Subplot 4")
    # axs[2, 0].set_title("Colorbar Subplot")

    cbar = fig.colorbar(s_map, cax = axs[2, 0], label = "Difference in diameter (mm)", shrink = 9, orientation = "horizontal")
    cax2 = axs[2, 0].twiny()
    
    # Define the size of each subplot
    width = (fig_width - 3*gap1)/2/fig_width
    height = width*fig_width/fig_height
    gap_width = gap1/fig_width
    gap_height = gap1/fig_height
    gap2_height = gap2/fig_height
    cbar_height = cbar_height/fig_width
    
    axs[0, 0].set_position([gap_width, 2*gap2_height+gap_height+cbar_height+height,
                            width, height])  # [left, bottom, width, height]
    axs[0, 1].set_position([2*gap_width+width, 2*gap2_height+gap_height+cbar_height+height,
                            width, height])
    axs[1, 0].set_position([gap_width, 2*gap2_height+cbar_height,
                            width, height])
    axs[1, 1].set_position([2*gap_width+width, 2*gap2_height+cbar_height,
                            width, height])
    axs[2, 0].set_position([gap_width, gap2_height,
                            1 - 2*gap_width, cbar_height-0.01])
    axs[2, 1].set_visible(False)
    
    subplots_order = ((0, 0), (0, 1), (1, 0), (1, 1)) #tuple to be called in the plotting loop below
    
    # axs[0, 1].set_frame_on(False)
    
    cax2.set_frame_on(False)
    cax2.xaxis.tick_top()
    cax2.xaxis.label_position = "top"
    cax2.set_xlim(norm_min[1], norm_max[1])
    cax2.set_xlabel("Difference in storage volume ($m^3$)")
    
    for i, pair in enumerate(designs):
        axs[subplots_order[i]].set_title(f"{pair[0]} - {pair[1]}")        
        for j, sub in enumerate(poly.index):
            x, y = zip(*poly.loc[sub].iloc[0])
            axs[subplots_order[i]].fill(x, y,
                                        color = "slategrey",
                                        alpha = 0.15,
                                        zorder = 1)
            axs[subplots_order[i]].plot(x, y, 
                                        color = "black",
                                        lw = 0.15,
                                        alpha = 0.15,
                                        zorder = 1)
        # axs[subplots_order[i]].set(xlim = (0, 10), ylim = (0, 10))
        
        d_storage = (self.storage[pair[0]] - self.storage[pair[1]])*norm_max[0]/norm_max[1]
        
        axs[subplots_order[i]].scatter(coords.loc["T1-002", "x"], coords.loc["T1-002", "y"],
                                       s = 250,
                                       color = s_map.to_rgba(d_storage),
                                       zorder = 3)
        
        for cond in self.conduits_diam.index:
            d_diam = (self.conduits_diam.loc[cond, pair[0]] - self.conduits_diam.loc[cond, pair[1]])*1000
            # print(d_diam)
            from_node = conds.loc[cond, "from_node"]
            to_node = conds.loc[cond, "to_node"]
            axs[subplots_order[i]].plot([coords.loc[from_node, "x"], coords.loc[to_node, "x"]], [coords.loc[from_node, "y"], coords.loc[to_node, "y"]], 
                                        color = s_map.to_rgba(d_diam),
                                        lw = 3,
                                        zorder = 2)
            