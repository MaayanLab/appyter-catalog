'''
slight edit of original code from ponnhide
https://github.com/ponnhide/pyCircos
'''

import math
import collections
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.path    as mpath
import matplotlib.patches as mpatches

class Gcircle(object):
    colors = ["#4E79A7","#F2BE2B","#E15759","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"]
    cmaps  = [plt.cm.Reds, plt.cm.Blues, plt.cm.Greens, plt.cm.Greys]  
    def __init__(self):
        self.locus_dict = collections.OrderedDict() 
        self.interspace = np.pi / 60 
        self.bottom      = 500 
        self.height      = 50 
        self.facecolor   = "#DDDDDD"
        self.edgecolor   = "#000000"
        self.linewidth   = 1.0 
        self.markersize  = 2.0 
        self.color_cycle = 0 
        self.cmap_cycle  = 0
    
    def get_dict(self):
        return self.locus_dict

    def add_locus(self, name, length, bottom=None, height=None,
                  facecolor=None, edgecolor=None, linewidth=None,
                  interspace=None):
        self.locus_dict[name]                 = {}
        self.locus_dict[name]["length"]       = length
        self.locus_dict[name]["features"] = [] 
        
        if bottom is None:
            self.locus_dict[name]["bottom"] = self.bottom
        else:
            self.locus_dict[name]["bottom"] = bottom
        
        if height is None:
            self.locus_dict[name]["height"] = self.height 
        else:
            self.locus_dict[name]["height"] = height

        if facecolor is None: 
            self.locus_dict[name]["facecolor"] = self.facecolor 
        else:
            self.locus_dict[name]["facecolor"] = facecolor

        if edgecolor is None:
            self.locus_dict[name]["edgecolor"] = self.edgecolor
        else:
            self.locus_dict[name]["edgecolor"] = edgecolor

        if interspace is None:
            self.locus_dict[name]["linewidth"] = self.linewidth
        else:
            self.locus_dict[name]["linewidth"] = linewidth

        if interspace is None:
            self.locus_dict[name]["interspace"] = self.interspace
        else:
            self.locus_dict[name]["interspace"] = interspace 


        sum_length       = sum(list(map(lambda x:  self.locus_dict[x]["length"], list(self.locus_dict.keys()))))
        sum_interspace   = sum(list(map(lambda x:  self.locus_dict[x]["interspace"], list(self.locus_dict.keys()))))
        self.theta_list  = np.linspace(0.0, 2 * np.pi - sum_interspace, sum_length, endpoint=True)
        s = 0
        sum_interspace = 0 
        for key in self.locus_dict.keys():
            self.locus_dict[key]["positions"] = sum_interspace + self.theta_list[s:s+self.locus_dict[key]["length"]+1]
            if s+self.locus_dict[key]["length"]+1 > len(self.theta_list):
                self.locus_dict[key]["positions"] = self.locus_dict[key]["positions"] + self.theta_list[:s+self.locus_dict[key]["length"] + 1- len(self.theta_list)]
            s = s + self.locus_dict[key]["length"]
            sum_interspace += self.locus_dict[key]["interspace"]
            
    def ax(self):
        return self.ax 
    
    
    def set_locus(self, figsize=(6, 6), lw=1): 
        self.figure = plt.figure(figsize=figsize)
        self.ax     = plt.subplot(111, polar=True)
        self.ax.set_theta_zero_location("N")
        self.ax.set_theta_direction(-1)
        self.ax.set_ylim(0,1000)
        self.ax.spines['polar'].set_visible(False)
        self.ax.xaxis.set_ticks([])
        self.ax.xaxis.set_ticklabels([])
        self.ax.yaxis.set_ticks([])
        self.ax.yaxis.set_ticklabels([])  
                
        pre_e = 0 
        for i, key in enumerate(self.locus_dict.keys()):
            pos       = self.locus_dict[key]["positions"][0] 
            width     = self.locus_dict[key]["positions"][-1] - self.locus_dict[key]["positions"][0]
            height    = self.locus_dict[key]["height"]
            bottom    = self.locus_dict[key]["bottom"]
            facecolor = self.locus_dict[key]["facecolor"]
            edgecolor = self.locus_dict[key]["edgecolor"]
            linewidth = self.locus_dict[key]["linewidth"]
            self.locus_dict[key]["bar"] = self.ax.bar([pos], [height], bottom=bottom,
                                                      width=width, facecolor=facecolor, 
                                                      linewidth=linewidth, edgecolor=edgecolor, align="edge")
    
    
    
    def chord_plot(self, start_list, end_list,  bottom=500, center=0, color="#1F77B4", alpha=0.5):
        #start_list and end_list is composed of "locus_id", "start", "end". 
        sstart = self.locus_dict[start_list[0]]["positions"][start_list[1]]
        send   = self.locus_dict[start_list[0]]["positions"][start_list[2]+1]   
        if len(start_list) == 4:
            stop = int(start_list[3]) 
        else:
            stop = bottom

        ostart = self.locus_dict[end_list[0]]["positions"][end_list[1]]
        oend   = self.locus_dict[end_list[0]]["positions"][end_list[2]+1] 
        if len(end_list) == 4:
            etop = int(end_list[3]) 
        else:
            etop = bottom

        z1 = stop - stop * math.cos(abs((send-sstart) * 0.5)) 
        z2 = etop - etop * math.cos(abs((oend-ostart) * 0.5)) 
        if sstart == ostart: 
            pass 
        else:
            Path      = mpath.Path
            path_data = [(Path.MOVETO,  (sstart, stop)),
                         (Path.CURVE3,  (sstart, center)),     
                         (Path.CURVE3,  (oend,   etop)),
                         (Path.CURVE3,  ((ostart+oend)*0.5, etop+z2)),
                         (Path.CURVE3,  (ostart, etop)),
                         (Path.CURVE3,  (ostart, center)),
                         (Path.CURVE3,  (send,   stop)),
                         (Path.CURVE3,  ((sstart+send)*0.5, stop+z1)),
                         (Path.CURVE3,  (sstart, stop)),
                        ]
            codes, verts = list(zip(*path_data)) 
            path  = mpath.Path(verts, codes)
            patch = mpatches.PathPatch(path, facecolor=color, alpha=alpha, linewidth=0, zorder=0)
            self.ax.add_patch(patch)
    
    def bar_plot(self, locus_name, data, bottom=None,
                 height=None, positions=None, facecolor=None,
                 linewidth=0.0, edgecolor="k"): 
        start = self.locus_dict[locus_name]["positions"][0] 
        end   = self.locus_dict[locus_name]["positions"][-1]
        if positions == None:
            positions = np.linspace(start, end, len(data), endpoint=False)
        else:
            pass 
        
        if bottom is None:
            bottom = self.bottom 
        
        if height is None:
            top = bottom + self.height
        else:
            top = bottom + height
        
        if facecolor is None:
            facecolor = Gcircle.colors[self.color_cycle % len(Gcircle.colors)] 
            self.color_cycle += 1
        
        width = positions[1] - positions[0] 
#         max_value = max(data) 
#         min_value = min(data)
#         data = np.array(data) - min_value
#         data = np.array(data * ((top - bottom) / (max_value - min_value)))
        self.ax.bar(positions, data, bottom=bottom,
                    facecolor=facecolor, width=width,
                    linewidth=linewidth, edgecolor=edgecolor, align="edge")
#         self.ax.set_ylim(bottom+ylim[0],bottom+ylim[1])