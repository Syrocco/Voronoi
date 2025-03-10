import numpy as np
import matplotlib.pyplot as plt
from parserL import Dump
import re
from scipy.spatial import Voronoi, voronoi_plot_2d
from cmap import Colormap
import os
import subprocess
from matplotlib.collections import PolyCollection

class Cell(Dump):
    def __init__(self, loc):
        self.loc = loc
        
        Dump.__init__(self, loc)
        self.extract_variables()
        self.L = np.sqrt(self.N)
        self.path = path = loc.split(".dump")[0]
        if loc[-1] == "o":
            path = loc
        else:
            path = self.path + ".thermo"
        with open(path, 'r') as file:
            first_line = file.readline().strip()
            column_names = first_line.split()
            try:
                data = np.loadtxt(path, skiprows = 1)
            except:
                data = np.loadtxt(path, skip_footer = 1, skiprows = 1)

            if data.ndim == 1:
                data = data.reshape(1, -1)

            for i, col_name in enumerate(column_names):
                setattr(self, col_name, data[:, i])
    def extract_variables(self):
        # Define the pattern to match variable names and values
        pattern = r'(\w+?)_(\d*\.?\d+)'

        # Find all matches
        matches = re.findall(pattern, self.loc)

        # Create a dictionary to store the variable names and values
        variables = {}

        # Iterate over matches and populate the dictionary
        for match in matches:
            variable_name = match[0]
            variable_value = float(match[1]) if '.' in match[1] else int(match[1])
            variables[variable_name] = variable_value

        # Assign values to variables with the same names
        for name, value in variables.items():
            setattr(self, name, value)
    def pbc(self,frame = -1, prop = None):
        if frame < 0:
            frame = self.nframes + frame
        
        self.jump_to_frame(frame)
        L = self.L
        dL = self.get_boxx()[0] - L
        
        size = 8 
        x = self.get_atompropf("x")
        y = self.get_atompropf("y")
        xadd=list(x[dL*y>L*(x-size)]+L)+list(x[y<size]+dL)+list(x[np.logical_and(dL*y>L*(x-size),y<size)]+L+dL)+list(x[dL*y<L*(x+size-L)]-L)+list(x[y>L-size]-dL)+list(x[np.logical_and(dL*y<L*(x+size-L),y>L-size)]-L-dL)+list(x[np.logical_and(dL*y>L*(x-size),y>L-size)]+L-dL)+list(x[np.logical_and(dL*y<L*(x+size-L),y<size)]-L+dL)
        yadd=list(y[dL*y>L*(x-size)])+list(y[y<size]+L)+list(y[np.logical_and(dL*y>L*(x-size),y<size)]+L)+list(y[dL*y<L*(x+size-L)])+list(y[y>L-size]-L)+list(y[np.logical_and(dL*y<L*(x+size-L),y>L-size)]-L)+list(y[np.logical_and(dL*y>L*(x-size),y>L-size)]-L)+list(y[np.logical_and(dL*y<L*(x+size-L),y<size)]+L)
        # xadd=list(x[x<size]+self.L)+list(x[y<size])+list(x[np.logical_and(x<size,y<size)]+self.L)+list(x[x>self.L-size]-self.L)+list(x[y>self.L-size])+list(x[np.logical_and(x>self.L-size,y>self.L-size)]-self.L)+list(x[np.logical_and(x<size,y>self.L-size)]+self.L)+list(x[np.logical_and(x>self.L-size,y<size)]-self.L)
        # yadd=list(y[x<size])+list(y[y<size]+self.L)+list(y[np.logical_and(x<size,y<size)]+self.L)+list(y[x>self.L-size])+list(y[y>self.L-size]-self.L)+list(y[np.logical_and(x>self.L-size,y>self.L-size)]-self.L)+list(y[np.logical_and(x<size,y>self.L-size)]-self.L)+list(y[np.logical_and(x>self.L-size,y<size)]+self.L)
        X=np.array(list(x)+xadd)
        Y=np.array(list(y)+yadd)
        if prop == None:
            return X, Y
        c= self.get_atompropf(prop)
        # cadd=list(c[x<size])+list(c[y<size])+list(c[np.logical_and(x<size,y<size)])+list(c[x>self.L-size])+list(c[y>self.L-size])+list(c[np.logical_and(x>self.L-size,y>self.L-size)])+list(c[np.logical_and(x<size,y>self.L-size)])+list(c[np.logical_and(x>self.L-size,y<size)])
        cadd=list(c[dL*y>L*(x-size)])+list(c[y<size])+list(c[np.logical_and(dL*y>L*(x-size),y<size)])+list(c[dL*y<L*(x+size-L)])+list(c[y>L-size])+list(c[np.logical_and(dL*y<L*(x+size-L),y>L-size)])+list(c[np.logical_and(dL*y>L*(x-size),y>L-size)])+list(c[np.logical_and(dL*y<L*(x+size-L),y<size)])
        C=np.array(list(c)+cadd)
        return X,Y,C
    def voronoi(self, frame = -1, prop = "area", cmap = "cmocean:delta", fig = None, ax = None):
        if frame < 0:
            frame = self.nframes + frame
        self.jump_to_frame(frame)
        X, Y, C = self.pbc(frame, prop)
        
        # Optimize colormap calculation
        C_norm = (C - np.min(C)) / (np.max(C) - np.min(C))
        colormap = Colormap(cmap)
        C_colors = colormap(C_norm)
        
        # Optimize point calculation
        points = np.vstack((X, Y)).T
        vor = Voronoi(points, qhull_options="Qc")
    
        
        # Pre-compute vertices for plotting
        regions = [vor.regions[vor.point_region[i]] for i in range(self.N)]
        valid_regions = [(i, region) for i, region in enumerate(regions) 
                         if region and -1 not in region]
        
        # Use collections for faster plotting
        polys = []
        colors = []
        
        for i, region in valid_regions:
            if i < len(C_colors):  # Ensure we have a color for this region
                polygon = vor.vertices[region]
                polys.append(polygon)
                colors.append(C_colors[i])
        
        if fig == None or ax == None:
            fig = plt.figure(figsize=(4*(1 + self.gamma), 4), layout='constrained')
            ax = plt.gca()
        plt.axis("square")
        plt.axis("off")
        
        pc = PolyCollection(polys, facecolors=colors, edgecolors='k', alpha=0.8)
        ax.add_collection(pc)
        
        # Plot points once
        plt.plot(X[:self.N], Y[:self.N], '.', c='k')
        plt.xlim(0, self.L*(1+self.gamma))
        plt.ylim(0, self.L)
        
    def batch(self,prop,cmap="cmasher:viola"):
        last = self.nframes
        os.makedirs(self.path+f'batch_{prop}', exist_ok = True)
        fig = plt.figure(figsize=(4*(1 + self.gamma), 4), layout='constrained')
        ax = plt.gca()
        for i in range(last):
            print(i, "/", last)
            self.voronoi(i,prop,cmap,fig = fig, ax = ax)
            plt.savefig(self.path+f'batch_{prop}/{i:04}.png')
            ax.clear()
        plt.close()
    def save(self, prop,frame_rate=10,cmap="cmasher:viola"):
        if not os.path.exists(self.path+f'batch_{prop}'):
            self.batch(prop,cmap=cmap)
        input_pattern = self.path+fr'batch_{prop}/%04d.png' 
        output_video = self.path+f'batch_{prop}.gif'

        ffmpeg_cmd = [
            "ffmpeg",
            "-y",
            "-framerate", str(frame_rate),
            "-i", input_pattern,
            output_video
        ]

        subprocess.run(ffmpeg_cmd, check=True)
        
        
        

a = Cell("/mnt/ssd/Documents/Voronoi/dump/N_4010Po_4.000000gamma_1.000000.dump")
a.save("area")