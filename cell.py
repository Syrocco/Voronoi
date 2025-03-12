import numpy as np
import matplotlib.pyplot as plt
from parserL import Dump
import re
from scipy.spatial import Voronoi, voronoi_plot_2d
from cmap import Colormap
import os
import subprocess
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1 import make_axes_locatable  # Add this import
from glob import glob
import imageio
cbar = None

def cellArray(loc, dump = False):
    return [Cell(data, dump) for data in sorted(glob(loc))]

class Cell(Dump):
    def __init__(self, loc, dump = True):
        self.loc = loc
        if dump:
            Dump.__init__(self, loc)
        self.extract_variables()
        self.L = np.sqrt(self.N)
        try:
            self.path = path = loc.split(".dump")[0]
        except:
            self.path = path = loc.split(".thermo")[0]
        if loc[-1] == "o":
            path = loc
        else:
            path = self.path + ".thermo"
        with open(path, 'r') as file:
            try:
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
            except:
                print("could not recover .thermo")
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
        dL = self.get_boxxy()
        if dL >= 0:
            dL = self.get_boxx()[0] - L
        
        size = 3 
        x = self.get_atompropf("x")
        y = self.get_atompropf("y")
        # Add points around x and y in a sheared box with periodic boundary conditions
        mask1 = dL * y > L * (x - size)
        mask2 = y < size
        mask3 = dL * y < L * (x + size - L)
        mask4 = y > L - size

        xadd = (
            list(x[mask1] + L) +
            list(x[mask2] + dL) +
            list(x[mask3] - L) +
            list(x[mask4] - dL) +
            list(x[mask1 & mask2] + L + dL) +
            list(x[mask3 & mask4] - L - dL) +
            list(x[mask1 & mask4] + L - dL) +
            list(x[mask3 & mask2] - L + dL)
        )
        yadd = (
            list(y[mask1]) +
            list(y[mask2] + L) +
            list(y[mask3]) +
            list(y[mask4] - L) +
            list(y[mask1 & mask2] + L) +
            list(y[mask3 & mask4] - L) +
            list(y[mask1 & mask4] - L) +
            list(y[mask3 & mask2] + L)
        )
        X=np.array(list(x)+xadd)
        Y=np.array(list(y)+yadd)
        
        if prop == None:
            return X, Y
        
        
        variables = re.findall(r'\b\w+\b', prop)
        local_vars = {var: self.get_atompropf(var) for var in variables}
    
        # Evaluate the property expression
        try:
            c = eval(prop, {"__builtins__": None, "np": np}, local_vars)
        except Exception as e:
            raise ValueError(f"Error evaluating property expression '{prop}': {e}")
        
        # cadd=list(c[x<size])+list(c[y<size])+list(c[np.logical_and(x<size,y<size)])+list(c[x>self.L-size])+list(c[y>self.L-size])+list(c[np.logical_and(x>self.L-size,y>self.L-size)])+list(c[np.logical_and(x<size,y>self.L-size)])+list(c[np.logical_and(x>self.L-size,y<size)])
        cadd=list(c[dL*y>L*(x-size)])+list(c[y<size])+list(c[np.logical_and(dL*y>L*(x-size),y<size)])+list(c[dL*y<L*(x+size-L)])+list(c[y>L-size])+list(c[np.logical_and(dL*y<L*(x+size-L),y>L-size)])+list(c[np.logical_and(dL*y>L*(x-size),y>L-size)])+list(c[np.logical_and(dL*y<L*(x+size-L),y<size)])
        C=np.array(list(c)+cadd)
        return X,Y,C
    def voronoi(self, frame = -1, prop = "area", cmap = "cmocean:delta", fig = None, ax = None, shear = True, bar = True):
        global cbar
        
        if frame < 0:
            frame = self.nframes + frame
        self.jump_to_frame(frame)
        X, Y, C = self.pbc(frame, prop)
        
        # Optimize colormap calculation
        C_norm = (C - np.min(C)) / (np.max(C) - np.min(C))
        try:
            colormap = Colormap(cmap).to_mpl()
            C_colors = colormap(C_norm)
        except:
            import matplotlib.cm as cm
            C_colors = cm.get_cmap(cmap)
        
        # Optimize point calculation
        points = np.vstack((X, Y)).T
        vor = Voronoi(points, qhull_options="Qc")
    
        
        # Pre-compute vertices for plotting
        if shear:
            lenght = self.N
        else:
            lenght = len(X)
            
        regions = [vor.regions[vor.point_region[i]] for i in range(lenght)]
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
        
        if shear:
            l = 1 + 2*self.gamma
        else:
            l = 1
        if fig == None or ax == None:
            fig = plt.figure(figsize=(4*l, 5.2), layout='constrained')
            ax = plt.gca()
        ax.axis("square")
        ax.axis("off")
        
        pc = PolyCollection(polys, facecolors=colors, edgecolors='k', alpha=0.8)
        ax.add_collection(pc)
        
        # Plot points once
        ax.plot(X[:self.N], Y[:self.N], '.', c='k')
        ax.set_xlim(-self.L*self.gamma, self.L*(1 + self.gamma))
        ax.set_ylim(0, self.L)

        if bar:
            if cbar is None:
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("bottom", size="5%", pad=0.05)
                cax.autoscale_view()
                cax.autoscale(False)
                norm = Normalize(vmin=np.min(C), vmax=np.max(C))
                sm = ScalarMappable(cmap=colormap, norm=norm)
                sm.set_array([])
                cbar = plt.colorbar(sm, cax=cax, orientation  = "horizontal")
                cbar.set_label(prop)
                fig.set_size_inches(4*l, 5.2)
            else:
                norm = Normalize(vmin=np.min(C), vmax=np.max(C))
                sm = ScalarMappable(cmap=colormap, norm=norm)
                sm.set_array([])
                
                cbar.update_normal(sm)
        
    def batch(self, prop, start, end, step, cmap="cmasher:viola", shear = True, bar = True):
        
        os.makedirs(self.path + f'batch_{prop}', exist_ok = True)
        if shear:
            l = 1 + 2*self.gamma
        else:
            l = 1
        fig = plt.figure(figsize = (4*l, 4), layout='constrained')
        ax = plt.gca()
        ax.autoscale_view()
        ax.autoscale(False)
        count = 0
        if end < 0:
            last = self.nframes + 1 + end
        else:
            last = end
        for i in range(start, last, step):
            print(i, "/", last)
            self.voronoi(i, prop, cmap, fig = fig, ax = ax, shear = shear, bar = bar)
            plt.savefig(self.path + f'batch_{prop}/{count:04}.png', pad_inches = 0.1, bbox_inches=None)
            ax.clear()
            count += 1
        plt.close()
        
    def save(self, prop, start = 0, end = -1, step = 1, frame_rate = 10, cmap="cmasher:viola", shear = True, bar = True, overide = False):
        if not os.path.exists(self.path + f'batch_{prop}') or overide:
            self.batch(prop, start, end, step, cmap = cmap, shear = shear, bar = bar)
        if True:
            inp = self.path + fr'batch_{prop}/*'
            output_video = self.path + f'batch_{prop}.gif'
    
            images = []
            for i in sorted(glob(inp)):
                images.append(imageio.imread(i))
            
            imageio.mimsave(output_video, images, fps=frame_rate)
        else:
            input_pattern = self.path + fr'batch_{prop}/%04d.png' 
            output_video = self.path + f'batch_{prop}.gif'
    
            ffmpeg_cmd = [
                "ffmpeg",
                "-y",
                "-framerate", str(frame_rate),
                "-i", input_pattern,
                output_video
            ]
            subprocess.run(ffmpeg_cmd, check=True)

if 0:      

    a = Cell("/mnt/ssd/Documents/Voronoi/dump/N_300qo_3.850000gamma_0.500000gammarate_0.100000Ka_0.000000v_3202.dump")
    #a.voronoi(frame = 15, shear = False, prop = "shear")
    'cmasher:pepper'
    a.save("shear", cmap ='jet', bar = True, shear = True, overide=False, frame_rate = 3)

if 0:
    files = glob("dump2/*.dump")
    cells = [Cell(file, False) for file in files]
    gammarate = np.array(sorted(list(set([cell.gammarate for cell in cells]))))
    gamma = np.array(sorted(list(set([cell.gamma for cell in cells]))))
    active = np.zeros((len(gamma), len(gammarate)))*np.nan
    time = np.zeros((len(gamma), len(gammarate)))*np.nan
    co = np.zeros((len(gamma), len(gammarate))).astype(str)
    marker = ["*", "s", "."]
    for cell in cells:
        g = np.where(cell.gamma == gamma)[0][0]
        G = np.where(cell.gammarate == gammarate)[0][0]
        active[g, G] = cell.frac_active[-1]/300
        time[g, G] = len(cell.frac_active)
        if active[g, G] > 0.5:
            co[g, G] = "red"
        else:
            co[g, G] = "black"
        
    plt.subplot(121)
    for i in range(len(gammarate)):
        plt.scatter(gamma, active[:, i], marker = marker[i], alpha = 0.7, label = gammarate[i])
    plt.legend(title = r"$\Delta \gamma$")
    plt.xlabel(r"$\gamma_{max}$")
    plt.ylabel(r"frac active")
    plt.subplot(122)
    for i in range(len(gammarate)):
        plt.scatter(gamma, time[:, i], marker = marker[i], alpha = 0.7, label = gammarate[i], c = co[:, i])
    plt.legend(title = r"$\Delta \gamma$")
    plt.xlabel(r"$\gamma_{max}$")
    plt.ylabel(r"number of cycles")

if 0:
    plt.plot(a.gamma_actual[1:], a.shear[1:])
    
    plt.xlabel(r"$\gamma$")
    plt.ylabel(r"$\sigma_{xy}$ (a.u)")