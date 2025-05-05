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
from matplotlib.collections import EllipseCollection
cbar = None
axval = None
from scipy.interpolate import interp1d

def polygon_area(vertices):
    """Calculate the area of a polygon using the Shoelace formula."""
    x = vertices[:, 0]
    y = vertices[:, 1]
    area = 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
    return area

def polygon_perimeter(vertices):
    """Calculate the perimeter of a polygon."""
    perimeter = 0
    for i in range(len(vertices)):
        # Calculate the distance between consecutive vertices
        p1 = vertices[i]
        p2 = vertices[(i + 1) % len(vertices)]  # Wrap around to the first vertex
        perimeter += np.linalg.norm(p2 - p1)
    return perimeter

def cellArray(loc, dump = False):
    return [Cell(data, dump) for data in sorted(glob(loc))]

class Cell(Dump):
    def __init__(self, loc, dump = True):
        self.loc = loc
        if dump and loc[-1] == "p":
            Dump.__init__(self, loc)
        self.extract_variables()
        self.L = np.sqrt(self.N)
        if "dump" in loc:
            self.path = path = loc.split(".dump")[0]
        else:
            self.path = path = loc.split(".thermo")[0]
        
        try:
            if loc[-1] == "o":
                path = loc
            else:
                path = self.path + ".thermo"
                
            self.stop = 0
            with open(path, 'r') as file:
                try:
                    first_line = file.readline().strip()
                    column_names = first_line.split()
                    try:
                        data = np.loadtxt(path, skiprows = 1)
                    except:
                        data = np.genfromtxt(path, skip_footer = 1, skip_header = 1)
                        self.stop = 1
        
                    if data.ndim == 1:
                        data = data.reshape(1, -1)
        
                    for i, col_name in enumerate(column_names):
                        setattr(self, col_name, data[:, i])
                except:
                    print("could not recover .thermo")
        except:
            print("could not recover .thermo")
                
        if loc[-1] == "b":
            path = loc
        else:
            try:
                path = self.path + ".strob"
                
    
                with open(path, 'r') as file:
                    try:
                        first_line = file.readline().strip()
                        column_names = first_line.split()
                        if column_names[0] == "i":
                            try:
                                data = np.loadtxt(path, skiprows = 1)
                            except:
                                try:
                                    data = np.loadtxt(path, skip_footer = 1, skiprows = 1)
                                except:
                                    try:
                                        data = np.loadtxt(path, skip_footer = 2, skiprows = 1)
                                    except:
                                        pass
                
                            if data.ndim == 1:
                                data = data.reshape(1, -1)
                
                            for i, col_name in enumerate(column_names):
                                setattr(self, "s_"+col_name, data[:, i])
                        else:
                            data = np.loadtxt(path)
                            self.s_E = data[:, 1]
                            self.s_i = data[:, 0]
                            self.s_shear = data[:, 4]
                            self.s_frac_active = data[:, 5]
                    except:
                        print("could not recover .strob")
            except:
                print("could not recover .strob")
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
    def unwrap(self, start = 0):
        X = []
        Y = []
        
        for i in range(start, self.nframes):
            self.jump_to_frame(i)
            X.append(self.get_atompropf("x"))
            Y.append(self.get_atompropf("y"))
        
        X = np.array(X)
        Y = np.array(Y)
        count_x = np.zeros_like(X)
        count_y = np.zeros_like(Y)
        DX = np.diff(X, axis = 0)
        DY = np.diff(Y, axis = 0)
        
        for i in range(1, len(X)):
            count_x[i] = count_x[i - 1] - (DX[i - 1] > self.L/2).astype(int) + (DX[i - 1] < -self.L/2).astype(int)
            count_y[i] = count_y[i - 1] - (DY[i - 1] > self.L/2).astype(int) + (DY[i - 1] < -self.L/2).astype(int)
        
        return X + self.L*count_x, Y + self.L*count_y
      
    def msd(self, start=0, max_lag=None, dim='both'):
        X, Y = self.unwrap(start)
        
        n_frames = X.shape[0]
        n_particles = X.shape[1]
        
        if max_lag is None or max_lag >= n_frames:
            max_lag = n_frames - 1
        
        # Prepare arrays
        lags = np.arange(1, max_lag + 1)
        msd_values = np.zeros(max_lag)
        
        # Calculate MSD for each lag
        for lag in lags:
            # Calculate displacements
            if dim == 'x':
                displacements = X[lag:] - X[:-lag]
                squared_displacements = displacements**2
            elif dim == 'y':
                displacements = Y[lag:] - Y[:-lag]
                squared_displacements = displacements**2
            else:  # 'both'
                displacements_x = X[lag:] - X[:-lag]
                displacements_y = Y[lag:] - Y[:-lag]
                squared_displacements = displacements_x**2 + displacements_y**2
            
            msd_values[lag-1] = np.mean(squared_displacements)
        
        return lags, msd_values
    
    def diffusivity(self, time, msd):
        where = len(time)//5
        sta = where
        end = 2*where
        return np.polyfit(time[sta:end], msd[sta:end], 1)[0]
        
    def pbc(self, frame = -1, prop = None):
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
        if "radius" in self.get_atomheader():
            r = self.get_atompropf("radius")
        else:
            r = np.ones(len(x))*0.05
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
        radd = (
        list(r[mask1]) +
        list(r[mask2]) +
        list(r[mask3]) +
        list(r[mask4]) +
        list(r[mask1 & mask2]) +
        list(r[mask3 & mask4]) +
        list(r[mask1 & mask4]) +
        list(r[mask3 & mask2])
        )
        X = np.array(list(x) + xadd)
        Y = np.array(list(y) + yadd)
        R = np.array(list(r) + radd)

        
        if prop == None:
            return X, Y, R
        
        
        variables = re.findall(r'\b\w+\b', prop)
        local_vars = {var: self.get_atompropf(var) for var in variables}
    
        try:
            c = eval(prop, {"__builtins__": None, "np": np}, local_vars)
        except Exception as e:
            raise ValueError(f"Error evaluating property expression '{prop}': {e}")
        
        cadd = list(c[dL*y>L*(x-size)])+list(c[y<size])+list(c[np.logical_and(dL*y>L*(x-size),y<size)])+list(c[dL*y<L*(x+size-L)])+list(c[y>L-size])+list(c[np.logical_and(dL*y<L*(x+size-L),y>L-size)])+list(c[np.logical_and(dL*y>L*(x-size),y>L-size)])+list(c[np.logical_and(dL*y<L*(x+size-L),y<size)])
        C = np.array(list(c) + cadd)
        return X, Y, R, C
    def voronoi(self, frame = -1, prop = "area", cmap = "cmocean:delta", fig = None, ax = None, shear = True, bar = True, quant = None):
        global cbar, axval
        
        if frame < 0:
            frame = self.nframes + frame
        self.jump_to_frame(frame)
        X, Y, R, C = self.pbc(frame, prop)
        
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
            fig = plt.figure(figsize=(4*l, 4), layout='constrained')
            ax = plt.gca()
        ax.axis("square")
        ax.axis("off")
        
        pc = PolyCollection(polys, facecolors=colors, edgecolors='k', alpha=0.8)
        ax.add_collection(pc)
        

        #ax.plot(X[:self.N], Y[:self.N], '.', c='k')

        widths = 2*R[:self.N] 
        heights = 2*R[:self.N]

        circles = EllipseCollection(
            widths, heights, 
            np.zeros_like(widths),
            offsets=np.column_stack((X[:self.N], Y[:self.N])),
            units='x',
            edgecolors='none',
            facecolors='k',
            linewidths=0.5,
            transOffset=ax.transData
        )
        ax.add_collection(circles)

        if shear:
            ax.set_xlim(-self.L*self.gamma, self.L*(1 + self.gamma))
        else:
            ax.set_xlim(0, self.L)
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
                fig.set_size_inches(4*l, 5.2*1.1)
            else:
                norm = Normalize(vmin=np.min(C), vmax=np.max(C))
                sm = ScalarMappable(cmap=colormap, norm=norm)
                sm.set_array([])
                
                cbar.update_normal(sm)
        if quant != None:
            if axval is None:
                divider = make_axes_locatable(ax)
                axval = divider.append_axes("top", size="10%", pad=0.05)
                # get fig actual size
                a, b = fig.get_size_inches()
                fig.set_size_inches(a, b*1.3)
                axval.set_xlim(0, len(self.i))
                axval.set_xticks([])
            if axval is not None:
                axval.clear()
                axval.plot(self.i[:frame], self.shear[:frame])
                axval.set_xlim(0, len(self.i))
                axval.set_ylim(np.min(self.shear), np.max(self.shear))
                axval.set_xticks([])
        
    def batch(self, prop, start, end, step, cmap="cmasher:viola", shear = True, bar = True, quant = None):
        
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
            self.voronoi(i, prop, cmap, fig = fig, ax = ax, shear = shear, bar = bar, quant = quant)
            plt.savefig(self.path + f'batch_{prop}/{count:04}.png', pad_inches = 0.1, bbox_inches=None)
            ax.clear()
            count += 1
        plt.close()
        
    def save(self, prop, start = 0, end = -1, step = 1, frame_rate = 10, cmap="cmasher:viola", shear = True, bar = True, overide = False, quant = None):
        if not os.path.exists(self.path + f'batch_{prop}') or overide:
            self.batch(prop, start, end, step, cmap = cmap, shear = shear, bar = bar, quant = quant)
        if True:
            inp = self.path + fr'batch_{prop}/*'
            output_video = self.path + f'batch_{prop}.gif'
    
            images = []
            for i in sorted(glob(inp)):
                images.append(imageio.imread(i))
            
            imageio.mimsave(output_video, images, fps=frame_rate, loop = 0)
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
    cells = cellArray("/mnt/ssd/Documents/Voronoi/data/simple shear/test/rk0003/*.thermo", False)
    for i in range(len(cells)):
        cell = cells[i]
        plt.plot(cell.gamma_actual, cell.shear)
    plt.ylabel(r"$\sigma_{xy}$")
    plt.xlabel(r"$\gamma$")

if 1:
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    plt.subplots_adjust(hspace=1)

    
    allData = glob("/mnt/ssd/Documents/Voronoi/data/soft core/fire 0.1/*.thermo")
    
    cells = [Cell(data, dump = False) for data in allData]
    X = np.zeros((len(allData), 61500))*np.nan
    g = []
    for i, cell in enumerate(cells):
        stress = cell.shear
        X[i, :len(stress)] = stress
        if len(cell.gamma_actual) > len(g):
            g = cell.gamma_actual
    ax.plot(g, np.nanmean(X, axis = 0)[:len(g)], label = r"fire $dt_{max} = 0.1$")
    ax.set_ylabel(r"$\sigma_{xy}$")
    ax.set_xlabel(r"$\gamma$")
    ax.set_title(r"Mean shear stress")
    
    
    
    allData = glob("/mnt/ssd/Documents/Voronoi/data/soft core/fire 0.01/*.thermo")
    
    cells = [Cell(data, dump = False) for data in allData]
    X = np.zeros((len(allData), 61500))*np.nan
    g = []
    for i, cell in enumerate(cells):
        stress = cell.shear
        X[i, :len(stress)] = stress
        if len(cell.gamma_actual) > len(g):
            g = cell.gamma_actual
    ax.plot(g, np.nanmean(X, axis = 0)[:len(g)], label = r"fire $dt_{max} = 0.01$")
    ax.set_ylabel(r"$\sigma_{xy}$")
    ax.set_xlabel(r"$\gamma$")
    ax.set_title(r"Mean shear stress")
    
    
    allData = glob("/mnt/ssd/Documents/Voronoi/data/soft core/rk/*.thermo")
    
    cells = [Cell(data, dump = False) for data in allData]
    
    all_times = [cell.gamma_actual for cell in cells]  
    common_time = np.linspace(0, max([times[-1] for times in all_times]), 1500) 
    interpolated_X = np.zeros((len(cells), len(common_time))) * np.nan
    for i, (cell, times) in enumerate(zip(cells, all_times)):
        stress = cell.shear
        interp_func = interp1d(times, stress, bounds_error=False, fill_value=np.nan)
        interpolated_X[i, :] = interp_func(common_time)
    
    mean_stress = np.nanmean(interpolated_X, axis=0)
    
    ax.plot(common_time, mean_stress, label = r"continuous time, rk45 adaptative,  $\gamma_{dot}=0.001$")
    ax.set_ylabel(r"$\sigma_{xy}$")
    ax.set_xlabel(r"$\gamma$")
    ax.set_title(r"Mean shear stress")
    """for i in range(2):
        ax[1].plot(cells[i].gamma_actual, cells[i].shear, alpha = 0.5)
    ax[1].set_ylabel(r"$\sigma_{xy}$")
    ax[1].set_xlabel(r"$\gamma$")"""
    ax.legend()
    
    
if 0:
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    plt.subplots_adjust(hspace=1)

    
    allData = glob("/mnt/ssd/Documents/Voronoi/data/simple shear/test/0.1/*.thermo")
    
    cells = [Cell(data, dump = False) for data in allData]
    X = np.zeros((len(allData), 61500))*np.nan
    for i, cell in enumerate(cells):
        stress = cell.shear
        X[i, :len(stress)] = stress
    g = np.linspace(0, 0.0025*len(X[0]), len(X[0]))
    ax.plot(g, np.nanmean(X, axis = 0), label = r"fire $dt_{max} = 0.1$")
    ax.set_ylabel(r"$\sigma_{xy}$")
    ax.set_xlabel(r"$\gamma$")
    ax.set_title(r"Mean shear stress")
    
    
    
    allData = glob("/mnt/ssd/Documents/Voronoi/data/simple shear/test/0.01/*.thermo")
    
    cells = [Cell(data, dump = False) for data in allData]
    X = np.zeros((len(allData), 5000))*np.nan
    for i, cell in enumerate(cells):
        stress = cell.shear
        X[i, :len(stress)] = stress
    g = np.linspace(0, 0.0025*len(X[0]), len(X[0]))
    ax.plot(g, np.nanmean(X, axis = 0), label = r"fire $dt_{max} = 0.01$")
    ax.set_ylabel(r"$\sigma_{xy}$")
    ax.set_xlabel(r"$\gamma$")
    ax.set_title(r"Mean shear stress")
    """for i in range(0, 20):
        ax[1].plot(g, X[i], alpha = 0.5)
    ax[1].set_ylabel(r"$\sigma_{xy}$")
    ax[1].set_xlabel(r"$\gamma$")
    ax[1].set_title(r"Shear stress for 4 different realizations")"""
    
    
    allData = glob("/mnt/ssd/Documents/Voronoi/data/simple shear/test/rk003/*.thermo")
    
    cells = [Cell(data, dump = False) for data in allData]
    
    all_times = [cell.gamma_actual for cell in cells]  
    common_time = np.linspace(0, max([times[-1] for times in all_times]), 1500) 
    interpolated_X = np.zeros((len(cells), len(common_time))) * np.nan
    for i, (cell, times) in enumerate(zip(cells, all_times)):
        stress = cell.shear
        interp_func = interp1d(times, stress, bounds_error=False, fill_value=np.nan)
        interpolated_X[i, :] = interp_func(common_time)
    
    mean_stress = np.nanmean(interpolated_X, axis=0)
    
    ax.plot(common_time, mean_stress, label = r"continuous time, rk45 adaptative,  $\gamma_{dot}=0.003$")
    ax.set_ylabel(r"$\sigma_{xy}$")
    ax.set_xlabel(r"$\gamma$")
    ax.set_title(r"Mean shear stress")
    """for i in range(2):
        ax[1].plot(cells[i].gamma_actual, cells[i].shear, alpha = 0.5)
    ax[1].set_ylabel(r"$\sigma_{xy}$")
    ax[1].set_xlabel(r"$\gamma$")"""
    ax.legend()
    plt.show()
    
    
    allData = glob("/mnt/ssd/Documents/Voronoi/data/simple shear/test/rk0003/*.thermo")
    
    cells = [Cell(data, dump = False) for data in allData]
    
    all_times = [cell.gamma_actual for cell in cells]  
    common_time = np.linspace(0, max([times[-1] for times in all_times]), 1500) 
    interpolated_X = np.zeros((len(cells), len(common_time))) * np.nan
    for i, (cell, times) in enumerate(zip(cells, all_times)):
        stress = cell.shear
        interp_func = interp1d(times, stress, bounds_error=False, fill_value=np.nan)
        interpolated_X[i, :] = interp_func(common_time)
    
    mean_stress = np.nanmean(interpolated_X, axis=0)
    
    ax.plot(common_time, mean_stress, label = r"continuous time, rk45 adaptative,  $\gamma_{dot}=0.0003$")
    ax.set_ylabel(r"$\sigma_{xy}$")
    ax.set_xlabel(r"$\gamma$")
    ax.set_title(r"Mean shear stress")
    """for i in range(2):
        ax[1].plot(cells[i].gamma_actual, cells[i].shear, alpha = 0.5)
    ax[1].set_ylabel(r"$\sigma_{xy}$")
    ax[1].set_xlabel(r"$\gamma$")"""
    ax.legend()
    plt.show()

if 0:
    cell = "1"
    x, y, a = cell.pbc(0, "perimeter")
    points = np.vstack((x, y)).T
    vor = Voronoi(points, qhull_options="Qc")
    
    region_idx = vor.point_region[0]
    
    # Get the region vertices indices
    region = vor.regions[region_idx]
    
    if not region or -1 in region:
        print("Particle 0 has an invalid Voronoi cell (possibly extending to infinity)")
    else:
        # Get the actual vertices coordinates
        vertices = vor.vertices[region]
        
        print("Voronoi vertices for particle 0:")
        for i, vertex in enumerate(vertices):
            print(f"Vertex {i}: ({vertex[0]:.16f}, {vertex[1]:.16f})")
            
        # Calculate cell properties
        area = polygon_area(vertices)
        perimeter = polygon_perimeter(vertices)
        
        print(f"\nVoronoi cell properties:")
        print(f"Area: {area:.16f}")
        print(f"Perimeter: {perimeter:.16f}")
        print(f"Number of vertices: {len(vertices)}")
if 0:
    # Calculate areas for each cell
    areas = np.zeros(len(points))
    for i, region in enumerate(vor.regions):
        if not region or -1 in region:
            continue
        # Find which point this region belongs to
        point_idx = None
        for j, region_idx in enumerate(vor.point_region):
            if region_idx == i and j < len(x):  # Ensure it's an original point
                point_idx = j
                break
        if point_idx is not None:
            polygon = vor.vertices[region]
            areas[point_idx] = polygon_perimeter(polygon)

if 0:    
    data = Cell("/mnt/ssd/Documents/Voronoi/dump/N_300qo_3.900000gamma_0.010000gammarate_0.010000Ka_0.000000v_1.dump")
    fx = data.get_atompropf("fx")
    fxa = np.loadtxt("/home/syrocco/Downloads/a.txt", skiprows = 1, delimiter=",")[:, 4]
    diff = (fx - fxa)/fx
    print(np.round(diff, 2))
if 0:
    x = np.genfromtxt("landscape.txt", dtype='str')
    X = x[:, 0].astype(np.float128)
    Y = x[:, 1].astype(np.float128)
    """ try:
        #Y = [float(a.split(".")[1][10:]) for a in x[:, 1]]
    except:"""
    Y = [float(a) for a in x[:, 1]]
    plt.scatter(X, Y)
    plt.yscale("log")
    plt.xlabel("distance to starting point")
    plt.ylabel("E")
if 0:      
    cell = Cell("/mnt/ssd/Documents/Voronoi/data/simple shear/test/0.01/N_300qo_3.900000dt_0.010000gamma_3.000000gammarate_-0.002500Ka_0.000000v_10.dump")
    #cell.voronoi(frame = -2, shear = True)
    'cmasher:pepper'
    cell.save("shear", cmap ='jet', bar = False, shear = True, overide=False, frame_rate = 2, quant = None)


if 0:
    cell = Cell("/mnt/ssd/Documents/Voronoi/dump/N_300qo_3.900000gamma_5.000000gammarate_0.010000Ka_0.000000v_1.dump", False)
    plt.plot(cell.gamma_actual[-10000:], cell.shear[-10000:])
    plt.xlabel(r"$\gamma$")
    plt.ylabel(r"$\sigma_{xy}$")
if 0:
    files = glob("/mnt/ssd/Documents/Voronoi/data/absorbing/CG local/*.dump")
    cells = [Cell(file) for file in files]
    gammarate = np.array(sorted(list(set([round(cell.gammarate, 3) for cell in cells]))))
    gamma = np.array(sorted(list(set([cell.gamma for cell in cells]))))
    active = np.zeros((len(gamma), len(gammarate)))*np.nan
    diffusivity = np.zeros((len(gamma), len(gammarate)))*np.nan
    time = np.zeros((len(gamma), len(gammarate)))*np.nan
    co = np.zeros((len(gamma), len(gammarate))).astype(str)
    marker = ["*", "s", "."]
    for i, cell in enumerate(cells):
        print(i, "/", len(cells))
        g = np.where(cell.gamma == gamma)[0][0]
        G = np.where(round(cell.gammarate, 3) == gammarate)[0][0]
        if np.isnan(cell.E[-1]):
            active[g, G] = np.nan
            time[g, G] = np.nan
            diffusivity[g, G] = np.nan
        else:
            a = cell.s_frac_active
            active[g, G] = np.mean(a[-1])
            time[g, G] = len(a[a>0.001])
            """t, msd = cell.msd(cell.nframes//4)
            diffusivity[g, G] = cell.diffusivity(t, msd)"""
        
    plt.subplot(131)
    for i in range(len(gammarate)):
        plt.scatter(gamma, active[:, i], marker = marker[i], alpha = 0.7, label = gammarate[i])
    plt.legend(title = r"$\Delta \gamma$")
    plt.xlabel(r"$\gamma_{max}$")
    plt.ylabel(r"frac active")
    plt.subplot(132)
    for i in range(len(gammarate)):
        plt.scatter(gamma, time[:, i], marker = marker[i], alpha = 0.7, label = gammarate[i], c = co[:, i])
    plt.legend(title = r"$\Delta \gamma$")
    plt.xlabel(r"$\gamma_{max}$")
    plt.ylabel(r"number of cycles")
    plt.subplot(133)
    for i in range(len(gammarate)):
        plt.scatter(gamma, diffusivity[:, i], marker = marker[i], alpha = 0.7, label = gammarate[i], c = co[:, i])
    plt.legend(title = r"$\Delta \gamma$")
    plt.xlabel(r"$\gamma_{max}$")
    plt.ylabel(r"diffusivity")

if 0:
    cell = Cell("/mnt/ssd/Documents/Voronoi/dump/N_300qo_3.550000gamma_0.200000gammarate_-0.000000Ka_3.000000v_2.dump", False)
    plt.plot(cell.E)
    plt.xscale("log")

    
