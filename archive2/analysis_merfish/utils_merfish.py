import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class RefLineSegs:
    """
    """
    def __init__(self, ps):
        """
        line segments defined by ps
        """
        self.ps = ps
        self.rps = ps[:-1] # reference points (exclude the last point)
        self.ns = len(ps)-1 # number of segments
        # get ts and uts
        self.calc_tan_vec()
        
        self.psn, self.pst = self.dists_to_qps(ps)
        self.fit_poly()
        
        return
    
    def fit_poly(self):
        """
        """
        poly_fit = np.poly1d(np.polyfit(self.ps[:,0], self.ps[:,1], 4))
        self.poly_fit = poly_fit
        
        return 
    
    def calc_tan_vec(self):
        """get tangent vectors
        """
        ps = self.ps 
        ts = ps[1:] - ps[:-1]
        lts = np.linalg.norm(ts, axis=1)
        uts = ts/lts.reshape(-1,1)
        nuts = uts.dot(np.array([[0,-1],[1,0]]).T)
        
        # np.power(uts[:,0],2)+np.power(uts[:,1],2) # check normed
        
        self.ts = ts # tangents
        self.uts = uts # unitary tangents
        self.nuts = nuts # norm to unitary tangents
        
        self.lts = lts # tangent lengths
        self.cumlts = np.hstack([0,np.cumsum(lts)])[:-1] # cumulative tangent lengths
        
        return 
    
    def calc_seg_len(self):
        """get the length of each segment
        """
        ts = self.ts
        return np.linalg.norm(ts, axis=1)
        
    def ndist_to_qp(self, query_point):
        """return the distance from a point to a set curve
        measure distance to each segment and take the min
        """
        query_vec = (query_point - self.rps)
        ndist = np.min(np.abs(query_vec[:,0]*self.nuts[:,0] + query_vec[:,1]*self.nuts[:,1]))
        
        return ndist
    
    def dists_to_qps(self, XY):
        """fixed a bug in the two functions below
        """
        query_mtx = np.repeat(XY[:,np.newaxis,:], self.ns, axis=1) - self.rps
        # tmp_dline = np.abs(query_mtx[:,:,0]*self.nuts[:,0] + query_mtx[:,:,1]*self.nuts[:,1]) # this calculate the distance to the lines - infinite - flawed two ends - ditched
        tmp_dsegm = np.linalg.norm(query_mtx, axis=2) # dist to segment
         
        # which segment - smallest distance to segments
        min_seg_idx = np.argmin(tmp_dsegm, axis=1)
        
        # smallest distance to lines - keep sign
        # ndists = np.min(tmp_dline, axis=1) # two ends cannot be too wild
        # ndists = np.min(tmp_dsegm, axis=1) # this can add a small correction - but meh
        # ndists = tmp_dline[np.arange(self.ns), min_seg_idx]
        
        matched_qry_mtx = XY - self.rps[min_seg_idx]
        matched_nrm_mtx = self.nuts[min_seg_idx]
        ndists = np.abs(matched_qry_mtx[:,0]*matched_nrm_mtx[:,0] + matched_qry_mtx[:,1]*matched_nrm_mtx[:,1]) # this calculate the distance to the lines
        
        signs = ((np.interp(XY[:,0], self.ps[:,0], self.ps[:,1]) - XY[:,1])>0).astype(int)
        signs[signs==0] = -1
        
        # cumsum up to the segment
        tdists = self.cumlts[min_seg_idx]
                
        uts = self.uts[min_seg_idx] # tangent vector for each ref point
        qs_vec = XY - self.rps[min_seg_idx] # 
        tdists_correction = qs_vec[:,0]*uts[:,0] + qs_vec[:,1]*uts[:,1]
        
        return signs*ndists, tdists+tdists_correction

    def ndist_to_qps(self, XY):
        """OLD - return the distance from a point to a set curve
        measure distance to each segment and take the min
        """
        query_mtx = np.repeat(XY[:,np.newaxis,:], self.ns, axis=1) - self.rps
        tmp = np.abs(query_mtx[:,:,0]*self.nuts[:,0] + query_mtx[:,:,1]*self.nuts[:,1])
        ndists = np.min(tmp, axis=1)
        
        return ndists
        
    def tdist_to_qps(self, XY):
        """OLD - return the distance from a point to a set curve
        adds up the distance of each segment
        """
        
        query_mtx = np.repeat(XY[:,np.newaxis,:], self.ns, axis=1) - self.rps
        tmp = np.abs(query_mtx[:,:,0]*self.nuts[:,0] + query_mtx[:,:,1]*self.nuts[:,1])
        
        # which segment
        min_seg_idx = np.argmin(tmp, axis=1)
        
        # cumsum up to the segment
        tdists = self.cumlts[min_seg_idx]
        
        uts = self.uts[min_seg_idx] # tangent vector for each ref point
        qs_vec = XY - self.rps[min_seg_idx] # 
        tdists_correction = qs_vec[:,0]*uts[:,0] + qs_vec[:,1]*uts[:,1]
        
        return tdists+tdists_correction
    
    # grid system
    def set_grid(self, t_interval=500, t_range=None):
        """
        """

        # t
        if t_range is None: # full range
            t = np.arange(0, np.max(self.pst), t_interval) # parameter runs tangentially along the curve
        else:
            l, r = t_range
            t = np.arange(l, r, t_interval) # parameter runs tangentially along the curve
            
        # tx, ty
        tx = np.interp(t, self.pst, self.ps[:,0]) # t -> x
        dx = 0.01

        ty = self.poly_fit(tx)
        dy = self.poly_fit(tx+dx) - ty

        s = np.sqrt(dx**2+dy**2)
        dx, dy = 1/s*dx, 1/s*dy

        self.grid = (t, tx, ty, dx, dy)
        
        return

    def plot_grid(self, ax, t_interval=500, v_length=1000, t_range=None):
        """
        """
        self.set_grid(t_interval=t_interval, t_range=t_range)
        
        t, tx, ty, dx, dy = self.grid
        x, y = self.ps[:,0], self.ps[:,1]

        ax.plot(x, y, color='k', linestyle='--', linewidth=1)
        for i in range(len(t)):
            _tx, _ty = tx[i], ty[i]
            _dx, _dy = v_length*dx[i], v_length*dy[i]

            ax.plot([_tx, _tx+_dy], 
                    [_ty, _ty-_dx], color='k', linestyle='--', linewidth=1) #y, c=ref_width, cmap='rainbow', s=1) # color='k', linestyle='--')
        return
    

def line_fitting(XY, x_interval=100):
    """
    """
    # sort points in increasing order
    XY = XY[np.argsort(XY[:,0])]
    x, y = XY[:,0], XY[:,1]

    # fit a line (x, y) -> (x_fit, y_fit)
    poly_fit = np.poly1d(np.polyfit(x, y, 4))
    x_fit = np.linspace(min(x), max(x), x_interval)
    y_fit = poly_fit(x_fit)

    # calc depth use x_fit, y_fit to get d, w for x, y
    XY_fit = np.vstack([x_fit, y_fit]).T
    XY_obj = RefLineSegs(XY_fit)
    depth, width = XY_obj.dists_to_qps(XY)
    DW = np.vstack([depth, width]).T

    return XY, DW, XY_fit, XY_obj, poly_fit 

def two_step_cut(x, y, xcut, ycut1, ycut2):
    """apply y1 left of x, y2 right of x 
    """
    cond = x < xcut
    overall_cond = np.logical_or(
        np.logical_and(cond,  y > ycut1), 
        np.logical_and(~cond, y > ycut2),)
    
    return overall_cond
    
# functions (the teacher wrote for you to use later)
def rot2d(x, y, theta, unit='degree'):
    """ rotate data points defined by `x` and `y` by `theta` degree
    """
    a = np.vstack([x,y]).T
    if unit == 'degree':
        theta = theta*np.pi/180 # convert to radian

    R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    ar = a.dot(R.T)
    return ar[:,0], ar[:,1]

def st_scatter(x, y, gexp=None, vmax_p=98, unit_norm=False, 
               title='', s=1, cbar_label='', output='', cmap='rocket_r', axis_off=True, 
               vmin=None, **cbar_kwargs):
    """customized scatter plot 
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
  
    fig, ax = plt.subplots(figsize=(10,8))
    if gexp is not None:
        vmax = np.percentile(gexp, vmax_p)
        if unit_norm:
            rgexp = gexp/vmax
            g = ax.scatter(x, y, c=rgexp, s=s, edgecolor='none', vmin=vmin, vmax=1, cmap=cmap, rasterized=True)
            fig.colorbar(g, label=cbar_label, shrink=0.3, **cbar_kwargs)
            title = title + f" (max {vmax:.2g} at {vmax_p:.2g} pctl)"
        else:
            g = ax.scatter(x, y, c=gexp, s=s, edgecolor='none', vmin=vmin, vmax=vmax, cmap=cmap, rasterized=True)
            fig.colorbar(g, label=cbar_label, shrink=0.3, **cbar_kwargs)
    else:
        g = ax.scatter(x, y, s=s, edgecolor='none', cmap=cmap, rasterized=True)
  
    if axis_off:
        ax.axis('off')
    ax.set_title(title)
    ax.set_aspect('equal')
  
    if output:
        powerplots.savefig_autodate(fig, output)
        
    return 

def st_scatter_ax(fig, ax, x, y, gexp=None, 
    vmin=None, vmax=None, vmin_p=2, vmax_p=98, unit_norm=False, 
    cmap='rocket_r', 
    alpha=1,
    title='', 
    s=1, 
    axis_off=True, 
    output='', 
    # cbar_label='', 
    # **cbar_kwargs,
    ):
    """customized scatter plot 
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
  
    if gexp is None:
        # do not color it 
        g = ax.scatter(x, y, s=s, edgecolor='none', cmap=cmap, rasterized=True)
    else: # color by gexp
        if vmax is None: vmax = np.percentile(gexp, vmax_p)
        if vmin is None: vmin = np.percentile(gexp, vmin_p)
        if unit_norm:
            gexp = gexp/vmax
            vmax = 1 # set to 1
  
        g = ax.scatter(x, y, c=gexp, s=s, edgecolor='none', vmin=vmin, vmax=vmax, cmap=cmap, rasterized=True, alpha=alpha)
        # fig.colorbar(g, label=cbar_label, shrink=0.3, **cbar_kwargs)

    ax.set_title(title)
    if axis_off:
        ax.axis('off')
    ax.set_aspect('equal')
        
    return g

# visualize clusters
def plot_cluster_simple_ax(fig, ax, clsts, x, y, s=1, axis_off=True, cmap=plt.cm.jet, suptitle=None, cbar=True, **kwargs):
    """this assumes `clsts` is a integer that starts from 0
    """
    from matplotlib import colors
  
    # unq_clsts, inv = np.unique(clsts, return_inverse=True)
    # n_unq = len(unq_clsts)
  
    norm = colors.BoundaryNorm(np.arange(-0.5, cmap.N, 1), cmap.N)

    g = ax.scatter(x, y, c=clsts, norm=norm, cmap=cmap, s=s, edgecolor='none', **kwargs)
    ax.set_aspect('equal')
    if axis_off:
        ax.axis('off')
    if cbar:
        fig.colorbar(g, ax=ax, label='clusters', ticks=[0,cmap.N-1], shrink=0.7) #  ticks=np.arange(n_unq), shrink=0.7)
    return fig

# visualize clusters
def plot_cluster(clsts, x, y, ux, uy, s=1, axis_off=True, cmap=plt.cm.jet, suptitle=None):
    """this assumes `clsts` is a integer that starts from 0
    """
    from matplotlib import colors
  
    unq_clsts, inv = np.unique(clsts, return_inverse=True)
    n_unq = len(unq_clsts)
    # colors = np.array(sns.color_palette('husl', n_unq))
    # c_vec = colors[inv]
  
    norm = colors.BoundaryNorm(np.arange(-0.5, n_unq, 1), cmap.N)
  
    fig, axs = plt.subplots(1, 2, figsize=(8*2,6))
    if suptitle is not None:
        fig.suptitle(suptitle)
    ax = axs[0]
    g = ax.scatter(x, y, norm=norm, cmap=cmap, c=clsts, s=s, edgecolor='none')
    ax.set_title('XY (spatial distribution)')
    ax.set_aspect('equal')
    if axis_off:
        ax.axis('off')
    
    ax = axs[1]
    ax.scatter(ux, uy, norm=norm, cmap=cmap, c=clsts, s=s, edgecolor='none')
    ax.set_title('UMAP (molecular similarity)')
    ax.set_aspect('equal')
    if axis_off:
        ax.axis('off')
  
    fig.colorbar(g, ax=ax, label='clusters', ticks=np.arange(n_unq), shrink=0.7)
    return fig, axs

def get_xyg(adata, gn, layer_name, x='x', y='y',):
    """get x,y, and a gene _
    log10(1+CP100)

    Args:
        adata (_type_): _description_
        gn (_type_): _description_
    """

    x = adata.obs['x'].values
    y = adata.obs['y'].values
    g = np.log10(np.ravel(1+np.array(adata[:,gn].layers[layer_name])))

    return x, y, g

def get_xygmean(adata, gns, layer_name, x='x', y='y',):
    """get x,y, and the mean of a few genes _

    first calculate: log10(1+CP100)

    then take the mean across genes on zscored values across cells
    Args:
        adata (_type_): _description_
        gn (_type_): _description_
    """
    from scipy.stats import zscore

    x = adata.obs['x'].values
    y = adata.obs['y'].values

    mat = np.log10(1+np.array(adata[:,gns].layers[layer_name]))
    g = np.mean(zscore(mat, axis=0), axis=1)

    return x, y, g

def binning(val, n):
    """
    """
    bins = np.linspace(np.min(val), np.max(val), n)
    binned = pd.cut(val, bins=bins)
    
    return bins, binned

def generate_discrete_cmap(n_colors_groups, keys=['tab20']):
    """
    """
    import matplotlib.colors as mcolors
    
    palette_list = []
    for n_colors, key in zip(n_colors_groups, keys):
        palette = sns.color_palette(key, n_colors)
        palette_list.append(palette)
    
    palette_list = np.vstack(palette_list)
    cmap = mcolors.ListedColormap(palette_list)
    
    return palette_list, cmap