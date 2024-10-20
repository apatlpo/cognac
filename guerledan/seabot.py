# seabot utils
import os
from glob import glob

from flatten_dict import flatten

import xarray as xr
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

import cartopy.crs as ccrs
crs = ccrs.PlateCarree()
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt

#root_path = "/Users/aponte/Current_Projects/ensta/guerledan_202310/data_seabot"
#root_dir = "/Users/aponte/Current_Projects/ensta/guerledan/"

# guerledan coords
g_lon, g_lat = 48.1964758, -3.0196383

piston_diameter = 45e-3 # mm
piston_area = np.pi*piston_diameter**2*0.25

def position2course(position):
    """ return piston course and volume displacement from position in pulses"""
    # 8192 pulse per turn / 1mm displacement per turn
    course = -position*1e-3/8192 # m, positive out
    volume = course*piston_area
    return course, volume

def load_rosbag(
    data_dir, 
    set_time=False, 
    verbose=False, 
    _browse=True,
):
    """ load rosbag data in a directory
    
    Parameters
    ----------
    data_dir: str
        path to data directory where rosbag files are
    verbose: boolean, optional
        activates verbosity
    set_time: boolean, optional
        set time from GPS data
    _browse: boolean, optional
        used for internal recursions when multiple rosbags are encountered

    Returns
    -------
    D: dict
        Dictionary containing data
    
    """

    # move deeper in dir tree
    if _browse:
        ros_dir = sorted(os.listdir(data_dir))
        ros_dir = [d for d in ros_dir if "." not in d]
        #assert len(ros_dir)==1, ("multiple rosbags?", data_dir)
        if len(ros_dir)>1:
            print("multiple rosbags encountered in "+data_dir)
            return {
                d: 
                load_rosbag(
                    os.path.join(data_dir, d), verbose=verbose, set_time=set_time, _browse=False
                )  
                for d in ros_dir
            }
        else:
            ros_dir = ros_dir[0]
    else:
        ros_dir = data_dir.split("/")[-1]
        data_dir = "/".join(data_dir.split("/")[:-1])
    
    # extract timestamp
    #print(ros_dir.replace('rosbag2_', ''))
    time_start = pd.to_datetime(ros_dir.replace('rosbag2_', ''), format="%Y_%m_%d-%H_%M_%S")
    ros_dir = os.path.join(data_dir, ros_dir)
    print(ros_dir)
    if not "data" in os.listdir(ros_dir):
        if verbose: print("no data directory in :"+ros_dir)
        return
    data_dir = os.path.join(ros_dir, "data")

    # start walking content
    content = sorted(os.listdir(data_dir))
    dirs = [c for c in content if os.path.isdir(os.path.join(data_dir, c))]
    files = [c for c in content if os.path.isfile(os.path.join(data_dir, c))]

    D = {}
    for d in dirs:

        if verbose: print(f"{d}") 

        subcontent = sorted(os.listdir(os.path.join(data_dir, d)))

        if not subcontent:
            if verbose: print("no data found in : ", data_dir, d)
            continue
            
        c = subcontent[0]
        for c in subcontent:
            if ".npz" in c:
                label = c.replace(".npz", "")
                file_path = os.path.join(data_dir, d, c)
                try:
                    with np.load(file_path, allow_pickle=False) as data:
                        df = pd.DataFrame({k: v for k, v in data.items()})
                    D[d+"_"+label] = df
                    if verbose: print(f"  {label}: "+" / ".join(c for c in df.columns))
                except Exception as err:
                    if verbose:
                        print(f"Unexpected {err=}, {type(err)=}")
                        print(file_path)
            else:
                if verbose: print(f"{c} skipped")
        
    # fix time:
    if set_time=="gps" and "driver_fix" in D:
        # get GPS time to set a proper clock
        # note that a clock drift could be computed
        df = D["driver_fix"]
        dt = np.mean( df["time_gnss"] - df["time"] )
        # could collect start/end positions
        for k, df in D.items():
            if "time" in df.columns:
                df["time_posix"] = df["time"] + dt
                df["time"] = pd.to_datetime(df["time_posix"], unit="s")
            D[k] = df.set_index("time").sort_index()
    else:
        if verbose: print("use rosbag filename to set time")
        for k, df in D.items():
            if "time" in df.columns:
                df["time_since_start"] = pd.to_timedelta(df["time"], unit="s")
                df["time"] = time_start + df["time_since_start"]
            D[k] = df.set_index("time").sort_index()
    
    # could collect general informations
    D["data_dir"] = data_dir
    
    return D

# with allow_pickle=False cannot load : observer/parameter.npz and driver/profile.npz
# driver_profile: distance / confidence / transmit_duration / ping_number / scan_start / scan_length / gain_setting / profile_data_length / profile_data / time_posix

def walk_load_repo(data_dir, top_dir=None, _top=True, **kwargs):
    """walk seabot data repository and load rosbags
    
    Parameters
    ----------
    data_dir: str
        Root directory to explore
    top_dir: str, optional
        Narrow walk to one sub directory
    """
    if not os.path.isdir(data_dir):
        return
    content = sorted(os.listdir(data_dir))
    content = [
        c for c in content if 
        not (
            c in [".DS_Store"]
            or ".txt" in c
        )
    ]
    if any(["rosbag" in c for c in content]):
        return load_rosbag(data_dir, **kwargs)
    else:
        if top_dir is not None and top_dir in content:
            content = [top_dir]
        D = {c: walk_load_repo(os.path.join(data_dir, c), _top=False, **kwargs) for c in content}
        if _top:
            D = flatten(D)
            D = {k:v for k, v in D.items() if v is not None} # filter out empty items
        return D

def dfilter(D, key):
    """select a key in a flatten seabot data dict"""
    Do = {}
    for k, v in D.items():
        if key in k:
            kn = list(k)
            kn.remove(key)
            Do[tuple(kn)] = v
    return Do

def combine(
    D, 
    rule, 
    op="median",
    **keys,
):
    """ resample and combine multiple keys
    The input dict must be flattened
    
    Dout = combine(D, "1s", dict(ka="observer_kalman", te="observer_temperature"))
    """
    _D = dfilter(D, next(iter(keys.values())))
    Dout = {}
    for k_out in list(_D):
        _D = []
        for suffix, k_in in keys.items():
            k = k_out+(k_in,)
            if k in D:
                df = D[k]
                df = df.rename(columns={c: suffix+"_"+c for c in df.columns})
                _D.append(df.resample("1s").apply(op))
        Dout[k_out] = pd.concat(_D, join="inner", axis=1)
    return Dout

def match_filter_with_cp(D, cp, filter=True):
    """ match data dict keys with campaign yaml file info
    filter based on start/end times

    Parameters
    ----------
    D: dict
        ensemble of dicts
    cp: pynsitu.events.campaign
    filter: boolean, optional
        activates filtering
    """
    Dd = dfilter(D, "data_dir")

    M = {}
    for plabel, p in cp.platforms.items():
        if "seabot" in plabel:
            for dlabel, d in p["deployments"].items():
                if "data_dir" in d:
                    data_dir = d["data_dir"]
                    #print(plabel, dlabel, data_dir)
                    for key, Ddata_dir in Dd.items():
                        #print(plabel, dlabel, Ddata_dir) 
                        if data_dir in Ddata_dir:
                            #print(plabel, dlabel, Ddata_dir) 
                            #break
                            M[(plabel, dlabel)] = key
    iM = {v: k for k, v in M.items()} # reverse dict
    # filter
    if filter:
        # filter with start/end time
        for p, k in M.items():
            start = cp[p[0]][p[1]].start.time
            end = cp[p[0]][p[1]].end.time
            for kk, df in D.items():
                if kk[0]==k[0] and kk[1]==k[1] and isinstance(df, pd.DataFrame):
                    D[kk] = df.loc[ (df.index>=start) & (df.index<=end) ]        
    return M, iM

def append_depth_filtered(df, tau, key="depth"):
    """ filter pressure, inplace """
    from scipy.signal import lfilter
    #dt, tau = 1, 2.5
    dt = 1
    alpha = dt/2/tau
    b, a = [alpha], [1, -(1-alpha)]
    df["depth_filtered"] = lfilter(b, a, df[key].bfill())
    dt = (df.reset_index()["time"].diff().bfill() / pd.Timedelta("1s")).values
    df["velocity_filtered"] = df["depth_filtered"].diff().bfill().values / dt

def plot_depth(
    D, 
    dkey="ka_depth", 
    legend=True, 
    thumbnail=None,
    colors=None,
    figsize=(15,4),
    time_slice=None,
    **kwargs,
):
    """ plot depth for all deployments 
    
    Parameter
    ---------
    D: seabot data dict
    dkey: depth variable key
    legend: boolean, optional
        add legend
    thumbnail:
        show as thumbnails
    """

    colors = _get_colors(colors, D)
    
    ax = None
    if thumbnail is None:
        fig, ax = plt.subplots(1,1, figsize=figsize)
    else:
        if isinstance(thumbnail, tuple):
            nx, ny = thumbnail
        else:
            nx, ny = 4, 1
    
    i=0
    for key, c in zip(D, colors):
        if isinstance(key, tuple):
            df = D[key]
        else:
            df = key
        if time_slice is not None:
            if pd.to_datetime(time_slice[0]) > df.index[-1] or pd.to_datetime(time_slice[1]) < df.index[0]:
                continue
        dkwargs = dict(color=c, label="/".join(key))
        dkwargs.update(**kwargs)
        if thumbnail is None:
            ax.plot(df.index, -df[dkey], **dkwargs)
        else:
            if i>=nx*ny:
                i=0
            if i==0:
                fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(ny, nx, i+1)
            #_ax = ax.flatten()[i]
            ax.plot(df.index, -df[dkey], **dkwargs)
            #ax.autofmt_xdate()
            plt.xticks(rotation=45)
            ax.set_title("/".join(key))
            ax.grid()
        i+=1
            
    if thumbnail is None:
        ax.grid()
        if legend:
            ax.legend()
    
    return fig, ax

def hv_plot(D, v, revert_yaxis=False, colors=None):
    
    #colors = get_cmap_colors(len(D))
    colors = _get_colors(colors, D)
    
    p = None
    for key, c in zip(D, colors):
        if isinstance(key, tuple):
            df = D[key]
        else:
            df = key
        if p is None:
            if v in df:
                p = df[v].hvplot(color=c, grid=True)
        else:
            if v in df:
                p = p * df[v].hvplot(color=c)
    
    if revert_yaxis:
        p = p.opts(invert_yaxis=True)
    
    return p

def _get_colors(colors, D):
    """ helper method to massage colors argument"""
    if colors is None:
        colors = get_cmap_colors(len(D))
    elif isinstance(colors, str):
        colors = [colors]*len(D)
    elif isinstance(colors, dict):
        colors = [colors[key] for key in D]
    return colors

def get_cmap_colors(Nc, cmap="plasma"):
    """load colors from a colormap to plot lines

    Parameters
    ----------
    Nc: int
        Number of colors to select
    cmap: str, optional
        Colormap to pick color from (default: 'plasma')
    """
    scalarMap = cmx.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=Nc), cmap=cmap)
    return [scalarMap.to_rgba(i) for i in range(Nc)]


def show_gps(Dp, Dk, i, start, full=False):
    """ show start/end gps tracks """

    key = list(Dp)[i]
    print(key)

    df = Dp[key]
    dfk = Dk[key]

    # extract start/end times
    depth_threshold = 0.5
    dt_threshold = pd.Timedelta("5min")

    dfk = Dk[key]
    start = dfk.loc[dfk.depth>depth_threshold].index[0]
    end = dfk.loc[dfk.depth>depth_threshold].index[-1]

    #dfs = df.loc[ (df.index<start) & (df.index>start-dt_threshold) ]
    dfs = df.loc[ (df.index<start) ]
    #dfs.plot(subplots=True, figsize=(10,10), layout=(5, 5));

    #dfe = df.loc[ (df.index>end) & (df.index<end+dt_threshold) ] 
    dfe = df.loc[ (df.index>end) ]
    #dfe.plot(subplots=True, figsize=(10,10), layout=(5, 5));

    if full:
        _df, _dfk = df, dfk
    elif start:
        _df, _dfk = dfs, dfk.loc[dfk.index<start]
    else:
        _df, _dfk = dfe, dfk.loc[dfk.index>end]

    return (
        _dfk["depth"].hvplot()
        +_df["mode"].hvplot()
        +_df["longitude"].hvplot() 
        +_df["latitude"].hvplot() 
    ).cols(2)

def load_gps(deployment, Dp):
    """ load gps positions based on deployment dict, inplace"""

    for key in Dp:
        
        df = Dp[key]
        
        df = df.loc[ df["mode"]==3 ]
        
        if df.index.size>0 and key in deployment:
            start =  df.loc[ df.index>=deployment[key]["start"] ].iloc[0]
            end =  df.loc[ df.index>=deployment[key]["end"] ].iloc[0]

            _dict = dict(
                start_lon=float(start.longitude),
                start_lat=float(start.latitude),
                end_lon=float(end.longitude),
                end_lat=float(end.latitude),
            )
        else:
            _dict = dict(
                start_lon=np.nan,
                start_lat=np.nan,
                end_lon=np.nan,
                end_lat=np.nan,
            )
        if key in deployment:
            deployment[key].update(**_dict)
            

def key2title(key, splitter):
    """ convert data dict key into a single string """
    key = [k.replace("test_guerledan_", "") for k in key]
    return splitter.join(key)

def pprint(D):
    """ pretty print data dict content"""
    fmt = "%Y/%m/%d %H:%M:%S"
    for d, df in D.items():
        start = df.index.min().strftime(fmt)
        end = df.index.max().strftime(fmt)
        print(" / ".join(d)+f" {start} to {end}")

# load data
def load_bathy(tiff_path):
    from pyproj import Transformer

    #if coverage=="east":
    #    tiff_file = "Guerledan_Feb19_0.5m.tiff"
    #tiff_path = os.path.join(root_dir, tiff_file)

    xds = (
        xr.open_dataset(tiff_path, engine="rasterio")
        .squeeze()
        .rename(band_data="height")
    )

    # compute lon/lat
    _, x, y = xr.broadcast(xds.height, xds.x, xds.y)
    transformer = Transformer.from_crs("EPSG:2154", "EPSG:4326")
    ll = transformer.transform(x, y)
    xds["lat"], xds["lon"] = ((x.dims, l) for l in ll)
    xds = xds.set_coords(["lon", "lat"])
    
    return xds

def plot_map(
    zoom="east",
    extent=None,
    bathy=None, 
    bathy_kwargs={},
    dx=10, 
    tile=True, tile_level=None,
    figsize=(5,5),
):
    """ plot bathymetry Guerledan
    
    Parameters
    ----------
    coverage
    """
    
    fig = plt.figure(figsize=figsize)

    if zoom=="full":
        _extent, level = [-3.1, -3., 48.19, 48.22], 12 # large
    elif zoom=="east":
        _extent, level = [-3.04, -3.01, 48.193, 48.205], 15 # east
    if tile_level is None:
        tile_level = level
    if extent is None:
        extent = _extent
        
    akwargs = dict()
    if tile:
        #tile = cimgt.GoogleTiles(style='terrain')
        #tile = cimgt.GoogleTiles(style='street')
        tile = cimgt.OSM()
        proj = tile.crs
    else:
        #lon, lat = float(ds.lon.mean()), float(ds.lat.mean())
        lon = (extent[0]+extent[1])*.5
        lat = (extent[2]+extent[3])*.5
        proj = ccrs.LambertAzimuthalEqualArea(central_longitude=lon, central_latitude=lat)
        akwargs.update(facecolor=cfeature.COLORS['land'])
    
    ax = fig.add_subplot(111, projection=proj, **akwargs)
    ax.set_extent(extent, crs=crs)

    if bathy is not None:
        ds = load_bathy(bathy).isel(x=slice(0,None,dx), y=slice(0,None,dx))
        height = ds.height - ds.height.min()
        hdl = ax.pcolormesh(
            ds.lon, ds.lat, height, transform=crs,
            **bathy_kwargs
            )

    if tile:
        ax.add_image(tile, tile_level)
    else:
        gl = ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
        gl.right_labels = False
        gl.top_labels = False
    if bathy:
        fig.colorbar(hdl)
    
    return fig, ax


