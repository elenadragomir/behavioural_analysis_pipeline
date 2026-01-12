
import numpy as np
import matplotlib.pyplot as plt
from copy import copy
import seaborn as sns
from bouter.utilities import extract_segments_above_threshold


def stims_separation (stim_logs, beh_logs):
# returns stims, a dict with separated behavior log for each fish (each key is a diff stim period)
        stims={}
        for i in range (len(stim_logs)-1):
            stims[i]=beh_logs.loc[np.logical_and(beh_logs['t']>=stim_logs[i]['t_start'], 
                                                        beh_logs['t']<stim_logs[i+1]['t_start'])]
        stims[i+1]=beh_logs.loc[beh_logs['t']>=stim_logs[i+1]['t_start']]
        
        return (stims)
    
    
    
def bouts_separation (stim_logs, bout_summary):
# separating bouts per stim/adapt period
    bouts_per={}
    for i in range (len(stim_logs)-1):
        bouts_per[i] = bout_summary.loc[np.logical_and(bout_summary['t_start']>=stim_logs[i]['t_start'], 
                                                       bout_summary['t_start']<stim_logs[i+1]['t_start'])]
    bouts_per[i+1]=bout_summary.loc[bout_summary['t_start']>=stim_logs[i+1]['t_start']]
    
    return (bouts_per)

from matplotlib.patches import Wedge

def dual_half_circle(center, radius, angle=0, ax=None, colors=('w','k'), alpha=0.2,
                     **kwargs):
    """
    Add two half circles to the axes *ax* (or the current axes) with the 
    specified facecolors *colors* rotated at *angle* (in degrees).
    """
    if ax is None:
        ax = plt.gca()
    theta1, theta2 = angle, angle + 180
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[1], alpha=alpha, **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[0], alpha=0, **kwargs)
    for wedge in [w1, w2]:
        ax.add_artist(wedge)
    return [w1, w2]



def plot_trajectories (n_fish, id_fish, beh_stims, div_pos, light_side_1, light_side_3): 
    
    #exp['stimulus']['display_params']['size'][0]/2
    #div_pos=exp['general']['program_version']['arguments']['camera']['roi'][2]/2 #division line between dark and light 
    circle_l = plt.Circle((392,392), 370, color='k', fill=False)
    circle_d = plt.Circle((392,392), 370, color='k', alpha=0.2)

    sns.set(style='ticks', palette='deep')
    custom_xlim = (-10, 810)
    custom_ylim = (-10, 810)


    fig, axs = plt.subplots(n_fish, len(beh_stims), figsize=(9,5)) #plotting traj for every fish (row), during each stim (column)
    plt.setp(axs, xlim=custom_xlim, ylim=custom_ylim)

    #just for getting proper time in colorbar
    sc1=axs[0].scatter(beh_stims[0].f0_x, beh_stims[0].f0_y, c=beh_stims[0].t, cmap='magma', s=0.1, zorder=10)

    for i_stims in range (len(beh_stims)):
        axs[i_stims].set_aspect(1)
        axs[i_stims].scatter(beh_stims[i_stims]['f0_x'], beh_stims[i_stims]['f0_y'], c=beh_stims[i_stims].t, cmap='magma', s=0.1, zorder=10)
        axs[i_stims].set_axis_off()


    circ_l=copy(circle_l)
    circ_d=copy(circle_d)
    axs[0].add_patch(circ_l)
    circ_l=copy(circle_l)
    axs[2].add_patch(circ_d)
    axs[2].add_patch(circ_l)
    circ_l=copy(circle_l)
    axs[1].add_patch(circ_l)
    circ_l=copy(circle_l)
    axs[3].add_patch(circ_l)
    
    if light_side_1=='right':
        dual_half_circle((div_pos, div_pos), radius=370, angle=90, ax=axs[1])
    else:
        dual_half_circle((div_pos, div_pos), radius=370, angle=270, ax=axs[1])
    
    if light_side_3=='right':
        dual_half_circle((div_pos, div_pos), radius=370, angle=90, ax=axs[3])
    else:
        dual_half_circle((div_pos, div_pos), radius=370, angle=270, ax=axs[3])

    axs[0].set_title('Light adaptation')
    axs[1].set_title('Photo after light')
    axs[2].set_title('Dark adaptation')
    axs[3].set_title('Photo after dark')

    fig.subplots_adjust(wspace=0.2, hspace=0, right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.3, 0.02, 0.4]) # numbers are x,y, width, height

    fig.colorbar(sc1, cax=cbar_ax, label= 'Time (s)')
    fig.suptitle(id_fish, fontsize=16)
    #plt.tight_layout()
    return axs

def classify_bouts(bout_summary):
# here i categorize into bouts- adjust bout summary theta abs threshold for more stringency
    bout_summary['bout_classifier']=np.nan
    bout_summary['bout_classifier'].loc[(bout_summary['theta_abs']>15) & (bout_summary['theta_abs']<100)]=1 # 1 is turn to the right
    bout_summary['bout_classifier'].loc[(bout_summary['theta_abs']<-15) & (bout_summary['theta_abs']>-100)]=-1 # -1 is turn to the left
    bout_summary['bout_classifier'].loc[(bout_summary['theta_abs']>=-15) & (bout_summary['theta_abs']<=15)]=0 #forward swims
    bout_summary['bout_classifier'].loc[np.isnan(bout_summary['bout_classifier'])]=0 # 0 is forward swim (less than 15 deg)
    
    return (bout_summary)


def bout_numbers (bouts_fish_stim):
#calculating n of left, right and forward swims per each stim period for each fish 
#(input is bouts df separated per stim and per fish)

    n_left= len (bouts_fish_stim.loc[bouts_fish_stim['bout_classifier']==-1])
    n_right=len (bouts_fish_stim.loc[bouts_fish_stim['bout_classifier']==1])
    n_forward=len (bouts_fish_stim.loc[bouts_fish_stim['bout_classifier']==0])
    
    return (n_left, n_right, n_forward)

def time_in_light (beh_stim_per, divider_position):
#figuring out time spent in each half field during all conditions
    if beh_stim_per.light_side.iloc[0]=='right':
        times_light=beh_stim_per['t'].loc[beh_stim_per['f0_x']>divider_position]
        times_dark=beh_stim_per['t'].loc[beh_stim_per['f0_x']<divider_position]
    else:
        times_light=beh_stim_per['t'].loc[beh_stim_per['f0_x']<divider_position]
        times_dark=beh_stim_per['t'].loc[beh_stim_per['f0_x']>divider_position]
    
    fr_light= len(times_light)/(len(times_light)+len(times_dark))
    
    return (fr_light)


def thigmo_time(beh_stim_per, center):
# time spent along the edges of the dish
    dist= np.sqrt((center-beh_stim_per['f0_x'])**2+
                  (center-beh_stim_per['f0_y'])**2)
    t_center=beh_stim_per.loc[dist<300]
    t_wall=beh_stim_per.loc[dist>300]
    thig_t=len(t_wall)/(len(t_wall+t_center))

    return (thig_t)

def time_out (beh_stim_per):   
# time spent out of range (fish not being detected either because of freezing or on the edge)
    t_out=len(beh_stim_per[beh_stim_per['f0_x'].isnull()])/len(beh_stim_per)
    return (t_out) 

def bout_duration(bouts_fish_stim):
#calculates the average bout duration, per stim period, per fish
#input: bouts df, separated before by fish and stim

    bout_durations_stim=[]
    for i in range (len(bouts_fish_stim)):
        bout_durations_stim.append(bouts_fish_stim.t_end.iloc[i]-bouts_fish_stim.t_start.iloc[i])
        
    return (np.nanmean(bout_durations_stim))


def bout_displacement(bouts_fish_stim):
#calculates the average bout duration, per stim period, per fish
#input: bouts df, separated before by fish and stim

    bout_disps_stim=[]
    for i in range (len(bouts_fish_stim)):
        bout_disps_stim.append(np.sqrt((bouts_fish_stim['x_end'].iloc[i]-bouts_fish_stim['x_start'].iloc[i])**2+
                                       (bouts_fish_stim['y_end'].iloc[i]-bouts_fish_stim['y_start'].iloc[i])**2))
        
    return (np.nanmean(bout_disps_stim))   


def mean_vel (vel_per):
# calculates mean vel and onset mean vel (first 15 bouts) per fish, per stim period
# input: velocities df separated into stim periods, and one fish 
    vel_bouts, cont1= extract_segments_above_threshold(vel_per.values, threshold=1, min_length=15, min_between=5)
    vel_mean=[]
    #this splices the vel df during only bout times
    for k in range(len(vel_bouts)):
        vel_mean.append(np.nanmean(vel_per.iloc[vel_bouts[k][0]:vel_bouts[k][1], :].values)) 
        
    return (np.nanmean(vel_mean), np.nanmean(vel_mean[:15]))


def light_frs_bins (beh_stim_per, div_pos, n_bins, stim_per_t_start, stim_bin_duration):
#figuring out time spent in each half field during all conditions divided by 1 min bins
    
    times_bins_light=[]
    times_bins_dark=[]
    fr_bins_light=[]

    for i in range (int(n_bins)):
        if beh_stim_per.light_side.iloc[0]=='right':
            times_bins_light.append (beh_stim_per['t'].loc[np.logical_and 
                                                           (beh_stim_per['t']>=stim_per_t_start+ (stim_bin_duration*(i)),
                                                            beh_stim_per['t']<=stim_per_t_start+ (stim_bin_duration*(i+1)))].loc[beh_stim_per['f0_x']>div_pos])
            times_bins_dark.append (beh_stim_per['t'].loc[np.logical_and 
                                                           (beh_stim_per['t']>=stim_per_t_start+ (stim_bin_duration*(i)),
                                                            beh_stim_per['t']<=stim_per_t_start+ (stim_bin_duration*(i+1)))].loc[beh_stim_per['f0_x']<div_pos])
        else:
            times_bins_light.append (beh_stim_per['t'].loc[np.logical_and 
                                                           (beh_stim_per['t']>=stim_per_t_start+ (stim_bin_duration*(i)),
                                                            beh_stim_per['t']<=stim_per_t_start+ (stim_bin_duration*(i+1)))].loc[beh_stim_per['f0_x']<div_pos])
            times_bins_dark.append (beh_stim_per['t'].loc[np.logical_and 
                                                           (beh_stim_per['t']>=stim_per_t_start+ (stim_bin_duration*(i)),
                                                            beh_stim_per['t']<=stim_per_t_start+ (stim_bin_duration*(i+1)))].loc[beh_stim_per['f0_x']>div_pos])
        try:
            fr_bins_light.append (len(times_bins_light[i])/(len(times_bins_light[i])+len(times_bins_dark[i])))
        except:
            fr_bins_light.append(np.NaN)
    
    return (fr_bins_light)


def total_dist(beh_stim_per, px_in_mm):
# calculates total distance traveled per stimulus period
    abs_pos=np.sqrt(beh_stim_per['f0_x']**2 + beh_stim_per['f0_y']**2)

    if np.isnan(abs_pos.iloc[0]):
        first_pos=0
    else:
        first_pos=abs_pos.iloc[0]

    dist_total= px_in_mm*(np.nancumsum(np.abs(np.diff(abs_pos.values, prepend=0)))-first_pos)
    
    return (dist_total[-1])


def freezing_times(bouts_fish_stim):
#calculates nr of freezing episodes per stim period
    n=0
    for i in range (len(bouts_fish_stim)-1):
        if bouts_fish_stim.t_start.iloc[i+1]-bouts_fish_stim.t_start.iloc[i]>5:
            n+=1
            
    return (n)


def interbout_int(bouts_fish_stim):
# calculates average interbout interval
    ibi=[]

    for i in range (len(bouts_fish_stim)-1):
        if bouts_fish_stim.t_start.iloc[i+1]-bouts_fish_stim.t_start.iloc[i]<10:
            ibi.append(bouts_fish_stim.t_start.iloc[i+1]-bouts_fish_stim.t_start.iloc[i])

    return (np.nanmean(ibi))