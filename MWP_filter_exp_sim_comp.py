# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 14:22:53 2021

@author: HS
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import find_peaks
plt.close('all')


class QDL_comb(object):
    
    def __init__(self, linewidth, fsr, bw_comb, n_comb, w, ws_profile):
        self.linewidth = linewidth
        self.fsr = fsr
        self.bw_comb = bw_comb        
        self.w = w
        self.ws_profile = ws_profile
        if n_comb>round(bw_comb/fsr):
            self.n_comb = round(bw_comb/fsr)
            print('Simulation bandwidth is too narrow, please increase the bandwidth')
        else:
            self.n_comb = n_comb
        
    def get_comb(self):
        center_w = (min(self.w)+max(self.w))/2
        p_combline = np.arange(center_w - self.n_comb/2 * self.fsr, center_w + self.n_comb/2 * self.fsr, self.fsr)
        c = self.linewidth / (2*np.sqrt(2*np.log(2)))  
        comb = 0
        for p in p_combline:
            comb += np.exp(-(self.w - p) ** 2 / (2 *c ** 2)) * self.ws_profile 
        return comb

class MWP_filter_ideal(object): # ideal modualtor and PD
    
    def __init__(self, dispersion, w, comb_spec, N_rf, Fs, FSR):
        self.dispersioin = dispersion
        self.comb_spec = comb_spec
        self.N_w = len(w)
        # freqeuncy to time mapping
        self.f_rf = np.linspace(-Fs/2, Fs/2, N_rf) # frequency vector
        delta_t = 1/Fs
        self.t_rf = np.arange(-delta_t*N_rf/2, delta_t*N_rf/2, delta_t) # time vector
        self.N_rf = N_rf
        self.t_w = np.linspace(min(dispersion * w), max(dispersion * w), self.N_w)

        self.t_interval = 1/FSR
    
    def f_response(self):
        
        impulse = np.interp(self.t_rf, self.t_w, self.comb_spec)
            
        f_filter = np.fft.fftshift(np.fft.fft(impulse))
        
        return f_filter, impulse
        
        

def rect(T):
    """create a centered rectangular pulse of width $T"""
    return lambda t: (-T/2 <= t) & (t < T/2)

def pulse_train(t, at, shape):
    """create a train of pulses over $t at times $at and shape $shape"""
    return np.sum(shape(t - at[:,np.newaxis]), axis=0)


def comb_exp_read(comb_name, w_fit):
# read the measured OFC data
# interp OFC data and conver to linear scale
    data = np.loadtxt(comb_name)
    comb_w = data[:,0]*1e-9
    comb_p = data[:,1]
    comb_p_fit = np.interp(w_fit, comb_w, comb_p)
    comb_p_linear = 10**((comb_p_fit)/10)  # linear scale


    return comb_p_linear, comb_w, comb_p

########### parameter setting for simulations ##################
N = 5000

linewidth = 0.5e-10 # unit: meter
fsr = 0.196e-9
bw_comb = 20e-9
n_comb_list = np.array([ 10, 20,40])
# n_comb_list = np.array([1, 2, 4, 6])
ws_profile = 1

n_comb = 50
# w = np.linspace(-bw_comb/2, bw_comb/2, N)
w = np.linspace(1527e-9, 1537e-9, N)
# disperion = 10*17e-3 # e-3: ps/nm to s/m
dispersion = 1655e-3 

N_rf = 50000 # points of RF domain
Fs = 100e9 # bw of RF domain
# coef_sinc = 10e9
# ws_profile = np.abs(np.sinc(coef_sinc*w)) # sinc shape

# bw_profile_list = np.array([2e-9, 1000e-9])
# bw_profile_list = 1e-9

# bw_profile_list = np.array([5,10])
# bw_profile = 10e-9
# step = 5e-9
# DC= 0.8

# coef_sinc = 0.5e9

# ws_profile = pulse_train(
#     t=w,              # time domain
#     at=np.arange(-bw_comb/2, bw_comb/2, step),  # times of pulses
#     shape=rect(DC*step)                 # shape of pulse
# )          

########### setting dir and measured OFC data name ##############

path = "F:\QDL\MWP_filter\Experiment\OSA\\2021-06-05\\"

##### unflatten #######
# comb_name_list = [path+"5combline_edfa.txt",
#                  path+"10combline_edfa.txt",
#                  path+"20combline_edfa.txt",
#                  path+"40combline_edfa.txt"]

###### flatten ##########
# comb_name_list = [
#                  path+"10comblines_flatten.txt",
#                  path+"20comblines_flatten.txt",
#                  path+"40comblines_flatten.txt"]

###### gaussian ##########
# comb_name_list = [path+"gs_0125T.txt",
    
#                   path+"gs_025T.txt",
                  
#                   path+"gs_05T.txt"
#                  ]
###### IPC #######

# comb_name_list = [path+"10comblines_flatten.txt",
    
#                   path+"20comblines_flatten.txt",
                  
#                   path+"gs_05T.txt"
#                  ]

# ###### gaussian WS after EDFA ########## path = "F:\QDL\MWP_filter\Experiment\OSA\\2021-05-19\\"
# comb_name_list = [path+"gs_125G.txt",
    
#                   path+"gs_250G.txt",
                  
#                   path+"gs_500G.txt"
#                   ]

###### sample WS after EDFA ########## path = "F:\QDL\MWP_filter\Experiment\OSA\\2021-05-24\\"
# comb_name_list = [path+"2channel_50dc.txt",
    
#                   path+"4channel_50dc.txt",
                  
#                   path+"8channel_50dc.txt"
#                   ]

###### flatten after EDFA ########## path = "F:\QDL\MWP_filter\Experiment\OSA\\2021-05-30\\"
# comb_name_list = [
    
#                   path+"10lines.txt",
                  
#                   path+"20lines.txt",
                  
#                   path+"40lines.txt"
                  
#                   ]

###### flatten after EDFA ########## path = "F:\QDL\MWP_filter\Experiment\OSA\\2021-06-05\\"
comb_name_list = [[path+"bandpass_06_09G_port1_V3.txt", path+"bandpass_06_09G_port2_V3.txt"],
    
                  [path+"bandpass_05_1G_port1_V3.txt", path+"bandpass_05_1G_port2_V3.txt"],
                  
                  [path+"bandpass_04_1_1G_port1_V3.txt", path+"bandpass_04_1_1G_port2_V3.txt"]
                  
                  ]

###### flatten after EDFA ########## path = "F:\QDL\MWP_filter\Experiment\OSA\\2021-06-06\\"
# comb_name_list = [[path+"bandpass_09_1.2G_port1.txt", path+"bandpass_09_1.2G_port2.txt"],
    
#                   [path+"bandpass_1.2_1.5G_port1.txt", path+"bandpass_1.2_1.5G_port2.txt"],
                  
#                   [path+"bandpass_overlap_port1.txt", path+"bandpass_overlap_port2.txt"]
                  
#                   ]

###### triangle tune bw pb ########## "F:\QDL\MWP_filter\Experiment\OSA\\2021-06-08\\"
# comb_name_list = [[path+"triagnle_bw0.3_fshift_0.75_port1.txt", path+"triagnle_bw0.3_fshift_0.75_port2.txt"],
    
#                   [path+"triagnle_bw0.5_fshift_0.75_port1.txt", path+"triagnle_bw0.5_fshift_0.75_port2.txt"],
                  
#                   [path+"triagnle_bw0.7_fshift_0.75_port1.txt", path+"triagnle_bw0.7_fshift_0.75_port2.txt"]
                  
#                   ]

###### triangle tune bw ########## "F:\QDL\MWP_filter\Experiment\OSA\\2021-06-08\\"
# comb_name_list = [[path+"triagnle_bw0.3_fshift_0_port1.txt"],
    
#                   [path+"triagnle_bw0.5_fshift_0_port1.txt"],
                  
#                   [path+"triagnle_bw0.7_fshift_0_port1.txt"]
                  
#                   ]

###### flat tune bw ########## "F:\QDL\MWP_filter\Experiment\OSA\\2021-06-08\\"
# comb_name_list = [[path+"flat0.3_fshift_0_port1.txt", path+"flat0.3_fshift_0_port2.txt"],
    
#                   [path+"flat0.5_fshift_0_port1.txt", path+"flat0.5_fshift_0_port2.txt"],
                  
#                   [path+"flat0.7_fshift_0_port1.txt", path+"flat0.7_fshift_0_port2.txt"]
                  
#                   ]

##########################################################################################################################

lw = 1.2
lw2 = 1.5
fs = 16
i=0
f_exp_list = [rf_f_5[k:len(rf_f_5)], rf_f_10[k:len(rf_f_5)], rf_f_20[k:len(rf_f_5)]]#, rf_f_40[k:len(rf_f_5)]]
fres_exp_list = [rf_p_5[k:len(rf_f_5)], rf_p_10[k:len(rf_f_5)], rf_p_20[k:len(rf_f_5)]]#,  rf_p_40[k:len(rf_f_5)]]

# f_exp_list = [rf_f_10[k:len(rf_f_5)], rf_f_20[k:len(rf_f_5)], rf_f_40[k:len(rf_f_5)]]
# fres_exp_list = [ rf_p_10[k:len(rf_f_5)], rf_p_20[k:len(rf_f_5)],  rf_p_40[k:len(rf_f_5)]]

# label = [['(a)' ,'(e)'],['(b)', '(f)'], ['(c)', '(g)'], ['(d)', '(h)']]

label = [['(a)' ,'(d)'],['(b)', '(e)'], ['(c)', '(f)']]


fig1, axs = plt.subplots(len(n_comb_list), 2, constrained_layout=True)
plt.subplots_adjust(wspace=0.45, hspace = 0.0)
# fontsize = 18


for n_comb in  n_comb_list:
    
    # comb1 = QDL_comb(linewidth, fsr, bw_comb, n_comb, w, ws_profile)
    # comb_spec = comb1.get_comb()
    
    # if len(comb_name_list[i])>1:
    comb_spec1,  w_origin, comb_spec_origin = comb_exp_read(comb_name_list[i][0], w)
    comb_spec2,  w_origin, comb_spec_origin = comb_exp_read(comb_name_list[i][1], w)
    # else:   
    comb_spec = comb_spec1 - comb_spec2
    #     comb_spec = comb_spec/max(np.abs(comb_spec))
    # else:
        
    # comb_spec, w_origin, comb_spec_origin = comb_exp_read(comb_name_list[i], w)  
    # index1, _ = find_peaks((comb_spec_origin), distance=70, height=max(comb_spec_origin)-30)
    # index2, _ = find_peaks(10*np.log10(comb_spec2), distance=70, height=max(10*np.log10(comb_spec2))-30)
    axs[i,0].plot(w_origin*1e9,comb_spec_origin, linewidth=lw, label = 'tap='+str(n_comb))#, label= 'FWHM='+str(bw_profile*1e9)+'nm') # 'tap='+str(n_comb)
    axs[i,0].plot(w*1e9,10*np.log10(comb_spec1), color='tab:blue', linewidth=lw2, label = 'FWHM='+str(n_comb)+'nm') # gaussian
    axs[i,0].plot(w*1e9,10*np.log10(comb_spec2), color='tab:red', linewidth=lw2, label = 'FWHM='+str(n_comb)+'nm') # gaussian
    # axs[i,0].plot(w_origin[index1]*1e9,(comb_spec_origin[index1]), '.',color='red', linewidth=lw) # gaussian
    # axs[i,0].plot(w[index2]*1e9,10*np.log10(comb_spec2[index2]), '.',color='darkblue', linewidth=lw) # gaussian
    # axs[i,0].legend()
    if i < len(n_comb_list)-1:    
        axs[i,0].get_xaxis().set_visible(False)
    # plt.tight_layout()
    # axs[i,0].legend(loc='upper left')
    plt.setp(axs[i,0], ylim=[-60, 1.4])
    plt.setp(axs[i,0], xlim=[1527, 1536])
    axs[i,0].text(1527.2, -10, label[i][0],fontsize=fs)
    if i == 1:
        axs[i,0].set(ylabel='Power (dBm)')
        
    if i == len(n_comb_list)-1:
        axs[i,0].set(xlabel='Wavelength (nm)')
    axs[i,0].spines['bottom'].set_linewidth(lw)#set the width of the bottom axis
    axs[i,0].spines['left'].set_linewidth(lw)
    axs[i,0].spines['right'].set_linewidth(lw)
    axs[i,0].spines['top'].set_linewidth(lw)
        
    # plt.rcParams.update({'font.size': 8})

    # plt.setp(axs[i], xlim=[-2.5, 2.5])
    
    #     axs[i].set(xlabel='wavelength (nm)', ylabel='Norm. Power(a.u.)')
    # # Hide x labels and tick labels for top plots and y ticks for right plots.
    #     axs[i].label_outer()
  

    Filter1 = MWP_filter_ideal(dispersion, w-(min(w)+max(w))/2, comb_spec, N_rf, Fs, 25e9)
    f = Filter1.f_rf
    f_resp, impulse = Filter1.f_response()
    
    # plt.tight_layout()
    index = 25010
    s = 1
    f = f[index:len(f_resp)]
    f_resp = f_resp[index:len(f_resp)]
    f_resp = abs(f_resp)/max(abs(f_resp))
    axs[i,1].plot(f*1e-9, 20*np.log10(f_resp),color='tab:green',  linewidth = lw2, label = 'Sim')
    axs[i,1].plot(f_exp_list[i], fres_exp_list[i]-max(fres_exp_list[i]), color='tab:orange', linewidth = lw2, label = 'Exp')#, label='FWHM='+str(bw_profile*1e9)+'nm')
    # axs[i,1].legend(loc='lower right')
    # plt.tight_layout()
    if i == 1:
        axs[i,1].set(ylabel='Norm. $|S_{21}|$ (dB)')
    if i == len(n_comb_list)-1: 
        axs[i,1].legend(loc='upper right', bbox_to_anchor=(1.03, 1.25), fontsize=12) #bbox_to_anchor=(0.8, 1.1)
    if i == len(n_comb_list)-1:
        axs[i,1].set(xlabel='Frequency (GHz)')
   
    plt.setp(axs[i,1], ylim=[-50, 10])
    plt.setp(axs[i,1], xlim=[0, 5])
    if i < len(n_comb_list)-1:    
        axs[i,1].get_xaxis().set_visible(False)
    axs[i,1].text(0.25, 0, label[i][1], fontsize=fs)

    axs[i,1].spines['bottom'].set_linewidth(lw)#set the width of the bottom axis
    axs[i,1].spines['left'].set_linewidth(lw)
    axs[i,1].spines['right'].set_linewidth(lw)
    axs[i,1].spines['top'].set_linewidth(lw)

    i+=1  
    

# 
fig1.set_size_inches(8,6)
# fig1.text(0.3, 0.01, 'wavelength (nm)', ha='center')
# fig1.text(0.005, 0.5, 'Norm. Power (linear scale)', va='center', rotation='vertical')
# fig1.tight_layout()
# fig1.text(0.74, 0.01, 'Frequency (GHz)', ha='center')
# fig1.text(0.49, 0.5, 'Norm. S21 (dBm)', va='center', rotation='vertical')
plt.rcParams.update({'font.size': fs})
ax=fig1.gca();
ax.spines['bottom'].set_linewidth(lw)#set the width of the bottom axis
ax.spines['left'].set_linewidth(lw)
ax.spines['right'].set_linewidth(lw)
ax.spines['top'].set_linewidth(lw)

# plt.tight_layout()
# plt.show
plt.savefig('F:\QDL\MWP_filter\Journal\Figure\\pb_bw_tune.eps', dpi = 600, bbox_inches='tight')
# fig1.set_size_inches(8,6)

# fig1.tight_layout() 

# plt.figure(99)

# plt.plot(Filter1.f_rf, impulse/max(impulse))
