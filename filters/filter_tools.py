from scipy import signal
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import traceback

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def print_warning(text):
    print(bcolors.WARNING + text + bcolors.ENDC)

def check(value, message=None):
    if not value:
        stack = traceback.extract_stack()
        frame = stack[-2]
        print_warning("")
        print_warning(f"WARNING ! FAIL TO CHECK :")
        print_warning(f"    {frame.line}")
        print_warning(f"    file : {frame.filename} - line : {frame.lineno}")
        print_warning("")
        if not( message is None ):
            print_warning("    " + message)
            print_warning("")

Hz = 1
kHz = 10**3
MHz = 10**6
s = 1
ms = 10**-3
us = 10**-6
ns = 10**-9
F = 1
uF = 10**-6
nF = 10**-9
pF = 10**-12
Ohms = 1
kOhms = 10**3
MOhms = 10**6
dB = 1

minimal_human_percption_time = 13 * ms

def scaleToDb(scale):
    return -20 * np.log10(scale)

def dbToScale(db):
    return 10**(-db/20)

resolution = 2**14

def compute_phases(angles):
    last = angles[0]
    res = [angles[0]]
    offset = 0
    for i in range(1, len(angles)):
        angle = angles[i]
        if(abs(angle-last) > 2*np.pi - np.pi/8):
            if angle > last:
                offset -= 1
            else:
                offset += 1
        elif(abs(angle-last) > np.pi - np.pi/8):
            if angle > last:
                offset -= 1/2
            else:
                offset += 1/2
        last = angle
        res.append( angle + 2*np.pi*offset )
    return np.array(res)


################################################################################
# Resistor series
################################################################################
E6_base = [1.0, 1.5, 2.2, 3.3, 4.7, 6.8]
E12_base = [1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2]
E24_base = [
    1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1
]
E48_base = [1.00, 1.05, 1.10, 1.15, 1.21, 1.27, 1.33, 1.40, 1.47, 1.54, 1.62, 1.69, 1.78, 1.87, 1.96, 2.05, 2.15, 2.26, 2.37, 2.49, 2.61, 2.74, 2.87, 3.01, 3.16, 3.32, 3.48, 3.65, 3.83, 4.02, 4.22, 4.42, 4.64, 4.87, 5.11, 5.36, 5.62, 5.90, 6.19, 6.49, 6.81, 7.15, 7.50, 7.87, 8.25, 8.66, 9.09, 9.53]
E96_base = [1.00, 1.02, 1.05, 1.07, 1.10, 1.13, 1.15, 1.18, 1.21, 1.24, 1.27, 1.30, 1.33, 1.37, 1.40, 1.43, 1.47, 1.50, 1.54, 1.58, 1.62, 1.65, 1.69, 1.74, 1.78, 1.82, 1.87, 1.91, 1.96, 2.00, 2.05, 2.10, 2.15, 2.21, 2.26, 2.32, 2.37, 2.43, 2.49, 2.55, 2.61, 2.67, 2.74, 2.80, 2.87, 2.94, 3.01, 3.09, 3.16, 3.24, 3.32, 3.40, 3.48, 3.57, 3.65, 3.74, 3.83, 3.92, 4.02, 4.12, 4.22, 4.32, 4.42, 4.53, 4.64, 4.75, 4.87, 4.99, 5.11, 5.23, 5.36, 5.49, 5.62, 5.76, 5.90, 6.04, 6.19, 6.34, 6.49, 6.65, 6.81, 6.98, 7.15, 7.32, 7.50, 7.68, 7.87, 8.06, 8.25, 8.45, 8.66, 8.87, 9.09, 9.31, 9.53, 9.76]
E24_minus_12_base = []
for R in E24_base:
    if R not in E12_base:
        E24_minus_12_base.append(R)

E6 = []
E12 = []
E24_minus_12 = []
E24 = []
E24 = []
E48 = []
E96 = []
for N in [0, 1, 2, 3, 4, 5, 6]:
    factor = 10**N
    for Exx, Exx_base in [
        (E6, E6_base), (E12, E12_base), (E24_minus_12, E24_minus_12_base), (E24, E24_base), (E48, E48_base), (E96, E96_base)
    ]:
        for R in Exx_base:
            Exx.append(R*factor)

E1_base = {1.0}

C1 = []
C6 = []
C12 = []
C24 = []
C48 = []
C96 = []
for Cxx, Exx_base in [
    (C1, E1_base), (C6, E6_base), (C12, E12_base), (C24, E24_base), (C48, E48_base), (C96, E96_base)
]:
    for C in Exx_base:
        for N in [1, 10, 100]:
            for unit in [pF, nF, uF]:
                factor = N * unit
                Cxx.append(C*factor)
        Cxx.append(C*1000*uF)


def compute_best_tension_divisor(ratio, resistor_series):
    resistances = list(filter( lambda x : 10**3<= x < 10**6, resistor_series ))
    r1, r2 = 1.0, 1.0
    def div(r1, r2):
        return (r1+r2)/r1
    sol = div(r1, r2)
    for res1 in resistances:
        for res2 in resistances:
            new_sol = div(res1, res2)
            if abs(sol-ratio) > abs(new_sol-ratio):
                sol = new_sol
                r1, r2 = res1, res2
    return r1, r2

low_pass_constant = {
    'butterworth'
}

def compute_sallen_key_cell(pole, resistors, capacitors):
    # pole = pole/np.abs(pole)
    l2 = (pole * pole.conjugate()).real
    l1 = (pole + pole.conjugate()).real
    l0 = 1
    polynome = [
        1, # X 1
        -l1/l2, # X j.w
        l0/l2, # X (j.w)**2
    ]
    # polynome[1] = R*C*(3-gain)
    # polynome[2] = (R*C)**2
    # print(polynome)
    error = None
    gain_optimal = 1
    R_optimal, C_optimal = 0, 0
    all_resistors = list(filter(lambda x: 10*kOhms <= x <= 1*MOhms, resistors))
    assert(len(all_resistors)>0)
    all_capacitors = list(filter(lambda x: 1*nF <= x <= 220 * nF, capacitors))
    assert(len(all_capacitors)>0)
    for R in all_resistors :
        for C in all_capacitors:
            gain = 3 - polynome[1]/(R*C)
            # Stability check
            if gain <= 0 or gain >= 3:
                continue
            new_error = np.abs(polynome[2] - (R*C)**2)
            if error is None or  new_error < error :
                R_optimal = R
                C_optimal = C
                gain_optimal = gain
                error = new_error
    return {'R':R_optimal, 'C':C_optimal, 'gain': gain_optimal, 'error':error}

def print_sallen_key_cells(sallen_key_cells, resistors, set_vin=True, set_vout=True):
    for i in range(len(sallen_key_cells)):
        cell = sallen_key_cells[i]
        C = cell['C']
        R = cell['R']

        print(f"Cellule {i+1} de Sallen-key : ")
        print("")
        print(f"C1 = C2 = {C/nF} nF")
        print(f"R1 = R2 = {R/kOhms} kOhms")
        K_sallen_key = cell['gain']
        r1, r2 = compute_best_tension_divisor(
            ratio=K_sallen_key, resistor_series = resistors
        )
        print(f"R : {r1/kOhms} kOhms")
        print(f"(K-1)R : {r2/kOhms} kOhms")
        print(f"Gain, K : {K_sallen_key}")
        if set_vin and i == 0:
            marker_in = "Vin"
        else:
            marker_in = f" V{i}"
        if set_vout and i == len(sallen_key_cells)-1:
            marker_out = "Vout"
        else:
            marker_out = f" V{i+1}"
        print(
            f"""
                      C1
                -------||-----------------------------
                |                       _______      |
    {marker_in} --[R1]--o--[R2]--o-------------| +     |     |
                         |             | .     |---o-o-- {marker_out}
                         = C2       |--| -     |   |
                         |          |  ---------   |
                         |          |              |
                         o-----[R]--o--[(K-1).R]----
                         |
                        gnd
    """
        )

################################################################################
# Analog filter
################################################################################

def compute_filter_cut_frequency_range(b, a, low_attenuation, high_attenuation, sample_freq):
    factor = 2
    w, h = signal.freqs(b, a, worN=2*np.pi*np.linspace(0,factor*sample_freq, resolution))
    while 20*np.log10(np.abs(h[-1])) > -high_attenuation :
        factor += 1
        w, h = signal.freqs(b, a, worN=2*np.pi*np.linspace(0,factor*sample_freq, resolution))
    margin_attenuation = 0.01 * dB
    for i in range(len(h)-1, -1, -1):
        if( 20*np.log10(np.abs(h[i])) > -high_attenuation + margin_attenuation ): 
            f_cut_high_attenuation = w[i] / (2*np.pi)
            break
    for i in range(len(h)):
        if( 20*np.log10(np.abs(h[i])) < -low_attenuation - margin_attenuation ): 
            f_cut_low_attenuation = w[i]/ (2*np.pi)
            break
    return f_cut_low_attenuation, f_cut_high_attenuation

def error(x, args):
    high_frequency = x[0]
    sample_freq = args['sample_freq']
    low_attenuation = args['low_attenuation']
    high_attenuation = args['high_attenuation']
    target_f_cut_low_attenuation = args['target_f_cut_low_attenuation']
    target_f_cut_high_attenuation = args['target_f_cut_high_attenuation']
    filter_ = args['filter']
    b, a = filter_(high_frequency, args)
    _, f_cut_high_attenuation = compute_filter_cut_frequency_range(b, a, low_attenuation, high_attenuation, sample_freq)
    return np.abs( 
        f_cut_high_attenuation - target_f_cut_high_attenuation
    )

def find_analog_filter(
    order, sample_freq, f_max,
    low_attenuation = 3*dB, high_attenuation = 40*dB, 
    filter_name = 'butterworth', margin = 1/10.0
):
    data = {
        'low_attenuation' : low_attenuation,
        'high_attenuation' : high_attenuation,
        'target_f_cut_low_attenuation' : f_max * (1-margin) + margin * sample_freq/2,
        'target_f_cut_high_attenuation' : f_max * margin + (1-margin)*sample_freq/2,
        'sample_freq' : sample_freq,
    }
    x0 = .5*f_max + .5*sample_freq/2
    if filter_name == 'butterworth':
        data['filter'] = lambda f_cut, _: signal.butter(order, 2*np.pi * f_cut, 'lowpass', analog=True,output='ba')
        data['filter_zpk'] = lambda f_cut, _: signal.butter(order, 2*np.pi * f_cut, 'lowpass', analog=True,output='zpk')
        data['type'] = 'butterworth',
        data['f_cut'] = None
        data['order'] = order
    elif filter_name == 'chebychev1':
        data['filter'] = lambda f_cut, _: signal.cheby1(order, low_attenuation, 2*np.pi*f_cut, 'lowpass', analog=True, output='ba')
        data['filter_zpk'] = lambda f_cut, _: signal.cheby1(order, low_attenuation, 2*np.pi*f_cut, 'lowpass', analog=True, output='zpk')
        data['type'] = 'chebychev1',
        data['f_cut'] = None
        data['order'] = order
        data['rp'] = low_attenuation,
    
    cons = (
        {'type': 'ineq', 'fun': lambda x:  x[0] - f_max},
        {'type': 'ineq', 'fun': lambda x:  sample_freq/2 - x[0]},
    )
    res = optimize.minimize(
        error, x0, args=data, 
        options= {'eps':1.0}, tol=0.2, constraints=cons
    )
    # print(res)
    data['f_cut'] = res.x[0]
    b, a = data['filter'](res.x[0], data)
    data['zeros'], data['poles'], data['k'] = data['filter_zpk'](res.x[0], data)

    data['b'], data['a'] = b, a 
    data['f_cut_low_attenuation'], data['f_cut_high_attenuation'] = compute_filter_cut_frequency_range(
        b, a, low_attenuation, high_attenuation, sample_freq
    )

    data['w'], data['h'] = signal.freqs(b, a, worN=2*np.pi*np.linspace(0,2*sample_freq, resolution))
    data['f'] = data['w']/(2*np.pi)

    # phases = np.angle(h) - 2*np.pi*(np.angle(h)>0)
    phases = compute_phases(np.angle(data['h']))
    delays = phases[1:]/(2*np.pi*data['f'][1:])
    signal_delays = delays * (data['f'][1:]<=f_max)
    max_signal_delay = np.max(-signal_delays)
    data['max_delay'] = max_signal_delay
    data['phases'] = phases
    data['delays'] = delays

    return data

def compute_sallen_key_cells(poles, resistors, capacitors):
    sallen_key_cells = []
    for i in range(len(poles)//2):
        sallen_key_cells.append(
            compute_sallen_key_cell(poles[i], resistors, capacitors)
        )
    return sallen_key_cells

def compute_RC_cell(pole, resistors, capacitors):
    best_value = 0
    R_optimal = 0
    C_optimal = 0
    error = 1/pole 

    all_resistors = list(filter( lambda x: 10*kOhms<= x < 1*MOhms, resistors ))
    assert(len(all_resistors)>0)
    all_capacitors = list(filter(lambda x: 1*nF <= x <= 220 * nF, capacitors))
    assert(len(all_capacitors)>0)

    for R in all_resistors:
        for C in all_capacitors:
            new_value = R * C
            new_error = np.abs(new_value - 1/pole)
            if error < new_error :
                new_error = error
                value = new_value
                R_optimal = R
                C_optimal = C
    return {'R': R_optimal, 'C': C_optimal, 'error' : error}

def print_RC_cell(RC_cell, marker):
    C = RC_cell['C']
    R = RC_cell['R']

    print(f"Cellule {marker} RC : ")
    print("")
    print(f"C = {C/nF} nF")
    print(f"R = {R/kOhms} kOhms")
    print(
        f"""
                      _______
V{marker} ---[R]----o--------| +     |
             |        |       |---o-- Vout
             = C   |--| -     |   |
             |     |  ---------   |
             |     |              |
             o     ----------------
             |
            gnd
"""
    )


def compute_divisor(value, resistors):
    assert(value <= 1)
    all_resistors = list(filter( lambda x : 1*kOhms<= x < 1*MOhms, resistors ))
    assert(len(all_resistors)>0)

    print("### OK ")
    print(all_resistors)
    optimal_value = 1/2.0
    error = np.abs(value - optimal_value)
    optimal_r1 = all_resistors[0]
    optimal_r2 = all_resistors[0]
    for r1 in all_resistors:
        for r2 in all_resistors:
            new_value = r1/(r1+r2)
            new_error = np.abs(value - new_value)
            if new_error < error :
                error = new_error
                optimal_r1 = r1
                optimal_r2 = r2
                optimal_value = new_value
    return {'R1':optimal_r1, 'R2':optimal_r2, 'error':error}


def print_divisor(cell, marker_in="Vin", marker_out="Vout"):
    R1 = cell['R1']
    R2 = cell['R2']
    print(f"Diviseur : ")
    print("")
    print(f"R1 = {R1/kOhms} kOhms")
    print(f"R2 = {R2/kOhms} kOhms")
    tab = ' ' * len(marker_in)
    print(
        f"""
{tab}                      _______
{marker_in} ---[R2]----o--------| +     |
{tab}            |        |       |---o-- {marker_out}
{tab}           [R1]   |--| -     |   |
{tab}            |     |  ---------   |
{tab}            |     |              |
{tab}            |     ----------------
{tab}            |
{tab}           gnd
"""
    )

def print_buffer(marker_in="Vin", marker_out="Vout"):
    print(f"Buffer : ")
    tab = ' ' * len(marker_in)
    print(
        f"""
{tab}        _______
{marker_in} ------| +     |
{tab}       |       |---o-- {marker_out}
{tab}    |--| -     |   |
{tab}    |  ---------   |
{tab}    |              |
{tab}    ----------------
"""
    )

def design_analog_filter(
    analog_filter_order, wanted_attenutation, sample_freq,
    max_signal_freq, resistors, capacitors, 
    #filter_name = 'butterworth'
    filter_name = 'chebychev1', wanted_low_attenuation = 3*dB
):
    print(
    """
##########################################################################################
# Analog filter design
##########################################################################################

"""
    )

    generated_filter = find_analog_filter(
        analog_filter_order, sample_freq, max_signal_freq,
        low_attenuation = wanted_low_attenuation, high_attenuation = wanted_attenutation,
        filter_name = filter_name, margin = 1/10.0
    )

    # The Bode diagram
    plt.plot(generated_filter['f'], 20*np.log10(np.abs(generated_filter['h'])))
    plt.title('Analog filter frequency response')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Amplitude [dB]')

    # Some usefull vertical lines
    plt.axvline(x=max_signal_freq, color="red", label='Signal (Joystick) frequency')
    plt.axvline(x=generated_filter['f_cut'], color="green", label='cut frequency')
    plt.axvline(x=sample_freq, color="orange", label='sample frequency')
    plt.axvline(x=sample_freq/2, color="blue", label='sample frequency/2')
    plt.axvline(x=generated_filter['f_cut_high_attenuation'], color="magenta", linestyle='dashed' , label=f"frequence at {int(wanted_attenutation)} dB attenuation")
    plt.xscale('log')

    # Some usefull horizontal lines
    plt.axhline(
        y=20*np.log10(.9), color="0.7", linestyle="dashed",
        label='attenuation X0.9'
    )
    plt.axhline(
        y=-3, color="0.5", linestyle="dashed",
        label='attenuation X0.7 (-3db)'
    )
    plt.axhline(
        y=20*np.log10(.1), color="0.35", linestyle="dashed",
        label='attenuation X0.1 (-20db)'
    )
    plt.axhline(
        y=20*np.log10(.01), color="black", linestyle="dashed",
        label='attenuation X0.01 (-40db)'
    )
    plt.legend()
    #plt.ylim((-200, 10))
    plt.show()


    plt.axvline(
        x=max_signal_freq, color="red", linestyle="dashed", 
        label='Signal (joystick) frequency'
    )
    plt.plot(generated_filter['f'][1:], generated_filter['delays'] )
    plt.xlabel("Frequence f (Hz)")
    plt.ylabel("filter delay (s)")
    plt.title("Filter delays of the analog filter")
    plt.legend()
    plt.show()

    print(
        "The maximum delay for the filtered signal is : %s ms"%(
            generated_filter['max_delay']/ms
        )
    )

    sallen_key_cells = compute_sallen_key_cells(generated_filter['poles'], resistors, capacitors)
    total_gain = 1
    print("###############")
    for cell in sallen_key_cells:
        print(cell["gain"])
        total_gain *= cell['gain']
    print("")
    print_buffer("Vin", "Vcopy")
    if len(generated_filter['poles']) > 1:
        divisor = compute_divisor(1/total_gain, resistors)
        print_divisor(divisor, "Vcopy", "V0")
        set_vin = False
    else:
        set_vin = True
    print_sallen_key_cells(
        sallen_key_cells, resistors, set_vin=set_vin,
        set_vout=len(generated_filter['poles']) % 2 == 0
    )
    if len(generated_filter['poles']) % 2 == 1 :
        RC_cell = compute_RC_cell(
            generated_filter['poles'][len(generated_filter['poles'])//2],
            resistors, capacitors
        )
        print_RC_cell(RC_cell, len(generated_filter['poles'])//2)
    return generated_filter

################################################################################
# Digital filter
################################################################################

def compute_frequences_of_digital_filter(
    b, a, f_sample, wanted_attenuation_at_fsys_over_2, low_attenuation
):
    w, h = signal.freqz(b, a, fs=f_sample, worN = resolution)
    
    f_cut_low_attenuation = 0
    for i in range(len(h)):
        if 20*np.log10(np.abs(h[i])) < -low_attenuation :
            f_cut_low_attenuation = w[i]
            break

    attenuation_margin = .01
    f_cut_high_attenuation = 0
    for i in range(len(h)-1, -1, -1):
        if 20*np.log10(np.abs(h[i])) > -wanted_attenuation_at_fsys_over_2 + attenuation_margin :
            f_cut_high_attenuation = w[i]
            break
    return w, h, f_cut_low_attenuation, f_cut_high_attenuation


"""
Parameter :
    wanted_attenuation_at_fsys_over_2 : .
                            This parameter is ignored for the butterworth filter

"""
def compute_digital_filter(
    order, f_max, f_sys, f_sample, filter_name='butterworth',
    wanted_attenuation_at_fsys_over_2 = 40*dB, float_type = np.float32, low_attenuation = 3*dB
):
    margin = 1/10.0
    f_cut_low_attenuation = (1-margin)*f_max + margin*f_sys/2
    f_cut_high_attenuation = margin*f_max + (1-margin)*f_sys/2
    res = { 'filter':filter_name }

    res['low_attenuation'] = low_attenuation
    res['high_attenuation'] = wanted_attenuation_at_fsys_over_2
    res['target_f_cut_high_attenuation'] = f_cut_high_attenuation
    res['f_sample'] = f_sample


    if filter_name == 'chebychev2':
        stop_band_attenuation = wanted_attenuation_at_fsys_over_2 # dB
        b_driver, a_driver = signal.cheby2(
            order, stop_band_attenuation, f_cut_high_attenuation, analog=False, output='ba',
            fs=f_sample
        )
        res['f_cut'] = f_cut_high_attenuation
        res['stop_band_attenuation'] = stop_band_attenuation
    elif filter_name == 'butterworth':
        def filter_(x):
            return signal.butter(
                order, x[0], analog=False, output='ba', fs=f_sample
            )
        def error(x):
            b_driver, a_driver = signal.butter(
                order, x[0], analog=False, output='ba', fs=f_sample
            )
            _, _, _, obtained_f_cut_high_attenuation = compute_frequences_of_digital_filter(
                b_driver, a_driver, f_sample, wanted_attenuation_at_fsys_over_2, low_attenuation
            )
            return np.abs(f_cut_high_attenuation - obtained_f_cut_high_attenuation)

        cons = (
            {'type': 'ineq', 'fun': lambda x:  x[0] - f_cut_low_attenuation},
            {'type': 'ineq', 'fun': lambda x:  f_cut_high_attenuation - x[0]},
        )
        x0 = np.array([.5 * f_cut_low_attenuation + .5 * f_cut_high_attenuation])
        res = optimize.minimize(
            error, x0, 
            options= {'eps':1.0}, tol=0.2, constraints=cons
        )
        b_driver, a_driver = filter_(res.x)
        res['f_cut'] = res.x[0]

    res['b'] = float_type(b_driver)
    res['a'] = float_type(a_driver)

    w, h, obtained_f_cut_low_attenuation, obtained_f_cut_high_attenuation = compute_frequences_of_digital_filter(
        res['b'], res['a'], f_sample, wanted_attenuation_at_fsys_over_2, low_attenuation
    )

    res['w'] = w
    res['h'] = h
    res['f_cut_low_attenuation'] = obtained_f_cut_low_attenuation
    res['f_cut_high_attenuation'] = obtained_f_cut_high_attenuation

    return res 

def test_fiter(b_driver, a_driver, f_max, f_sample):
    f_signal_1 = f_max/4 * Hz
    f_signal_2 = f_max/3 * Hz
    f_noise_1 = 2*f_max * Hz
    f_noise_2 = 10*f_max * Hz
    t = np.arange(0, 3/f_signal_1, 1/f_sample)
    x = 3*np.sin(2*np.pi*f_signal_1*t) + 4*np.cos(2*np.pi*f_signal_2*t)
    x_noise = x + .5*np.cos(2*np.pi*f_noise_1*t) + .3*np.cos(2*np.pi*f_noise_2*t)
    y = signal.lfilter(b_driver, a_driver ,x_noise)
    plt.plot(t, x_noise, label='signal with noise')
    plt.plot(t, x, color="green", label='originel signal')
    plt.plot(t, y, label='filtered signal')
    plt.xlabel('time [s]')
    plt.ylabel('signal amplitude')
    plt.legend()
    plt.show()

def check_filter_stability(
    b_driver, a_driver, float_type, max_delay_in_sample_number,
    error_margin = 0.005
):
    # We check filter can be computed without overflow with 32 Bit float.
    mean = 1.0
    deviation = 0.01
    N_sample = 100000
    sig = np.random.normal(mean, deviation, N_sample)
    sig_32 = float_type( sig )
    b_32 = float_type(b_driver)
    a_32 = float_type(a_driver)
    fil_sig = [ float_type(0.0) for i in range(len(a_32)) ]
    for i in range( len(sig) ):
        fil_sig.append(
            sum([ b_32[j] * sig_32[i-j]  for j in range(len(b_32)) ])
            -
            sum([ a_32[j] * fil_sig[-j] for j in range(1, len(a_32)) ])
        )
    plt.plot( sig, label="input signal" )
    plt.plot( fil_sig[len(a_32):], label="filtered signal" )
    plt.legend()
    plt.title(f"Test filter using a constant input signal with a gaussian noise (N({mean},{deviation}))")
    plt.xlabel('time [s]')
    plt.ylabel('signal amplitude')
    plt.show()

    mean_filter = np.average( fil_sig[len(a_32)+3*max_delay_in_sample_number:] )
    error_filter = np.abs(mean - mean_filter)

    
    print("Test filter with a constant input signal with a gaussian noise : ")
    print("")
    print("   mean : " + str(mean))
    print("   filtered mean : " + str(mean_filter))
    print("   error : " + str(error_filter))

    check(
        error_filter < error_margin, 
        "Be carfull, the filter is not Stable or have a bias due to float approximations !" +
        "You need to reduce the order of the filter or to increase " + 
        "the number of bits of the real type." 
    )



"""
filter_name in ['chebychev2', 'butterworth']
"""
def design_digital_filter(
    order_digital, wanted_attenuation_at_fsys_over_2, max_signal_freq, sample_freq, 
    output_frequency_signal, float_type=np.float32, filter_name='butterworth',
    low_attenuation = 3*dB
):
    print(
    """
##########################################################################################
# Digital filter design
##########################################################################################

"""
    )
    print(f"float type : {float_type}")
    print("")

    digital_filter = compute_digital_filter(
        order_digital, max_signal_freq, output_frequency_signal, sample_freq, filter_name,
        wanted_attenuation_at_fsys_over_2, float_type, low_attenuation
    )
    b_driver = digital_filter['b']
    a_driver = digital_filter['a']
    w = digital_filter['w']
    h = digital_filter['h']

    # The Bode diagram
    plt.plot(w, 20*np.log10(np.abs(h)))
    #plt.plot(w,np.abs(h))
    plt.title('Digital filter frequency response')

    # Some usefull vertical lines
    plt.axvline(x=max_signal_freq, color="red", label='maximal signal frequency')
    plt.axvline(x=output_frequency_signal/2, color="orange", label='output frequency/2')
    plt.axvline(x=sample_freq/2, color="pink", label='sample frequency/2')
    plt.axvline(x=digital_filter['f_cut'], linestyle="dashed", color="green", label='filter :cut frequency')
    plt.axvline(x=digital_filter['f_cut_high_attenuation'], linestyle="dashed", color="magenta", label='cut frequency at hight attenuation')
    plt.xscale('log')

    # Some usefull horizontal lines
    plt.axhline(
        y=20*np.log10(.9), color="0.7", linestyle="dashed",
        label='attenuation X0.9'
    )
    plt.axhline(
        y=-3, color="0.5", linestyle="dashed",
        label='attenuation X0.7 (-3db)'
    )
    plt.axhline(
        y=20*np.log10(.1), color="0.35", linestyle="dashed",
        label='attenuation X0.1 (-20db)'
    )
    plt.axhline(
        y=20*np.log10(.01), color="black", linestyle="dashed",
        label='attenuation X0.01 (-40db)'
    )
    plt.axhline(
        y=20*np.log10(.001), color="black", linestyle="dashed",
        label='attenuation X0.001 (-60db)'
    )
    plt.legend()
    min_y = -60
    if -wanted_attenuation_at_fsys_over_2 < -60:
        min_y = -int(wanted_attenuation_at_fsys_over_2//10)*10 - 20
    plt.ylim((min_y, 10))

    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Amplitude [dB]')

    plt.show()

    phases = compute_phases(np.angle(h))
    delays = phases[1:]/(2*np.pi*w[1:])
    plt.axvline(x=max_signal_freq, color="red", linestyle="dashed", label='Signal frequency')
    plt.plot(w[1:], delays) #, linestyle='', marker="." )
    plt.xlabel("Frequence f (Hz)")
    plt.ylabel("filter delay (s)")
    plt.title("Filter delays of the digital filter")
    plt.legend()
    plt.show()

    delays_for_joystick = delays * (w[1:]<=max_signal_freq)
    max_delay_for_joystick = np.max(-delays_for_joystick)
    max_delay_for_joystick_in_sample = int(max_delay_for_joystick * sample_freq)
    print("")
    print(
        "The maximum delay for the filtered signal is : %s seconds (%s Hz) ~= %s samples at %s Hz"%(
            max_delay_for_joystick, 1/max_delay_for_joystick, max_delay_for_joystick_in_sample, sample_freq
        )
    )

    b_string = ""
    for val in b_driver:
        b_string += repr(val) + ", "

    a_string = ""
    for val in a_driver:
        a_string += repr(val) + ", "

    print("")
    print("The filter order is : " + str(order_digital))

    print("")
    print("The filter coefficients are :")
    print("")
    print(f"b : {b_driver}")
    print(f"a : {a_driver}")
    print("")
    print(
    """
We recall that :

a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[M]*x[n-M]
                      - a[1]*y[n-1] - ... - a[N]*y[n-N]

where X is the input signal and Y is the filtered signal.
"""
    )

    test_fiter(b_driver, a_driver, max_signal_freq, sample_freq)
    check_filter_stability(b_driver, a_driver, float_type, max_delay_for_joystick_in_sample)

    if float_type is np.float32 :
        type_string = "float"
    elif float_type is np.float64 :
        type_string = "double"
    else:
        raise ValueError("float type not supported.")

    print("The C++ code is :")
    if filter_name == 'chebychev2':
        print(
        f"""
#include "filter.hpp"
typedef {type_string} REAL_TYPE;
const unsigned int sample_frequence = {sample_freq}; // In Hz
const unsigned int output_frequence = {output_frequency_signal}; // In hz
const unsigned int decimation_factor = sample_frequence / output_frequence;
// Chebychev2 filter :
//    Parameter :
//        order : {order_digital}
//        cut frequency : {digital_filter['f_cut']} Hz,
//        stop band attenuation : {digital_filter['stop_band_attenuation']} dB,
//        sample frequency : {sample_freq/kHz} kHz
//    Analysis :
//        Maximal expected delay : {max_delay_for_joystick/ms:.4} ms
const int order = {order_digital};
RII_filter<REAL_TYPE, order> filter{'{'}
    {'{'}{b_string}{'}'}, // b : Numerator
    {'{'}{a_string}{'}'} // a : Denominator
{'}'};
"""
        )
    elif filter_name == 'butterworth':
        print(
        f"""
#include "filter.hpp"
typedef {type_string} REAL_TYPE;
const unsigned int sample_frequence = {sample_freq}; // In Hz
const unsigned int output_frequence = {output_frequency_signal}; // In hz
const unsigned int decimation_factor = sample_frequence / output_frequence;
// Butterworth filter :
//    Parameter :
//        order : {order_digital}
//        cut frequency : {digital_filter['f_cut']} Hz,
//        sample frequency : {sample_freq/kHz} kHz
//    Analisys :
//        Maximal expected delay : {max_delay_for_joystick/ms:.4} ms
const int order = {order_digital};
RII_filter<REAL_TYPE, order> filter{'{'}
    {'{'}{b_string}{'}'}, // b : Numerator
    {'{'}{a_string}{'}'} // a : Denominator
{'}'};
    """
        )
    digital_filter['max_delay'] = max_delay_for_joystick
    return digital_filter 


