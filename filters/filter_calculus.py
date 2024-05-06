from scipy import signal
import numpy as np
from matplotlib import pyplot as plt
import sys


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


def test_fiter(b_driver, a_driver, f_max, f_sample):
    f_signal_1 = f_max * Hz
    f_signal_2 = f_max*10 * Hz
    f_noise_1 = 20*f_max * Hz
    f_noise_2 = 50*f_max * Hz
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
    plt.figure()
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

    if not(error_filter < error_margin):
        print(
            "Be carfull, the filter is not Stable or have a bias due to float approximations !" +
            "You need to reduce the order of the filter or to increase " + 
            "the number of bits of the real type." )
    



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


resolution = 2**14
print(
"""
################################################################################
# Information
################################################################################
"""
)

max_signal_freq =  0.01* Hz # Useless in our case 

sample_freq =  1800 * Hz
niquist_sample_freq = sample_freq/2

print(f"Signal frequency contain : {max_signal_freq} Hz")
print(f"Sample frequence : {sample_freq/kHz} kHz")

print(
"""
################################################################################
# Digital filter
################################################################################
"""
)
float_type = np.float64
#float_type = np.float32
output_frequency_signal = sample_freq # we send to computer signal at 100 Hz
niquist_output_frequency_signal = output_frequency_signal/2
digital_cut_freq = .5* Hz
assert(digital_cut_freq <= niquist_output_frequency_signal )
# assert(max_signal_freq <= digital_cut_freq)

print(f"float type : {float_type}")
print("")

# Design of a butterworth filter :
# Order of the filter 
order_digital = 4

b_driver, a_driver = signal.butter(
    order_digital, digital_cut_freq, analog=False, output='ba', fs=sample_freq
)

# stop_band_attenuation = 10 # dB
# b_driver, a_driver = float_type(
#    signal.cheby2(
#        order_digital, stop_band_attenuation, .36*Hz, #digital_cut_freq,
#        analog=False, output='ba', fs=sample_freq
#    )
# )


# rp_attenuation = .1 # dB
# b_driver, a_driver = float_type(
#    signal.cheby1(
#        order_digital, rp_attenuation, 1*Hz, #digital_cut_freq,
#        analog=False, output='ba', fs=sample_freq
#    )
# )

# We draw the bode diagram ot the filter
w, h = signal.freqz(b_driver, a_driver, fs=sample_freq, worN = resolution)
plt.figure(0)
# The Bode diagram
plt.plot(w, 20*np.log10(np.abs(h)))
#plt.plot(w,np.abs(h))
plt.title('Digital filter frequency response')

# Some usefull vertical lines
plt.axvline(x=max_signal_freq, color="red", label='signal (Joystick) frequency')
plt.axvline(x=digital_cut_freq, color="green", label='cut frequency')
plt.axvline(x=niquist_output_frequency_signal, color="orange", label='output frequency/2')
plt.axvline(x=niquist_sample_freq, color="pink", label='sample frequency/2')
plt.xscale('log')

# Some usefull horizontal lines
plt.axhline(
    y=20*np.log10(.9), color="0.7", linestyle="dashed",
    label='attenuation X0.9'
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
plt.ylim((-80, 10))

plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude [dB]')

plt.figure(1)

#phases = np.angle(h) - 2*np.pi*(np.angle(h)>0)
phases = compute_phases(np.angle(h))
delays = phases[1:]/(2*np.pi*w[1:])
plt.axvline(x=max_signal_freq, color="red", linestyle="dashed", label='Signal frequency')
plt.plot(w[1:], delays) #, linestyle='', marker="." )
plt.xlabel("Frequence f (Hz)")
plt.ylabel("filter delay (s)")
plt.title("Filter delays of the filter")
plt.legend()

delays_for_joystick = delays * (w[1:]<=max_signal_freq)
max_delay_for_joystick = np.max(-delays_for_joystick)
max_delay_for_joystick_in_sample = int(max_delay_for_joystick * sample_freq)
print("")
print(
    "The maximum delay for a filtered signal is : %s seconds = %s samples (%s Hz)"%(
        max_delay_for_joystick, max_delay_for_joystick_in_sample, 1/max_delay_for_joystick
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
print(f"b : [{b_string}]")
print(f"a : [{a_string}]")
print("")
print(
"""
We recall that :

a[0]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ... + b[M]*x[n-M]
                      - a[1]*y[n-1] - ... - a[N]*y[n-N]

where X is the input signal and Y is the filtered signal.
"""
)

if float_type is np.float32 :
    type_string = "float"
elif float_type is np.float64 :
    type_string = "double"
else:
    raise ValueError("float type not supported.")
 
print("The C++ code is :")
'''
print(
f"""
// Chebychev2 filter :
//    Parameter :
//        cut frequency : {digital_cut_freq} Hz,
//        stop band attenuation : {stop_band_attenuation} dB,
//        sample frequency : {sample_freq/kHz} kHz
//    Analisys :
//        Maximal expected delay : {max_delay_for_joystick/ms:.4} ms
const int order = {order_digital};
RII_filter<{type_string}, order> filter{'{'}
    {'{'}{b_string}{'}'}, // b : Numerator
    {'{'}{a_string}{'}'} // a : Denominator
{'}'};
"""
)
'''
print(
f"""
// Butterworth filter :
//    Parameter :
//        cut frequency : {digital_cut_freq} Hz,
//        sample frequency : {sample_freq/kHz} kHz
//    Analisys :
//        Maximal expected delay : {max_delay_for_joystick/ms:.4} ms
const int order = {order_digital};
RII_filter<FLOATING_TYPE, order> filter{'{'}
    {'{'}{b_string}{'}'}, // b : Numerator
    {'{'}{a_string}{'}'} // a : Denominator
{'}'};
"""
)

if (len(sys.argv) > 1 and sys.argv[1] == "y"):
    with open("arduino/src/filter_in_use.h", "w") as f:
        f.write(f"""#define FLOATING_TYPE double \
                \
                \n\n#define SAMPLE_FREQUENCY 1800 \
                \n#define NB_MOTOR_PERIODS_PER_REVOLUTION 400 \
                \n\n#define FILTER_ORDER {order_digital}\
                \n#define FILTER_NUMERATOR """ + "{" + f"{b_string}" + "}" +  
                """\n#define FILTER_DENOMINATOR """ + "{" + f"{a_string}" + "}\n" )
else:

    check_filter_stability(
        b_driver, a_driver, float_type, max_delay_for_joystick_in_sample,
        error_margin = 0.005
    )

    test_fiter(b_driver, a_driver, digital_cut_freq, sample_freq)

    plt.show()
