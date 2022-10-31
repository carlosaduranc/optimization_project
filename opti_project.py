from matplotlib import pyplot as plt
from pandas import *
import numpy as np

# Model identification and heating strategy optimization for a small, single zone RC-lumped model
# The room in question is a single room with a single window facing south (maximum sun exposure throughout the day)
# The room is equipped with a single radiator with a max power output of 1kW which can be controlled with a valve


def create_model_onoff(gA, C, R, plot=False):
    data = read_csv("leuven_october2022_16-22.csv")

    # Time dependent parameters
    time = data['time'].tolist()
    temp = data['temp'].tolist()  # deg C
    for i in range(len(temp)):
        temp[i] = temp[i] + 273.15  # deg K

    Qg = np.zeros(len(time))
    Qh = np.zeros(len(time))
    for i in range(5):
        Qg[0 + 24 * (i + 1):7 + 24 * (i + 1)] = 100  # W. Human heat output
        Qg[18 + 24 * (i + 1):24 + 24 * (i + 1)] = 100  # W. For one person

        Qh[0 + 24 * (i + 1):6 + 24 * (i + 1)] = 1000  # W. Heater
        Qh[19 + 24 * (i + 1):24 + 24 * (i + 1)] = 1000  # W. Max power on/off strategy

    Qg[0:24] = 100
    Qg[144:] = 100
    Qh[0:24] = 1000
    Qh[144:] = 1000

    Qsun = data['solrad'].tolist()  # W/m2

    Tz = np.zeros(len(time))
    Tz[0] = 293.15  # initial zone temp

    delta_t = 3600  # s

    # calculating the temperature
    for i in range(len(Tz) - 1):
        Tz[i + 1] = delta_t * ((Qsun[i] * gA + Qh[i] + Qg[i]) / C + (temp[i] - Tz[i]) / (R * C)) + Tz[i]

    # calculating the running cost for energy use
    cost = np.zeros(len(time))
    rate = 0.410  # euro/kWh in Belgium

    for i in range(len(time) - 1):
        cost[i + 1] = cost[i] + rate * Qg[i] * 0.001

    if plot:
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7.5, 7.5))
        ax1.plot(time, temp, label=r'$T_{amb}$')
        ax1.plot(time, Tz, label=r'$T_{z}$')
        for i in [293.15, 298.15]:
            ax1.hlines(i, xmin=time[0], xmax=time[-1], lw=0.7, ls='-', color='k', zorder=1)
        ax2.plot(time, Qsun, label=r'$Q_{sun} [W/m^2]$')
        ax2.plot(time, Qg, label=r'$Q_g [W]$')
        ax2.plot(time, Qh, label=r'$Q_h [W]$')
        ax3.plot(time, cost)
        ax3.hlines(cost[-1], xmin=time[0], xmax=time[-1], lw=0.7, ls='--', color='b', zorder=1)

        for i in range(6):
            ax1.vlines(time[24 * (i + 1)], ymin=-100, ymax=1000, lw=0.7, ls='--', color='k', zorder=1)
            ax2.vlines(time[24 * (i + 1)], ymin=-100, ymax=2000, lw=0.7, ls='--', color='k', zorder=1)
            ax3.vlines(time[24 * (i + 1)], ymin=-100, ymax=2000, lw=0.7, ls='--', color='k', zorder=1)

        ax1.set_xticks([])
        ax2.set_xticks([])

        ax1.set_xlim([time[0], time[-1]])
        ax2.set_xlim([time[0], time[-1]])
        ax3.set_xlim([time[0], time[-1]])

        ax1.set_ylim([280, 305])
        ax2.set_ylim([0, 1300])
        ax3.set_ylim([0, 5])

        ax1.set_ylabel('Temperature [K]')
        ax2.set_ylabel('Heat')
        ax3.set_ylabel('Running cost [Euro]')

        ax3.set_xticks([time[24 * (i + 1)] for i in range(6)])
        ax3.set_xticklabels(['17/10', '18/10', '19/10', '20/10', '21/10', '22/10'])

        ax1.legend(loc='upper right')
        ax2.legend(loc='upper right')

        ax1.set_title("traditional operation for one radiator")

        plt.show()

    return Tz, time


def compare_models(Tz_a, Tz_b, time, plot=False):
    # calculating l-2 norm for each
    e = np.zeros(len(Tz_a))
    for i in range(len(Tz_a)):
        e[i] = 1 / 2 * (Tz_a[i] - Tz_b[i]) ** 2

    # calculating RMSE
    rmse = 1/len(e) * sum(e)

    if plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 7.5))

        ax1.plot(time, Tz_a, label='Actual')
        ax1.plot(time, Tz_b, label='Model')
        ax2.plot(time, e, label='Error')
        ax2.hlines(rmse, xmin=time[0], xmax=time[-1], lw=1, ls='-', color='k', zorder=1, label='RMSE')

        for i in range(6):
            ax1.vlines(time[24 * (i + 1)], ymin=-100, ymax=1000, lw=0.7, ls='--', color='k', zorder=1)
            ax2.vlines(time[24 * (i + 1)], ymin=-100, ymax=2000, lw=0.7, ls='--', color='k', zorder=1)

        ax1.set_xticks([])

        ax1.set_xlim([time[0], time[-1]])
        ax2.set_xlim([time[0], time[-1]])

        ax1.set_ylim([285, 315])
        ax2.set_ylim([0, 45])

        ax1.set_ylabel('Temperature [K]')
        ax2.set_ylabel('Error [$K^2$]')

        ax2.set_xticks([time[24 * (i + 1)] for i in range(6)])
        ax2.set_xticklabels(['17/10', '18/10', '19/10', '20/10', '21/10', '22/10'])

        ax1.legend(loc='upper right')
        ax2.legend(loc='upper right')

        ax1.set_title("actual vs model behavior")

        plt.show()

    return e, time


# Assume correct values
gA_true = 0.45 * 2.10  # glazing factor * window area. m2
C_true = 5000000  # J/K
R_true = 0.01  # K/W

# Guesses for the parameters
gA_guess = 0.1 * 0.9  # m2
C_guess = 1000000  # J/K
R_guess = 0.015  # K/W

# Loading true zone temperature and model zone temperature
[Tz_true, time_true] = create_model_onoff(gA=gA_true, C=C_true, R=R_true, plot=True)
[Tz_guess, time_guess] = create_model_onoff(gA=gA_guess, C=C_guess, R=R_guess, plot=False)

# Comparing models. Extracting error function to be minimized
[error, time_error] = compare_models(Tz_a=Tz_true, Tz_b=Tz_guess, time=time_true, plot=True)

# Minimizing error function
print("HelloWorld")