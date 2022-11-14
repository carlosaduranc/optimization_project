from matplotlib import pyplot as plt
from pandas import *
from casadi import *


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
        cost[i + 1] = cost[i] + rate * Qh[i] * 0.001

    if plot:
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7.5, 7.5), sharex=True)
        ax1.plot(time, temp, label=r'$T_{a}$')
        ax1.plot(time, Tz, label=r'$T_{z}$')
        for i in [293.15, 298.15]:
            ax1.hlines(i, xmin=time[0], xmax=time[-1], lw=0.7, ls='-', color='k', zorder=1)
        ax2.plot(time, Qsun, label=r'$\dot{Q}_{sun} [W/m^2]$')
        ax2.plot(time, Qg, label=r'$\dot{Q}_{g} [W]$')
        ax2.plot(time, Qh, label=r'$\dot{Q}_{h} [W]$')
        ax3.plot(time, cost)
        ax3.hlines(cost[-1], xmin=time[0], xmax=time[-1], lw=0.7, ls='--', color='b', zorder=1)
        ax3.text(x=time[1], y=cost[-1]+1.02, s=str(round(cost[-1], 2)) + '€')

        for i in range(6):
            ax1.vlines(time[24 * (i + 1)], ymin=-100, ymax=1000, lw=0.7, ls='--', color='k', zorder=1)
            ax2.vlines(time[24 * (i + 1)], ymin=-100, ymax=2000, lw=0.7, ls='--', color='k', zorder=1)
            ax3.vlines(time[24 * (i + 1)], ymin=-100, ymax=2000, lw=0.7, ls='--', color='k', zorder=1)
        for i in range(5):
            ax1.fill_between(time[7 + 24 * (i + 1):18 + 24 * (i + 1)], -100, 5000, color='c', alpha=0.15)
            ax2.fill_between(time[7 + 24 * (i + 1):18 + 24 * (i + 1)], -100, 5000, color='c', alpha=0.15)
            ax3.fill_between(time[7 + 24 * (i + 1):18 + 24 * (i + 1)], -100, 5000, color='c', alpha=0.15)

        ax1.set_xticks([])
        ax2.set_xticks([])

        ax1.set_xlim([time[0], time[-1]])
        ax2.set_xlim([time[0], time[-1]])
        ax3.set_xlim([time[0], time[-1]])

        ax1.set_ylim([280, 305])
        ax2.set_ylim([0, 1300])
        ax3.set_ylim([0, 50])

        ax1.set_ylabel('Temperature [K]')
        ax2.set_ylabel('Heat')
        ax3.set_ylabel(r'Running cost [€]')

        ax3.set_xticks([time[24 * (i + 1)] for i in range(6)])
        ax3.set_xticklabels(['17/10', '18/10', '19/10', '20/10', '21/10', '22/10'])

        ax1.legend(loc='upper right')
        ax2.legend(loc='upper right')

        plt.show()

    return Tz, time


def compare_models(Tz_a, Tz_b, time, plot=False):
    # calculating l-2 norm for each
    e = np.zeros(len(Tz_a))
    for i in range(len(Tz_a)):
        e[i] = 1 / 2 * (Tz_a[i] - Tz_b[i]) ** 2

    # calculating RMSE
    rmse = 1 / len(e) * sum(e)

    if plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 5), sharex=True)

        ax1.plot(time, Tz_a, label='Actual')
        ax1.plot(time, Tz_b, label='Model')
        ax2.plot(time, e, label='Error')
        ax2.hlines(rmse, xmin=time[0], xmax=time[-1], lw=1, ls='-', color='k', zorder=1, label='RMSE')
        ax2.text(x=time[1], y=rmse*1.01, s='RMSE: ' + str(format(rmse, '.2E')) + r'$K^2$')

        for i in range(6):
            ax1.vlines(time[24 * (i + 1)], ymin=-100, ymax=1000, lw=0.7, ls='--', color='k', zorder=1)
            ax2.vlines(time[24 * (i + 1)], ymin=-100, ymax=50000, lw=0.7, ls='--', color='k', zorder=1)

        ax1.set_xticks([])

        ax1.set_xlim([time[0], time[-1]])
        ax2.set_xlim([time[0], time[-1]])

        ax1.set_ylim([285, max(Tz_b)*1.01])
        ax2.set_ylim([0, max(e)*1.01])

        ax1.set_ylabel('Temperature [K]')
        ax2.set_ylabel('Error [$K^2$]')

        ax2.set_xticks([time[24 * (i + 1)] for i in range(6)])
        ax2.set_xticklabels(['17/10', '18/10', '19/10', '20/10', '21/10', '22/10'])

        ax1.legend(loc='upper right')
        ax2.legend(loc='upper right')

        ax1.set_title("actual vs model behavior")

        plt.show()

    return e, time


def minimize_function(Tz_a):
    opti = Opti()
    N = len(Tz_a)

    delta_t = 3600  # s

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

    # Initiating optimization variables
    R = opti.variable()
    C = opti.variable()
    gA = opti.variable()

    #  parameter vector
    p = vertcat(R, C, gA)

    # Initial temperature guess
    Tz = 293.15

    # Error function: initial value
    f = 0

    for i in range(N):
        f = f + (Tz - Tz_a[i]) ** 2

        Tz_next = delta_t * ((Qsun[i] * gA + Qh[i] + Qg[i]) / C + (temp[i] - Tz) / (R * C)) + Tz
        Tz = Tz_next

    f = f + (Tz - Tz_a[N - 1]) ** 2

    opti.minimize(f)

    opti.solver('ipopt')

    # filling initial guess parameter vector
    p_hat = vertcat(R_guess, C_guess, gA_guess)

    opti.set_initial(p, p_hat)

    sol = opti.solve()
    return sol.value(p)

# Assume correct values
gA_true = 0.45 * 2.10  # glazing factor * window area. m2
C_true = 5000000  # J/K
R_true = 0.01  # K/W

# Guesses for the parameters
gA_guess = 0.1 * 1  # m2
C_guess = 1000000  # J/K
R_guess = 1  # K/W

# Loading true zone temperature and model zone temperature
[Tz_true, time_true] = create_model_onoff(gA=gA_true, C=C_true, R=R_true, plot=False)
[Tz_guess, time_guess] = create_model_onoff(gA=gA_guess, C=C_guess, R=R_guess, plot=False)

# Comparing models. Extracting error function to be minimized
[error, time_error] = compare_models(Tz_a=Tz_true, Tz_b=Tz_guess, time=time_true, plot=False)

# Minimizing error function. Output optimal parameters
[R_opt, C_opt, gA_opt] = minimize_function(Tz_a=Tz_true)
print(R_opt, C_opt, gA_opt)

# Loading model with calculated variables and comparing to true model
[Tz_opt, time_opt] = create_model_onoff(gA=gA_opt, C=C_opt, R=R_opt, plot=False)
[error_opt, time_error_opt] = compare_models(Tz_a=Tz_true, Tz_b=Tz_opt, time=time_true, plot=True)

print('\nRMSE: ', 1 / len(error_opt) * np.sum(error_opt), '\n\n\n\n')

a = np.array(np.array([[0, 5.58e6], [1, 5.55e6], [2, 5.45e6], [3, 5.09e6], [4, 7.29e5], [5, 7.33e4], [6, 1.12e4],
                       [7, 1.89e3], [8, 5.75e2], [9, 1.55e2], [10, 2.58e1],[11, 1.32e0],[12, 7.15e-3],[13, 2.18e-3],
                       [14, 2.18e-3]]))

fig, ax1 = plt.subplots(1, 1, figsize=(7.5, 4))
ax1.plot(a[:, 0], a[:, 1], label='Error')
ax1.set_yscale('log')
ax1.set_xlabel('Iteration')
ax1.set_ylabel('Error [$K^2$]')
ax1.grid()
plt.show()