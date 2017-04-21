import numpy as np
import matplotlib.pyplot as plt

def Plot():
    # To simplify loog at Plot Once
    filename = "/Users/clementine.fourrier/build/AwesomeNetwork"

    x_vec        =np.array([])
    noise        =np.array([])
    g_vec        =np.array([])
    tau_var_vec  =np.array([])
    tau_mean_vec =np.array([])
    ksi_var_vec  =np.array([])
    ksi_mean_vec =np.array([])

    for line in open(filename, 'r'):
        line_vec = line.split()
        x_vec = np.append(x_vec, line_vec[0])
        noise = np.append(noise, line_vec[1])
        g_vec = np.append(g_vec, line_vec[2])
        tau_mean_vec = np.append(tau_mean_vec, line_vec[3])
        tau_var_vec = np.append(tau_var_vec, line_vec[4])
        ksi_mean_vec = np.append(ksi_mean_vec, line_vec[5])
        ksi_var_vec = np.append(ksi_var_vec, line_vec[6])

    plt.figure(1)
    plt.subplot(211)
    plt.plot(x_vec, noise, 'r')
    plt.title("Noise")

    plt.subplot(212)
    plt.plot(x_vec, g_vec, 'r')
    plt.title("G")

    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.10, right=0.95, hspace=0.25,
                        wspace=0.35)
    plt.figure(2)
    plt.subplot(221)
    plt.plot(x_vec, tau_mean_vec, 'b')
    plt.title("Mean Tau")

    plt.subplot(222)
    plt.plot(x_vec, ksi_mean_vec, 'g')
    plt.title("Mean Ksi")

    plt.subplot(223)
    plt.plot(x_vec, tau_var_vec, 'b')
    plt.title("Variance Tau")

    plt.subplot(224)
    plt.plot(x_vec, ksi_var_vec, 'g')
    plt.title("Variance Ksi")

    plt.subplots_adjust(top=0.95, bottom=0.05, left=0.10, right=0.95, hspace=0.25,
                        wspace=0.35)

    plt.show()

    print "end"