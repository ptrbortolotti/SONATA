#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 14:49:51 2018

@author: gu32kij
"""
# Third party modules
import matplotlib.pyplot as plt
import numpy as np

# First party modules

# from matplotlib2tikz import save as tikz_save


def plot_histogram_2Ddata(data, **kwargs):
    """
    simple plotting procedure to illustrate multiple datasets of a 2D array in 
    in a histogram.
    
    Parameters
    ----------
    data : ndarray
        array.shape = (sample Nb, i, j)
        
    ref : ndarray, optional
    upper_tri : bool
    title : str
    ylabel : list of str
        
    """

    # Define function Args and Defaults for kwargs
    shape = data.shape[1:]
    ref = None
    upper_tri = True
    title = "No Title"
    ylabel = [[("k%s%s" % (i + 1, j + 1)) for j in range(shape[0])] for i in range(shape[1])]

    if "ref" in kwargs:
        if isinstance(kwargs["ref"], (np.ndarray)) and kwargs["ref"].shape == data.shape[1:]:
            ref = kwargs["ref"]
        else:
            print("ref must provide a ndarray with shape %s" % str(data.shape[1:]))

    if "upper_tri" in kwargs:
        if type(kwargs["upper_tri"]) == bool:
            upper_tri = kwargs["upper_tri"]
        else:
            print("upper_tri must provide a bool value")

    if "title" in kwargs:
        if type(kwargs["title"]) == str:
            title = kwargs["title"]
        else:
            print("title must provide as string")

    if "ylabel" in kwargs:
        # list of strings
        ylabel = kwargs["ylabel"]

    # init figure
    shape = data.shape[1:]
    fig, ax = plt.subplots(shape[0], shape[1], sharex=False)
    fig.suptitle(title, fontsize=14)
    fig.subplots_adjust(wspace=0.3, hspace=0.3)

    # upper triangle
    if upper_tri:
        ut = np.zeros((shape[0], shape[1]))
        ut[np.triu_indices(shape[0])] = 1
    else:
        ut = np.ones((shape[0], shape[1]))

    for i in range(shape[0]):
        for j in range(shape[1]):
            print(i, j)
            # if ut[i,j] == 1 and ref[i,j] != 0:
            if ut[i, j] == 1:
                if isinstance(ref, (np.ndarray)):
                    if ref[i, j] != 0:
                        arr = (data[:, i, j] - ref[i, j]) / ref[i, j] * 100
                        sigma = arr.std()
                        mu = arr.mean()
                        string = r"$\sigma$ = %.2f $\%%$" % sigma
                        ax[i][j].text(mu, 0, string)
                        ax[i][j].set_xlim(-20, 20)

                    else:
                        arr = data[:, i, j]

                    count, bins, ignored = ax[i][j].hist(arr, 30, density=True, alpha=0.5, label="histogram")
                    # print(count,bins,ignored)
                    sigma = arr.std()
                    mu = arr.mean()
                    if sigma > 0:
                        ax[i][j].plot(bins, 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-((bins - mu) ** 2) / (2 * sigma ** 2)), linewidth=2, color="r", linestyle="--", alpha=0.5, label="ML approx.")
                else:
                    arr = data[:, i, j]
                    count, bins, ignored = ax[i][j].hist(arr, 30, density=True, alpha=0.5, label="histogram")

            elif ut[i, j] == 0:
                fig.delaxes(ax[i][j])

            if i == j:
                ax[i][j].set_xlabel(r"$\%%$ deviation from baseline value")
                if j == shape[1] - 1:
                    ax[i][j].legend()

            ax[i][j].set_ylabel(ylabel[i][j])

    return None


def plot_2dhist(data, labels, ref=None, title=None, xlim=(-20, 20), bins=21, upper_tri=False, ptype="hist", MLE=False, **kwargs):
    """
    
    
    """

    if len(data.shape) == 2:
        data = np.expand_dims(data, 2)
        labels = np.expand_dims(labels, 1)
    imax = data.shape[1]
    jmax = data.shape[2]

    fig, ax = plt.subplots(imax, jmax)
    fig.suptitle(title, fontsize=16)
    fig.subplots_adjust(wspace=0.3, hspace=0.3)

    # upper triangle
    if upper_tri:
        ut = np.zeros((imax, jmax))
        ut[np.triu_indices(imax)] = 1
    else:
        ut = np.ones((imax, jmax))

    for i in range(imax):
        for j in range(jmax):
            if imax == 1 and jmax == 1:
                axh = ax
            elif imax > 1 and jmax == 1:
                axh = ax[i]
            else:
                axh = ax[i][j]

            if ut[i, j] == 1:
                arr = data[:, i, j]

                if ptype == "hist":
                    if ref is not None:
                        arr = (data[:, i, j] - ref[i, j]) / ref[i, j] * 100
                        axh.set_xlabel((labels[i, j] + r" (\% deviation from baseline value)"))
                        axh.axvline(c="k")
                    else:
                        arr = (data[:, i, j] - np.mean(data[:, i, j], axis=0)) / np.mean(data[:, i, j], axis=0) * 100
                        axh.set_xlabel((labels[i, j] + r" (\% deviation from mean value)"))

                    count, bins, ignored = axh.hist(arr, density=True, color="black", rwidth=0.8)
                    sigma = arr.std()
                    mu = arr.mean()
                    if MLE and sigma > 0:
                        string = r"MLE $\sigma$ = %.1f $\%%$" % sigma
                        axh.plot(bins, 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-((bins - mu) ** 2) / (2 * sigma ** 2)), linewidth=2, color="r", linestyle="--", label=string)
                        axh.legend()
                        # axh.text(mu,0,string, color = 'r')
                    ylim = np.asarray(axh.get_ylim())
                    # print(ylim*1)
                    axh.set_ylim(ylim * 1.2)
                    axh.set_xlim(xlim)
                    axh.set_ylabel("Probability of occurance")

                elif ptype == "scatter":
                    stdvs = []
                    x = np.arange(2, arr.shape[0])
                    for r in x:
                        stdvs.append(np.std(arr[:r], axis=0))
                    stdvs = np.asarray(stdvs)
                    axh.scatter(x, stdvs, s=0.5, c="k", marker=".")
                    axh.set_xlabel(r"samples")
                    # print((r'SD of '+labels[i,j]))
                    axh.set_ylabel((r"SD of " + labels[i, j]))

            elif ut[i, j] == 0:
                # fig.delaxes(ax[i][j])
                pass
    plt.legend()
    plt.show()


def plot_beam_properties(data, sigma=None, ref=None, x_offset=0, description=True):
    """
    generates a plot of the beamproperties using the 
    
    Parameters
    ----------
    Data : np.ndarray 
        Massterms(6), Stiffness(21), damping(1), and coordinate(1) stacked 
        horizontally for the use in DYMORE/PYMORE/MARC  
    sigma: np.ndarray
        standart deviation of Data
    
    ref: np.ndarray
        of reference to select = 'all', 'Massterms', 'Stiffness'
    
    x_offset : float 
        x offset if curvilinear coordinate doen't start at the rotational 
        center of the rotordefault UH-60A value = 0.8178698
            
    """

    # Inertia Properties:
    ylabel = [(r"m_{00}", "kg/m"), (r"X_{m2}", "m"), (r"X_{m3}", "m"), (r"m_{33}", "kg m"), (r"m_{23}", "kg m"), (r"m_{22}", "kg m")]
    # x = data[:,-1] + x_offset
    x = data[:, -1] + x_offset
    x = x / x[-1]
    data = np.copy(data)
    data[:, 1] = data[:, 1] / data[:, 0]
    data[:, 2] = data[:, 2] / data[:, 0]
    if isinstance(ref, (np.ndarray)):
        ref = np.copy(ref)
        ref_x = ref[:, -1] + x_offset
        ref_x = ref_x / ref_x[-1]
        ref[:, 1] = ref[:, 1] / ref[:, 0]
        ref[:, 2] = ref[:, 2] / ref[:, 0]
    fig1, ax1 = plt.subplots(3, 2, sharex=True)
    fig1.subplots_adjust(wspace=0.3, hspace=0.3)
    # fig1.suptitle('Inertia Properties', fontsize=14)

    c = 0
    for i in range(3):
        for j in range(2):
            if isinstance(ref, (np.ndarray)):
                (ref_line,) = ax1[i][j].plot(ref_x, ref[:, c], "-.b")
            if isinstance(sigma, (np.ndarray)):
                ax1[i][j].fill_between(
                    x, data[:, c] - sigma[:, c], data[:, c] + sigma[:, c], alpha=0.6, edgecolor="r", linestyle=":", facecolor="r", antialiased=True,
                )
                # ax1[i][j].errorbar(x, data[:,c], yerr=sigma[:,c])
            (line,) = ax1[i][j].plot(x, data[:, c], "--k.", markersize=2)
            ax1[i][j].ticklabel_format(axis="y", style="sci")
            tmp = r"$%s \quad [%s]$" % (ylabel[c][0], ylabel[c][1])
            ax1[i][j].set_ylabel(tmp)

            if i == 2:
                ax1[i][j].set_xlabel(r"radius $\quad [1/R]$")

            c += 1

    # ax1[2,1].legend([ref_line, line], ['UH-60A Reference', 'Inertial Properties'])
    ax1[2, 1].set_ylim(0, 0.15)

    # 6x6 Stiffness Properties:
    TSunits = np.array(
        [
            ["N", "N", "N", "Nm", "Nm", "Nm"],
            ["N", "N", "N", "Nm", "Nm", "Nm"],
            ["N", "N", "N", "Nm", "Nm", "Nm"],
            ["Nm", "Nm", "Nm", "Nm^2", "Nm^2", "Nm^2"],
            ["Nm", "Nm", "Nm", "Nm^2", "Nm^2", "Nm^2"],
            ["Nm", "Nm", "Nm", "Nm^2", "Nm^2", "Nm^2"],
        ]
    )

    ut = np.zeros((6, 6))
    ut[np.triu_indices(6)] = 1

    fig2, ax2 = plt.subplots(6, 6)
    fig2.subplots_adjust(wspace=0.5, hspace=0.5)
    fig2.suptitle("6x6 Stiffness Matrix", fontsize=14)
    c = 6
    for j in range(6):
        for i in range(6):
            if ut[i, j] == 1:
                if isinstance(ref, (np.ndarray)):
                    ax2[i][j].plot(ref_x, ref[:, c], "--b")
                if isinstance(sigma, (np.ndarray)):
                    ax2[i][j].fill_between(
                        x, data[:, c] - sigma[:, c], data[:, c] + sigma[:, c], alpha=0.6, edgecolor="r", linestyle=":", facecolor="r", antialiased=True,
                    )
                ax2[i][j].plot(x, data[:, c], "--k.")
                ylabel = r"$k_{%s%s} \quad [%s]$" % (i + 1, j + 1, TSunits[i, j])
                ax2[i][j].set_ylabel(ylabel)
                # ax2[i][j].set_ylim(0)
                # ax2[i][j].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

                if i == j:
                    ax2[i][j].set_xlabel(r"radius $\quad [1/R]$")

                c += 1
            else:
                fig2.delaxes(ax2[i][j])
                # fig2.text(0,0,r'$A=B\cdotC$')
                # fig2.text(0,0,r'$\left(\begin{matrix}k11&k12&k13&k14&k15&k16\\k12&k22&k23&k24&k25&k26\\k13&k23&k33&k34&k35&k36\\k14&k24&k34&k44&k45&k46\\k15&k25&k35&k45&k55&k56\\k16&k26&k36&k46&k56&k66\end{matrix}\right)$')

    desc = (
        r"\begin{minipage}[b]{10cm} "
        r"\underline{\textbf{Description:}} \\"
        r"The 6x6 sectional stiffness matrix, TS (Timoshenko Stiffness Matrix) "
        r"(1-extension; 2,3-shear, 4-twist; 5,6-bending) relates the sectional axial "
        r"strain, $\epsilon_1$, transverse shearing strains, $\epsilon_2$ and $\epsilon_3$, "
        r"twisting curvatures, $\kappa_1$ and two bending curvatures, $\kappa_2$ and $\kappa_3$, "
        r"to the axial force, $F_1$, transverse shear forces, $F_2$ and $F_3$, twisting "
        r"moment, $M_1$, and two bending moments, $M_2$ and $M_3$. The relationship between "
        r"these sectional strains and sectional stress resultants takes the form of a"
        r"symmetric, 6x6 matrix: \\"
        r"$$ \left( \begin{matrix} F_{1} \\ F_{2} \\ F_{3} \\ M_{1} \\ M_{2} \\ M_{3} \end{matrix} \right) = "
        r"\left( \begin{matrix} k_{11} & k_{12} & k_{13} & k_{14} & k_{15} & k_{16} \\ "
        r"k_{12} & k_{22} & k_{23} & k_{24} & k_{25} & k_{26} \\ "
        r"k_{13} & k_{23} & k_{33} & k_{34} & k_{35} & k_{36} \\ "
        r"k_{14} & k_{24} & k_{34} & k_{44} & k_{45} & k_{46} \\ "
        r"k_{15} & k_{25} & k_{35} & k_{45} & k_{55} & k_{56} \\ "
        r"k_{16} & k_{26} & k_{36} & k_{46} & k_{56} & k_{66} \end{matrix} \right) "
        r"\cdot \left( \begin{matrix} \epsilon_{1} \\ \epsilon_{2} \\ \epsilon_{3} \\ \kappa_{1} \\ \kappa_{2} \\ \kappa_{3} \end{matrix} \right) $$ \end{minipage} "
    )

    if description == True:
        from matplotlib import rcParams

        plt.figtext(0.05, 0.05, desc, usetex=True, wrap=True, bbox=dict(ec=(1.0, 0.5, 0.5), fc=(1.0, 0.8, 0.8)))

    # ============ ONLY FOR PAPER!
    ut = np.zeros((6, 6))
    ut[np.triu_indices(6)] = 1

    fig3, ax3 = plt.subplots(3, 2)
    fig3.subplots_adjust(wspace=0.5, hspace=0.3)
    # fig3.suptitle('6x6 Stiffness Matrix', fontsize=14)
    c = 6
    p = 0
    q = 0
    for j in range(6):
        for i in range(6):
            if ut[i, j] == 1:
                if i == j:
                    if p == 3:
                        p = 0
                        q = 1

                    if isinstance(ref, (np.ndarray)):
                        ax3[p][q].plot(ref_x, ref[:, c], "-.b")
                    if isinstance(sigma, (np.ndarray)):
                        ax3[p][q].fill_between(
                            x, data[:, c] - sigma[:, c], data[:, c] + sigma[:, c], alpha=0.6, edgecolor="r", linestyle=":", facecolor="r", antialiased=True,
                        )
                    ax3[p][q].plot(x, data[:, c], "--k.", markersize=2)
                    ylabel = r"$k_{%s%s} \quad [%s]$" % (i + 1, j + 1, TSunits[i, j])
                    ax3[p][q].set_ylabel(ylabel)
                    # ax2[i][j].set_ylim(0)
                    # ax2[i][j].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
                    ax3[p][q].set_ylim(0)
                    if p == 2:
                        ax3[p][q].set_xlabel(r"radius $\quad [1/R]$")
                    p += 1
                #                else:
                #                    fig3.delaxes(ax3[i][j])

                c += 1

    # ax3[2,1].legend([ref_line, line], ['UH-60A Reference', 'Inertial Properties'], loc='upper center')
    ax3[2, 1].set_ylim(0, 1.4e7)

    plt.show()


def sim_plot(fq_dym, RPM_vec, target_dir):
    # plt.figure(figsize=(12,8))
    plt.figure(figsize=(12, 8))
    array_per = RPM_vec
    plt.plot(array_per, fq_dym[:, 1], "k*")
    plt.plot(array_per, fq_dym[:, 2], "k*")
    plt.plot(array_per, fq_dym[:, 3], "k*")
    plt.plot(array_per, fq_dym[:, 6], "k*")
    plt.plot(array_per, fq_dym[:, 7], "k*")
    plt.plot(array_per, fq_dym[:, 10], "k*")
    plt.plot(array_per, fq_dym[:, 0], "r*")
    plt.plot(array_per, fq_dym[:, 4], "r*")
    plt.plot(array_per, fq_dym[:, 5], "g*")
    plt.plot(array_per, fq_dym[:, 9], "k*")

    plt.xlabel(r"RPM / RPM_{ref}")
    plt.ylabel(r"Freq / RPM_{ref}")
    plt.legend()
    # tikz_save(target_dir+'F2.tikz',figurewidth='\\figurewidth', figureheight='\\figureheight')
    plt.savefig(target_dir + "F2.png")
    # plt.plot(np.array([0, 100]),np.array([0, 100]),'k-')
