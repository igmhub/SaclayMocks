import os, sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
import fitsio
import picca.wedgize
import argparse
from SaclayMocks import powerspectrum


parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, nargs="*", help="ex: Out/v2.7.1/from_transmission")
parser.add_argument("--to-do", type=str, nargs="*", help="ex: cf xcf")
parser.add_argument("--rt-min", type=float, default=0.)
parser.add_argument("--rt-max", type=float, default=200.)
parser.add_argument("--nrt", type=int, default=50)
parser.add_argument("--rp-min", type=float, default=0.)
parser.add_argument("--rp-max", type=float, default=200.)
parser.add_argument("--nrp", type=int, default=50)
parser.add_argument("--mu-min", type=float, default=-1.)
parser.add_argument("--mu-max", type=float, default=1.)
parser.add_argument("--r-pow", type=int, default=2)
parser.add_argument("--title", default="")
parser.add_argument("--pred", action='store_true')
parser.add_argument("--z-bin", nargs='*', default=None)
parser.add_argument("--save-cf", default=False)

args = parser.parse_args()

if os.path.isdir(args.i[0]):
    dir_option = True
    indir = args.i[0]
else:  # If files specify, put fit in first, and then cf and xcf 
    dir_option = False
    if "cf" in args.to_do or "co" in args.to_do:
        fitcf_file = args.i[0]
        cf_file = args.i[1]
    if "xcf" in args.to_do:
        fitxcf_file = args.i[0]
        xcf_file = args.i[1]
    if "cf" in args.to_do and "xcf" in args.to_do:
        fitcf_file = args.i[0]
        fitxcf_file = args.i[1]
        cf_file = args.i[2]
        xcf_file = args.i[3]

rpmin = args.rp_min
rpmax = args.rp_max
nrp = args.nrp
rtmin = args.rt_min
rtmax = args.rt_max
nrt = args.nrt
mumin = args.mu_min
mumax = args.mu_max
r_pow = args.r_pow
title=args.title
if dir_option:
    print("Indir is :{}".format(indir))

# Treat CF
if "cf" in args.to_do:
    print("Starting plotting for CF...")
    # Read data
    if dir_option:
        filename = indir+"/Correlations/e_cf"
        if args.z_bin:
            filename += "_z_"+args.z_bin[0]+"_"+args.z_bin[1]
        filename += ".fits"
        data = fitsio.FITS(filename)
    else:
        data = fitsio.FITS(cf_file)
    da = data[1].read()['DA']
    co = data[1].read()['CO']
    data.close()

    # Read fit
    if dir_option:
        filename = indir+"/Fit/result_cf"
        if args.z_bin:
            filename += "_z_"+args.z_bin[0]+"_"+args.z_bin[1]
        filename += ".h5"
        ff = h5py.File(filename)
    else:
        ff = h5py.File(fitcf_file)

    # # Read xi pred
    # data_pred = np.loadtxt("data/xi_pred.txt")
    # r_pred = data_pred[:,0]

    # Plot X0
    w = picca.wedgize.wedge(mumin=0.,mumax=1., rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)
    data_wedge = w.wedge(da,co)
    try:
        r,f,_ = w.wedge(ff["LYA(LYA)xLYA(LYA)/fit"][...],co)
    except KeyError:
        print("Can't find LYA(LYA)xLYA(LYA)/fit")
        try:
            i = int(fitcf_file.find(".h5"))
            j = int(fitcf_file.rfind("/"))+1
            r,f,_ = w.wedge(ff[fitcf_file[j:i]+"/fit"][...],co)
        except KeyError:
            print("Can't find {}".format(fitcf_file[j:i]+"/fit"))
            try:
                i = int(cf_file.find(".fits"))
                j = int(cf_file.rfind("/"))+1
                r,f,_ = w.wedge(ff[cf_file[j:i]+"/fit"][...],co)
            except KeyError:
                print("Can't find {}".format(cf_file[j:i]+"/fit"))
                try:
                    r, f, _ =w.wedge(ff["cf_z_0_10/fit"][...],co)
                except KeyError:
                    print("Can't find {}".format("cf_z_0_10/fit"))
                    try:
                        r,f,_ = w.wedge(ff['QSOxQSO/fit'][...],co)
                    except:
                        print("Can't find QSOxQSO/fit")
                        try:
                            r,f,_ = w.wedge(ff['HCDxHCD/fit'][...],co)
                        except:
                            print("Can't find HCDxHCD/fit")
                            print("Exit!")
                            sys.exit(1)
    coef = data_wedge[0]**r_pow
    if args.pred:
        data = fitsio.FITS(indir+"/Correlations/e_cf_pred.fits")
        da_pred = data[1].read()['DA']
        co_pred = data[1].read()['CO']
        data.close()
        data_wedge_pred = w.wedge(da_pred, co_pred)

    fig, ax = plt.subplots()
    ax.errorbar(data_wedge[0],coef*data_wedge[1],yerr=coef*np.sqrt(np.diag(data_wedge[2])),fmt='+', label='mock')
    ax.plot(r, f*r**r_pow, label='fit')
    if args.pred:
        # ax.plot(r_pred, data_pred[:,1]*r_pred**r_pow, label='pred', linestyle=':')
        ax.plot(data_wedge[0],coef*data_wedge_pred[1],linestyle=':', label='pred')
    ax.grid()
    ax.legend()
    ax.set_title("CF - {}".format(title))
    ax.set_xlabel(r"$r \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 2:
        ax.set_ylabel(r"$r^{2}\xi(r) \; [(h^{-1}\mathrm{Mpc})^2]$")
    if r_pow == 1:
        ax.set_ylabel(r"$r\xi(r) \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 0:
        ax.set_ylabel(r"$\xi(r)$")
    plt.tight_layout()

    # Plot wedges
    mu0, mu1, mu2, mu3, mu4 = 0, 0.5, 0.8, 0.95, 1
    w1 = picca.wedgize.wedge(mumin=mu0,mumax=mu1, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)
    w2 = picca.wedgize.wedge(mumin=mu1,mumax=mu2, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)
    w3 = picca.wedgize.wedge(mumin=mu2,mumax=mu3, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)
    w4 = picca.wedgize.wedge(mumin=mu3,mumax=mu4, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)
    data_wedge1 = w1.wedge(da,co)
    coef1 = data_wedge1[0]**r_pow
    data_wedge2 = w2.wedge(da,co)
    coef2 = data_wedge2[0]**r_pow
    data_wedge3 = w3.wedge(da,co)
    coef3 = data_wedge3[0]**r_pow
    data_wedge4 = w4.wedge(da,co)
    coef4 = data_wedge4[0]**r_pow

    if args.save_cf:
        print("Saving cf in txt file...")
        table0 = [data_wedge[0], data_wedge[1], np.sqrt(np.diag(data_wedge[2]))]
        table1 = [data_wedge1[0], data_wedge1[1], np.sqrt(np.diag(data_wedge1[2]))]
        table2 = [data_wedge2[0], data_wedge2[1], np.sqrt(np.diag(data_wedge2[2]))]
        table3 = [data_wedge3[0], data_wedge3[1], np.sqrt(np.diag(data_wedge3[2]))]
        table4 = [data_wedge4[0], data_wedge4[1], np.sqrt(np.diag(data_wedge4[2]))]
        table0 = np.array(table0).T
        table1 = np.array(table1).T
        table2 = np.array(table2).T
        table3 = np.array(table3).T
        table4 = np.array(table4).T
        np.savetxt(indir+"/Correlations/cf.txt", table0)
        np.savetxt(indir+"/Correlations/cf_{}_{}.txt".format(mu0, mu1), table1)
        np.savetxt(indir+"/Correlations/cf_{}_{}.txt".format(mu1, mu2), table2)
        np.savetxt(indir+"/Correlations/cf_{}_{}.txt".format(mu2, mu3), table3)
        np.savetxt(indir+"/Correlations/cf_{}_{}.txt".format(mu3, mu4), table4)
        print("Done.")

    try:
        r1,f1,_ = w1.wedge(ff["LYA(LYA)xLYA(LYA)/fit"][...],co)
        r2,f2,_ = w2.wedge(ff["LYA(LYA)xLYA(LYA)/fit"][...],co)
        r3,f3,_ = w3.wedge(ff["LYA(LYA)xLYA(LYA)/fit"][...],co)
        r4,f4,_ = w4.wedge(ff["LYA(LYA)xLYA(LYA)/fit"][...],co)
    except KeyError:
        print("Can't find LYA(LYA)xLYA(LYA)/fit")
        try:
            i = int(fitcf_file.find(".h5"))
            j = int(fitcf_file.rfind("/"))+1
            r1,f1,_ = w1.wedge(ff[fitcf_file[j:i]+"/fit"][...],co)
            r2,f2,_ = w2.wedge(ff[fitcf_file[j:i]+"/fit"][...],co)
            r3,f3,_ = w3.wedge(ff[fitcf_file[j:i]+"/fit"][...],co)
            r4,f4,_ = w4.wedge(ff[fitcf_file[j:i]+"/fit"][...],co)
        except:
            print("Can't find {}".format(fitcf_file[j:i]+"/fit"))
            try:
                i = int(cf_file.find(".fits"))
                j = int(cf_file.rfind("/"))+1
                r1,f1,_ = w1.wedge(ff[cf_file[j:i]+"/fit"][...],co)
                r2,f2,_ = w2.wedge(ff[cf_file[j:i]+"/fit"][...],co)
                r3,f3,_ = w3.wedge(ff[cf_file[j:i]+"/fit"][...],co)
                r4,f4,_ = w4.wedge(ff[cf_file[j:i]+"/fit"][...],co)
            except:
                print("Can't find {}".format(cf_file[j:i]+"/fit"))
                try:
                    r1, f1, _ =w1.wedge(ff["cf_z_0_10/fit"][...],co)
                    r2, f2, _ =w2.wedge(ff["cf_z_0_10/fit"][...],co)
                    r3, f3, _ =w3.wedge(ff["cf_z_0_10/fit"][...],co)
                    r4, f4, _ =w4.wedge(ff["cf_z_0_10/fit"][...],co)
                except KeyError:
                    try:
                        r1,f1,_ = w1.wedge(ff['QSOxQSO/fit'][...],co)
                        r2,f2,_ = w2.wedge(ff['QSOxQSO/fit'][...],co)
                        r3,f3,_ = w3.wedge(ff['QSOxQSO/fit'][...],co)
                        r4,f4,_ = w4.wedge(ff['QSOxQSO/fit'][...],co)
                    except:
                        print("Can't find QSOxQSO/fit")
                        try:
                            r1,f1,_ = w1.wedge(ff['HCDxHCD/fit'][...],co)
                            r2,f2,_ = w2.wedge(ff['HCDxHCD/fit'][...],co)
                            r3,f3,_ = w3.wedge(ff['HCDxHCD/fit'][...],co)
                            r4,f4,_ = w4.wedge(ff['HCDxHCD/fit'][...],co)
                        except:
                            print("Can't find HCDxHCD/fit")
                            print("Exit!")
                            sys.exit(1)

    if args.pred:
        data_wedge1_pred = w1.wedge(da_pred,co_pred)
        data_wedge2_pred = w2.wedge(da_pred,co_pred)
        data_wedge3_pred = w3.wedge(da_pred,co_pred)
        data_wedge4_pred = w4.wedge(da_pred,co_pred)
    
    fig, ax = plt.subplots()
    ax.errorbar(data_wedge1[0],coef1*data_wedge1[1],yerr=coef1*np.sqrt(np.diag(data_wedge1[2])),fmt='+', label=r"${}<\mu<{}$".format(mu0, mu1), color='b')
    ax.errorbar(data_wedge2[0],coef2*data_wedge2[1],yerr=coef2*np.sqrt(np.diag(data_wedge2[2])),fmt='+', label=r"${}<\mu<{}$".format(mu1, mu2), color='g')
    ax.errorbar(data_wedge3[0],coef3*data_wedge3[1],yerr=coef3*np.sqrt(np.diag(data_wedge3[2])),fmt='+', label=r"${}<\mu<{}$".format(mu2, mu3), color='orange')
    ax.errorbar(data_wedge4[0],coef4*data_wedge4[1],yerr=coef4*np.sqrt(np.diag(data_wedge4[2])),fmt='+', label=r"${}<\mu<{}$".format(mu3, mu4), color='red')
    ax.plot(r1, f1*r1**r_pow, color='b')
    ax.plot(r2, f2*r2**r_pow, color='g')
    ax.plot(r3, f3*r3**r_pow, color='orange')
    ax.plot(r4, f4*r4**r_pow, color='red')
    if args.pred:
        # ax.plot(r_pred, data_pred[:,2]*r_pred**r_pow, color='b', linestyle=':')
        # ax.plot(r_pred, data_pred[:,3]*r_pred**r_pow, color='g', linestyle=':')
        # ax.plot(r_pred, data_pred[:,4]*r_pred**r_pow, color='r', linestyle=':')
        ax.plot(data_wedge1[0],coef1*data_wedge1_pred[1], linestyle=':', color='b')
        ax.plot(data_wedge2[0],coef2*data_wedge2_pred[1], linestyle=':', color='g')
        ax.plot(data_wedge3[0],coef3*data_wedge3_pred[1], linestyle=':', color='orange')
        ax.plot(data_wedge4[0],coef4*data_wedge4_pred[1], linestyle=':', color='red')
        
    ax.grid()
    ax.set_title("CF - {}".format(title))
    ax.set_xlabel(r"$r \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 2:
        ax.set_ylabel(r"$r^{2}\xi(r) \; [(h^{-1}\mathrm{Mpc})^2]$")
    if r_pow == 1:
        ax.set_ylabel(r"$r\xi(r) \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 0:
        ax.set_ylabel(r"$\xi(r)$")

    ax.legend()
    plt.tight_layout()

# Treat XCF
if "xcf" in args.to_do:
    print("Starting plotting for XCF...")
    # Read data
    if dir_option:
        filename = indir+"/Correlations/e_xcf"
        if args.z_bin:
            filename += "_z_"+args.z_bin[0]+"_"+args.z_bin[1]
        filename += ".fits"
        data = fitsio.FITS(filename)
    else:
        data = fitsio.FITS(xcf_file)
    da = data[1].read()['DA']
    co = data[1].read()['CO']
    data.close()

    # Read fit
    if dir_option:
        filename = indir+"/Fit/result_xcf"
        if args.z_bin:
            filename += "_z_"+args.z_bin[0]+"_"+args.z_bin[1]
        filename += ".h5"
        ff = h5py.File(filename)
    else:
        ff = h5py.File(fitxcf_file)

    # Plot
    w = picca.wedgize.wedge(rpmin=rpmin,rpmax=rpmax,nrp=nrp,rtmin=rtmin,rtmax=rtmax,nrt=nrt,mumin=mumin,mumax=mumax,absoluteMu=True)
    # w = picca.wedgize.wedge(rpmin=rpmin,rpmax=rpmax,nrp=nrp,rtmin=rtmin,rtmax=rtmax,nrt=nrt,mumin=0,mumax=1,absoluteMu=False)
    # w = picca.wedgize.wedge(rpmin=rpmin,rpmax=rpmax,nrp=nrp,rtmin=rtmin,rtmax=rtmax,nrt=nrt,mumin=-1,mumax=0, absoluteMu=False)
    data_wedge = w.wedge(da,co)
    try:
        r,f,_ = w.wedge(ff["LYA(LYA)xQSO/fit"][...],co)
    except KeyError:
        print("Can't find LYA(LYA)xQSO/fit")
        try:
            r,f,_ = w.wedge(ff["xcf_z_0_10-exp/fit"][...],co)
        except KeyError:
            print("Can't find xcf_z_0_10-exp/fit")
            try:
                i = int(fitxcf_file.find(".h5"))
                j = int(fitxcf_file.rfind("/"))+1
                r,f,_ = w.wedge(ff[fitxcf_file[j:i]+"/fit"][...],co)
            except:
                print("Can't find {}".format(fitxcf_file[j:i]+"/fit"))
                try:
                    i = int(xcf_file.find(".fits"))
                    j = int(xcf_file.rfind("/"))+1
                    r,f,_ = w.wedge(ff[xcf_file[j:i]+"/fit"][...],co)
                except:
                    print("Can't find {}".format(xcf_file[j:i]+"/fit"))
                    sys.exit(1)

    coef = data_wedge[0]**r_pow

    fig, ax = plt.subplots()
    ax.errorbar(data_wedge[0],coef*data_wedge[1],yerr=coef*np.sqrt(np.diag(data_wedge[2])),fmt='+', label='mock')
    ax.plot(r, f*r**r_pow, label='fit')
    ax.grid()
    ax.set_title("XCF - {}".format(title))
    ax.set_xlabel(r"$r \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 2:
        ax.set_ylabel(r"$r^{2}\xi(r) \; [(h^{-1}\mathrm{Mpc})^2]$")
    if r_pow == 1:
        ax.set_ylabel(r"$r\xi(r) \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 0:
        ax.set_ylabel(r"$\xi(r)$")
    
    # Plot wedges
    mu0, mu1, mu2, mu3, mu4 = 0, 0.5, 0.8, 0.95, 1
    # mu0, mu1, mu2, mu3, mu4 = -1, -0.95, -0.8, -0.5, 0
    w1 = picca.wedgize.wedge(rpmin=rpmin,rpmax=rpmax,nrp=nrp,rtmin=rtmin,rtmax=rtmax,nrt=nrt,mumin=mu0,mumax=mu1,absoluteMu=True)
    w2 = picca.wedgize.wedge(rpmin=rpmin,rpmax=rpmax,nrp=nrp,rtmin=rtmin,rtmax=rtmax,nrt=nrt,mumin=mu1,mumax=mu2,absoluteMu=True)
    w3 = picca.wedgize.wedge(rpmin=rpmin,rpmax=rpmax,nrp=nrp,rtmin=rtmin,rtmax=rtmax,nrt=nrt,mumin=mu2,mumax=mu3,absoluteMu=True)
    w4 = picca.wedgize.wedge(rpmin=rpmin,rpmax=rpmax,nrp=nrp,rtmin=rtmin,rtmax=rtmax,nrt=nrt,mumin=mu3,mumax=mu4,absoluteMu=True)
    data_wedge1 = w1.wedge(da,co)
    coef1 = data_wedge1[0]**r_pow
    data_wedge2 = w2.wedge(da,co)
    coef2 = data_wedge2[0]**r_pow
    data_wedge3 = w3.wedge(da,co)
    coef3 = data_wedge3[0]**r_pow
    data_wedge4 = w4.wedge(da,co)
    coef4 = data_wedge4[0]**r_pow
    try:
        r1,f1,_ = w1.wedge(ff["LYA(LYA)xQSO/fit"][...],co)
        r2,f2,_ = w2.wedge(ff["LYA(LYA)xQSO/fit"][...],co)
        r3,f3,_ = w3.wedge(ff["LYA(LYA)xQSO/fit"][...],co)
        r4,f4,_ = w4.wedge(ff["LYA(LYA)xQSO/fit"][...],co)
    except KeyError:
        print("Can't find LYA(LYA)xQSO/fit")
        try:
            r1,f1,_ = w1.wedge(ff["xcf_z_0_10-exp/fit"][...],co)
            r2,f2,_ = w2.wedge(ff["xcf_z_0_10-exp/fit"][...],co)
            r3,f3,_ = w3.wedge(ff["xcf_z_0_10-exp/fit"][...],co)
            r4,f4,_ = w4.wedge(ff["xcf_z_0_10-exp/fit"][...],co)
        except KeyError:
            print("Can't find xcf_z_0_10-exp/fit")
            try:
                i = int(fitxcf_file.find(".h5"))
                j = int(fitxcf_file.rfind("/"))+1
                r1,f1,_ = w1.wedge(ff[fitxcf_file[j:i]+"/fit"][...],co)
                r2,f2,_ = w2.wedge(ff[fitxcf_file[j:i]+"/fit"][...],co)
                r3,f3,_ = w3.wedge(ff[fitxcf_file[j:i]+"/fit"][...],co)
                r4,f4,_ = w4.wedge(ff[fitxcf_file[j:i]+"/fit"][...],co)
            except:
                print("Can't find {}".format(fitxcf_file[j:i]+"/fit"))
                try:
                    i = int(xcf_file.find(".fits"))
                    j = int(xcf_file.rfind("/"))+1
                    r1,f1,_ = w1.wedge(ff[xcf_file[j:i]+"/fit"][...],co)
                    r2,f2,_ = w2.wedge(ff[xcf_file[j:i]+"/fit"][...],co)
                    r3,f3,_ = w3.wedge(ff[xcf_file[j:i]+"/fit"][...],co)
                    r4,f4,_ = w4.wedge(ff[xcf_file[j:i]+"/fit"][...],co)
                except:
                    print("Can't find {}".format(xcf_file[j:i]+"/fit"))

    fig, ax = plt.subplots()
    ax.errorbar(data_wedge1[0],coef1*data_wedge1[1],yerr=coef1*np.sqrt(np.diag(data_wedge1[2])),fmt='+', label=r"${}<\mu<{}$".format(mu0, mu1), color='b')
    ax.plot(r1, f1*r1**r_pow, color='b')
    ax.errorbar(data_wedge2[0],coef2*data_wedge2[1],yerr=coef2*np.sqrt(np.diag(data_wedge2[2])),fmt='+', label=r"${}<\mu<{}$".format(mu1, mu2), color='g')
    ax.plot(r2, f2*r2**r_pow, color='g')
    ax.errorbar(data_wedge3[0],coef3*data_wedge3[1],yerr=coef3*np.sqrt(np.diag(data_wedge3[2])),fmt='+', label=r"${}<\mu<{}$".format(mu2, mu3), color='orange')
    ax.plot(r3, f3*r3**r_pow, color='orange')
    ax.errorbar(data_wedge4[0],coef4*data_wedge4[1],yerr=coef4*np.sqrt(np.diag(data_wedge4[2])),fmt='+', label=r"${}<\mu<{}$".format(mu3, mu4), color='r')
    ax.plot(r4, f4*r4**r_pow, color='r')

    ax.grid()
    ax.set_title("XCF - {}".format(title))
    ax.set_xlabel(r"$r \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 2:
        ax.set_ylabel(r"$r^{2}\xi(r) \; [(h^{-1}\mathrm{Mpc})^2]$")
    if r_pow == 1:
        ax.set_ylabel(r"$r\xi(r) \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 0:
        ax.set_ylabel(r"$\xi(r)$")

    ax.legend()
    plt.tight_layout()

# Treat CO
if "co" in args.to_do:
    print("Starting plotting for CO...")
    # Read data
    if dir_option:
        filename = indir+"/Correlations/e_co_{}".format(args.to_do[1])
        if args.z_bin:
            filename += "_z_"+args.z_bin[0]+"_"+args.z_bin[1]
        filename += ".fits"
        data = fitsio.FITS(filename)
    else:
        data = fitsio.FITS(cf_file)
    da = data[1].read()['DA']
    co = data[1].read()['CO']
    data.close()

    # Read fit
    if dir_option:
        filename = indir+"/Fit/result_co_{}".format(args.to_do[1])
        if args.z_bin:
            filename += "_z_"+args.z_bin[0]+"_"+args.z_bin[1]
        filename += ".h5"
        ff = h5py.File(filename)
    else:
        ff = h5py.File(fitcf_file)

    # # Read xi pred
    # data_pred = np.loadtxt("data/xi_pred.txt")
    # r_pred = data_pred[:,0]

    # Plot X0
    w = picca.wedgize.wedge(mumin=0.,mumax=1., rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp)
    data_wedge = w.wedge(da,co)
    try:
        r_fit,cf_fit,_ = w.wedge(ff["QSOxQSO/fit"][...],co)
    except KeyError:
        print("Can't find QSOxQSO/fit")
        try:
            r_fit,cf_fit,_ = w.wedge(ff["HCDxHCD/fit"][...],co)
        except KeyError:
            print("Can't find HCDxHCD/fit")
            try:
                i = int(fitcf_file.find(".h5"))
                j = int(fitcf_file.rfind("/"))+1
                r_fit,cf_fit,_ = w.wedge(ff[fitcf_file[j:i]+"/fit"][...],co)
            except:
                print("Can't find {}".format(fitcf_file[j:i]+"/fit"))

    coef = data_wedge[0]**r_pow
    # if args.pred:
    #     data = fitsio.FITS(indir+"/Correlations/e_cf_pred.fits")
    #     da_pred = data[1].read()['DA']
    #     co_pred = data[1].read()['CO']
    #     data.close()
    #     data_wedge_pred = w.wedge(da_pred, co_pred)

    if args.pred:
        #.............................  xi_CAMB
        Pcamb = powerspectrum.P_0()
        k = np.linspace(0,10,10000)
        P = Pcamb.P(k)
        r,xi = powerspectrum.xi_from_pk(k,P)
        cut=np.where(r<300)
        r=r[cut]
        xi=xi[cut]
        if args.to_do[1] == "QSO":
            xx = 1.87  # bias and growth factor for QSO: D(z=2.4)**2 * 3.7**2
        if args.to_do[1] == "DLA":
            xx = 0.58  # bias and growth factor for DLA: 0.38**2 * 2**2
        xi *= xx

        #.........................  xibar, xibarbar, xi_0. xi_2, xi4 from Hamilton
        xi_Ham = powerspectrum.xi_Hamilton(r,xi,300)
        if args.to_do[1] == "QSO":
            bias = 3.7  # according to Constant.py
            beta = 0.96 / bias
        if args.to_do[1] == "DLA":
            bias = 2.  # default value in DLA_Saclay.py
            beta = 0.96 / bias
        xi0 = xi_Ham.xi0(beta,r)

        #............................  aproximate wedges
        xi0_5 = np.zeros(len(r))
        for mu in np.linspace(0,0.5,10) :
            xi0_5 += xi_Ham.xi(beta,r,mu)
        xi0_5 /= 10
        xi5_8 = np.zeros(len(r))
        for mu in np.linspace(0.5,0.8,10) :
            xi5_8 += xi_Ham.xi(beta,r,mu)
        xi5_8 /= 10
        xi8_10 = np.zeros(len(r))
        for mu in np.linspace(0.8,1,10) :
            xi8_10 += xi_Ham.xi(beta,r,mu)
        xi8_10 /= 10

    # fig, ax = plt.subplots(figsize=(12,8))
    fig, ax = plt.subplots()
    ax.errorbar(data_wedge[0],coef*data_wedge[1],yerr=coef*np.sqrt(np.diag(data_wedge[2])),fmt='+', label='mock')
    ax.plot(r_fit, cf_fit*r_fit**r_pow, label='fit')
    # if args.pred:
    #     # ax.plot(r_pred, data_pred[:,1]*r_pred**r_pow, label='pred', linestyle=':')
    #     ax.plot(data_wedge[0],coef*data_wedge_pred[1],linestyle=':', label='pred')
    if args.pred:
        ax.plot(r, xi0 * r**r_pow, label='pred', linestyle=':')
    ax.grid()
    ax.legend()
    # ax.set_title("CF - {}".format(title), fontsize=20)
    ax.set_title(title)
    ax.set_xlabel(r"$r \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 2:
        ax.set_ylabel(r"$r^{2}\xi(r) \; [(h^{-1}\mathrm{Mpc})^2]$")
    if r_pow == 1:
        ax.set_ylabel(r"$r\xi(r) \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 0:
        ax.set_ylabel(r"$\xi(r)$")

    # Plot wedges
    mu0, mu1, mu2, mu3 = 0, 0.2, 0.5, 1
    w1 = picca.wedgize.wedge(mumin=mu0,mumax=mu1, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp)
    w2 = picca.wedgize.wedge(mumin=mu1,mumax=mu2, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp)
    w3 = picca.wedgize.wedge(mumin=mu2,mumax=mu3, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp)
    data_wedge1 = w1.wedge(da,co)
    coef1 = data_wedge1[0]**r_pow
    data_wedge2 = w2.wedge(da,co)
    coef2 = data_wedge2[0]**r_pow
    data_wedge3 = w3.wedge(da,co)
    coef3 = data_wedge3[0]**r_pow
    try:
        r1,f1,_ = w1.wedge(ff["QSOxQSO/fit"][...],co)
        r2,f2,_ = w2.wedge(ff["QSOxQSO/fit"][...],co)
        r3,f3,_ = w3.wedge(ff["QSOxQSO/fit"][...],co)
    except KeyError:
        print("Can't find QSOxQSO/fit")
        try:
            r1,f1,_ = w1.wedge(ff["HCDxHCD/fit"][...],co)
            r2,f2,_ = w2.wedge(ff["HCDxHCD/fit"][...],co)
            r3,f3,_ = w3.wedge(ff["HCDxHCD/fit"][...],co)
        except KeyError:
            print("Can't find HCDxHCD/fit")
            try:
                i = int(fitcf_file.find(".h5"))
                j = int(fitcf_file.rfind("/"))+1
                r1,f1,_ = w1.wedge(ff[fitcf_file[j:i]+"/fit"][...],co)
                r2,f2,_ = w2.wedge(ff[fitcf_file[j:i]+"/fit"][...],co)
                r3,f3,_ = w3.wedge(ff[fitcf_file[j:i]+"/fit"][...],co)
            except:
                print("Can't find {}".format(fitcf_file[j:i]+"/fit"))

    # if args.pred:
    #     data_wedge1_pred = w1.wedge(da_pred,co_pred)
    #     data_wedge2_pred = w2.wedge(da_pred,co_pred)
    #     data_wedge3_pred = w3.wedge(da_pred,co_pred)
    
    # fig, ax = plt.subplots(figsize=(12,8))
    fig, ax = plt.subplots()
    ax.errorbar(data_wedge1[0],coef1*data_wedge1[1],yerr=coef1*np.sqrt(np.diag(data_wedge1[2])),fmt='+', label=r"${}<\mu<{}$".format(mu0, mu1), color='b')
    ax.errorbar(data_wedge2[0],coef2*data_wedge2[1],yerr=coef2*np.sqrt(np.diag(data_wedge2[2])),fmt='+', label=r"${}<\mu<{}$".format(mu1, mu2), color='g')
    ax.errorbar(data_wedge3[0],coef3*data_wedge3[1],yerr=coef3*np.sqrt(np.diag(data_wedge3[2])),fmt='+', label=r"${}<\mu<{}$".format(mu2, mu3), color='r')
    ax.plot(r1, f1*r1**r_pow, color='b')
    ax.plot(r2, f2*r2**r_pow, color='g')
    ax.plot(r3, f3*r3**r_pow, color='r')
    if args.pred:
        # ax.plot(r_pred, data_pred[:,2]*r_pred**r_pow, color='b', linestyle=':')
        # ax.plot(r_pred, data_pred[:,3]*r_pred**r_pow, color='g', linestyle=':')
        # ax.plot(r_pred, data_pred[:,4]*r_pred**r_pow, color='r', linestyle=':')
        # ax.plot(data_wedge1[0],coef1*data_wedge1_pred[1], linestyle=':', color='b')
        # ax.plot(data_wedge2[0],coef2*data_wedge2_pred[1], linestyle=':', color='g')
        # ax.plot(data_wedge3[0],coef3*data_wedge3_pred[1], linestyle=':', color='r')
        ax.plot(r, xi0_5 * r**r_pow, color='b', linestyle=':')
        ax.plot(r, xi5_8 * r**r_pow, color='g', linestyle=':')
        ax.plot(r, xi8_10 * r**r_pow, color='r', linestyle=':')
    ax.grid()
    # ax.set_title("CF - {}".format(title), fontsize=20)
    ax.set_title(title)
    ax.set_xlabel(r"$r \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 2:
        ax.set_ylabel(r"$r^{2}\xi(r) \; [(h^{-1}\mathrm{Mpc})^2]$")
    if r_pow == 1:
        ax.set_ylabel(r"$r\xi(r) \; [h^{-1}\mathrm{Mpc}]$")
    if r_pow == 0:
        ax.set_ylabel(r"$\xi(r)$")

    ax.legend()
    plt.tight_layout()

plt.show()
