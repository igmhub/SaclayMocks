import numpy as np
from numpy import ma
from SaclayMocks import constant, transmissions
from matplotlib import pyplot as plt
import argparse
import glob


parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="ex: DesiMocks/.../output/")
parser.add_argument("-o", type=str, help="ex: Out/v2.7.1/mean_trans/")
parser.add_argument("--to-do", type=str, nargs="*", help="ex: compute plot")
parser.add_argument("--nfiles", default=None)
parser.add_argument("--title", default="")
args = parser.parse_args()

outdir = args.o

if "compute" in args.to_do:
    print("Computing mean transmission vs z...")
    lmin = constant.lylimit
    lmax = constant.lya
    lmin = [constant.lylimit, 900, 950, 1000, 1050, 1100, 1120, 1140, 1160, 1180, 1200]
    lmax = [constant.lya, 950, 1000, 1050, 1100, 1120, 1140, 1160, 1180, 1200, constant.lya]

    trans = transmissions.ReadTransmission(args.i, read_dla=False, nfiles=args.nfiles)
    lambda_rf = trans.wavelength / (1+trans.metadata['Z'].reshape(-1,1))
    redshift = trans.wavelength / constant.lya - 1
    np.save(outdir+"/redshift.npy", redshift)
    
    for i in range(len(lmin)):
        msk = (lambda_rf < lmax[i]) & (lambda_rf > lmin[i])
        transmissions = ma.array(trans.transmission, mask=~msk)
        mean_trans = ma.mean(transmissions, axis=0)
        np.save(outdir+"/mean_trans_{}_{}.npy".format(lmin[i], lmax[i]), mean_trans.data)
        print("iteration {} done".format(i))

if "plot" in args.to_do:
    files = np.sort(glob.glob(outdir+"/mean_trans_*"))
    names = [files[i][files[i].find('trans_')+6:files[i].rfind('.')] for i in range(len(files))]
    z = np.load(outdir+"/redshift.npy")
    f, ax = plt.subplots()
    for i in range(len(files)):
        t = np.load(files[i])
        ax.plot(z, t, label=names[i])

    ax.legend()
    ax.set_xlabel('z')
    ax.set_ylabel('mean transmission')
    ax.grid()
    ax.set_title(args.title)
    plt.show()
