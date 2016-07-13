import sys

if len(sys.argv)!=3:
    print "Usage: %s directory basename" % sys.argv[0]
    sys.exit(0)

directory = sys.argv[1]
basename = sys.argv[2]

from txbr import TxBR_Reconstruction

reconstruction = TxBR_Reconstruction(directory,basename)
XYZ = reconstruction.XYZ

def chi_boundaries(coeff):
    chi = 0;
    for xyz in XYZ:
        f1 = coeff[0]*xyz[0] + coeff[1]*xyz[1] - xyz[2] + coeff[2]
        f2 = coeff[0]*xyz[0] + coeff[1]*xyz[1] - xyz[2] + coeff[3]
        chi += pow(f1,2)*pow(f2,2)
    return chi

from scipy.optimize import fmin

x0 = [0, 0, -10, 10 ]
xopt = fmin(chi_boundaries, x0, xtol=1e-6, ftol=1e-6, maxiter=None, maxfun=10000,full_output=0)

def save_boundaries():
    filename = reconstruction.basename + ".pitch.txt"
    from os.path import join
    print filename
    filename = join(reconstruction.directory, "txbr-setup", filename)
    f = open(filename, 'w')
    f.write("%f    %f    %f\n" % (xopt[2],xopt[0],xopt[1]))
    f.write("%f    %f    %f\n" % (xopt[3],xopt[0],xopt[1]))
    f.close()

save_boundaries()

print xopt[2],xopt[0],xopt[1]
print xopt[3],xopt[0],xopt[1]