import os, sys, shutil
import numpy as np
import segyio
from segyio import BinField
from mpi4py import MPI

# Ensure all processes abort if one does
def mpiAbort(msg):
     print(msg)
     sys.stdout.flush()
     MPI.COMM_WORLD.Abort(-1)
# Additional printing (dev only)
DEBUG = False

# hard-link to Madagascar installation so user doesn't have to do much
os.environ["PATH"] += os.pathsep + "/apps/madacascar/2.0-1/bin/"
os.environ["LD_LIBRARY_PATH"] += os.pathsep + "/apps/madacascar/2.0-1/lib/"

# MPI basic information
size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name() 
#srnk = MPI.COMM_WORLD.Split(color=int(name[-2:]), key=rank).Get_rank() # unimplemented Split_type
srnk = MPI.COMM_WORLD.Split_type(MPI.COMM_TYPE_SHARED).Get_rank()
if DEBUG: print("Rank %d of %d, node %s" % (rank, size, name))

# Parse arguments ################################################################################# 
if len(sys.argv) == 10:                                                                           #
    VTI = True                                                                                    #
    if rank == 0:                                                                                 #
        print("Isotropy: VTI")                                                                    #
else:                                                                                             #
    VTI = False                                                                                   #
    if rank == 0:                                                                                 #
        if len(sys.argv) != 8:                                                                    #
            mpiAbort("Arguments: VM.sgy dx dy xOffset yOffset xSkip ySkip [eps.sgy del.sgy]")     #
        print("Isotropic")                                                                        #
try:                                                                                              #
    vmFname = sys.argv[1]                                                                         #
    dx = max(float(sys.argv[2]), 0.1)                                                             #
    dy = max(float(sys.argv[3]), 0.1)                                                             #
    rx = max(float(sys.argv[4]), 0.01)                                                            #
    ry = max(float(sys.argv[5]), 0.01)                                                            #
    xSkip = max(int(sys.argv[6]), 1)                                                              #
    ySkip = max(int(sys.argv[7]), 1)                                                              #
except ValueError:                                                                                #
    if rank == 0:                                                                                 #
        mpiAbort("Could not parse arguments (file, 4 floats, 2 ints, [2 more files])")            #
if rank == 0:                                                                                     #
    if not os.path.exists(vmFname):                                                               #
        mpiAbort('Model file not found, argument 1')                                              #
    if DEBUG:                                                                                     #
        print("Found VM file %s" % (os.path.basename(vmFname)))                                   #
vmFname = os.path.abspath(vmFname)                                                                #
SHAREDPATH = os.path.dirname(vmFname) + "/"                                                       #
vmRsf = SHAREDPATH + "vm.rsf"                                                                     #
if VTI:                                                                                           #
    epsFname = sys.argv[8]                                                                        #
    delFname = sys.argv[9]                                                                        #
    vxRsf = SHAREDPATH + "vx.rsf"                                                                 #
    etaRsf = SHAREDPATH + "eta.rsf"                                                               #
    if rank == 0:                                                                                 #
        if not os.path.exists(epsFname):                                                          #
            mpiAbort('Epsilon file not found, argument 8')                                        #
        else:                                                                                     #
            if os.path.getsize(epsFname) != os.path.getsize(vmFname):                             #
                mpiAbort('Epsilon file size must match VM, argument 8')                           #
        if not os.path.exists(delFname):                                                          #
            mpiAbort('Delta file not found, argument 9')                                          #
        else:                                                                                     #
            if os.path.getsize(delFname) != os.path.getsize(vmFname):                             #
                mpiAbort('Delta file size must match VM, argument 9')                             #
    epsFname = os.path.abspath(epsFname)                                                          #
    delFname = os.path.abspath(delFname)                                                          #
# End parse arguments #############################################################################

# Prepare scratch directory for Madagascar, and 'cd' there
DATAPATH = "/dev/shm/TTtmp/"
if srnk == 0:
    shutil.rmtree(DATAPATH, ignore_errors=True)
    os.makedirs(DATAPATH)
MPI.COMM_WORLD.barrier()
os.environ["DATAPATH"] = DATAPATH
os.chdir(DATAPATH)

# Establish model geometry and convert to Madagascar format (1 process, then spread)
nx = ny = nz = 0
dz = 0.
if rank == 0:
    with segyio.open(vmFname, "r", ignore_geometry=True) as vmSgy:
        vmSgy.mmap()
        ilines = vmSgy.attributes(segyio.TraceField.INLINE_3D)[:]
        if ilines[-1] == ilines[0]:
            nx = 1
            ny = len(ilines)
            if dy == 0.1 or ry == 0.01:
                mpiAbort('Model is 2D (y,z) but y spacing or radius is too small.')
        else:
            decr = np.less(ilines[1:], ilines[:-1])
            if not any(decr):
                nx = len(ilines)
                ny = 1
                if dx == 0.1 or rx == 0.01:
                    mpiAbort('Model is 2D (x,z) but x spacing or radius is too small.')
            else:
                nx = list(decr).index(True) + 1
                ny = len(ilines) // nx
        dz = vmSgy.bin[BinField.Interval] / 1000
        nz = vmSgy.bin[BinField.Samples]
        print("Model Dimensions (il, xl, z): ", nx, ny, nz)
        print("Model Spacing [m] (il, xl, z): ", dx, dy, dz)
        print("Traveltime radii [m] (il, xl)", rx, ry)
        sys.stdout.flush()
    cmd = "< %s sfsegyread read=d" % (vmFname)
    cmd += " | sfput n2=%d n3=%d d1=%f d2=%f d3=%f o3=0" % (ny, nx, dz, dy, dx)
    cmd += " --out=stdout > %s" % (vmRsf)
    os.system(cmd)
    if VTI:
        cmd = "< %s sfsegyread read=d" % (epsFname)
        cmd += " | sfput n2=%d n3=%d d1=%f d2=%f d3=%f o3=0" % (ny, nx, dz, dy, dx)
        cmd += ' | sfmath v=%s output="v * (1 + input)"' % (vmRsf)
        cmd += " --out=stdout > %s" % (vxRsf)
        os.system(cmd)
        cmd = "< %s sfsegyread read=d" % (delFname)
        cmd += " | sfput n2=%d n3=%d d1=%f d2=%f d3=%f o3=0" % (ny, nx, dz, dy, dx)
        cmd += ' | sfmath v=%s x=%s output="(x/v - 1 - input) / (1 + 2*input)"' % (vmRsf, vxRsf)
        cmd += " --out=stdout > %s" % (etaRsf)
        os.system(cmd)
data = np.array([nx, ny, nz, dz])
MPI.COMM_WORLD.Bcast(data, root=0)
nx = int(data[0])
ny = int(data[1])
nz = int(data[2])
dz = data[3]
# All processes obtain a copy of the model in scratch
if srnk == 0:
    os.system("cp %s %s" % (vmRsf, DATAPATH))
    if VTI:
        os.system("cp %s %s" % (vxRsf, DATAPATH))
        os.system("cp %s %s" % (etaRsf, DATAPATH))
MPI.COMM_WORLD.barrier()
# Clean up the copy on shared storage
if rank == 0:
    os.remove(vmRsf)
    if VTI:
        os.remove(vxRsf)
        os.remove(etaRsf)
vmRsf = os.path.basename(vmRsf)
if VTI:
    vxRsf = os.path.basename(vxRsf)
    etaRsf = os.path.basename(etaRsf)
MPI.COMM_WORLD.barrier()

# MPI domain decomposition of traveltime centres
nxCen = (nx - 1) // xSkip + 1
nyCen = (ny - 1) // ySkip + 1
N = nxCen * nyCen
myN = N // size
myN0 = myN * rank
if rank < N % size:
    myN = myN + 1
    myN0 += rank
else:
    myN0 += N % size

# The main loop
bs = 0 # Will store size of a single centred block of traveltimes
tmpRsf = "tmpRsf_%d" % (rank)   # scratch file for 1 traveltime block
headRsf = "headRsf_%d" % (rank) # scratch file for 1 block's segy headers
tmpSegy = "tmpSegy_%d" % (rank) # scratch file accumulating local blocks to segy
for n in range(myN0, myN0 + myN): # loop over local traveltime centres / blocks
    g = (n // nyCen) * xSkip + 1  # Convert MPI index to model grid 
    h = (n % nyCen) * ySkip + 1
    if rank == 0 and n%10 == 1:   # print progress
        print("%d%% done (%d of %d)" % (n / myN * 100, n, myN))
        sys.stdout.flush()

    # Calculate parameters for Madagascar commands for this block
    xc = (g - 1) * dx
    x0 = max(xc - rx, 0)
    x1 = min(xc + rx, (nx - 1) * dx)
    pad0x = max(int(rx/dx) - (g - 1), 0)
    pad1x = max(int(rx/dx) - (nx - g), 0)
    yc = (h - 1) * dy
    y0 = max(yc - ry, 0)
    y1 = min(yc + ry, (ny - 1) * dy)
    pad0y = max(int(ry/dy) - (h - 1), 0)
    pad1y = max(int(ry/dy) - (ny - h), 0)

    # Run Madagascar: Extract VM subset, run Eikonal, convert units, and pad if needed
    if VTI:
        os.system("< %s sfwindow min2=%f max2=%f min3=%f max3=%f squeeze=n > vx%d" 
                % (vxRsf, y0, y1, x0, x1, rank))
        os.system("< %s sfwindow min2=%f max2=%f min3=%f max3=%f squeeze=n > et%d" 
                % (etaRsf, y0, y1, x0, x1, rank))
    cmd = "< %s sfwindow min2=%f max2=%f min3=%f max3=%f squeeze=n" % (vmRsf, y0, y1, x0, x1)
    if VTI:
        cmd += (" | sfeikonalvti vx=vx%d eta=et%d yshot=%f xshot=%f 2>/dev/null" 
                % (rank, rank, yc, xc))
    else:
        cmd += " | sfeikonal yshot=%f xshot=%f 2>/dev/null" % (yc, xc)
    cmd += ' | sfmath output="input*1e3"'
    cmd += " | sfpad beg2=%d end2=%d beg3=%d end3=%d" % (pad0y, pad1y, pad0x, pad1x)
    cmd += " > %s" % (tmpRsf)
    if DEBUG: print(cmd)
    os.system(cmd)

    # Create segy headers for this block
    nj = int(ry / dy) * 2 + 1
    ni0 = g - int(rx / dx)
    nj0 = h - int(ry / dy)
    cmd = "< %s sfsegyheader" % (tmpRsf)
    cmd += " | sfheadermath key=iline output=%d" % (g)
    cmd += " | sfheadermath key=xline output=%d" % (h)
    cmd += ' | sfheadermath key=swdep output="%d+N/%d"' % (ni0, nj)
    cmd += ' | sfheadermath key=gwdep output="%d+N%%%d"' % (nj0, nj)
    cmd += " > %s" % (headRsf)
    if DEBUG: print(cmd)
    os.system(cmd)

    # Create segy data (headers and traces) for this block
    cmd = "< %s sfsegywrite tfile=%s | tail -c +3601 >> %s" % (tmpRsf, headRsf, tmpSegy)
    if DEBUG: print(cmd)
    os.system(cmd)
   
    # Store the size of a single traveltime block in the segy file
    if (n == myN0):
        bs = os.path.getsize(tmpSegy)
    
# All blocks done, so prepare a segy text and reel header on disk
finalSegy = SHAREDPATH + "TT.segy"
if rank == 0:
    cmd = "< %s sfsegywrite tfile=%s > %s" % (tmpRsf, headRsf, finalSegy)
    if DEBUG: print(cmd)
    os.system(cmd)
    cmd = "truncate -s 3600 %s" % (finalSegy)
    if DEBUG: print(cmd)
    os.system(cmd)
MPI.COMM_WORLD.barrier()

# Write the local segy scratch data to the shared global file
offset = 3600 + bs * myN0
outmode = MPI.MODE_WRONLY | MPI.MODE_APPEND
myOut = MPI.File.Open(MPI.COMM_WORLD, finalSegy, outmode)
CHUNK = 1<<30
if myN > 0:
    myIn = open(tmpSegy, "rb")
    bytesIn = myIn.read(CHUNK)
    while bytesIn:
        myOut.Write_at(offset, bytesIn)
        offset += len(bytesIn)
        bytesIn = myIn.read(CHUNK)
    myIn.close()
myOut.Close()
MPI.COMM_WORLD.barrier()

# Clean up the scratch data
if srnk == 0:
    shutil.rmtree(DATAPATH, ignore_errors=True)
