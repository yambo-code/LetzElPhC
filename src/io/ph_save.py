#// Creates a ph_elph.save that is used by the main code. 
# It creates a netcdf file which contain eigen vectors, q points, their symmetries 
# and dVscfs for all q points 



import numpy as np
import xml.dom.minidom
from netCDF4 import Dataset 




single_prec = True
FFT_GRID = [45, 45,225] ### Dense grid, but local, wfcs are in smooth grid
prefix = 'wse2'
nmodes = 3*3
nspin_mag = 1
ph_dir = '.'

prefix = prefix.strip()
dvscf_file = ph_dir+"/_ph0/"+prefix+".dvscf1_1"

yambo_gkpp = 'elph_dir_full/elph_yambo_s.dbph_000001.h5'



patterns_file = ph_dir+'/_ph0/'+prefix+'.phsave/patterns.1.xml'
print("Reading Patterns ...")
### Read Patterns
pattern_vec = np.zeros((nmodes,nmodes,2))

pattern_xml = xml.dom.minidom.parse(patterns_file)

root = pattern_xml.getElementsByTagName('Root')[0]

irrep_info = root.getElementsByTagName('IRREPS_INFO')[0]

Q_point = int(irrep_info.getElementsByTagName('QPOINT_NUMBER')[0].firstChild.data)

nsymq = int(irrep_info.getElementsByTagName('QPOINT_GROUP_RANK')[0].firstChild.data)

minus_q = irrep_info.getElementsByTagName('MINUS_Q_SYM')[0].firstChild.data

nirr = int(irrep_info.getElementsByTagName('NUMBER_IRR_REP')[0].firstChild.data)

imode0 = 0

#####
for i in range(nirr):

    rep_temp = irrep_info.getElementsByTagName('REPRESENTION.'+str(i+1))[0]
    
    irrep = int(rep_temp.getElementsByTagName('NUMBER_OF_PERTURBATIONS')[0].firstChild.data)
    
    for j in range(irrep):
    
        jpert = rep_temp.getElementsByTagName('PERTURBATION.'+str(j+1))[0]
    
        pattern_str = jpert.getElementsByTagName('DISPLACEMENT_PATTERN')[0].firstChild.data
    
        pattern_vec[imode0,:,:] = np.fromstring(pattern_str, sep = ' ' ).reshape(-1,2)
    
        imode0 += 1

if imode0 != nmodes :
    exit("Total pertubations doesnot match with total modes")

if (minus_q.lower() == "true"):
    minus_q = True
else:
    minus_q = False
####
pattern_vec = pattern_vec[...,0] + 1j*pattern_vec[...,1]

pattern_vec = pattern_vec.reshape(nmodes,-1,3) #(modes,natom,3,)

#//pattern_vec = pattern_vec.transpose(0,2,1)
print(pattern_vec.shape)
## Read dVscf ######
print("Reading dVscf ...")

f = open(dvscf_file,'rb')

Record = np.dtype((np.double, 2*nspin_mag*np.product(FFT_GRID)))

dVscf = np.fromfile(f,dtype=Record, count=nmodes).astype(np.double).reshape(nmodes,nspin_mag,FFT_GRID[2],FFT_GRID[1],FFT_GRID[0],2) # Real array

dVscf = dVscf[...,0] + 1j*dVscf[...,1] #(nmodes,nspin_mag,Fz,Fy,Fx)

dVscf = dVscf.transpose(0,1,4,3,2) ## in pattern basis ( now the shape is (nmodes,nspin_mag,Fx,Fy,Fz))

# convert to carsian basis
dVscf = np.einsum('vsijk, vax-> sijkax',dVscf,np.conj(pattern_vec),optimize=True) ##(nspin_mag,Fx,Fy,Fz,natom,3) ## screened Electric field of 
print(dVscf.shape)

# dVscf = (\epsilon - 1 )dVbare

def get_nc(file,varname,complx):
    aa = Dataset(file,mode='r')
    aa = aa[varname][...]
    if complx:
        aa = aa[...,0]+1j*aa[...,1]
    return aa

pol_vecs = get_nc(yambo_gkpp,"POLARIZATION_VECTORS",True).data.transpose(2,1,0) ##(pol,atom,modes) -> (modes,atom,pol)

print(pol_vecs.shape)


### Write to netCDF file 

ncfile = Dataset('nc.dVscf_new',mode='w',format='NETCDF4_CLASSIC')

# ---
nq = ncfile.createDimension('nqpts', 1) # // FIX ME, this needs to be changed to nqpts
nmag = ncfile.createDimension('nmag', dVscf.shape[0])
Nx = ncfile.createDimension('Nx', FFT_GRID[0])
Ny = ncfile.createDimension('Ny', FFT_GRID[1])
Nz = ncfile.createDimension('Nz', FFT_GRID[2])
natom = ncfile.createDimension('natom', dVscf.shape[4])
xcart = ncfile.createDimension('xcart', dVscf.shape[5])
re_im = ncfile.createDimension('re_im', 2)
nmode= ncfile.createDimension('nmodes', pol_vecs.shape[0])
const_vals = ncfile.createDimension('single_val', 1)

write_type = np.single
if (not single_prec):
    write_type = np.double


dvscf = ncfile.createVariable('dVscfs', write_type, ('nqpts','nmodes','nmag','Nx','Ny','Nz','re_im'))
ph_modes = ncfile.createVariable('ph_pol_vec', write_type, ('nqpts', 'nmodes', 'natom','xcart','re_im' ))

fft_dims = ncfile.createVariable('fft_dims',np.intc, ('xcart'))
nq_pts   = ncfile.createVariable('qpts',np.intc, ('single_val'))

dVscf = np.einsum('sijkax,vax->vsijk',dVscf,pol_vecs)#.reshape(nmodes,nspin_mag,-1)
dvscf[0,...,0] = dVscf.real
dvscf[0,...,1] = dVscf.imag
fft_dims[:] = np.array(FFT_GRID)
nq_pts[0] = 1 ##// FIX ME, this needs to be changed to nqpts

ph_modes[0,...,0] = pol_vecs.real
ph_modes[0,...,1] = pol_vecs.imag
## ----



ncfile.close()


#np.save('Dvscf',dVscf)

#print('Done .')
