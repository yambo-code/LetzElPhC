from netCDF4 import Dataset
import numpy as np
from time import time

np.set_printoptions(suppress=False)

file_path_original = 'nc.temp'

file_path_test = 'nc.temp_2s'

file_path_test_sub = 'nc.temp_2ssub'

def get_data(filename):
    db = Dataset(filename)
    db.set_always_mask(False)
    db = db['exc_elph'][...]
    return db[...,:]


time1 = time()

original = get_data(file_path_original)
test = get_data(file_path_test)

test_sub =get_data(file_path_test_sub)
original_sub = original[:,1:3,173:356:7, 324:894:9,:]#[:, 1:2:5, 13:23:3, 17: 64:1]
#print(original_sub.shape)

print('\n ********* Time : ',time()-time1, '*********** \n')
print(np.allclose(test,original))
print(np.allclose(original_sub,test_sub))

time1 = time()
#time1 = time()
print(test.ravel().shape[0])
print("{:e}".format(test[0,2,999,999,0]))
print("{:e}".format(test[0,1,348,749,1]))

print('\n ********* Time : ',time()-time1, '*********** \n')
time1 = time()

reshaped = test.reshape(6,500,4,125,4)
print('Reshape: {:e}'.format(reshaped[3,127,1,73,2]))

transposed = reshaped.transpose(4,2,0,1,3)
print('Transpose: {:e}'.format(transposed[2,1,3,273,39]))


sliced_arr = transposed[2:4,0,3:5,179:500:3,53:117:5]
print('Sliced: {:e}'.format(sliced_arr[1,0,37,7]))

stripped_arr = sliced_arr[1,0,:,:]
print('Stripped: {:e}'.format(stripped_arr[13,5]))

print('\n ********* Time : ',time()-time1, '*********** \n')
time1 = time()

matmul = test[0,1,300,:,:]@test[0,2,579,:,:].T
einsum= np.einsum('ijkkp,ijlkp->ijl',test,test)

#print('Sliced: {:e}'.format())

print('matmul: ',format(matmul[239,739]))
print('einsum: ',format(einsum[0,2,473]))
print('einsum 2: ',format(np.einsum('ijkkp,ijlkp->',test,test)))
print('\n ********* Time : ',time()-time1, '*********** \n')
