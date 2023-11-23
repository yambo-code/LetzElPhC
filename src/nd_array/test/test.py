from netCDF4 import Dataset
import numpy as np
from time import time

np.set_printoptions(suppress=False)

file_path_original = 'nc.temp'

file_path_test = 'nc.temp2'

file_path_test_sub = 'nc.temp_sub'

def get_data(filename):
    db = Dataset(filename)
    db.set_always_mask(False)
    db = db['exc_elph'][...]
    return db[...,0] + 1j*db[...,1]


time1 = time()

original = get_data(file_path_original)
test = get_data(file_path_test)

test_sub =get_data(file_path_test_sub)
original_sub = original[:,1:3,173:356:7, 324:894:9]#[:, 1:2:5, 13:23:3, 17: 64:1]
#print(original_sub.shape)

print('\n ********* Time : ',time()-time1, '*********** \n')
print(np.allclose(test,original))
print(np.allclose(original_sub,test_sub))

time1 = time()
#time1 = time()
print(test.ravel().shape[0])
print("{:e}".format(test[0,2,999,999]))
print("{:e}".format(test[0,1,348,749]))

print('\n ********* Time : ',time()-time1, '*********** \n')
time1 = time()

reshaped = test.reshape(6,500,2,125,4)
print('Reshape: {:e}'.format(reshaped[3,127,1,73,2]))

transposed = reshaped.transpose(4,2,0,1,3)
print('Transpose: {:e}'.format(transposed[2,1,3,273,39]))


sliced_arr = transposed[2:4,0,3:5,179:500:3,53:117:5]
print('Sliced: {:e}'.format(sliced_arr[1,0,37,7]))

stripped_arr = sliced_arr[1,0,:,:]
print('Stripped: {:e}'.format(stripped_arr[13,5]))

print('\n ********* Time : ',time()-time1, '*********** \n')
time1 = time()

matmul = test[0,1,:,:].T@np.conj(test[0,2,:,:].T)
einsum= np.einsum("ijkk,ijlk->ijl",test,test)

#print('Sliced: {:e}'.format())

print('matmul: {:e}'.format(matmul[239,739]))
print('einsum: {:e}'.format(einsum[0,2,473]))
print('\n ********* Time : ',time()-time1, '*********** \n')
