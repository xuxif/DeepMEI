import numpy as np
np_del=[0,1,2]
pad_range=[20,21,19,22,18,23,17,24,16,25,15,26,14,27,13,28,12,29,11,30,10,31,9,32,8,33,7,34,6,35,5,36,4,37,3,38,2,39,1]
print(pad_range)
pad_range=np.delete(pad_range,np_del,0)
print(pad_range)
