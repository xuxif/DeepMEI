import numpy as np;
import random
def sample_read(img_row_num_ref,img_row_num_r,img_row_num_l,shuffle_depth):
    img_row_num_ref_range=np.array(range(5,5+img_row_num_ref))
    img_row_num_r_range=np.array(range(img_row_num_ref_range[-1],img_row_num_ref_range[-1]+img_row_num_r))
    img_row_num_l_range=np.array(range(img_row_num_r_range[-1],img_row_num_r_range[-1]+img_row_num_l))

    while (len(img_row_num_ref_range)+len(img_row_num_r_range)+len(img_row_num_l_range)) >shuffle_depth:
      if len(img_row_num_ref_range) > (len(img_row_num_r_range)+len(img_row_num_l_range)):
        img_row_num_ref_range=np.delete(img_row_num_ref_range,0,0)
        img_row_num_ref_range=np.delete(img_row_num_ref_range,-1,0)
      else:
          if len(img_row_num_r_range)>len(img_row_num_l_range):
            img_row_num_r_range=np.delete(img_row_num_r_range,0,0)
            img_row_num_r_range=np.delete(img_row_num_r_range,-1,0)
          else:
            img_row_num_l_range=np.delete(img_row_num_l_range,0,0)
            img_row_num_l_range=np.delete(img_row_num_l_range,-1,0)
    rs_t=np.concatenate((img_row_num_ref_range,img_row_num_r_range,img_row_num_l_range),axis=0)
    return rs_t
rs=sample_read(70,30,40,95);
rs_t= np.array(random.sample(range(5, 150), 95))
print(rs)
print(rs_t)
