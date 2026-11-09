import numpy as np
from common import *
di,dj,dist=downstream(); order,rec=topo_order(di,dj,land)
assert len(order)==land.sum()
nr,nc=land.shape; outlet=np.full(nr*nc,-1)
for i in order[::-1]: outlet[i]= i if rec[i]<0 else outlet[rec[i]]
outlet=outlet.reshape(land.shape); outlet[~land]=-1
np.save(f'{OUT}/outlet.npy',outlet); print('land ha %.1f'%(land.sum()*cell*cell/1e4))
