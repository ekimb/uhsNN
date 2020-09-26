import numpy as np
import pandas as pd
import itertools
import os
mink=5
maxk=9
maxmaxk=15

total_labels=np.empty(0,dtype=np.int8)
total_lst=np.empty([0,2*maxmaxk+2],dtype=np.int8)

for k in range (mink,maxk+1,1):
	print(k)
	decycling=pd.read_csv("decyc"+str(k)+".int")
	lst = np.array(list(itertools.product([0, 1], repeat=2*k)))
	for L in range(20,201,1):
		print(L)
		lst = np.array(list(itertools.product([0, 1], repeat=2*k)))
		padlst = np.pad(lst, ((0, 0),(0,maxmaxk*2-lst.shape[1])), 'constant', constant_values=((2, 2),(2,2)))
		lstL = np.append(np.full(lst.shape[0], L).reshape(lst.shape[0],1), padlst, axis=1)
		lstkL = np.append(np.full(lst.shape[0], k).reshape(lstL.shape[0],1), lstL, axis=1)
		labels=np.zeros(np.power(4,k), dtype=np.int8)
		labels[decycling]=2
		print("PDOCKS"+str(k)+str(L)+".int")
		print(os.path.isfile("PDOCKS"+str(k)+str(L)+".int"))
		if (os.path.isfile("PDOCKS"+str(k)+str(L)+".int")):
			docks=np.array(pd.read_csv("PDOCKS"+str(k)+str(L)+".int", header=None))
			labels[docks]=1
		combined = np.append(lstkL, labels.reshape(labels.shape[0],1), axis=1)
		combined = combined[np.logical_not(combined[:,2*maxmaxk+2]==2)]
		total_labels=np.append(total_labels, combined[:,2*maxmaxk+2],axis=0)
		total_lst=np.append(total_lst, combined[:,0:2*maxmaxk+2],axis=0)

print(total_labels.shape)
print(total_lst.shape)
from keras.models import Sequential
from keras.layers import Dense,Masking,LSTM,GRU
from keras import metrics

model = Sequential()
model.add(Masking(mask_value=2., input_shape=(maxmaxk*2+2, 1)))
model.add(GRU(100, input_dim=1))
#model.add(Dense(100, input_dim=2*k+1, activation='relu'))
model.add(Dense(50, activation='relu'))
model.add(Dense(1, activation='sigmoid'))
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=[metrics.AUC(),'accuracy'])
model.fit(total_lst, total_labels, epochs=10, batch_size=1024)

model.save("PASHAmodel_uhs_L_k.h5")
