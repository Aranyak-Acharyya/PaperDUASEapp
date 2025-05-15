#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
data = np.load('X1_39.npz')['arr_0'] #code to record 768-d vector answers for all LLMs, at all times, per replicate


# In[2]:


v=data.shape
print(v)


# In[3]:


n_time=v[0]
n_model=v[1]
n_question=v[2]
n_component=v[3]


# In[222]:


D=np.zeros((1,n_component))
for time in range(n_time):
    Dt=data[time]
    Dtt=np.concatenate(Dt,axis=0)
    D=np.concatenate((D,Dtt),axis=0)
D=np.delete(D,0,axis=0)
D.shape


# In[223]:


import pandas as pd
df=pd.DataFrame(D)
df.to_csv('LLMdnew1_39.csv')


# In[ ]:





# In[79]:


S=np.concatenate((S,ll),axis=1)


# In[ ]:


print(S)


# In[ ]:




