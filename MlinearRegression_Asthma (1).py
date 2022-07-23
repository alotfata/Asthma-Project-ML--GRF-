#!/usr/bin/env python
# coding: utf-8

# In[26]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# In[27]:


dataset = pd.read_csv('/Users/aynaz/Desktop/FinalAstma.csv')


# In[21]:


X = dataset.iloc[:, 1:11].values
y = dataset.iloc[:, -3].values


# In[22]:


from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 0)


# In[23]:


from sklearn.linear_model import LinearRegression
regressor = LinearRegression()
regressor.fit(X_train, y_train)


# In[24]:


y_pred = regressor.predict(X_test)
np.set_printoptions(precision=2)
print(np.concatenate((y_pred.reshape(len(y_pred),1), y_test.reshape(len(y_test),1)),1))


# In[25]:


from sklearn.metrics import r2_score
r2_score(y_test, y_pred)


# In[ ]:





# In[ ]:




