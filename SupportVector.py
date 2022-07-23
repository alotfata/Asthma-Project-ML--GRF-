#!/usr/bin/env python
# coding: utf-8

# In[1]: Import libraries


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# In[2]:  Upload dataset


dataset = pd.read_csv('/Users/aynaz/Desktop/FinalAstma.csv')


# In[27]: Independent and dependent variables


X = dataset.iloc[:, 1:11].values
y = dataset.iloc[:, -3].values


# In[28]: reshape indepednt variables 


y = y.reshape(len(y),1)


# In[29]:


y


# In[30]: split data to test and train


from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 0)


# In[31]:


from sklearn.preprocessing import StandardScaler
sc_X = StandardScaler()
sc_y = StandardScaler()
X_train = sc_X.fit_transform(X_train)
y_train = sc_y.fit_transform(y_train)


# In[ ]:





# In[38]:


from sklearn.svm import SVR
regressor = SVR(kernel = 'rbf')
regressor.fit(X_train, y_train)


# In[39]: prediction on test data 


y_pred = sc_y.inverse_transform(regressor.predict(sc_X.transform(X_test)))
np.set_printoptions(precision=2)
print(np.concatenate((y_pred.reshape(len(y_pred),1), y_test.reshape(len(y_test),1)),1))


# In[40]:


from sklearn.metrics import r2_score
r2_score(y_test, y_pred)


# In[ ]:




