#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# In[3]:


dataset = pd.read_csv('/Users/aynaz/Desktop/FinalAstma.csv')


# In[5]:


X = dataset.iloc[:, 1:11].values
y = dataset.iloc[:, -3].values


# In[6]:


from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 0)


# In[19]:


from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
poly_reg = PolynomialFeatures(degree = 2)
X_poly = poly_reg.fit_transform(X_train)
regressor = LinearRegression()
regressor.fit(X_poly, y_train)


# In[20]:


y_pred = regressor.predict(poly_reg.transform(X_test))
np.set_printoptions(precision=2)
print(np.concatenate((y_pred.reshape(len(y_pred),1), y_test.reshape(len(y_test),1)),1))


# In[21]:


from sklearn.metrics import r2_score
r2_score(y_test, y_pred)


# In[ ]:




