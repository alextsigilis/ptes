# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 09:34:09 2023

@author: sted
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, KFold
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler, OrdinalEncoder, OneHotEncoder, LabelEncoder
from sklearn.metrics import cohen_kappa_score, confusion_matrix, ConfusionMatrixDisplay
from sklearn.decomposition import PCA
 

# Data Loading
data = pd.read_csv('features.csv');
data = data.dropna()


# !!!! Used for Classification N1 vs. N2
#data = data[(data.annotation=='Sleep stage N1') & (data.annotation=='Sleep stage N2')]

# Spliting Features from class label
X = data.loc[:, data.columns != 'annotation']
y = data.loc[:, data.columns == 'annotation']

# Scale & Encode the data
le = LabelEncoder()
ohe = OneHotEncoder(sparse=False)
X = StandardScaler().fit_transform(X)
y = le.fit_transform(y).reshape((y.shape[0], ))

# Apply PCA
pca = PCA(n_components=5)
X = pca.fit_transform(X)

# Split training and testing data
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.15, stratify=y)

# Instatiate the MLP
mlp = MLPClassifier(hidden_layer_sizes=[64, 64, 64, 16, 8, ],
                    solver='sgd',
                    learning_rate='adaptive',
                    max_iter=300,
                    warm_start=True,
                    verbose=True)


# Do k-Fold 
kf = KFold(n_splits=2)
train_loss = []
train_kappa = []
validation_kappa = []

for i, (train, valid) in enumerate(kf.split(X_train, y_train)):
    print(f"|=====> ITERATION {i+1} <=====|")
    mlp.fit(X_train[train, :], y_train[train])
    train_loss.append(mlp.loss_)
    y_pred_train = mlp.predict(X_train[train, :])
    y_pred = mlp.predict(X_train[valid,:])
    train_kappa.append(cohen_kappa_score(y_train[train], y_pred_train))
    validation_kappa.append(cohen_kappa_score(y_train[valid], y_pred))
    print(f"Training kappa: {train_kappa[-1]}")
    print(f"Validation kappa: {validation_kappa[-1]}")
    print("=========================================\n")
    
# Compute the Confusion matrix
y_pred = mlp.predict(X_test)
C = confusion_matrix(y_test, y_pred)
cmd = ConfusionMatrixDisplay(C, display_labels=le.inverse_transform([0,1,2,3,4]))
cmd.plot()
plt.plot()