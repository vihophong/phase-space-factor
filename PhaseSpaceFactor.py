
# coding: utf-8

# In[1]:


import ROOT
import matplotlib.pyplot as plt
import numpy as np
import sys


# ### Read csv file with parameters and put it into an array of dictionary

# In[2]:


tree = ROOT.TTree("tree","tree");
tree.ReadFile("f_table_wilkinson74.csv");
Z = np.zeros(1,dtype=int)
A = np.zeros(1,dtype=int)
An = np.zeros(20)
tree.SetBranchAddress("Z",Z);
tree.SetBranchAddress("A",A);
tree.SetBranchAddress("An",An);
ftable = []
for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    An_c=[]
    Z_c=Z[0]
    A_c=A[0]
    for j in An:
        An_c.append(j)
    ftable.append({"Z": Z_c,"A": A_c,"An": An_c})
#np.save("f_table_wilkinson74.npy",ftable)


# ### Load .npy file and create a function

# In[3]:


#ftable  = np.load("f_table_wilkinson74.npy",allow_pickle='TRUE')


# In[7]:


def f_func(x,p):
    Z = int(round(p[0]))
    E = p[1]-x[0]
    W = E/0.510998918+1 #unit of mec2
    fZ0 = 1./60.*(2*W**4-9*W**2-8)*np.sqrt(W**2-1) +1/4*W*np.log(W+np.sqrt(W**2-1))
    An = ftable[Z-2]["An"]
    An_i=[]
    if (E>0 and E<0.079):
        An_i.append(An[0])
        An_i.append(An[1])
        An_i.append(An[2])
        An_i.append(An[3])
    elif (E>=0.079 and E<0.501):
        An_i.append(An[4])
        An_i.append(An[5])
        An_i.append(An[6])
        An_i.append(An[7])
    elif (E>=0.501 and E<3.162):
        An_i.append(An[8])
        An_i.append(An[9])
        An_i.append(An[10])
        An_i.append(An[11])
    elif (E>=3.162 and E<12.589):
        An_i.append(An[12])
        An_i.append(An[13])
        An_i.append(An[14])
        An_i.append(An[15])
    elif (E>=12.589 and E<25.044):
        An_i.append(An[16])
        An_i.append(An[17])
        An_i.append(An[18])
        An_i.append(An[19])
    #print(An_i)
    S = np.exp(An_i[0]+An_i[1]*np.log(E)+An_i[2]*(np.log(E)**2)+An_i[3]*(np.log(E)**3))
    return fZ0*S          


# In[17]:


def f_func_num(x,p):
    return f_func(x,p)*(x[0]-p[2])


# In[24]:


Qb =14.464
Sn = 3.63119
ZZ = float(sys.argv[1])
print(sys.argv[0])
f_d = ROOT.TF1("f_func",f_func,0.,Qb,2);
f_d.SetParameter(0,ZZ)
f_d.SetParameter(1,Qb)

f_n = ROOT.TF1("f_func_num",f_func_num,0.,Qb,3);
f_n.SetParameter(0,ZZ)
f_n.SetParameter(1,Qb)
f_n.SetParameter(2,Sn)
print(f_n.Integral(Sn,Qb)/f_d.Integral(Sn,Qb))

