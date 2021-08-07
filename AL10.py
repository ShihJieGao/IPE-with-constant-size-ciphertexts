# -*- coding: utf-8 -*-
"""
Created on Sun May 16 05:37:17 2021

@author: j23793276
"""
"""
from timeit import default_timer as timer
"""
from charm.toolbox.pairinggroup import PairingGroup,ZR,G1,GT,pair
from charm.toolbox.ABEnc import ABEnc
import time

pk_t = { 'g':G1, 'g^alpha0':G1, 'g^alphavec':G1, 'e_gg^alpha':GT }
mk_t = {'g^alpha':G1 }
sk_t = { 'D_0':G1, 'D_1':G1, 'K_i':G1,}
ct_t = { 'E_0':GT, 'E_1':G1, 'E_2':G1 }
X = []
Y = []

class IPE_AL10(ABEnc):
    def __init__(self, groupObj, vecSize, verbose=False):
        ABEnc.__init__(self)
        global group, vec_size
        group = groupObj
        vec_size = vecSize
    
    
    def setup(self):
        g = group.random(G1)
        alpha, alpha0 =  group.random(ZR), group.random(ZR)
        g.initPP()

        A = []
        for i in range(vec_size):
            A.append(g ** group.random(ZR))
        g_alpha0 = g ** alpha0
        e_gg_alpha = pair(g, g ** alpha)
        g_alpha = g ** alpha
        
        pk = {'g':g, 'g^alpha0':g_alpha0, 'g^A':A, 'e_gg^alpha':e_gg_alpha}
        mk = {'g^alpha':g_alpha}
        return (pk, mk)
        
    def keygen(self, pk, mk, X):
            
        if X[0]==0: 
            print("X[0] can not be 0")
            return
        t =  group.random(ZR)
        D_0 = pk['g'] ** t 
        D_1 = mk['g^alpha'] * (pk['g^alpha0'] ** t)
        K = []
        for i in range(vec_size - 1):
            
            tmp0 = (X[i+1]/X[0])*(-1)
            tmp1 = pk['g^A'][0] ** tmp0
            tmp2 = tmp1 * pk['g^A'][i+1]
            tmp = tmp2 ** t
            K.append(tmp) 
            
        sk = {'D_0':D_0, 'D_1':D_1,'K':K}
        return sk
    
    def encrypt(self, pk , M, Y):
       
        s =  group.random(ZR)
        E_0 = M * (pk['e_gg^alpha'] ** s)
        prod = 1
        for i in range(vec_size):
            prod = prod * (pk['g^A'][i] ** Y[i])
        E_1 = (pk['g^alpha0'] * prod) ** s
        E_2 = pk['g'] ** s
        
        ct = {'E_0':E_0, 'E_1':E_1, 'E_2':E_2}
        return ct
        
    def decrypt(self, pk, sk, ct, Y):
        prodOfK = 1
        for i in range(vec_size - 1):
            prodOfK = prodOfK * (sk['K'][i] ** Y[i+1])
            
        numerator = pair(sk['D_1'] * prodOfK, ct['E_2'])
        denominator = pair(ct['E_1'], sk['D_0'])
        
        maskTerm = numerator/denominator
        
        return ct['E_0']/maskTerm
        
def dot(V1,V2):
    if len(V1) != len(V2):
        return -1
    return sum(i[0] * i[1] for i in zip(V1,V2))        
        
def main(): 
    groupObj = PairingGroup('SS512')
    vecSize = 4
    ipe = IPE_AL10(groupObj, vecSize)
    
    startSetup = time.time()
    (pk, mk) = ipe.setup()
    endSetup = time.time()
    print("\npk :=>", pk) 
    print("\nmk :=>", mk)
    
    for i in range(vec_size):
            X.append(group.random(ZR))
    
    startKeygen = time.time()    
    sk = ipe.keygen(pk,mk,X)
    endKeygen = time.time()
    print("\nsk :=>", sk)
    
    sum = 0
    for i in range(len(X)-1):
        tmp = group.random(ZR)
        sum = sum + X[i] * tmp
        Y.append(tmp)
    Y.append(-sum/X[len(X)-1])
    #print(dot(X,Y))    
    
    M = groupObj.random(GT)
    print("\nM :=>", M)
    
    startEnc = time.time()
    ct = ipe.encrypt(pk,M,Y)
    endEnc = time.time()
    print("\nct :=>", ct)
    
    startDec = time.time()
    rec_M = ipe.decrypt(pk,sk,ct,Y)
    endDec = time.time()
    print("\nRec_M =>", rec_M)
    
    print("\nSetupTime:",endSetup - startSetup)
    print("KeygenTime:",endKeygen - startKeygen)
    print("EncTime:",endEnc - startEnc)
    print("DecTime:",endDec - startDec)    
    
if __name__ == "__main__":
    main()