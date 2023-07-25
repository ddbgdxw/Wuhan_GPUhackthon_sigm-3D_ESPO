<<<<<<< HEAD
import scipy.io as io
import numpy as np
import cupy as cp
import scipy
import time

data=io.loadmat('./data/TSPolInCM_filter_1.mat')
# data=io.loadmat('./data/TSPolInCM_filter_3.mat')
# data=io.loadmat('./data/TSPolInCM_filter_9.mat')
# data=io.loadmat('./data/TSPolInCM_filter_15.mat')
TSPolInCM_filter=data['data']
print(TSPolInCM_filter.shape)
# 数据格式为：6*6*m*n*N
# print(type(TSPolInCM_filter))
TSPolInCM_filter = cp.asarray(TSPolInCM_filter)
# print(type(TSPolInCM_filter))

def favecoh_fullpol(x,pTST11,pTST12,pTST22,nintf):
    # 计算三个ω值，x为求解的4个参数值
    _omegaconj=np.zeros((3,1),dtype='complex_')
    [_omegaconj[0],_omegaconj[1],_omegaconj[2]]=[np.conjugate(np.cos(x[0])),np.conjugate(np.sin(x[0])*np.cos(x[1])*np.exp(x[2]*1j)),np.conjugate(np.sin(x[0])*np.sin(x[1])*np.exp(x[3]*1j))]
    _omegaconj=_omegaconj.T
    
    [_npol,qq,_nintf]=pTST11.shape
    _omegaSet=np.kron(np.ones((1,nintf)),np.conjugate(_omegaconj))

    pTST11opt=np.multiply(np.dot(_omegaconj,np.reshape(pTST11,[_npol,_npol*_nintf],order='F')),_omegaSet)
    pTST12opt=np.multiply(np.dot(_omegaconj,np.reshape(pTST12,[_npol,_npol*_nintf],order='F')),_omegaSet)
    pTST22opt=np.multiply(np.dot(_omegaconj,np.reshape(pTST22,[_npol,_npol*_nintf],order='F')),_omegaSet)
    pTST11opt=np.sum(np.reshape(pTST11opt,[pTST11opt.shape[0],_nintf,_npol],order='C'),axis=2)
    pTST12opt=np.sum(np.reshape(pTST12opt,[pTST11opt.shape[0],_nintf,_npol],order='C'),axis=2)
    pTST22opt=np.sum(np.reshape(pTST22opt,[pTST11opt.shape[0],_nintf,_npol],order='C'),axis=2)
    TSCoh=np.divide(pTST12opt,np.sqrt(np.multiply(pTST11opt,pTST22opt)))
    
    meanTScoh=np.mean(np.abs(TSCoh),axis=1)
    return -meanTScoh

def ESPO_TSPolDS_1(TSPolInCM,interdeg=10):
    [q,qq,nrow,ncol,nintf]=TSPolInCM.shape
    # step1
    # 计算枚举4个参数的步长,并枚举4个参数（枚举需要包含右端点）
    npol=q//2
    num_1=90//interdeg
    num_2=360//interdeg
    if 90%interdeg==0:
        num_1+=1
        num_2+=1
    num_1=num_1*1j
    num_2=num_2*1j

    a=cp.mgrid[0:90:num_1,0:90:num_1,-180:180:num_2,-180:180:num_2]
    [aaa,bbb,ccc,ddd,eee]=a.shape
    sum=bbb*ccc*ddd*eee
    dff=cp.zeros((4,sum))
    b=cp.ravel(a[0],order='F')
    c=cp.ravel(a[1],order='F')
    d=cp.ravel(a[2],order='F')
    e=cp.ravel(a[3],order='F')
    dff[0]=b
    dff[1]=c
    dff[2]=d
    dff[3]=e
    para=dff*cp.pi/180

    # 使用4个参数计算ω值
    omegaconj=cp.zeros((3,sum),dtype='complex_')
    [omegaconj[0],omegaconj[1],omegaconj[2]]=[cp.conjugate(cp.cos(para[0,:])),cp.conjugate(cp.sin(para[0,:])*cp.cos(para[1,:])*cp.exp(para[2,:]*1j)),cp.conjugate(cp.sin(para[0,:])*cp.sin(para[1,:])*cp.exp(para[3,:]*1j))]
    omegaconj=omegaconj.T
    para=para.T
    # print(omegaconj.shape)
    # ω的数据格式：枚举数目×3
    omegaSet=cp.kron(cp.ones((1,nintf)),cp.conjugate(omegaconj))
    # print(omegaSet.shape)
    # ω的数据格式：枚举数目×（3*N）
    optpara=cp.zeros((nrow,ncol,4))

    TSPolInCMT11=TSPolInCM[0:npol,0:npol,:,:,:]
    TSPolInCMT12=TSPolInCM[0:npol,npol:npol*2,:,:,:]
    TSPolInCMT22=TSPolInCM[npol:npol*2,npol:npol*2,:,:,:]    
    # T矩阵数据格式：3×3×n×m×N
    # step2
    # 开始对每个像素进行处理
    for ii in range(nrow):
        for jj in range(ncol):

            pTST11=cp.squeeze(TSPolInCMT11[:,:,ii,jj,:])
            pTST12=cp.squeeze(TSPolInCMT12[:,:,ii,jj,:])
            pTST22=cp.squeeze(TSPolInCMT22[:,:,ii,jj,:])
            
            pTST11opt=cp.multiply(cp.dot(omegaconj,cp.reshape(pTST11,[npol,npol*nintf],order='F')),omegaSet)
            pTST12opt=cp.multiply(cp.dot(omegaconj,cp.reshape(pTST12,[npol,npol*nintf],order='F')),omegaSet)
            pTST22opt=cp.multiply(cp.dot(omegaconj,cp.reshape(pTST22,[npol,npol*nintf],order='F')),omegaSet)
            
            pTST11opt=cp.sum(cp.reshape(pTST11opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST12opt=cp.sum(cp.reshape(pTST12opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST22opt=cp.sum(cp.reshape(pTST22opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)

            TSCoh=cp.divide(pTST12opt,cp.sqrt(cp.multiply(pTST11opt,pTST22opt)))
            meanTScoh=cp.mean(cp.abs(TSCoh),axis=1)
            maxInd=cp.argmax(meanTScoh, axis=0)
            optpara[ii,jj,:]=para[maxInd,:]


    return optpara

def ESPO_TSPolDS_2(TSPolInCM,optparaDS0_lag3,interdeg):
    # optparaDS0_lag3是step2求解出的最优参数
    [q,qq,nrow,ncol,nintf]=TSPolInCM.shape
    npol=q//2
    optomegaconj=np.zeros((nrow,ncol,npol),dtype='complex_')
    optmeanCoh=np.zeros((nrow,ncol)) 
    optTSCoh=np.zeros((nrow,ncol,nintf),'complex_')
    TSPolInCMT11=TSPolInCM[0:npol,0:npol,:,:,:]
    TSPolInCMT12=TSPolInCM[0:npol,npol:npol*2,:,:,:]
    TSPolInCMT22=TSPolInCM[npol:npol*2,npol:npol*2,:,:,:]

    # step3
    # 非线性优化部分
    for ii in range(nrow):
        for jj in range(ncol):
            pTST11=np.squeeze(TSPolInCMT11[:,:,ii,jj,:])
            pTST12=np.squeeze(TSPolInCMT12[:,:,ii,jj,:])
            pTST22=np.squeeze(TSPolInCMT22[:,:,ii,jj,:])
            x0=optparaDS0_lag3[ii,jj,:].T
            lb0 = [0,0,-np.pi,-np.pi]
            ub0 = [np.pi/2,np.pi/2,np.pi,np.pi]
            lb1 = [x0[0],x0[1],x0[2],x0[3]]-np.deg2rad(interdeg)
            ub1 = [x0[0],x0[1],x0[2],x0[3]]+np.deg2rad(interdeg)
            lb=np.max([lb0,lb1],axis=0)
            ub=np.min([ub0,ub1],axis=0)
            bounds=((lb[0],ub[0]),(lb[1],ub[1]),(lb[2],ub[2]),(lb[3],ub[3]))
            fval=scipy.optimize.minimize(favecoh_fullpol,x0,args=(pTST11,pTST12,pTST22,nintf),bounds=bounds)
            x0=fval.x
            optmeanCoh[ii,jj]=-fval.fun
            fomegaconj=np.zeros((3,1),dtype='complex_')
            [fomegaconj[0],fomegaconj[1],fomegaconj[2]]=[np.conjugate(np.cos(x0[0])),np.conjugate(np.sin(x0[0])*np.cos(x0[1])*np.exp(x0[2]*1j)),np.conjugate(np.sin(x0[0])*np.sin(x0[1])*np.exp(x0[3]*1j))]
            fomegaconj=fomegaconj.T
            optomegaconj[ii,jj,:]=fomegaconj
            fomegaSet=np.kron(np.ones((1,nintf)),np.conjugate(fomegaconj))
            pTST11opt=np.multiply(np.dot(fomegaconj,np.reshape(pTST11,[npol,npol*nintf],order='F')),fomegaSet)
            pTST12opt=np.multiply(np.dot(fomegaconj,np.reshape(pTST12,[npol,npol*nintf],order='F')),fomegaSet)
            pTST22opt=np.multiply(np.dot(fomegaconj,np.reshape(pTST22,[npol,npol*nintf],order='F')),fomegaSet)
            pTST11opt=np.sum(np.reshape(pTST11opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST12opt=np.sum(np.reshape(pTST12opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST22opt=np.sum(np.reshape(pTST22opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            TSCoh=np.divide(pTST12opt,np.sqrt(np.multiply(pTST11opt,pTST22opt)))
            optTSCoh[ii,jj,:]=TSCoh

    s3=time.time()
    return optomegaconj,optTSCoh,optmeanCoh

s1=time.time()
optparaDS0=ESPO_TSPolDS_1(TSPolInCM_filter,30)
s2=time.time()

TSPolInCM_filter=cp.asnumpy(TSPolInCM_filter)
optparaDS0=cp.asnumpy(optparaDS0)


omegaconjDS,optTSCoh,optmeanCoh=ESPO_TSPolDS_2(TSPolInCM_filter,optparaDS0,30)
s3=time.time()
print('step1&step2 time:')
print(s2-s1,'s')
print('step3 time:')
print(s3-s2,'s')


def TSPolInCM2TSInCMSet(TSPolInCM,omegaconj):
    [q,qq,nrow,ncol,nintf]=TSPolInCM.shape
    optTSInCMSet=np.zeros((2,2,nrow,ncol,nintf),dtype='complex_')
    for ii in range(nrow):
        for jj in range(nrow):
            for kk in range(nintf):
                pomegaconj=np.kron(np.identity(2),(np.squeeze(omegaconj[ii,jj,:])).T)
                optTSInCMSet[:,:,ii,jj,kk]=np.dot(np.dot(pomegaconj,TSPolInCM[:,:,ii,jj,kk]),np.conjugate(pomegaconj.T))
    return optTSInCMSet
    
def Normalize_TSInCM(TSInCM):
    ndim=TSInCM.ndim
    if ndim>=4:
        [q,qq,nrow,ncol,npol]=TSInCM.shape
        norTSInCM=np.zeros((q,q,nrow,ncol,npol),dtype='complex_')
        for ii in range(npol):
            for kk in range(nrow):
                for ll in range(ncol):
                    TSIntensity=np.diag(TSInCM[:,:,kk,ll,ii])
                    TSIntensity=TSIntensity.reshape(q,1)
                    TSIntCM1=np.dot(TSIntensity,np.ones((1,q)))
                    TSIntCM2=np.ones((q,1))*TSIntensity.T
                    norTSInCM[:,:,kk,ll,ii]=TSInCM[:,:,kk,ll,ii]/np.sqrt(TSIntCM1*TSIntCM2)
    return norTSInCM


TSInCMSet_ESPODS=Normalize_TSInCM(TSPolInCM2TSInCMSet(TSPolInCM_filter,omegaconjDS))#ESPO


optCor=np.squeeze((TSInCMSet_ESPODS[0,1,:,:,]))
=======
<<<<<<< HEAD
import scipy.io as io
import numpy as np
import cupy as cp
import scipy
import time

data=io.loadmat('./data/TSPolInCM_filter_1.mat')
# data=io.loadmat('./data/TSPolInCM_filter_3.mat')
# data=io.loadmat('./data/TSPolInCM_filter_9.mat')
# data=io.loadmat('./data/TSPolInCM_filter_15.mat')
TSPolInCM_filter=data['data']
print(TSPolInCM_filter.shape)
# 数据格式为：6*6*m*n*N
# print(type(TSPolInCM_filter))
TSPolInCM_filter = cp.asarray(TSPolInCM_filter)
# print(type(TSPolInCM_filter))

def favecoh_fullpol(x,pTST11,pTST12,pTST22,nintf):
    # 计算三个ω值，x为求解的4个参数值
    _omegaconj=np.zeros((3,1),dtype='complex_')
    [_omegaconj[0],_omegaconj[1],_omegaconj[2]]=[np.conjugate(np.cos(x[0])),np.conjugate(np.sin(x[0])*np.cos(x[1])*np.exp(x[2]*1j)),np.conjugate(np.sin(x[0])*np.sin(x[1])*np.exp(x[3]*1j))]
    _omegaconj=_omegaconj.T
    
    [_npol,qq,_nintf]=pTST11.shape
    _omegaSet=np.kron(np.ones((1,nintf)),np.conjugate(_omegaconj))

    pTST11opt=np.multiply(np.dot(_omegaconj,np.reshape(pTST11,[_npol,_npol*_nintf],order='F')),_omegaSet)
    pTST12opt=np.multiply(np.dot(_omegaconj,np.reshape(pTST12,[_npol,_npol*_nintf],order='F')),_omegaSet)
    pTST22opt=np.multiply(np.dot(_omegaconj,np.reshape(pTST22,[_npol,_npol*_nintf],order='F')),_omegaSet)
    pTST11opt=np.sum(np.reshape(pTST11opt,[pTST11opt.shape[0],_nintf,_npol],order='C'),axis=2)
    pTST12opt=np.sum(np.reshape(pTST12opt,[pTST11opt.shape[0],_nintf,_npol],order='C'),axis=2)
    pTST22opt=np.sum(np.reshape(pTST22opt,[pTST11opt.shape[0],_nintf,_npol],order='C'),axis=2)
    TSCoh=np.divide(pTST12opt,np.sqrt(np.multiply(pTST11opt,pTST22opt)))
    
    meanTScoh=np.mean(np.abs(TSCoh),axis=1)
    return -meanTScoh

def ESPO_TSPolDS_1(TSPolInCM,interdeg=10):
    [q,qq,nrow,ncol,nintf]=TSPolInCM.shape
    # step1
    # 计算枚举4个参数的步长,并枚举4个参数（枚举需要包含右端点）
    npol=q//2
    num_1=90//interdeg
    num_2=360//interdeg
    if 90%interdeg==0:
        num_1+=1
        num_2+=1
    num_1=num_1*1j
    num_2=num_2*1j

    a=cp.mgrid[0:90:num_1,0:90:num_1,-180:180:num_2,-180:180:num_2]
    [aaa,bbb,ccc,ddd,eee]=a.shape
    sum=bbb*ccc*ddd*eee
    dff=cp.zeros((4,sum))
    b=cp.ravel(a[0],order='F')
    c=cp.ravel(a[1],order='F')
    d=cp.ravel(a[2],order='F')
    e=cp.ravel(a[3],order='F')
    dff[0]=b
    dff[1]=c
    dff[2]=d
    dff[3]=e
    para=dff*cp.pi/180

    # 使用4个参数计算ω值
    omegaconj=cp.zeros((3,sum),dtype='complex_')
    [omegaconj[0],omegaconj[1],omegaconj[2]]=[cp.conjugate(cp.cos(para[0,:])),cp.conjugate(cp.sin(para[0,:])*cp.cos(para[1,:])*cp.exp(para[2,:]*1j)),cp.conjugate(cp.sin(para[0,:])*cp.sin(para[1,:])*cp.exp(para[3,:]*1j))]
    omegaconj=omegaconj.T
    para=para.T
    # print(omegaconj.shape)
    # ω的数据格式：枚举数目×3
    omegaSet=cp.kron(cp.ones((1,nintf)),cp.conjugate(omegaconj))
    # print(omegaSet.shape)
    # ω的数据格式：枚举数目×（3*N）
    optpara=cp.zeros((nrow,ncol,4))

    TSPolInCMT11=TSPolInCM[0:npol,0:npol,:,:,:]
    TSPolInCMT12=TSPolInCM[0:npol,npol:npol*2,:,:,:]
    TSPolInCMT22=TSPolInCM[npol:npol*2,npol:npol*2,:,:,:]    
    # T矩阵数据格式：3×3×n×m×N
    # step2
    # 开始对每个像素进行处理
    for ii in range(nrow):
        for jj in range(ncol):

            pTST11=cp.squeeze(TSPolInCMT11[:,:,ii,jj,:])
            pTST12=cp.squeeze(TSPolInCMT12[:,:,ii,jj,:])
            pTST22=cp.squeeze(TSPolInCMT22[:,:,ii,jj,:])
            
            pTST11opt=cp.multiply(cp.dot(omegaconj,cp.reshape(pTST11,[npol,npol*nintf],order='F')),omegaSet)
            pTST12opt=cp.multiply(cp.dot(omegaconj,cp.reshape(pTST12,[npol,npol*nintf],order='F')),omegaSet)
            pTST22opt=cp.multiply(cp.dot(omegaconj,cp.reshape(pTST22,[npol,npol*nintf],order='F')),omegaSet)
            
            pTST11opt=cp.sum(cp.reshape(pTST11opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST12opt=cp.sum(cp.reshape(pTST12opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST22opt=cp.sum(cp.reshape(pTST22opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)

            TSCoh=cp.divide(pTST12opt,cp.sqrt(cp.multiply(pTST11opt,pTST22opt)))
            meanTScoh=cp.mean(cp.abs(TSCoh),axis=1)
            maxInd=cp.argmax(meanTScoh, axis=0)
            optpara[ii,jj,:]=para[maxInd,:]


    return optpara

def ESPO_TSPolDS_2(TSPolInCM,optparaDS0_lag3,interdeg):
    # optparaDS0_lag3是step2求解出的最优参数
    [q,qq,nrow,ncol,nintf]=TSPolInCM.shape
    npol=q//2
    optomegaconj=np.zeros((nrow,ncol,npol),dtype='complex_')
    optmeanCoh=np.zeros((nrow,ncol)) 
    optTSCoh=np.zeros((nrow,ncol,nintf),'complex_')
    TSPolInCMT11=TSPolInCM[0:npol,0:npol,:,:,:]
    TSPolInCMT12=TSPolInCM[0:npol,npol:npol*2,:,:,:]
    TSPolInCMT22=TSPolInCM[npol:npol*2,npol:npol*2,:,:,:]

    # step3
    # 非线性优化部分
    for ii in range(nrow):
        for jj in range(ncol):
            pTST11=np.squeeze(TSPolInCMT11[:,:,ii,jj,:])
            pTST12=np.squeeze(TSPolInCMT12[:,:,ii,jj,:])
            pTST22=np.squeeze(TSPolInCMT22[:,:,ii,jj,:])
            x0=optparaDS0_lag3[ii,jj,:].T
            lb0 = [0,0,-np.pi,-np.pi]
            ub0 = [np.pi/2,np.pi/2,np.pi,np.pi]
            lb1 = [x0[0],x0[1],x0[2],x0[3]]-np.deg2rad(interdeg)
            ub1 = [x0[0],x0[1],x0[2],x0[3]]+np.deg2rad(interdeg)
            lb=np.max([lb0,lb1],axis=0)
            ub=np.min([ub0,ub1],axis=0)
            bounds=((lb[0],ub[0]),(lb[1],ub[1]),(lb[2],ub[2]),(lb[3],ub[3]))
            fval=scipy.optimize.minimize(favecoh_fullpol,x0,args=(pTST11,pTST12,pTST22,nintf),bounds=bounds)
            x0=fval.x
            optmeanCoh[ii,jj]=-fval.fun
            fomegaconj=np.zeros((3,1),dtype='complex_')
            [fomegaconj[0],fomegaconj[1],fomegaconj[2]]=[np.conjugate(np.cos(x0[0])),np.conjugate(np.sin(x0[0])*np.cos(x0[1])*np.exp(x0[2]*1j)),np.conjugate(np.sin(x0[0])*np.sin(x0[1])*np.exp(x0[3]*1j))]
            fomegaconj=fomegaconj.T
            optomegaconj[ii,jj,:]=fomegaconj
            fomegaSet=np.kron(np.ones((1,nintf)),np.conjugate(fomegaconj))
            pTST11opt=np.multiply(np.dot(fomegaconj,np.reshape(pTST11,[npol,npol*nintf],order='F')),fomegaSet)
            pTST12opt=np.multiply(np.dot(fomegaconj,np.reshape(pTST12,[npol,npol*nintf],order='F')),fomegaSet)
            pTST22opt=np.multiply(np.dot(fomegaconj,np.reshape(pTST22,[npol,npol*nintf],order='F')),fomegaSet)
            pTST11opt=np.sum(np.reshape(pTST11opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST12opt=np.sum(np.reshape(pTST12opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST22opt=np.sum(np.reshape(pTST22opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            TSCoh=np.divide(pTST12opt,np.sqrt(np.multiply(pTST11opt,pTST22opt)))
            optTSCoh[ii,jj,:]=TSCoh

    s3=time.time()
    return optomegaconj,optTSCoh,optmeanCoh

s1=time.time()
optparaDS0=ESPO_TSPolDS_1(TSPolInCM_filter,30)
s2=time.time()

TSPolInCM_filter=cp.asnumpy(TSPolInCM_filter)
optparaDS0=cp.asnumpy(optparaDS0)


omegaconjDS,optTSCoh,optmeanCoh=ESPO_TSPolDS_2(TSPolInCM_filter,optparaDS0,30)
s3=time.time()
print('step1&step2 time:')
print(s2-s1,'s')
print('step3 time:')
print(s3-s2,'s')


def TSPolInCM2TSInCMSet(TSPolInCM,omegaconj):
    [q,qq,nrow,ncol,nintf]=TSPolInCM.shape
    optTSInCMSet=np.zeros((2,2,nrow,ncol,nintf),dtype='complex_')
    for ii in range(nrow):
        for jj in range(nrow):
            for kk in range(nintf):
                pomegaconj=np.kron(np.identity(2),(np.squeeze(omegaconj[ii,jj,:])).T)
                optTSInCMSet[:,:,ii,jj,kk]=np.dot(np.dot(pomegaconj,TSPolInCM[:,:,ii,jj,kk]),np.conjugate(pomegaconj.T))
    return optTSInCMSet
    
def Normalize_TSInCM(TSInCM):
    ndim=TSInCM.ndim
    if ndim>=4:
        [q,qq,nrow,ncol,npol]=TSInCM.shape
        norTSInCM=np.zeros((q,q,nrow,ncol,npol),dtype='complex_')
        for ii in range(npol):
            for kk in range(nrow):
                for ll in range(ncol):
                    TSIntensity=np.diag(TSInCM[:,:,kk,ll,ii])
                    TSIntensity=TSIntensity.reshape(q,1)
                    TSIntCM1=np.dot(TSIntensity,np.ones((1,q)))
                    TSIntCM2=np.ones((q,1))*TSIntensity.T
                    norTSInCM[:,:,kk,ll,ii]=TSInCM[:,:,kk,ll,ii]/np.sqrt(TSIntCM1*TSIntCM2)
    return norTSInCM


TSInCMSet_ESPODS=Normalize_TSInCM(TSPolInCM2TSInCMSet(TSPolInCM_filter,omegaconjDS))#ESPO


optCor=np.squeeze((TSInCMSet_ESPODS[0,1,:,:,]))
=======
import scipy.io as io
import numpy as np
import cupy as cp
import scipy
import time

data=io.loadmat('./data/TSPolInCM_filter_1.mat')
# data=io.loadmat('./data/TSPolInCM_filter_3.mat')
# data=io.loadmat('./data/TSPolInCM_filter_9.mat')
# data=io.loadmat('./data/TSPolInCM_filter_15.mat')
TSPolInCM_filter=data['data']
print(TSPolInCM_filter.shape)
# 数据格式为：6*6*m*n*N
# print(type(TSPolInCM_filter))
TSPolInCM_filter = cp.asarray(TSPolInCM_filter)
# print(type(TSPolInCM_filter))

def favecoh_fullpol(x,pTST11,pTST12,pTST22,nintf):
    # 计算三个ω值，x为求解的4个参数值
    _omegaconj=np.zeros((3,1),dtype='complex_')
    [_omegaconj[0],_omegaconj[1],_omegaconj[2]]=[np.conjugate(np.cos(x[0])),np.conjugate(np.sin(x[0])*np.cos(x[1])*np.exp(x[2]*1j)),np.conjugate(np.sin(x[0])*np.sin(x[1])*np.exp(x[3]*1j))]
    _omegaconj=_omegaconj.T
    
    [_npol,qq,_nintf]=pTST11.shape
    _omegaSet=np.kron(np.ones((1,nintf)),np.conjugate(_omegaconj))

    pTST11opt=np.multiply(np.dot(_omegaconj,np.reshape(pTST11,[_npol,_npol*_nintf],order='F')),_omegaSet)
    pTST12opt=np.multiply(np.dot(_omegaconj,np.reshape(pTST12,[_npol,_npol*_nintf],order='F')),_omegaSet)
    pTST22opt=np.multiply(np.dot(_omegaconj,np.reshape(pTST22,[_npol,_npol*_nintf],order='F')),_omegaSet)
    pTST11opt=np.sum(np.reshape(pTST11opt,[pTST11opt.shape[0],_nintf,_npol],order='C'),axis=2)
    pTST12opt=np.sum(np.reshape(pTST12opt,[pTST11opt.shape[0],_nintf,_npol],order='C'),axis=2)
    pTST22opt=np.sum(np.reshape(pTST22opt,[pTST11opt.shape[0],_nintf,_npol],order='C'),axis=2)
    TSCoh=np.divide(pTST12opt,np.sqrt(np.multiply(pTST11opt,pTST22opt)))
    
    meanTScoh=np.mean(np.abs(TSCoh),axis=1)
    return -meanTScoh

def ESPO_TSPolDS_1(TSPolInCM,interdeg=10):
    [q,qq,nrow,ncol,nintf]=TSPolInCM.shape
    # step1
    # 计算枚举4个参数的步长,并枚举4个参数（枚举需要包含右端点）
    npol=q//2
    num_1=90//interdeg
    num_2=360//interdeg
    if 90%interdeg==0:
        num_1+=1
        num_2+=1
    num_1=num_1*1j
    num_2=num_2*1j

    a=cp.mgrid[0:90:num_1,0:90:num_1,-180:180:num_2,-180:180:num_2]
    [aaa,bbb,ccc,ddd,eee]=a.shape
    sum=bbb*ccc*ddd*eee
    dff=cp.zeros((4,sum))
    b=cp.ravel(a[0],order='F')
    c=cp.ravel(a[1],order='F')
    d=cp.ravel(a[2],order='F')
    e=cp.ravel(a[3],order='F')
    dff[0]=b
    dff[1]=c
    dff[2]=d
    dff[3]=e
    para=dff*cp.pi/180

    # 使用4个参数计算ω值
    omegaconj=cp.zeros((3,sum),dtype='complex_')
    [omegaconj[0],omegaconj[1],omegaconj[2]]=[cp.conjugate(cp.cos(para[0,:])),cp.conjugate(cp.sin(para[0,:])*cp.cos(para[1,:])*cp.exp(para[2,:]*1j)),cp.conjugate(cp.sin(para[0,:])*cp.sin(para[1,:])*cp.exp(para[3,:]*1j))]
    omegaconj=omegaconj.T
    para=para.T
    # print(omegaconj.shape)
    # ω的数据格式：枚举数目×3
    omegaSet=cp.kron(cp.ones((1,nintf)),cp.conjugate(omegaconj))
    # print(omegaSet.shape)
    # ω的数据格式：枚举数目×（3*N）
    optpara=cp.zeros((nrow,ncol,4))

    TSPolInCMT11=TSPolInCM[0:npol,0:npol,:,:,:]
    TSPolInCMT12=TSPolInCM[0:npol,npol:npol*2,:,:,:]
    TSPolInCMT22=TSPolInCM[npol:npol*2,npol:npol*2,:,:,:]    
    # T矩阵数据格式：3×3×n×m×N
    # step2
    # 开始对每个像素进行处理
    for ii in range(nrow):
        for jj in range(ncol):

            pTST11=cp.squeeze(TSPolInCMT11[:,:,ii,jj,:])
            pTST12=cp.squeeze(TSPolInCMT12[:,:,ii,jj,:])
            pTST22=cp.squeeze(TSPolInCMT22[:,:,ii,jj,:])
            
            pTST11opt=cp.multiply(cp.dot(omegaconj,cp.reshape(pTST11,[npol,npol*nintf],order='F')),omegaSet)
            pTST12opt=cp.multiply(cp.dot(omegaconj,cp.reshape(pTST12,[npol,npol*nintf],order='F')),omegaSet)
            pTST22opt=cp.multiply(cp.dot(omegaconj,cp.reshape(pTST22,[npol,npol*nintf],order='F')),omegaSet)
            
            pTST11opt=cp.sum(cp.reshape(pTST11opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST12opt=cp.sum(cp.reshape(pTST12opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST22opt=cp.sum(cp.reshape(pTST22opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)

            TSCoh=cp.divide(pTST12opt,cp.sqrt(cp.multiply(pTST11opt,pTST22opt)))
            meanTScoh=cp.mean(cp.abs(TSCoh),axis=1)
            maxInd=cp.argmax(meanTScoh, axis=0)
            optpara[ii,jj,:]=para[maxInd,:]


    return optpara

def ESPO_TSPolDS_2(TSPolInCM,optparaDS0_lag3,interdeg):
    # optparaDS0_lag3是step2求解出的最优参数
    [q,qq,nrow,ncol,nintf]=TSPolInCM.shape
    npol=q//2
    optomegaconj=np.zeros((nrow,ncol,npol),dtype='complex_')
    optmeanCoh=np.zeros((nrow,ncol)) 
    optTSCoh=np.zeros((nrow,ncol,nintf),'complex_')
    TSPolInCMT11=TSPolInCM[0:npol,0:npol,:,:,:]
    TSPolInCMT12=TSPolInCM[0:npol,npol:npol*2,:,:,:]
    TSPolInCMT22=TSPolInCM[npol:npol*2,npol:npol*2,:,:,:]

    # step3
    # 非线性优化部分
    for ii in range(nrow):
        for jj in range(ncol):
            pTST11=np.squeeze(TSPolInCMT11[:,:,ii,jj,:])
            pTST12=np.squeeze(TSPolInCMT12[:,:,ii,jj,:])
            pTST22=np.squeeze(TSPolInCMT22[:,:,ii,jj,:])
            x0=optparaDS0_lag3[ii,jj,:].T
            lb0 = [0,0,-np.pi,-np.pi]
            ub0 = [np.pi/2,np.pi/2,np.pi,np.pi]
            lb1 = [x0[0],x0[1],x0[2],x0[3]]-np.deg2rad(interdeg)
            ub1 = [x0[0],x0[1],x0[2],x0[3]]+np.deg2rad(interdeg)
            lb=np.max([lb0,lb1],axis=0)
            ub=np.min([ub0,ub1],axis=0)
            bounds=((lb[0],ub[0]),(lb[1],ub[1]),(lb[2],ub[2]),(lb[3],ub[3]))
            fval=scipy.optimize.minimize(favecoh_fullpol,x0,args=(pTST11,pTST12,pTST22,nintf),bounds=bounds)
            x0=fval.x
            optmeanCoh[ii,jj]=-fval.fun
            fomegaconj=np.zeros((3,1),dtype='complex_')
            [fomegaconj[0],fomegaconj[1],fomegaconj[2]]=[np.conjugate(np.cos(x0[0])),np.conjugate(np.sin(x0[0])*np.cos(x0[1])*np.exp(x0[2]*1j)),np.conjugate(np.sin(x0[0])*np.sin(x0[1])*np.exp(x0[3]*1j))]
            fomegaconj=fomegaconj.T
            optomegaconj[ii,jj,:]=fomegaconj
            fomegaSet=np.kron(np.ones((1,nintf)),np.conjugate(fomegaconj))
            pTST11opt=np.multiply(np.dot(fomegaconj,np.reshape(pTST11,[npol,npol*nintf],order='F')),fomegaSet)
            pTST12opt=np.multiply(np.dot(fomegaconj,np.reshape(pTST12,[npol,npol*nintf],order='F')),fomegaSet)
            pTST22opt=np.multiply(np.dot(fomegaconj,np.reshape(pTST22,[npol,npol*nintf],order='F')),fomegaSet)
            pTST11opt=np.sum(np.reshape(pTST11opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST12opt=np.sum(np.reshape(pTST12opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            pTST22opt=np.sum(np.reshape(pTST22opt,[pTST11opt.shape[0],nintf,npol],order='C'),axis=2)
            TSCoh=np.divide(pTST12opt,np.sqrt(np.multiply(pTST11opt,pTST22opt)))
            optTSCoh[ii,jj,:]=TSCoh

    s3=time.time()
    return optomegaconj,optTSCoh,optmeanCoh

s1=time.time()
optparaDS0=ESPO_TSPolDS_1(TSPolInCM_filter,30)
s2=time.time()

TSPolInCM_filter=cp.asnumpy(TSPolInCM_filter)
optparaDS0=cp.asnumpy(optparaDS0)


omegaconjDS,optTSCoh,optmeanCoh=ESPO_TSPolDS_2(TSPolInCM_filter,optparaDS0,30)
s3=time.time()
print('step1&step2 time:')
print(s2-s1,'s')
print('step3 time:')
print(s3-s2,'s')


def TSPolInCM2TSInCMSet(TSPolInCM,omegaconj):
    [q,qq,nrow,ncol,nintf]=TSPolInCM.shape
    optTSInCMSet=np.zeros((2,2,nrow,ncol,nintf),dtype='complex_')
    for ii in range(nrow):
        for jj in range(nrow):
            for kk in range(nintf):
                pomegaconj=np.kron(np.identity(2),(np.squeeze(omegaconj[ii,jj,:])).T)
                optTSInCMSet[:,:,ii,jj,kk]=np.dot(np.dot(pomegaconj,TSPolInCM[:,:,ii,jj,kk]),np.conjugate(pomegaconj.T))
    return optTSInCMSet
    
def Normalize_TSInCM(TSInCM):
    ndim=TSInCM.ndim
    if ndim>=4:
        [q,qq,nrow,ncol,npol]=TSInCM.shape
        norTSInCM=np.zeros((q,q,nrow,ncol,npol),dtype='complex_')
        for ii in range(npol):
            for kk in range(nrow):
                for ll in range(ncol):
                    TSIntensity=np.diag(TSInCM[:,:,kk,ll,ii])
                    TSIntensity=TSIntensity.reshape(q,1)
                    TSIntCM1=np.dot(TSIntensity,np.ones((1,q)))
                    TSIntCM2=np.ones((q,1))*TSIntensity.T
                    norTSInCM[:,:,kk,ll,ii]=TSInCM[:,:,kk,ll,ii]/np.sqrt(TSIntCM1*TSIntCM2)
    return norTSInCM


TSInCMSet_ESPODS=Normalize_TSInCM(TSPolInCM2TSInCMSet(TSPolInCM_filter,omegaconjDS))#ESPO


optCor=np.squeeze((TSInCMSet_ESPODS[0,1,:,:,]))
>>>>>>> 75e45bb (first commit)
>>>>>>> 1b4de05 (first commit)
print(np.mean(np.abs(optCor)))