from scipy.fftpack import fft
import matplotlib.pyplot as plt
import numpy as np
import math


plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False   #用来正常显示负号


'''
将txt文件导入
'''
dataair=np.loadtxt("G:\medcine\\20191022\\air.txt")
x_air=(dataair[:,0])
y_air=list(dataair[:,1])

dataref=np.loadtxt("G:\medcine\\20191022\\G-1-0.txt") #双\\可以添加绝对路径
x_sam1=list(dataref[:,0])
y_sam1=list(dataref[:,1])


y_sam=[]
for i in range(len(x_air)):
    y_sam.append(y_sam1[i])


'''
对数据进行FFT变换，并取模
'''
yy_air=fft(y_air)
yy_sam=fft(y_sam)
yf_sam=abs(fft(y_sam))
yf_air=abs(fft(y_air))

'''
将时间轴转化为频率轴
'''
m=int((len(x_air)-1)/2)
d=1750  #sample厚度um
f=1/(3.3*10**(-14))/10**12
F=f/(len(x_air)-1)
xf=np.arange(0,f+F,F)
shuzi=np.arange(1,len(x_air)+1,1)
W=shuzi*F*2*math.pi           #求解透射率需要的角频率
x_air=[i*10**12 for i in x_air]
'''
做图
'''
graph_one=plt.figure(1)

plt.subplot(211)
plt.plot(x_air, y_sam,color='green',label='sam Original wave')
plt.plot(x_air,y_air,color='red',label='air Original wave')
plt.grid(linestyle='-.')
plt.xlabel('Time (ps)')
plt.ylabel('Amplitude (a.u)')
plt.legend()


plt.subplot(212)
plt.plot(xf,yf_sam,color='green',label='sam Original wave')
plt.plot(xf,yf_air,color='red',label='air Original wave')
# plt.title('Spectral amplitude in frequency domain for 100 ms')
plt.grid(linestyle='-.')
plt.xlabel('Frequency (THz)')
plt.ylabel('Amplitude (a.u)')
plt.xlim((0, 3.8))
x=np.arange(0,3.8,0.2)
plt.xticks(x)
plt.legend()

graph_two=plt.figure(2)
plt.subplot(211)
pw=yf_sam/yf_air  #复投射系数幅度
xf1=list(xf)
pw1=list(pw)
plt.plot(xf1[0:m],pw1[0:m])
plt.grid(linestyle='-.')
plt.xlabel('Frequence(THz)')
plt.ylabel('Transmission (%)')
plt.title('Transmission coefficient (%)')


plt.subplot(212)
transmission=yy_sam/yy_air
angle_y=np.angle(transmission)  #计算模角
transmission_phase=-(np.unwrap(angle_y))
plt.plot(xf,transmission_phase)
plt.grid(linestyle='-.')
plt.xlabel('Frequency (THz)')
plt.ylabel('Phase (rad)')
plt.title('Phase in frequency domain')



graph_three=plt.figure(3)
plt.subplot(211)  #透射率谱
n=(transmission_phase*3*10**8/(d*10**(-6)*W*10**12))+1
plt.plot(xf,n)
plt.grid(linestyle='-.')
plt.xlabel('Frequency (THz)')
plt.title('refractive index (n)')
plt.xlim((0, 3))
x=np.arange(0,3,0.2)
plt.xticks(x)


plt.subplot(212)
n=np.array(n)
b=np.array(4*n/(pw*(n+1)**2))
a=np.array([math.log(x) for x in b])
A=a*2/(d*10**(-4))
plt.plot(xf,A,color='red',label='A')
# plt.plot(xf,ca,color='green',label='ca')   #tianjia
plt.grid(linestyle='-.')
# plt.legend()
plt.xlabel('Frequency (THz)')
plt.title('Absorption coefficient (cm-1)')
plt.xlim((0, 3))
x=np.arange(0,3,0.2)
plt.xticks(x)

plt.show()




