import math
#from decimal import Decimal
import numpy as np

def processLine(vehLine):
    vehVals = vehLine.split(' ')
    vehXyz = Xv(vehVals)
    Tv = float(vehVals[0])
    B = HorizonCheck(vehXyz, Tv)
    for i in range(0,len(B)):
       if B[i]:

           term1 = (Xs(i, Tv)[0] - vehXyz[0])**2
           term2 = (Xs(i, Tv)[1] - vehXyz[1])**2
           term3 = (Xs(i, Tv)[2] - vehXyz[2])**2
           insideNorm = (term1 + term2 + term3)
           t0 = Tv - np.sqrt(insideNorm)/c

           timeofSat = Ts(i, Tv, t0, 0, vehXyz)
           posOfSat = Xs(i,timeofSat)
           ret = str(i)+" "+str(timeofSat)+" "+str(posOfSat[0])+" "+str(posOfSat[1])+" "+str(posOfSat[2])
           print(ret)
    return


# Xv position in cartesian
def Xv(vals):
    t = int(float(vals[0]))
    longarray = np.array([vals[1],vals[2],vals[3]])
    theta = int(vals[4])*(DegreesToRadians(longarray))
    latarray = np.array([vals[5],vals[6],vals[7]])
    phi = int(vals[8])*(DegreesToRadians(latarray))
    h = int(float(vals[9]))
    x = (R+h) * (math.cos(theta))*(math.cos(phi))
    y = (R+h) * (math.cos(theta))*(math.sin(phi))
    z = (R+h) * (math.sin(theta))
    alpha = (2*pi*t)/s

    return np.array([math.cos(alpha)*x-math.sin(alpha)*y,math.sin(alpha)*x+math.cos(alpha)*y,z])

#Checks to see if the satellite can be seen from Xv
def HorizonCheck(vehXyz, Tv):
    retval = np.zeros(24)
    for i in range (0,24):
        satXyz = Xs(i,Tv)
        if (np.dot(satXyz - vehXyz, vehXyz, out=None))>=0:
            retval[i] = True

    return retval


def f(satNum, t, Xv):

    return np.linalg.norm(np.subtract(Xs(satNum,t), Xv)) - c*(Tv-t)


#Completing fixed point iteration to solve for Ts.
#Once we have Ts we can plug it in to Xs(t) to solve for our position
def Ts(satNum, Tv, tk, depth, Xv):
   posofSatUnknownT = Xs(satNum, tk)
   TkPlusOne = Tv - (np.sqrt(((posofSatUnknownT[0] - Xv[0])**2) + ((posofSatUnknownT[1] - Xv[1])**2) + ((posofSatUnknownT[2] - Xv[2])**2) ))/c
   if abs(TkPlusOne-tk) < 0.01/c:
        return TkPlusOne
   elif depth>=9:
        return -1
   else:
       return Ts(satNum,Tv, TkPlusOne,depth+1, Xv)



#Given degrees, minutes, and seconds, we can find the angle in radians
def DegreesToRadians(st):
    return 2*float(pi)*((int(st[0])/360) + (int(st[1])/(360*60)) + (float(st[2])/(360*60*60)))

#Orbitals of satellites at some time t (computed from parent function in XsSingle)
def Xs(satNum,t):
    x = XsSingle(satNum, t, 0, 3)
    y = XsSingle(satNum, t, 1, 4)
    z = XsSingle(satNum, t, 2, 5)
    return np.array([x,y,z])


# This is computing the cartesian coordinates of the satellite in orbit (Values sent in to Xs to get final coordinates)
def XsSingle(satNum, t, u, v):
   rPlusH = R+V[satNum][7]
   uCos = V[satNum][u]*math.cos((2*pi*t)/V[satNum][6]+V[satNum][8])
   vSin = V[satNum][v]*math.sin((2*pi*t)/V[satNum][6]+V[satNum][8])

   return rPlusH*(uCos+vSin)


file = open("data.dat", "rt")
r = (file.readline())

pi = float((r[0:27]).strip())
r = (file.readline())
c = float((r[0:27]).strip())
r = (file.readline())
R = float((r[0:27]).strip())
r = (file.readline())
s = float((r[0:27]).strip())
Tv = 0

#Populate V with the data points (9) of the 24 satellites
V = np.empty((24,9),float)
for i in range(0,24):
    for j in range(0,9):
        d = file.readline()[0:27].strip()
        V[i,j] = float(d)

file.close()

defveh = []
flag = True

while flag:
    try:
        vehicle = input()
        vehicle.rstrip()
        defveh.append(vehicle)
    except EOFError:
        flag = False
    processLine(vehicle)
