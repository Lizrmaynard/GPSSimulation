import math
import numpy as np
import sys

def ReadData():
    bounds = CalculateBounds()
    retval = np.empty(shape=[bounds,4])
    global aI

    for j in range(0,bounds):
        retval[j][0] = vals[aI+1]
        retval[j][1] = vals[aI+2]
        retval[j][2] = vals[aI+3]
        retval[j][3] = vals[aI+4]
        aI = aI + 5
    global M
    M = bounds
    return retval


def CalculateBounds():
    m = 1

    while aI+1+5*m < len(vals) and abs(float(vals[aI+1+5*(m)])-float(vals[aI+1+5*(m-1)])) <0.5:
        m = m + 1
    return m

def ComputeLocation():
    XV = NewtonsMethod(X0,0)
    TV = TimeAt(XV)
    val = [TV,XV[0],XV[1],XV[2]]
    retval = PositionInGeographic(val)
    ret = ""
    for i in range(0,10):
        ret = ret + str(retval[i]) + " "
    print(ret)

def NewtonsMethod(x,d):

    global Sk

    jacx = Jacobian(x)
    gradx = gradf(x)

    Sk = Solve3x3(Jacobian(x), gradf(x))
    retval = x[:]

    for i in range(0,3):
        retval[i] -= Sk[i]

    if abs(diff(retval,x)[0])<0.0000001 and abs(diff(retval,x)[1])<0.0000001 and abs(diff(retval,x)[2])<0.0000001:
        return retval
    elif d>9:
        return None
    else:
        d=d+1
        return NewtonsMethod(retval,d)

def Solve3x3(A,b):
    x = (A[0][1]*A[1][2]*b[2] - A[0][1]*A[2][2]*b[1] - A[0][2]*A[1][1]*b[2] + A[0][2]*A[2][1]*b[1] + A[1][1]*A[2][2]*b[0] - A[1][2]*A[2][1]*b[0])/(A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0])
    y = -1 * (A[0][0]*A[1][2]*b[2] - A[0][0]*A[2][2]*b[1] - A[0][2]*A[1][0]*b[2] + A[0][2]*A[2][0]*b[1] + A[1][0]*A[2][2]*b[0] - A[1][2]*A[2][0]*b[0])/(A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0])
    z = (A[0][0]*A[1][1]*b[2] - A[0][0]*A[2][1]*b[1] - A[0][1]*A[1][0]*b[2] + A[0][1]*A[2][0]*b[1] + A[1][0]*A[2][1]*b[0] - A[1][1]*A[2][0]*b[0])/(A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1] - A[0][2]*A[1][1]*A[2][0])
    return [x,y,z]


def Jacobian(x):

    retval = np.empty(shape=[3,3])
    for i in range(0,3):
        for j in range(0,3):
            retval[i][j] = J_ij(i,j,x)

    return retval

def diff(a,b):

    if len(a) == len(b):
        retval = np.empty(shape=[len(a)])
        for i in range(0,len(a)):
            retval[i] = a[i]-b[i]
        return retval

    return None

def twonorm(a):
    sum = 0.0
    for d in a:
        sum = sum + d*d
    return math.sqrt(sum)

def gradf(x):

    diffs = np.empty(shape=[M,3])
    for j in range(0,M):
        diffs[j] = diff([data[j][1],data[j][2],data[j][3]],x)

    N = np.empty(shape=[M])
    for j in range(0,M):
        N[j] = twonorm(diffs[j])

    A = np.empty(shape=[M-1])
    for j in range(0,M-1):
        A[j] = N[j+1]-N[j]-c*(data[j][0]-data[j+1][0])

    XYZ = np.empty(shape=[3,M-1])
    for i in range(0,3):
        for j in range(0,M-1):
            XYZ[i][j] = diffs[j][i]/N[j] - diffs[j+1][i]/N[j+1]

    retval = np.empty(shape=[3])
    for j in range(0,3):
        retval[j] = 0.0
        for i in range(0,M-1):
            retval[j] += A[i]*XYZ[j][i]
        retval[j] = 2.0*retval[j]

    return retval

def gradF_ij(i,j,x):

    retval = 0.5 * (2.0 * x[j] - 2.0 * data[i][j+1]) / (math.sqrt(pow((x[0]-data[i][1]),2) + pow((x[1]-data[i][2]),2) + pow((x[2]-data[i][3]),2)));
    retval -= 0.5 * (2.0 * x[j] - 2.0 * data[i+1][j+1]) / (math.sqrt(pow((x[0]-data[i+1][1]),2) + pow((x[1]-data[i+1][2]),2) + pow((x[2]-data[i+1][3]),2)));
    return retval

def J_ij(i,j,x):
    sum = 0.0
    for k in range(0,M-1):
        sum += gradF_ij(k,i,x) * gradF_ij(k,j,x)

    return 2.0 * sum

def R3(alpha,x):
    retval = [math.cos(alpha)*x[0] - math.sin(alpha)*x[1],math.sin(alpha)*x[0] + math.cos(alpha)*x[1],x[2]]
    return retval

def RadiansToDegrees(a):
    b = a*180/math.pi
    if a<0:
        b = -1*b
    d = math.floor(b)
    m = math.floor(60*(b-d))
    s = 60*(60*(b-d)-m)
    if a<0:
        return [str(d),str(m),str(s), str(-1)]
    return [str(d),str(m),str(s), str(1)]

def PositionInGeographic(a):

    t = a[0]
    x1 = a[1]
    y1 = a[2]
    z1 = a[3]

    xyz = R3(-2*math.pi*t/S, [x1,y1,z1])
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]

    global psi
    if x*x+y*y == 0:
        if z>=0:
            psi = math.pi/2
        else:
            psi = -1*math.pi/2
    else:
        psi = math.atan2(z,math.sqrt(x*x+y*y))

    global lambd

    if x>0 and y>0:
        lambd = math.atan2(y,x)
    elif x < 0:
        lambd = math.pi + math.atan2(y,x)
    else:
        lambd = 2*math.pi + math.atan2(y,x)

    lambd-=math.pi;
    Psi = RadiansToDegrees(psi);
    Lambd = RadiansToDegrees(lambd);
    h = math.sqrt(x*x + y*y + z*z) - R
    return [t, Psi[0], Psi[1], Psi[2], Psi[3], Lambd[0], Lambd[1], Lambd[2], Lambd[3], h]

def TimeAt(Xv):
    retval = data[0][0] + (1/c) * math.sqrt(pow(Xv[0]-data[0][1],2)+pow(Xv[1]-data[0][2],2)+pow(Xv[2]-data[0][3],2))
    return retval


global aI
aI = 0
c = float("2.997924580000000000E+08")
R = float("6.367444500000000000E+06")
S = float("8.616408999999999651E+04")
X0 = [0,0,0]
global vals
global data
vals = [0]

while 1:
    data = (sys.stdin.readline()).strip()
    if len(data) < 1:
        break
    tmp = data.split(" ")
    for x in tmp:
        vals.append(x)
vals.pop(1)

while aI<len(vals):
    data = ReadData()
    #print(data)
    ComputeLocation()