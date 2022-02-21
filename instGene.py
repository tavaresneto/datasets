import math
import numpy
import os


##This function format a value into a GAMS style. Used by the formatGAMSTable
##@param value the value to be formatted
##@param size the lenght of the returning string
def formatGAMSCol(value, size=9):
    sp = " "*size
    sptrim = sp[len(str(value)):]
    return (str(value)+sptrim)

##This function format a vector given in "data" to a GAMS-formatted 
##parameter named as defined in "label"
def formatGAMSParameter(data, label):
    ret = label + """
/
0\t0
"""
    #for k,v in data.iteritems():
    #    ret = ret + str(k) + "\t" + str(v) +"""
    for i in range(data.shape[0]):
        ret = ret + str(i+1) + "\t" + str(data[i]) + """
"""
    ret = ret + "/"
    return(ret)

##This function format a dictionary given in "data" to a GAMS-formatted 
##table named as defined in "label"
def formatGAMS2DTable(data, label):
    ret = label + """
""" + formatGAMSCol("")
    for i in range(data.shape[0]):
        ret = ret + formatGAMSCol(i)
    
    for i in range(data.shape[0]):
        ret = ret + "\n" + formatGAMSCol(i)
        for j in range(data.shape[1]):
            ret = ret + formatGAMSCol(data[i][j])
    ret = ret + "\n;\n\n"
    return ret

def getRhoi(n):
    ret = numpy.random.randint(low=1,high=100,size=n)
    return ret

def getGamma(n):
    ret = numpy.random.randint(low=1,high=100,size=n)
    return ret

def getEpisilon(n):
    ret = numpy.random.randint(low=1,high=10,size=n)
    return ret

def getSigmai(n):
    ret = numpy.random.randint(low=1,high=10,size=n)
    return ret

def getKappaij(n, theta, rhoi):
    size = theta * (n+1)
    pos = numpy.zeros((n+1,2), dtype=numpy.int32)
    pos[0][0] = numpy.random.randint(low=0, high=size)
    pos[0][1] = numpy.random.randint(low=0, high=size)
    for i in range(1,n+1):
        npos = [numpy.random.randint(low=0, high=size), numpy.random.randint(low=0, high=size)]
        while npos in pos:
            npos = [numpy.random.randint(low=0, high=size), numpy.random.randint(low=0, high=size)]
        pos[i][0] = npos[0]
        pos[i][1] = npos[1]
    
    ret = numpy.zeros((n+1,n+1), dtype=numpy.int32)
    for i in range(n+1):
        for j in range(n+1):                
            ret[i][j] = long(math.sqrt((pos[i][0] - pos[j][0])**2 + (pos[i][1] - pos[j][1])**2))

    return(ret)

"""
    ret = numpy.zeros((n+1,n+1), dtype=numpy.int8)
    for i in range(n+1):
        for j in range(n+1):
            if i == j:
                ret[i][j] = 0
            else:
                tmp = (0.1 - theta) * numpy.random.random() + theta
                r1, r2 = 0,0
                if i > 0:
                    r1 = rhoi[i-1]
                if j > 0:
                    r2 = rhoi[j-1]
                ret[i][j] = int(tmp * (r1 + r2))
    return ret
"""
def getDelta(n, theta, g):
    size = theta * (n+1)
    pos = numpy.zeros((n+1,2), dtype=numpy.int32)
    pos[0][0] = numpy.random.randint(low=0, high=size)
    pos[0][1] = numpy.random.randint(low=0, high=size)
    for i in range(1,n+1):
        npos = [numpy.random.randint(low=0, high=size), numpy.random.randint(low=0, high=size)]
        while npos in pos:
            npos = [numpy.random.randint(low=0, high=size), numpy.random.randint(low=0, high=size)]
        pos[i][0] = npos[0]
        pos[i][1] = npos[1]
    
    ret = numpy.zeros((n+1,n+1), dtype=numpy.int32)
    for i in range(n+1):
        for j in range(n+1):
            ret[i][j] = long(math.sqrt((pos[i][0] - pos[j][0])**2 + (pos[i][1] - pos[j][1])**2))
            if i==0 or j==0:
                ret[i][j] = ret[i][j] * g
    return(ret)
    
def getPsi(sigma, thetaSigma):
    s = max(sigma)
    return numpy.random.randint(low=s, high=s*thetaSigma)

dir0 = os.getcwd()
j = 0
for n in [5, 10, 20, 40, 80]:
    for thetaS in [10, 50, 100]:
        for thetaD in [5, 10, 30]:
            for thetaG in [1, 10, 30]:
                for thetaSigma in [5, 10, 20]:
                    dirname = dir0+"/"+str(n) + "_" + str(thetaS) + "_" + str(thetaD) + "_" + str(thetaG) + "_" + str(thetaSigma)
                    if not os.path.isdir(dirname):
                        os.makedirs(dirname)
                        os.chdir(dirname)
                        print dirname

                        i = 0
                        while i < 30:
                            rho = getRhoi(n)
                            gamma = getGamma(n)
                            sigma = getSigmai(n)
                            episilon = getEpisilon(n)
                            kappa = getKappaij(n, thetaS, rho)
                            delta = getDelta(n, thetaD, thetaG)
                            psi = getPsi(sigma, thetaSigma)

                            ret = """

sets
i        orders          /0*"""
                            ret = ret + str(n)+ """/
k        routes          /1*"""
                            ret = ret + str(n)+ """/
p        positions       /0*"""
                            ret = ret + str(n)+ """/
;

alias(i,j)
alias(i,h)
alias(i,hh)
alias(i,jj)
alias(k,kk)
                            """
                            ret = ret + "Parameters\n\n" + formatGAMSParameter(rho, "rho(i)\tProcessing time")
                            ret = ret + "\n\n" + formatGAMSParameter(gamma, "gamma(i)\tDuedate")
                            ret = ret + "\n\n" + formatGAMSParameter(sigma, "sigma(i)\tSize")
                            ret = ret + "\n\n" + formatGAMSParameter(episilon, "epsilon(i)\tPriority")
                            ret = ret + "\n\n" + "psi\tVehicle Capacity\n/ "+str(psi)+" /\n;"
                            ret = ret + "\n\n" + formatGAMS2DTable(kappa, "Table kappa(i,j)                 Setup")
                            ret = ret + "\n\n" + formatGAMS2DTable(delta, "Table delta(i,j)                 Distance")
                            
                            filename = "rand_"+str(n) + "_" + str(thetaS) + "_" + str(thetaD)  + "_" + str(thetaG) + "_" + str(thetaSigma)+ "_" + str(i) + ".txt"
                            f = open(filename, "w")
                            f.write(ret)
                            f.flush()
                            f.close()
                            
                            #print filename+ "\tOK"
                            i = i + 1
                            j = j + 1
                    os.chdir(dir0)
                    print dirname
print str(j) + " files generated"
