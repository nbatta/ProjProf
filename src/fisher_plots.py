import numpy as np
import matplotlib.pyplot as plt
import itertools

def marginal_pars(pars,fishmats,fishlab=None):
    #print out marginalized constraints
    diag = np.diagonal(fishmats)
    if (fishlab):
        for i in xrange(len(diag)):
            print fishlab[i], pars[i],np.sqrt(diag[i])
    else:
        for i in xrange(len(diag)):
            print pars[i],np.sqrt(diag[i])
    return
    
def make_circle(numpnts=360):
    xx = np.array(np.arange(numpnts) / 180. * np.pi)
    circl = np.array([np.cos(xx),np.sin(xx)])
    return circl

def make_elsp(chi2,chifac = 2.3):
    
    circ = make_circle()
    Lmat = np.linalg.cholesky(chi2)

    ans1sig = np.dot(np.sqrt(chifac)*Lmat,circ)
    ans2sig = np.dot(np.sqrt(2.0*chifac)*Lmat,circ)
    return ans1sig, ans2sig

def tri_plot(pars,fishmats,setNames,cols=itertools.repeat(None),lss=itertools.repeat(None),parlabs=itertools.repeat(None),saveFile="default.png",thk=3):
    #
    numpars = len(pars)
    print numpars
    plt.figure(figsize=(4*numpars,4*numpars))
    plt.rc('axes', linewidth=thk)
    
    print np.shape(fishmats)
    bigcount = 0
    for setName,col,ls,lab in zip(setNames,cols,lss,parlabs):
        if (len(setName) > 1):
            covmat = np.linalg.inv(fishmats[:][bigcount])
        else:
            covmat = np.linalg.inv(fishmats)
        for i in xrange(0,numpars):
            for j in xrange(i+1,numpars):
                count = 1+(j-1)*(numpars-1) + i
                chi2 = np.array([[covmat[i,i],covmat[i,j]],[covmat[j,i],covmat[j,j]]])
                ansout, ansout2 = make_elsp(chi2)
                plt.subplot(numpars-1,numpars-1,count)
                if (bigcount == 0): 
                    plt.tick_params(size=14,width=thk,labelsize = 16)
                    plt.plot(0,0,'xk',mew=thk)
                    if (count ==1):
                        plt.ylabel(lab[1][:], fontsize=32,weight='bold')
                    if (count == 3):
                        plt.ylabel(lab[2][:], fontsize=32,weight='bold')
                        plt.xlabel(lab[0][:], fontsize=32,weight='bold')
                    if (count == 4):
                        plt.xlabel(lab[1][:], fontsize=32,weight='bold')
    
                plt.plot(ansout[0,:]/pars[i]*100. , ansout[1,:]/pars[j]*100.,ansout2[0,:]/pars[i]*100.,ansout2[1,:]/pars[j]*100.,'--',color=col,linewidth=thk)
#                plt.plot(ansout2[0,:]/pars[i]*100. , ansout2[1,:]/pars[j]*100.,'--',color=col,linewidth=thk)
        bigcount = bigcount + 1
        
    plt.savefig(saveFile, bbox_inches='tight',format='png')

    return

