import numpy as np
from scipy import special
from scipy import stats
from scipy.stats import beta
from scipy.stats import binom
from scipy.stats import gamma
import pandas as pd
from numpy.linalg import eig
from numpy.linalg import multi_dot
import math 
import pickle


def fr(r,pos,neg,u,v):
    '''
    fr(r,pos,neg,u,v)
        Computes the numerator of the binomial posterior.
    
    r - seroprevalence [0,1]
    pos - number of positive serological results
    neg - number of negative serological results
    u - false positive rate (1-specificity)
    v - false negative rate (1-sensitivity)

    for log version see log_fr
    '''
    t = u+(1-u-v)*r
    return (t**pos) * ((1-t)**neg)

def log_fr(r,pos,neg,u,v):
    '''
    log_fr(r,pos,neg,u,v)
        Computes the log of the numerator of the binomial posterior.
    
    r - seroprevalence [0,1]
    pos - number of positive serological results
    neg - number of negative serological results
    u - false positive rate (1-specificity)
    v - false negative rate (1-sensitivity)

    for non-log version see fr
    '''
    t = u+(1-u-v)*r
    return pos*np.log(t) + neg*np.log(1-t)

def sample_post_r_log(pos,neg,u,v,size=1):
    '''
    sample_post_r_log(pos,neg,u,v,size=1)
        Samples from the posterior distribution of seroprevalence 
        under a uniform prior, given data and parameters below.
        (default 1 sample)
        Uses the log of the posterior numerator to avoid underflow issues.
        After max_failures of the accept-reject, this code returns -1.
        This is due to our lazy proposal distribution: unif[0,1]. We could
        do better, but didn't optimize right away. 

    pos - number of positive serological results
    neg - number of negative serological results
    u - false positive rate (1-specificity)
    v - false negative rate (1-sensitivity)
    '''
    rm = (pos/(pos+neg)-u)/(1-u-v)
    
    if (rm>0) and (rm < 1):
        log_M = np.array([log_fr(0,pos,neg,u,v),log_fr(1,pos,neg,u,v),
            log_fr(rm,pos,neg,u,v)]).max()
    else:
        log_M = np.array([log_fr(0,pos,neg,u,v),log_fr(1,pos,neg,u,v)]).max()
    # uniform accept/reject algorithm
    samples = []
    n_accepted = 0
    max_failures = 1e6
    failures = 0
    while n_accepted < size:
        r_proposal = np.random.rand()
        if (np.log(np.random.rand()) < (log_fr(r_proposal,pos,neg,u,v)-log_M)):
            samples.append(r_proposal)
            n_accepted += 1
            failures = 0
        else:
            failures += 1
        if failures == max_failures:
            return [-1]
    if size==1:
        return samples[0]
    else:
        return np.array(samples)


def sample_post_r(pos,neg,u,v,size=1):
    '''
    sample_post_r(pos,neg,u,v,size=1)
        Samples from the posterior distribution of seroprevalence 
        under a uniform prior, given data and parameters below.
        (default 1 sample)
        Uses the posterior numerator without logs, so may run into underflow issues.
        After max_failures of the accept-reject, this code returns -1.
        This is due to our lazy proposal distribution: unif[0,1]. We could
        do better, but didn't optimize right away. 

    pos - number of positive serological results
    neg - number of negative serological results
    u - false positive rate (1-specificity)
    v - false negative rate (1-sensitivity)
    '''
    rm = (pos/(pos+neg)-u)/(1-u-v)
    
    if (rm>0) and (rm < 1):
        M = np.array([fr(0,pos,neg,u,v),fr(1,pos,neg,u,v),
            fr(rm,pos,neg,u,v)]).max()
    else:
        M = np.array([fr(0,pos,neg,u,v),fr(1,pos,neg,u,v)]).max()
    # uniform accept/reject algorithm
    samples = []
    n_accepted = 0
    max_failures = 1e6
    failures = 0
    while n_accepted < size:
        r_proposal = np.random.rand()
        if (np.random.rand() < (fr(r_proposal,pos,neg,u,v)/M)):
            samples.append(r_proposal)
            n_accepted += 1
            failures = 0
        else:
            failures += 1
        if failures == max_failures:
            return [-1]
    if size==1:
        return samples[0]
    else:
        return np.array(samples)

def interpolate_cos(x, x0, y0, x1, y1):
    '''
    interpolate_cos(x, x0, y0, x1, y1)
        The paper by Davies et al uses cosine interpolation. 
        We didn't see this built in anywhere, so we coded it up from docs.
        
    Interpolate between points (x0, y0) and (x1, y1) using cosine interpolation. 
    For x < x0, returns y0; for x > x1, returns y1; 
    For x0 < x < x1, returns the cosine interpolation between y0 and y1.
    '''
    if x < x0:
        y = y0
    elif x > x1:
        y = y1
    else:
        y = y0 + (y1 - y0) * (0.5 - 0.5 * math.cos(math.pi * (x - x0) / (x1 - x0)))
    return y

def get_clinical_fraction(age_bin, age_cp=[14,55,64], y_cp=[0.056,0.49,0.74]):
    '''
    get_clinical_fraction(age_bin, age_cp=[14,55,64], y_cp=[0.056,0.49,0.74]
        The paper by Davies et al uses cosine interpolation with control points
        to smooth their estimates of clinical cases by age.

    Given a set of points age_cp and y_cp, return an interpolated value for a particular age_bin.
    Control points default to values from Davies paper.
    Age bins are hardcoded assuming age groups by 5, max age = 80 (i.e. last age bin is 75-80)
    '''
    ages =  np.linspace(2.5, 77.5, 16, endpoint = True) 
    mat = np.zeros(len(ages))
    for j in range(len(ages)):
        if ages[j] < age_cp[1]:
            mat[j] = interpolate_cos(ages[j], age_cp[0], y_cp[0], age_cp[1], y_cp[1])
        else:
            mat[j] = interpolate_cos(ages[j], age_cp[1], y_cp[1], age_cp[2], y_cp[2])
    y_i = mat[age_bin]
    return(mat, y_i)

def get_reff(my_N,my_r=np.zeros(16)):
    '''
    get_reff(N,r=np.zeros(16))
        Compute R_effective from a next-generation matrix a vector of 
        seropositivity estimates. In the absence of a vector, returns R0.

    N - next generation matrix as numpy array
    r - seroprevalence estimates [0,1] for subpopulations corresp. to N.
    '''
    X = np.matmul(np.diag(1-my_r),my_N)
    evals,evecs = eig(X)
    evals = np.abs(evals)
    return np.max(evals)

def get_evec(N):
    _,evecs = eig(N)
    evec = evecs[:,0]
    evec = np.array(list(np.array(np.abs(evec).transpose())[0]))
    return evec/np.sum(evec)

def sample_posterior_r_mcmc_hyper(samps,posi,ni,se,sp,gam0):
    '''
    sample_posterior_r_mcmc_hyper(samps,posi,ni,se,sp,gam0)
        Sample from the seroprevalence posterior.
        Uses the Bayesian hierarchy described in the manuscript.
        NOTE: this function is not used due to the fact that, for whatever
        reason, the R code that Bailey wrote is actually faster, even though
        this is a line for line implementation of that R code. Unclear why,
        though we suspect it has to do with scipy's implementations of logprobs.

    samps - number of samples to obtain
    posi - array of positive serological tests in each age bin
    ni - array of total tests in each age bin
    se - sensitivity
    sp - specificity
    gam0 - hyperprior variance parameter
    
    returns a LIST of arrays of samples. 
    '''
    fn = 1-se
    fp = 1-sp
    posi = np.array(posi)
    ni = np.array(ni)

    nu = 1
    if np.any(ni==0)==True:
        ri = (posi+1)/(ni+2)
    else:    
        ri = (posi+1)/(ni+1)
    r = np.mean(ri)
    gam = r*(1-r)/np.max([np.var(ri),.001]) - 1
    if gam <= 0:
        gam = gam0
    
    ri_posterior = []
    r_posterior = []
    gam_posterior = []
    
    # MCMC params
    delta_r = 200
    delta_ri = 100
    delta_gam = 10
    thin = 50
    burn_in = 2*thin
    
    for s in range(samps*thin+burn_in):  
        # MH step for gamma
        # propose gam_prop | gamma ~ gamma(delta_gam,scale=gamma/delta_gam)
        gam_prop = gamma.rvs(delta_gam,scale=gam/delta_gam)
        ar_gam = (
            np.sum(beta.logpdf(ri,r*gam_prop,(1-r)*gam_prop)) -
            np.sum(beta.logpdf(ri,r*gam,(1-r)*gam)) +
            gamma.logpdf(gam_prop, nu,       scale=gam0/nu) -
            gamma.logpdf(gam,      nu,       scale=gam0/nu) +
            gamma.logpdf(gam,      delta_gam,scale=gam_prop/delta_gam) -
            gamma.logpdf(gam_prop, delta_gam,scale=gam/delta_gam)
        )
        if(np.log(np.random.rand()) <ar_gam):
            gam = gam_prop
        
        # MH step to update r
        # propose r_prop | r ~ Beta(r*delta_r,(1-r)*delta_r)
        r_prop = beta.rvs(r*delta_r,(1-r)*delta_r)
        ar_r = (
            np.sum(beta.logpdf(ri,r_prop*gam,(1-r_prop)*gam)) - 
            np.sum(beta.logpdf(ri,r*gam,(1-r)*gam)) +
            beta.logpdf(r,r_prop*delta_r,(1-r_prop)*delta_r) - 
            beta.logpdf(r_prop,r*delta_r,(1-r)*delta_r)
        )
        if np.log(np.random.rand()) < ar_r:
            r = r_prop
        #MH step to update each ri
        #propose ri_prop | ri ~ B(ri*delta_ri,(1-ri)*delta_ri)
        for k in range(len(ri)):
            ri_prop = beta.rvs(ri[k]*delta_ri,(1-ri[k])*delta_ri)
            ar_ri = (
                binom.logpmf(posi[k],ni[k],ri_prop*(1-fn)+(1-ri_prop)*fp) -
                binom.logpmf(posi[k],ni[k],ri[k]*(1-fn)+(1-ri[k])*fp)+
                beta.logpdf(ri_prop,r*gam,(1-r)*gam) -
                beta.logpdf(ri[k],r*gam,(1-r)*gam) +
                beta.logpdf(ri[k],ri_prop*delta_ri,(1-ri_prop)*delta_ri) -
                beta.logpdf(ri_prop,ri[k]*delta_ri,(1-ri[k])*delta_ri))

            if np.log(np.random.rand())<ar_ri:
                ri[k] = ri_prop
        if (s%thin==0) and (s>=burn_in):
            r_posterior.append(r.copy())
            ri_posterior.append(ri.copy())
            gam_posterior.append(gam)
        
    return ri_posterior #r_posterior,gam_posterior

def NGM(Cmatrix,uvec,yvec,f,EdS,EdP,EdC):
    '''
    NGM(Cmatrix,uvec,yvec,f,EdS,EdP,EdC)
        Computes the next generation matrix from Davies et al 2020

    Cmatrix - 16x16 contact matrix (POLYMOD)
    uvec - suceptibility (16 array)
    yvec - probability of clinical infection (16 array)
    f - relative infectiousness of subclinical cases
    EdS - expected duration subclinical infectiousness
    EdP - expected duration preclinical infectiousness
    EdC - expected duration clinical infectiousness
    See Davies et al Table in Methods for values and original references

    returns the next generation matrix N
    '''
    a = EdP+EdC-f*EdS
    b = f*EdS
    Du = np.diag(uvec)
    Dq = np.diag(a*yvec+b)
    N = multi_dot([Du,Cmatrix,Dq])
    return N

def p_seropositive_r(r,sensitivity,specificity):
    '''
    p_seropositive_r(r,sensitivity,specificity)
        Computes the probability that a test comes back seropositive, given
        the true seroprevalence in the population and test specs

    r - true seroprevalence
    sensitivity - sensitivity
    specificity - specificity
    '''
    return r*sensitivity+(1-r)*(1-specificity)

def simulate_serology(ri,ni,sensitivity=1,specificity=1):
    '''
    simulate_serology(ri,ni,sensitivity=1,specificity=1)
        Simulates the serological sampling process for ni samples from  
        a set of subpopulations with true seroprevalences ri, 
        and given test characteristics

    ri - true seroprevalences array
    ni - number of samples in each subpopulation array
    sensitivity - sensitivity
    specificity - specificity
    '''
    pi = [p_seropositive_r(r,sensitivity,specificity) for r in ri]
    return [np.random.binomial(ni[idx],pi[idx]) for idx in range(len(pi))]


def load_contact_data():
    '''
    load_contact_data()
        Loads the contact data from Prem et al, which 
        has been preprocessed and pickled already. 

    It returns a big dictionary which can be called by country and then by contact type. 
    '''
    with open('../contact_matrices_premetal2017/all_data.pickle','rb') as file:
        x = pickle.load(file)
    return x

def get_population_demographics(country):
    '''
    get_population_demographics(country)
        Loads the population demographics in 5-year age bins from UN Data.
        This is a convenience function that wraps a Pandas call, which reads Excel.

    country - name of country to be loaded. Open the xlsx file to see all names available.
    returns a vector, which sums to 1, of the population distribution.
    '''
    pop_file_name = '../population_data_UNWPP/WPP2019_POP_F07_1_POPULATION_BY_AGE_BOTH_SEXES.xlsx'
    sheet_name = 'ESTIMATES'
    df = pd.read_excel(pop_file_name,sheet_name=sheet_name)
    popdata = df.loc[((df['Country']==country) & (df['Year']==2020))]
    n = [
        popdata['0-4'].values[0],
        popdata['5-9'].values[0],
        popdata['10-14'].values[0],
        popdata['15-19'].values[0],
        popdata['20-24'].values[0],
        popdata['25-29'].values[0],
        popdata['20-24'].values[0],
        popdata['35-39'].values[0],
        popdata['40-44'].values[0],
        popdata['45-49'].values[0],
        popdata['50-54'].values[0],
        popdata['55-59'].values[0],
        popdata['60-64'].values[0],
        popdata['65-69'].values[0],
        popdata['70-74'].values[0],
        popdata['75-79'].values[0],
    ]
    n = np.array(n)/np.sum(n)
    return n

def mle_reff(posi,ni,fp,fn,N):
    '''
    This returns the maximum likelihood Reffective value, given serological data
    across subpopulations of a next-generation matrix. 
    Note that the maximum likelihood may occur at the edge of the interval, 0 or 1.

    posi - array of positive serological sample counts
    ni - array of total serological sample counts
    fp - false positive rate (1-specificity)
    fn - false negative rate (1-sensitivity)
    N - next generation matrix

    returns the MLE R_effective
    '''
    dirty = (posi/ni - fp)/(1-fp-fn)
    cleanleft = [np.max([0,x]) for x in dirty]
    cleanright = [np.max([0,x]) for x in cleanleft]
    return get_reff(N,np.array(cleanright))





##### DEPRECATED #####
class MyRV(stats.rv_continuous):
    '''
    This CLASS is a subclass of the scipy continuous random variable class.
    As you can see below, we're just going to define a PDF for the posterior distribution.
    Then, we can call it just like we would any other scipy distribution. 

    Unfortunately, I don't think we use this much because it turned out to be slower than
    the accept-reject approach. Keeping in the codebase because I'd rather not have to
    relearn it, if I need it. ;)
    '''
    def _pdf(self,x,pos,neg,sensitivity,specificity):
        fn = (1-sensitivity)
        fp = (1-specificity)
        p = x*(1-fn)+(1-x)*fp
        num = (p**pos)*((1-p)**neg)
        denom = (special.betainc(pos+1,neg+1,1-fn)-special.betainc(pos+1,neg+1,fp))*special.beta(pos+1,neg+1) /(1-fp-fn)
        return num/denom
    def _argcheck(self, *args):
        """Default check for correct values on args and keywords.
        Returns condition array of 1's where arguments are correct and
         0's where they are not.
        """
        cond = 1
        for arg in args:
            cond = np.logical_and(cond, (np.asarray(arg) > -np.inf))
        return cond
