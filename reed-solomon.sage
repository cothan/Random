

def toBin(coordinate):
    # Convert variable to number 
    temp = bin(int( coordinate._int_repr() )).lstrip('0b')
    temp = '0'*(m - len(temp)) + temp
    return map(int, temp)

def toInt(value):
    return value._int_repr()

def myprint(vector):
    # print vector to number
    return map(lambda x: x.integer_representation(), nonzeroelements)



debug = True; verbosim = True; decimal = True;

m = 4; 
N = 2^m-1;  # with this setting the gRS code will be primitive
K = 12;            # dimension of the gRS code
D = N-K+1;         # gRS codes are MDS, so D = N-K+1.
# 
# checked m=4, K=13 with goppapol X^2 (k=11)
# checked m=4, K=12 with goppapol X^3 (k=7), X^3+X+1 (k=3)
# checked m=4, K=11 with goppapol X^4 (k=7), X^5 (k=5), X^4+X^2+X (k=1)
# checked m=4, K=10 with goppapol X^5 (k=5), X^5+X^2+1 (k=1)
# checked m=4, K= 9 with goppapol X^6+X+1 (k=0), X^6 (k=5), X^5 (k=4)
# checked m=4, K= 8 with goppapol X^7 (k=1), X^7+X+1 (k=0), X^6 (k=4), X^5 (k=0)
# checked m=4, K= 7 with goppapol X^6 and X^5 (k=0), X^4 and X^3 (k=2)
#
print 'N =',N,'K =',K,'D =',D;

F = GF(2); 
Phi.<x> = GF(2^m); 
if verbosim: 
    print 'F is',F,'with modulus',F.modulus();
    print 'Phi is',Phi,'with modulus',Phi.modulus();



def goppapolynomial(F,z): # return a Goppa polynomial in z over ring or field F
    X = PolynomialRing(F,repr(z)).gen();
    return X;



nonzeroelements = set([(x^i) for i in range(N)]);
if verbosim:
    print 'non zero elements are', map(lambda x: x.integer_representation(), nonzeroelements)
if len(nonzeroelements)<>N:
    print ('Alarm');



def printmatrix(nam,mat,decimal=True):
    #print mat.ncols() , '========', type(mat.nrows())
    #print mat.nrows(), '========', type(mat.nrows())
    if verbosim:
        print nam,'is a matrix of',mat.parent();
        if decimal:
            matint = matrix(ZZ,[[toInt(mat[i,j]) for j in xrange(mat.ncols())] \
                                                for i in xrange(mat.nrows())]);
            pretty_print(matint); 
        else:
            pretty_print(mat); 
    return


def printmatrix1(nam,mat,decimal=True):

    if verbosim:
        print nam,'is a matrix of',mat.parent();
        if decimal:
            matint = matrix(ZZ,[[int(mat[i,j]) for j in xrange(mat.ncols())] \
                                                for i in xrange(mat.nrows())]);
            pretty_print(matint); 
        else:
            pretty_print(mat); 
    return



z = var('z'); 
g = goppapolynomial(Phi,z);

while True:
   
    if g.is_irreducible():

        codelocators = [x^i for i in range(N)]; 
        
        if any( g(codelocators[i])==Phi(0) for i in range(N)):
            continue
        else:
            break





columnmultipliers = [1/g(codelocators[i]) for i in range(N)]; 
if verbosim:
    print 'codelocators are', myprint(codelocators);
    print 'columnmultipliers are', myprint(columnmultipliers);


print '========= REED SOLOMON'

H_gRS = matrix([[codelocators[j]^(i) for j in range(N)] for i in range(N-K)]);

H_gRS = H_gRS*diagonal_matrix([1/g(codelocators[i]) for i in range(N)]);

printmatrix('H_gRS',H_gRS);

print '========= H_gRS above'

####################################

A = matrix([[codelocators[j]^(i) for j in range(N)]for i in range(N-1)]);
A = A*diagonal_matrix([1/g(codelocators[i]) for i in range(N)]);
decimal = False;
#printmatrix('A',A);
B = matrix(Phi,N,N,0); B[:N-1,:]=A; B[N-1,N-1] = Phi(1);
#printmatrix('B',B);
if rank(B)==N:
    if verbosim:
        print 'Coefficient matrix has full rank!';
else:
    print ('Alarm');
    


rhs = matrix(Phi,N,1); 
rhs[N-1]=1; 
vdash = B\rhs; 
printmatrix("vdash'",vdash.transpose())

####################################

print '============= G_gRS'

G_gRS = matrix([[codelocators[j]^(i) for j in range(N)] for i in range(K)]);
# printmatrix('G_gRS',G_gRS);
G_gRS = G_gRS * diagonal_matrix(Phi,[vdash[i,0] for i in range(N)]);
printmatrix('G_gRS',G_gRS);


print '============= G_gRS above'

#print G_gRS


# printmatrix('check G_gRS*transpose(H_gRS)',G_gRS*transpose(H_gRS));
if G_gRS*transpose(H_gRS) == matrix(Phi,K,N-K):
    if verbosim:
        print 'G_gRS*transpose(H_gRS) is the', str(K)+'x'+str(N-K), 'zero matrix!';
else:
    print ('Alarm');

## HAVE G_GRS, I CAN FIND BACK H_GRS, VICE VERSA


if (N==G_gRS.nrows()+H_gRS.nrows()):
    if verbosim:
        print 'N == G_gRS.nrows()+H_gRS.nrows()'
else:
    print ('Alarm');


u = vector([x^randint(0,N) for _ in range(K)]); 
c = u*G_gRS; 

c = matrix(c)

if verbosim:
    print 'info word u =';
    print map(lambda x: x.integer_representation(), u)
    #pretty_print(u);
    print 'code word c ='; 
    print map(lambda x: x.integer_representation(), vector(c))
    #pretty_print(c);


if H_gRS*c.transpose()<>matrix(Phi,N-K,1):
    print "(H_gRS*(c'))' =",transpose(H_gRS*transpose(c));
    print ('Alarm');




H_Goppa = matrix(F,m*H_gRS.nrows(),H_gRS.ncols()); 
for i in range(H_gRS.nrows()):
    for j in range(H_gRS.ncols()):                 # transform each element of H_gRS
        
        H_Goppa[m*i:m*(i+1),j] = vector(toBin(H_gRS[i,j]));   # to a binary vector


r = rank(H_Goppa); print 'rank(H_Goppa) =',r;
H_Goppa = H_Goppa.rref()[:r,:]; # reduced row echelon form without all zero rows

printmatrix1('H_Goppa',H_Goppa);

if rank(H_Goppa)<>r:
    print 'removing zero rows from H_Goppa.rref() varies rank!'; 
    print ('Alarm');


n = H_Goppa.ncols(); 
k = n-r; 
print 'n =',n,'k =',k;

Krnl = H_Goppa.right_kernel(); 
G_Goppa = Krnl.basis_matrix();

print '========= G_GOPPA'
printmatrix1('G_Goppa',G_Goppa);
print '========= G_GOPPA above'



if H_Goppa*transpose(G_Goppa) == matrix(F,n-k,k):
    if verbosim:
        print H_Goppa*transpose(G_Goppa)
        print 'H_Goppa*transpose(G_Goppa) is the', str(n-k)+'x'+str(k),'zero matrix!';
else:
    print ('Alarm');




## HAVE H_GOPPA CAN FIND BACK G_GOPPA, VICE VERSA


if k>=n-(D-1)*m:
    if verbosim:
        print 'rank(G_Goppa) =',rank(G_Goppa),'>= N-(D-1)*m =',N-(D-1)*m;
else:
    print 'rank(G_Goppa) =',rank(G_Goppa),'must be >= N-(D-1)*m =',N-(D-1)*m;
    print ('Alarm');


R.<X> = PolynomialRing(Phi);
RmodGoppa = R.quotient(goppapolynomial(Phi,X),'X'); 
if verbosim:
    print 'R is',R;
    print 'RmodGoppa is',RmodGoppa;
    print 'g =',g,'with g.coeffs() =',g.coefficients()


def theta(j): 
    # returns (X-codelocators[j])^(-1) in RmodGoppa
    degGoppa = g.degree();
    sum = RmodGoppa(0);
    for ell in range(degGoppa):
        sc = Phi(0);
        for i in range(ell+1,degGoppa+1):
            sc += Phi(g.coeffs()[i])*(codelocators[j])^(i-ell-1);
        sum += sc*(X^ell); 
    return -sum/g(codelocators[j]);

if debug:
    OK = True; 
    for j in range(N):
        OK = OK and (X-codelocators[j])*theta(j) == RmodGoppa(1);
    if OK:
        if verbosim:
            print OK
    else:
        print ('Alarm');


def encode(u):
    return u*G_Goppa;



z = var('z'); PR = PolynomialRing(Phi,'z'); 
# init integer bigN
bigN = D-1;
# declare sigma's (sigma_{-1} ... sigma_bigN)
# declare omega's (omega_{-1} ... omega_bigN) 
sigma = vector(PR,bigN+2); omega = vector(PR,bigN+2); 
delta = vector(Phi,bigN+2);
# init sigma_{-1} and sigma_0 as well as omega_{-1} and omega_0 
sigma[-1+1] = PR(0); sigma[0+1] = PR(1);
flag = 2*bigN; # z^flag represents the rational function 1/z
omega[-1+1] = z^(flag); omega[0+1] = PR(0);
# init mu and delta
mu = -1; delta[-1+1] = 1;


def decode(y):
    y = matrix(y)
    s = H_gRS*y.transpose();
    if s==matrix(Phi,H_gRS.nrows(),1):
        print 'NO ERROR, RETURN DIRECTLY'
        return y;
    if verbosim:
        print "syndrome s' =",s.transpose();
        print "syndrome s' =",[int(s[_,0]) for _ in range(s.nrows())];
    b = PR([s[_,0] for _ in range(s.nrows())]);
    # init sigma_{-1} and sigma_0 as well as omega_{-1} and omega_0 
    sigma[-1+1] = PR(0); sigma[0+1] = PR(1);
    flag = 2*bigN; # z^flag represents the rational function 1/z
    omega[-1+1] = z^(flag); omega[0+1] = PR(0);
    # init mu and delta
    mu = -1; delta[-1+1] = 1;
    for i in range(bigN):
        delta[i+1] = (sigma[i+1]*b).coeffs()[i];
        sigma[i+1+1] = sigma[i+1](z)\
                       -(delta[i+1]/delta[mu+1])*z^(i-mu)*sigma[mu+1](z);
        if (omega[mu+1].degree()==flag):
            omega[i+1+1] = omega[i+1](z)\
                       -(delta[i+1]/delta[mu+1])*z^(i-mu-1);
        else:
            omega[i+1+1] = omega[i+1](z)\
                       -(delta[i+1]/delta[mu+1])*z^(i-mu)*omega[mu+1](z);
        rord = max(sigma[i+1].degree(),1+omega[i+1].degree()); # recurrence order
        if (delta[i+1]<>0)and(2*rord<=i):
            mu = i;
        if verbosim:
            print 'after step i =',i,'we have mu =',mu;
            print 'new delta[i] =',delta[i+1],'delta[mu] =',delta[mu+1];
            print 'new sigma[i+1] =',sigma[i+2],'sigma[mu] =',sigma[mu+1];
            print 'new omega[i+1] =',omega[i+2],'omega[mu] =',omega[mu+1];
    ELP = sigma[bigN+1]; # the Error Locator Polynomial
    # compute zeroes of ELP, compute error positions, compute error vector ee
    ee = vector(F,[0 for _ in range(n)]);
    for i in range(N):
        if (ELP(x^i)==Phi(0)):
            if verbosim:
                print 'x^'+str(i),'=',x^i,'is a root of ELP:', \
                      'an error occured in position', N-i;
            ee[mod(N-i,N)] += 1;
    cc = vector(y)+ee;
    if verbosim:
        print 'computed error =',ee; 
        print 'corrected y+ee =',y+ee;
    return cc;



tau = floor((D-1)/2); print 'D =',D,'tau =',tau,'k =',k; cntr = 0; verbosim = True;
wt = walltime(); ntrials =1;
for trials in range(ntrials):
    u = vector(F,[randint(0,1) for _ in range(k)]); 
    c = encode(u); 
    if verbosim:
        print 'info word u =',u; print 'code word c =', c;
    e = vector(F,[0 for _ in range(n)]); 
    for trial in range(tau):
        j = randint(0,n-1); e[j] += 1; 
    y = c+e; 
    if verbosim:
        print 'error vec e =', e; print 'received  y =',y;
    verbosim = False; cc = decode(y); verbosim = True;
    if (cc == c):
        cntr += 1;
        if verbosim:
            print 'corrected received == sent';
    else:
        print 'info word u =',u; print 'code word c =',c;
        print 'error vec e =',e; print 'recieved  y =',y;
        print 'corrected y =',cc;
        print html('<font color ="red">Alarm</font>');
print 'to execute decode takes about',(walltime()-wt)/cntr,'seconds on average'