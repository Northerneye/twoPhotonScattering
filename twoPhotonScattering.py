import sympy as sy
from itertools import permutations
import dill
import pickle
import time
total_start = time.time()
#
#
#Two Photon Scattering Diagram
#
#                K-Q2-Q1
#P1 ~~~~~~~~|----<-----|~~~~~~~~ Q1
#           |          |
#           V K-P2     ^ K-Q2
#           |          |
#P2 ~~~~~~~~|---->-----|~~~~~~~~ Q2
#                K
#
#
#TWO PHOTON EQUATION COMES IN THREE PARTS
#  M_a, M_b, M_c
#
#  UNPOLARIZED VVV
#  M = integral_K TRACE((kslash + m)(kslash + q2slash + m)(kslash + q2slash + q1slash + m)(kslash + p2slash + m))/((k^2 - m^2)((k+q2)^2 - m^2)((k+q2+q1)^2 - m^2)((k+p2)^2 - m^2))
#  
#  POLARIZED VVV
#  M = 
sy.var('identity4 overlap F a epsilon g m holder1x holder1y holder1z holder1u holder1d holder2x holder2y holder2z holder2u holder2d holder3x holder3y holder3z holder3u holder3d holder4x holder4y holder4z holder4u holder4d P1 p1t p1x p1y p1z P2 p2t p2x p2y p2z P3 p3t p3x p3y p3z P4 p4t p4x p4y p4z K1 k1t k1x k1y k1z K2 k2t k2x k2y k2z K3 k3t k3x k3y k3z Q1 q1t q1x q1y q1z Q2 q2t q2x q2y q2z s1 s1u s1d s2 s2u s2d s3 s3u s3d s4 s4u s4d z1 z1u z1d z2 z2u z2d UP1 UP2 VbarP3 UP4 UbarQ1 UbarQ2 Gamma Pauli y y24')

if(True): #Creates new overlap |integral(<q1,q2|F(a)|p1,p2,p3,p4>)|^2
    s1 = sy.Matrix([s1u, s1d]) #change for different spin values
    s2 = sy.Matrix([s2u, s2d])
    z1 = sy.Matrix([z1u, z1d])
    z2 = sy.Matrix([z2u, z2d])

    identity4 = sy.Matrix([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])
    g = sy.Matrix([[
        [-1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ]])
    Gamma = [
        sy.Matrix(
        [[1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, -1, 0],
        [0, 0, 0, -1]]),
        sy.Matrix(
        [[0, 0, 0, 1],
        [0, 0, 1, 0],
        [0, -1, 0, 0],
        [-1, 0, 0, 0]]),
        sy.Matrix(
        [[0, 0, 0, -sy.I],
        [0, 0, sy.I, 0],
        [0, sy.I, 0, 0],
        [-sy.I, 0, 0, 0]]),
        sy.Matrix(
        [[0, 0, 1, 0],
        [0, 0, 0, -1],
        [-1, 0, 0, 0],
        [0, 1, 0, 0]])
    ]
    Pauli = sy.Matrix([
        [[1, 0],
        [0, 1]],
        [[0, 1],
        [1, 0]],
        [[0, -sy.I],
        [sy.I, 0]],
        [[1, 0],
        [0, -1]]
    ])
    P1 = sy.Matrix([p1t, p1x, p1y, p1z])
    P2 = sy.Matrix([p2t, p2x, p2y, p2z])
    Q1 = sy.Matrix([q1t, q1x, q1y, q1z])
    Q2 = sy.Matrix([q2t, q2x, q2y, q2z])

    K = sy.Matrix([k1t, k1x, k1y, k1z])

    kslash = K[0]*Gamma[0]+K[1]*Gamma[1]+K[2]*Gamma[2]+K[3]*Gamma[3]+m*identity4
    q1slash = Q1[0]*Gamma[0]+Q1[1]*Gamma[1]+Q1[2]*Gamma[2]+Q1[3]*Gamma[3]+m*identity4
    q2slash = Q2[0]*Gamma[0]+Q2[1]*Gamma[1]+Q2[2]*Gamma[2]+Q2[3]*Gamma[3]+m*identity4
    p2slash = P2[0]*Gamma[0]+P2[1]*Gamma[1]+P2[2]*Gamma[2]+P2[3]*Gamma[3]+m*identity4
    

    print("Creating Simple Momentum State Amplitude (y)...")
    start = time.time()
    y=0 
    myi = 0
    for mu in range(4):#Summing over all einstein indicies
        for nu in range(4):
            for sigma in range(4):
                for tau in range(4):
                    #y = y+sy.simplify((UbarQ1*Gamma[mu]*UP1)[0]*(VbarP3*Gamma[nu]*k2slash*Gamma[sigma]*UP2)[0]*(UbarQ2*Gamma[tau]*UP4)[0]*g[mu][nu]*g[sigma][tau])
                    y = y+(UbarQ1*Gamma[mu]*UP1)[0]*(VbarP3*Gamma[nu]*k2slash*Gamma[sigma]*UP2)[0]*(UbarQ2*Gamma[tau]*UP4)[0]*g[mu][nu]*g[sigma][tau]
    print(str(time.time()-start)+" seconds")

    y = y/(K1.T*K1)[0]/(K3.T*K3)[0]/((K2.T*K2)[0]-m*m+sy.I*epsilon)
    y = y.subs([(p1t, sy.sqrt(p1x*p1x + p1y*p1y + p1z*p1z + m*m)), (p2t, sy.sqrt(p2x*p2x + p2y*p2y + p2z*p2z + m*m)), (p3t, sy.sqrt(p3x*p3x + p3y*p3y + p3z*p3z + m*m)), (p4t, sy.sqrt(p4x*p4x + p4y*p4y + p4z*p4z + m*m)), (q1t, sy.sqrt(q1x*q1x + q1y*q1y + q1z*q1z + m*m)), (q2t, sy.sqrt(q2x*q2x + q2y*q2y + q2z*q2z + m*m))])
    #substitutes actual values for the variables
    #spin of electron 1 down, (up+down)^2 !must! = 1

    #now for all permutations of the particles
    y24 = 0
    for p in permutations([[p1x,p1y,p1z,s1u,s1d], [p2x,p2y,p2z,s2u,s2d], [p3x,p3y,p3z,s3u,s3d], [p4x,p4y,p4z,s4u,s4d]]):
        print("\nAdding permutation: "+str(p))
        start = time.time()
        y24 = y24+y.subs([(p1x, holder1x), (p1y, holder1y), (p1z, holder1z), (s1u, holder1u), (s1d, holder1d), (p2x, holder2x), (p2y, holder2y), (p2z, holder2z), (s2u, holder2u), (s2d, holder2d), (p3x, holder3x), (p3y, holder3y), (p3z, holder3z), (s3u, holder3u), (s3d, holder3d), (p4x, holder4x), (p4y, holder4y), (p4z, holder4z), (s4u, holder4u), (s4d, holder4d)]).subs([(holder1x, p[0][0]), (holder1y, p[0][1]), (holder1z, p[0][2]), (holder1u, p[0][3]), (holder1d, p[0][4]), (holder2x, p[1][0]), (holder2y, p[1][1]), (holder2z, p[1][2]), (holder2u, p[1][3]), (holder2d, p[1][4]), (holder3x, p[2][0]), (holder3y, p[2][1]), (holder3z, p[2][2]), (holder3u, p[2][3]), (holder3d, p[2][4]), (holder4x, p[3][0]), (holder4y, p[3][1]), (holder4z, p[3][2]), (holder4u, p[3][3]), (holder4d, p[3][4])])
        print(str(time.time()-start)+" seconds")
    print("Adding Permutation of Outgoing Particles...")
    start = time.time()
    y24 = y24 + y24.subs([(q1x, holder1x), (q1y, holder1y), (q1z, holder1z), (s1u, holder1u), (s1d, holder1d), (q2x, holder2x), (q2y, holder2y), (q2z, holder2z), (s2u, holder2u), (s2d, holder2d)]).subs([(holder1x, q2x), (holder1y, q2y), (holder1z, q2z), (holder1u, s2u), (holder1d, s2d), (holder2x, q1x), (holder2y, q1y), (holder2z, q1z), (holder2u, s1u), (holder2d, s1d)])
    print(str(time.time()-start)+" seconds")


    print("\nConstructing Arbitrary Wavefunction F(a)...")
    #NEED TO ADD ARBITRARY WAVEFUNCTION F(a,p) => integral(F*y24)
    start = time.time()
    F=sy.exp(-a*P1.T*P1)*sy.exp(-a*P2.T*P2)*sy.exp(-a*P3.T*P3)*sy.exp(-a*P4.T*P4)#Input States
    F= F*sy.exp(-a*Q1.T*Q1)*sy.exp(-a*Q2.T*Q2)#Output States, NEED TO FIND A BETTER OUTPUT THAN A GAUSSIAN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print(str(time.time()-start)+" seconds")

    print("\nJoining Functions, (F*y24)")
    start = time.time()
    y24 = F*y24
    print(str(time.time()-start)+" seconds")

    print("\nEnforcing Momentum Conservation...")
    start = time.time()
    y24 = y24.subs([(q2x, p1x+p2x+p3x+p4x-q1x), (q2y, p1y+p2y+p3y+p4y-q1y), (q2z, p1z+p2z+p3z+p4z-q1z)])#MOMENTUM CONSERVATION
    print(str(time.time()-start)+" seconds")

    """
    print("\nSimplifying Amplitude...")
    start = time.time()
    y24 = sy.simplify(y24)
    print(str(time.time()-start)+" seconds")
    #"""

    print("\nIntegrating over Momentum...")
    start = time.time()
    overlap = sy.integrate(y24, (p1x, -sy.oo, sy.oo), (p1y, -sy.oo, sy.oo), (p1z, -sy.oo, sy.oo), (p2x, -sy.oo, sy.oo), (p2y, -sy.oo, sy.oo), (p2z, -sy.oo, sy.oo), (p3x, -sy.oo, sy.oo), (p3y, -sy.oo, sy.oo), (p3z, -sy.oo, sy.oo), (p4x, -sy.oo, sy.oo), (p4y, -sy.oo, sy.oo), (p4z, -sy.oo, sy.oo), (q1x, -sy.oo, sy.oo), (q1y, -sy.oo, sy.oo), (q1z, -sy.oo, sy.oo))
    #Still dependent on spins, m, epsilon, and input/output state parameter a
    print(str(time.time()-start)+" seconds")
    print(overlap)

    #Lambdafication and Saving
    print("\nLambdifying Expression for Overlap...")
    start = time.time()
    lambda_overlap = sy.lambdify((a, s1u, s1d, s2u, s2d, s3u, s3d, s4u, s4d, z1u, z1d, z2u, z2d, m, epsilon), overlap, modules="numpy")
    print(str(time.time()-start)+" seconds")
    

    print("\nSAVING...")
    start = time.time()
    dill.settings['recurse'] = True
    dill.dumps(lambda_overlap)
    print(str(time.time()-start)+" seconds")

    """#For Momentum Eigenstate y24, (without F).  Allows you to investigate momentum eigenstate scattering
    print("Inputting Desired Particle Momentum...")
    result = result.subs([(p1x, 1), (p1y, 1), (p1z, 0)])#momentum vector of electron 1
    result = result.subs([(p2x, 1), (p2y, 0), (p2z, 2)])
    result = result.subs([(p3x, 1), (p3y, 0), (p3z, 1)])#incoming electron 3 is positron
    result = result.subs([(p4x, 1), (p4y, 3), (p4z, 0)])

    result = result.subs([(q1x, 0), (q1y, 1), (q1z, 0)])
    #result = result.subs([(q2x, 0), (q2y, 0), (q2z, 0)])#q2z is automatically determined due to momentum conservation
    """

    print("Inputting Desired Electron Spins...")
    result = y24.subs([(s1u, 1),(s1d, 0)])#spin of electron 1 up
    result = result.subs([(s2u, 1),(s2d, 0)])
    result = result.subs([(s3u, 1),(s3d, 0)])
    result = result.subs([(s4u, 1),(s4d, 0)])
    result = result.subs([(z1u, 1),(z1d, 0)])
    result = result.subs([(z2u, 1),(z2d, 0)])
    
    print("Inputting Electron Mass...")
    result = result.subs(m, 1)

    print("Approximating Epsilon...")
    result = result.subs(epsilon, 0.0001)

    print("Inputting Test a...")
    result = result.subs(a, 1)

    print("\n")
    print("TEST RESULT IS...")
    print(result) #zoo = complex infinity
    print("\n")
    print("OR APPROXIMATELY:")
    result = result.evalf()
    print(result)
    print("\n")
    print("|<q1,q2|F(a)>|^2:")#NEED TO FIND A BETTER OUTPUT STATE
    result = (sy.Pow(sy.Abs(result), 2))
    print(result)
    print("ENTIRE PROCESS TOOK "+str(total_start - time.time()))


#Now the Fun Begins
# time to Test over values of F(a), and use gradient descent to find the local max efficiency

