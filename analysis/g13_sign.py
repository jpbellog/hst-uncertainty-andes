import sys, numpy as np
def R(lat, lon, ang):
    la, lo, a = np.radians([lat, lon, ang])
    k = np.array([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])
    K = np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0]])
    return np.eye(3) + np.sin(a)*K + (1-np.cos(a))*K@K
def pole(M):
    a = np.degrees(np.arccos(np.clip((np.trace(M)-1)/2, -1, 1)))
    v = np.array([M[2,1]-M[1,2], M[0,2]-M[2,0], M[1,0]-M[0,1]]); v /= np.linalg.norm(v)
    if v[2] > 0: v, a = -v, -a   # pole in the southern hemisphere, as in Table 1
    return np.degrees(np.arcsin(v[2])), np.degrees(np.arctan2(v[1], v[0])) % 360, a
T1 = {"12o": (-84.78,190.99,-1.33), "13o": (-81.09,200.37,-3.17),
      "16y": (-70.72,338.78,-1.38), "18o": (-85.87,220.49,-4.48)}
S = {}
for l in open(sys.argv[1]).read().splitlines()[1:]:
    t = l.split()
    S[(t[0], t[1])] = R(*map(float, t[2:5]))
for an in T1:
    AE, AW = S[(an,"A.-E.A.")], S[(an,"A.-W.A.")]
    print(an, "Table 1 as printed:", T1[an])
    for nom, M in [("AE·AWT", AE@AW.T), ("AWT·AE", AW.T@AE), ("AET·AW", AE.T@AW), ("AW·AET", AW@AE.T)]:
        p = pole(M); print("   %-7s lat %7.2f lon %7.2f ang %6.2f" % (nom, *p))
