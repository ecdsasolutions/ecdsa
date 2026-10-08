import sys, json, math
d = json.load(sys.stdin)
p=int(d["p"]);a=int(d["a"]);n=int(d["n"])
G=(int(d["Gx"]),int(d["Gy"]));P=(int(d["Px"]),int(d["Py"]))
def add(A,B):
    if A is None: return B
    if B is None: return A
    if A[0]==B[0]:
        if (A[1]+B[1])%p==0: return None
        l=(3*A[0]*A[0]+a)*pow(2*A[1],-1,p)%p
    else:
        l=(B[1]-A[1])*pow(B[0]-A[0],-1,p)%p
    x=(l*l-A[0]-B[0])%p
    return (x,(l*(A[0]-x)-A[1])%p)
def neg(A): return None if A is None else (A[0],(-A[1])%p)
m=math.isqrt(n)+1
baby={}
R=None
for j in range(m):
    baby.setdefault(R,j) if R is None else baby.setdefault(R,j)
    R=add(R,G)
S=R  # m*G
Q=P
for i in range(m+2):
    if Q in baby:
        k=(i*m+baby[Q])%n
        if k and (lambda: True)():
            print("k=%d"%k);break
    Q=add(Q,neg(S))
