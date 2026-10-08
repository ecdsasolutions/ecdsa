import sys,json,math
D=json.load(sys.stdin)
p,a,b,Gx,Gy,n,Px,Py=(int(D[k]) for k in ('p','a','b','Gx','Gy','n','Px','Py'))
G=(Gx,Gy); P=(Px,Py)
def add(U,V):
    if U is None:return V
    if V is None:return U
    x,y=U; X,Y=V
    if x==X:
        if (y+Y)%p==0:return None
        s=(3*x*x+a)*pow(2*y,-1,p)%p
    else:s=(Y-y)*pow(X-x,-1,p)%p
    z=(s*s-x-X)%p
    return z,(s*(x-z)-y)%p
def mul(k,U):
    R=None
    while k:
        if k&1:R=add(R,U)
        U=add(U,U);k>>=1
    return R

def points(start,step,count,B=32):
    lanes=[start]
    for _ in range(1,min(B,count)):lanes.append(add(lanes[-1],step))
    jump=mul(len(lanes),step)
    X,Y=jump
    off=0
    while off<count:
        for R in lanes[:min(len(lanes),count-off)]:yield R
        off+=len(lanes)
        if off>=count:break
        ds=[]; prefixes=[]; acc=1
        for R in lanes:
            d=(X-R[0])%p if R is not None else 0
            ds.append(d); prefixes.append(acc)
            if d:acc=acc*d%p
        inv=pow(acc,-1,p)
        for h in range(len(lanes)-1,-1,-1):
            d=ds[h]
            if not d:lanes[h]=add(lanes[h],jump);continue
            invd=inv*prefixes[h]%p
            inv=inv*d%p
            x,y=lanes[h]
            s=(Y-y)*invd%p
            z=(s*s-x-X)%p
            lanes[h]=(z,(s*(x-z)-y)%p)

m=math.isqrt(n)//2+1
baby={}
for j,R in enumerate(points(G,G,m),1):baby[R[0]]=(j<<1)|(R[1]&1)
s=2*m+1
step=mul(s,G); step=(step[0],(-step[1])%p)
for i,R in enumerate(points(P,step,(n+s-1)//s+1)):
    if R is None:
        print('k='+str(i*s%n));break
    match=baby.get(R[0])
    if match is not None:
        j=match>>1
        if (match&1)!=(R[1]&1):j=-j
        print('k='+str((i*s+j)%n));break
