import sys,json,math
from array import array
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

def points(start,step,count,B=256):
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

def solve(n,G,P,points,mul):
    m=math.isqrt(n)//2+1
    ibits=(2*m+1).bit_length()
    fmask=(1<<(64-ibits))-1
    imask=(1<<ibits)-1
    size=1<<((3*m//2).bit_length())
    mask=size-1
    baby=array('Q',[0])*size
    for j,R in enumerate(points(G,G,m),1):
        x,y=R
        f=x&fmask
        h=f&mask
        while baby[h]:h=(h+1)&mask
        baby[h]=(f<<ibits)|(j<<1)|(y&1)
    s=2*m+1
    step=mul(s,G); step=(step[0],(-step[1])%p)
    for i,R in enumerate(points(P,step,(n+s-1)//s+1)):
        if R is None:return i*s%n
        f=R[0]&fmask
        h=f&mask
        while baby[h]:
            v=baby[h]
            if v>>ibits==f:
                match=v&imask
                j=match>>1
                if (match&1)!=(R[1]&1):j=-j
                k=(i*s+j)%n
                if mul(k,G)==P:return k
            h=(h+1)&mask
print('k='+str(solve(n,G,P,points,mul)))
