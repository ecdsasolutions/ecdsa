import sys, json, random
def main():
    d = json.load(sys.stdin)
    p=int(d['p']);a=int(d['a']);n=int(d['n'])
    G=(int(d['Gx']),int(d['Gy']));P=(int(d['Px']),int(d['Py']))
    inv=lambda x:pow(x,-1,p)
    def add(A,B):
        if A is None: return B
        if B is None: return A
        if A[0]==B[0]:
            if (A[1]+B[1])%p==0: return None
            l=(3*A[0]*A[0]+a)*inv(2*A[1])%p
        else:
            l=(B[1]-A[1])*inv(B[0]-A[0])%p
        x=(l*l-A[0]-B[0])%p
        return (x,(l*(A[0]-x)-A[1])%p)
    def mul(k,A):
        R=None
        while k:
            if k&1:R=add(R,A)
            A=add(A,A);k>>=1
        return R
    # Pollard rho, distinguished points, negation-free, batch parallel walkers
    S=32
    ws=[]
    for _ in range(S):
        u=random.randrange(1,n);v=random.randrange(1,n)
        ws.append((u,v))
    ta=[random.randrange(1,n) for _ in range(S)];tb=[random.randrange(1,n) for _ in range(S)]
    tp=[add(mul(ta[i],G),mul(tb[i],P)) for i in range(S)]
    # note tp[i] could be None; ignore (negligible)
    L=len(ta)
    # use separate table
    tA=[random.randrange(1,n) for _ in range(32)];tB=[random.randrange(1,n) for _ in range(32)]
    T=[add(mul(tA[i],G),mul(tB[i],P)) for i in range(32)]
    # steps of the walk
    W=256
    pa=[random.randrange(1,n) for _ in range(W)];pb=[random.randrange(1,n) for _ in range(W)]
    pts=[add(mul(pa[i],G),mul(pb[i],P)) for i in range(W)]
    mask=(1<<12)-1  # distinguished: x & mask ==0 ; expected sqrt(n)~1.2e7*1.25 steps
    dp={}
    while True:
        # batch inversion
        den=[];
        for i in range(W):
            x,y=pts[i];j=x&31
            q=T[j]
            if q[0]==x: den.append(2*y%p or 1)
            else: den.append((q[0]-x)%p)
        pre=[1]*(W+1)
        for i in range(W): pre[i+1]=pre[i]*den[i]%p
        iv=pow(pre[W],-1,p)
        for i in range(W-1,-1,-1):
            dinv=iv*pre[i]%p; iv=iv*den[i]%p
            x,y=pts[i];j=x&31;q=T[j]
            if q[0]==x:
                if q[1]==y: l=(3*x*x+a)*dinv%p
                else: pa[i]=(pa[i]+tA[j])%n;pb[i]=(pb[i]+tB[j])%n; pts[i]=(x,y);continue
            else: l=(q[1]-y)*dinv%p
            nx=(l*l-x-q[0])%p
            pts[i]=(nx,(l*(x-nx)-y)%p)
            pa[i]=(pa[i]+tA[j])%n;pb[i]=(pb[i]+tB[j])%n
            if nx&mask==0:
                key=nx
                if key in dp:
                    a2,b2,y2=dp[key]
                    if y2==pts[i][1]:
                        da=(pa[i]-a2)%n;db=(b2-pb[i])%n
                        if db:
                            k=da*pow(db,-1,n)%n
                            if mul(k,G)==P:
                                print("k=%d"%k);return
                    else:
                        da=(pa[i]+a2)%n;db=(-b2-pb[i])%n
                        if db:
                            k=da*pow(db,-1,n)%n
                            if mul(k,G)==P:
                                print("k=%d"%k);return
                else: dp[key]=(pa[i],pb[i],pts[i][1])
                # restart walker
                pa[i]=random.randrange(1,n);pb[i]=random.randrange(1,n)
                pts[i]=add(mul(pa[i],G),mul(pb[i],P))
main()
