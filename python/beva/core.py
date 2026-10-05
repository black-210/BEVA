import ctypes as C, os, platform
from pathlib import Path

class _Vector(C.Structure):
    _fields_=[("n",C.c_size_t),("data",C.POINTER(C.c_double))]
class _Matrix(C.Structure):
    _fields_=[("rows",C.c_size_t),("cols",C.c_size_t),("data",C.POINTER(C.c_double))]
class _Lattice(C.Structure):
    _fields_=[("rows",C.c_size_t),("cols",C.c_size_t),("data",C.POINTER(C.c_int64))]
class _Stats(C.Structure):
    _fields_=[("dimension",C.c_size_t),("swaps",C.c_size_t),("size_reductions",C.c_size_t),("iterations",C.c_size_t),("converged",C.c_int)]
class _Options(C.Structure):
    _fields_=[("dimension",C.c_size_t),("delta",C.c_double),("eta",C.c_double),("max_iterations",C.c_size_t),("use_siegel",C.c_int),("track_transform",C.c_int)]

def _load():
    names=[]
    if os.getenv("BEVA_LIBRARY"): names.append(os.getenv("BEVA_LIBRARY"))
    if platform.system()=="Windows": names += ["beva.dll"]
    elif platform.system()=="Darwin": names += ["libbeva.dylib","libbeva.so"]
    else: names += ["libbeva.so","libbeva.so.1"]
    roots=[Path(__file__).resolve().parent,Path.cwd(),Path(__file__).resolve().parents[2]/"build"]
    for n in names:
        for r in roots:
            p=Path(n) if Path(n).is_absolute() else r/n
            if p.exists():
                return C.CDLL(str(p))
    try:return C.CDLL(names[0])
    except OSError as e: raise RuntimeError("BEVA shared library not found; build it with CMake and set BEVA_LIBRARY.") from e
_lib=_load()
# lifecycle
for fn in ["beva_vector_init","beva_vector_free","beva_matrix_init","beva_matrix_free","beva_lattice_init","beva_lattice_free"]:
    getattr(_lib,fn).restype=C.c_int if fn.endswith("_init") else None
_lib.beva_lattice_lll.argtypes=[C.POINTER(_Lattice),C.POINTER(_Options),C.POINTER(_Stats),C.c_void_p];_lib.beva_lattice_lll.restype=C.c_int
_lib.beva_lattice_gso.argtypes=[C.POINTER(_Lattice),C.POINTER(_Matrix),C.POINTER(_Matrix),C.POINTER(C.c_double)];_lib.beva_lattice_gso.restype=C.c_int
_lib.beva_matrix_det.argtypes=[C.POINTER(_Matrix)];_lib.beva_matrix_det.restype=C.c_double
_lib.beva_matrix_rank.argtypes=[C.POINTER(_Matrix),C.c_double];_lib.beva_matrix_rank.restype=C.c_size_t
_lib.beva_matrix_solve.argtypes=[C.POINTER(_Matrix),C.POINTER(_Vector),C.POINTER(_Vector)];_lib.beva_matrix_solve.restype=C.c_int
_lib.beva_matrix_inverse.argtypes=[C.POINTER(_Matrix),C.POINTER(_Matrix)];_lib.beva_matrix_inverse.restype=C.c_int
_lib.beva_number_theory_dummy = None

class Vector:
    def __init__(self, values):
        self._v=_Vector(); self._owner=[]
        vals=[float(x) for x in values]; self._v.n=len(vals)
        arr=(C.c_double*len(vals))(*vals) if vals else None; self._owner.append(arr); self._v.data=arr
    @property
    def values(self): return [self._v.data[i] for i in range(self._v.n)]
    def __len__(self): return self._v.n

class Matrix:
    def __init__(self, rows):
        vals=[list(map(float,r)) for r in rows]; cols=len(vals[0]) if vals else 0
        if any(len(r)!=cols for r in vals): raise ValueError("ragged matrix")
        self._m=_Matrix(len(vals),cols,None); arr=(C.c_double*(len(vals)*cols))(*[x for r in vals for x in r]) if vals and cols else None
        self._owner=arr; self._m.data=arr
    @property
    def shape(self): return (self._m.rows,self._m.cols)
    def to_list(self): return [[self._m.data[i*self._m.cols+j] for j in range(self._m.cols)] for i in range(self._m.rows)]

class Lattice:
    def __init__(self, rows):
        vals=[list(map(int,r)) for r in rows]; cols=len(vals[0]) if vals else 0
        if any(len(r)!=cols for r in vals): raise ValueError("ragged lattice")
        self._l=_Lattice(len(vals),cols,None); arr=(C.c_int64*(len(vals)*cols))(*[x for r in vals for x in r]) if vals and cols else None
        self._owner=arr;self._l.data=arr
    @property
    def shape(self): return (self._l.rows,self._l.cols)
    def to_list(self): return [[int(self._l.data[i*self._l.cols+j]) for j in range(self._l.cols)] for i in range(self._l.rows)]
    def lll(self, delta=.99, eta=.501, max_iterations=100000, siegel=False):
        o=_Options(self._l.rows,delta,eta,max_iterations,int(siegel),0);s=_Stats()
        rc=_lib.beva_lattice_lll(C.byref(self._l),C.byref(o),C.byref(s),None)
        if rc: raise RuntimeError(f"BEVA LLL error {rc}")
        return LLLStats(s),self
    def gso(self):
        m=self._l.rows; bs=_Matrix();mu=_Matrix(); arr=(C.c_double*m)(); rc=_lib.beva_lattice_gso(C.byref(self._l),C.byref(bs),C.byref(mu),arr)
        if rc: raise RuntimeError(f"BEVA GSO error {rc}")
        return ([[bs.data[i*self._l.cols+j] for j in range(self._l.cols)] for i in range(m)],
                [[mu.data[i*m+j] for j in range(m)] for i in range(m)],
                [arr[i] for i in range(m)])

class LLLStats:
    def __init__(self,s): self.dimension=s.dimension;self.swaps=s.swaps;self.size_reductions=s.size_reductions;self.iterations=s.iterations;self.converged=bool(s.converged)
    def __repr__(self): return f"LLLStats(iterations={self.iterations}, swaps={self.swaps}, reductions={self.size_reductions})"

def lll_reduce(rows, **kwargs): return Lattice(rows).lll(**kwargs)
def gram_schmidt(rows): return Lattice(rows).gso()
def determinant(rows):
    m=Matrix(rows);return _lib.beva_matrix_det(C.byref(m._m))
def rank(rows,tol=0.0):
    m=Matrix(rows);return int(_lib.beva_matrix_rank(C.byref(m._m),tol))
def solve(a,b):
    m=Matrix(a);v=Vector(b);out=_Vector();rc=_lib.beva_matrix_solve(C.byref(m._m),C.byref(v._v),C.byref(out))
    if rc: raise RuntimeError(f"solve error {rc}")
    return [out.data[i] for i in range(out.n)]
def inverse(a):
    m=Matrix(a);out=_Matrix();rc=_lib.beva_matrix_inverse(C.byref(m._m),C.byref(out))
    if rc: raise RuntimeError(f"inverse error {rc}")
    return [[out.data[i*out.cols+j] for j in range(out.cols)] for i in range(out.rows)]
def gcd(a,b):
    # portable Python-side convenience
    while b:a,b=b,a%b
    return abs(a)
def mod_inverse(a,m):
    a%=m
    for x in range(1,m):
        if (a*x)%m==1:return x
    raise ValueError("inverse does not exist")
def cholesky(a):
    import math
    n=len(a);L=[[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            s=a[i][j]-sum(L[i][k]*L[j][k] for k in range(j))
            if i==j:
                if s<=0: raise ValueError("matrix is not positive definite")
                L[i][j]=math.sqrt(s)
            else:L[i][j]=s/L[j][j]
    return L
def babai(rows,target,delta=.99):
    lat=Lattice(rows); lat.lll(delta=delta)
    gs,mu,n2=lat.gso();res=[float(x) for x in target];m=len(rows);coeff=[0]*m
    for ii in range(m):
        i=m-1-ii; c=round(sum(res[j]*gs[i][j] for j in range(len(target)))/n2[i]);coeff[i]=c
        for j in range(len(target)):res[j]-=c*rows[i][j]
    return coeff,res
