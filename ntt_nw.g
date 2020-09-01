Class(TField, AtomicTyp, rec(
    hash := (val,size) -> 1+(10047871*val mod size),
    check := v -> Cond(IsInt(v), v, Error("<v> must be an int")),
    realType := self >> TInt,
    complexType := self >> self,
    zero := self >> self.value(0),
    one := self >> self.value(1),
));

Class(NWLin, DiagFunc, rec(
    checkParams := (self, params) >> Checked(Length(params)=3,
    IsPosInt0Sym(params[1]), IsType(params[3]), params),
    lambda := self >> let(j:=Ind(self.params[1]), d := self.params[2],
        Lambda(j, 1-d+2*j).setRange(self.params[3])),
    range := self >> self.params[3],
    domain := self >> self.params[1]
));

Class(T_Field, T_Type, rec(
    hash:=TField.hash,
    realType := self >> self.params[1],
    complexType := self >> self,
#    product := self >> (v1, v2) -> (v1*v2) mod self.params[2],
#    sum := self >> (v1,v2) -> (v1+v2) mod self.params[2],
));

Class(modOmega, AutoFoldExp, rec(
	ev:= self >> (self.args[1].ev() ^ self.args[2].ev()) mod self.args[3],
	computeType := self >> TInt));


# Replaces dOmega with root in finite field
# defines f: TInt -> TInt : i -> w_n^(k*i) mod p
# params( <w_n> , <k> , <p>)
Class(MyOmega, DiagFunc, rec(
	checkParams := (self,params) >> Checked(Length(params)=3, params),
	lambda := self >> let(i:=Ind(), Lambda(i,modOmega(self.params[1],self.params[2]*i,self.params[3]))),
	range := self >> TInt,
	domain := self >> TInt
));

# T_i,j = r^((1-d)i+2ij)
ModTw1 := (n,d,k,m,r) -> Checked(
	IsPosIntSym(n), IsPosIntSym(d), IsIntSym(k), IsIntSym(m), IsIntSym(r),
	fCompose(MyOmega(r,k,m),diagTensor(dLin(div(n,d),1,0,TInt),NWLin(d,d,TInt))));

# NTT_NW(<size>, <field>, <root>)
# 2*size | field-1
# root has order 2*size in field
# size is a power of 2
Class(NTT_NW, TaggedNonTerminal, rec(
	abbrevs := [ (n,p,psi) -> Checked(IsInt(n), n>0, IsInt(p), p>0, [n,p,psi]) ],

	dims := self >> [ self.params[1],self.params[1] ],
	
	terminate := self >> let(N := self.params[1], P:=self.params[2], Psi:=self.params[3],
		Mat(List([0..N-1], r -> List([0..N-1], c -> (Psi^(c+2*r*c) mod P))))),

	isReal:=self >> true,
    SmallRandom := () -> Random([2..16]),
    LargeRandom := () -> 2^Random([6..15]),
    TType := T_Field(TUnknown)
));


# Note that NTT_NW_CT applies the Cooley-Tukey algorithm for radix-m
NewRulesFor(NTT_NW,rec(
    NTT_NW_CT := rec(
	    maxSize:=false,
	    forcePrimeFactor := false,
	    applicable := (self,nt) >> nt.params[1]>2,
	    children := nt -> Map2(DivisorPairs(nt.params[1]),
		    (m,n) -> [ NTT_NW(m,nt.params[2],nt.params[3]^n mod nt.params[2]),
			    		NTT_NW(n,nt.params[2],nt.params[3]^m mod nt.params[2]) ]),
	    apply := (nt,C,cnt) -> let(mn := nt.params[1], p:=nt.params[2], r:=nt.params[3], m:= Rows(C[1]), n:=Rows(C[2]),
		    Tensor(C[1],I(n)) * 
		    Diag(fPrecompute(ModTw1(mn,n,1,p,r)))*
		    Tensor(I(m),C[2])*
		    L(mn,m))
    ),

	NTT_Base := rec(
		forTransposition := false,
		applicable := nt->nt.params[1]=2,
		apply := (nt,C,cnt) -> F(2)
	)
));

#	Transforms
#       Notice the roots (third argument) now have order 2N
#t2:=NTT(2,13,8);
#t4:=NTT(4,17,9);
t8:=NTT_NW(8,17,3);

# n=8 transform matrix and negative wrapped twiddle matrix
#  t8.terminate();
#  MatSPL(Diag(fPrecompute(ModTw1(8,4,1,17,3))));

#	Initialize
  opts:=CopyFields(SpiralDefaults,rec(breakdownRules:=rec(NTT_NW:=[NTT_Base, NTT_NW_CT])));
  rt:=RandomRuleTree(t8,opts);

#	Intermediate Steps
#SPL formula
  spl:=SPLRuleTree(rt);
  ss1:=spl.sums();
#Sigma-SPL
  ss:=SumsRuleTree(rt,opts);
#Sigma-SPL -> code
  c1:=CodeSums(ss,opts);

#	End Result
  c:=CodeRuleTree(rt,opts);
  PrintCode("ntt"::StringInt(8),c,opts);
  PrintTo("GC.c",PrintCode("ntt"::StringInt(8),c,opts));


