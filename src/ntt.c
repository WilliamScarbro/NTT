#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ntt.h"

/*
* Computes NTT using the FFT Algorithm
* 
*/
void fntt(long* vec, long* dest, long* temp, long N, long P, long root, long stride, long start, long depth){
	//printf("stride: %ld, start: %ld, depth: %ld\n",stride,start,depth);
	//base case
	if (depth==0){
		fntt_help(vec,dest,N,P,root,stride,start);
		return;
	}
	//recursive calls (notice temp and dest are reversed)
	fntt(vec,temp,dest,N,P,root,stride*2,start,depth-1);
	fntt(vec,temp,dest,N,P,root,stride*2,start+stride,depth-1);
	//butterfly
	long i,r2i;
	long ti;
	long di;
	for (i=0; i<N/2; i+=stride){
		r2i=power_mod(root,i,P);
		ti=2*i+start;
		di=i+start;
		dest[di]=(temp[ti]+r2i*temp[ti+stride])%P;
		dest[di+N/2]=(temp[ti]-r2i*temp[ti+stride])%P;
	}
	#ifdef DEBUG
	printf("fntt: Intermediate Result: \n");
	printPoly(dest,N);
	#endif
}

/*
* NTT base case for fntt
* The effective length of the polynomial is N/stride
* root must be adjusted to correct for shorter N
* vec and dest are only accessed on stride indexes (to represent shorter vec)
*/
void fntt_help(long* vec, long* dest, long N, long P, long root, long stride, long start){
	long i,j,sum;
	#ifdef DEBUG
	printf("fntt_help\nInput:\n  ");
	for (i=0; i*stride+start<N; i++)
		printf("%ld*x%ld + ",vec[start+i*stride],i);
	#endif
	
	//adjust root for smaller N
	root=power_mod(root,stride,P);
	for (i=0; i<N/stride; i++){
		sum=0;
		for (j=0; j*stride+start<N; j++){
			sum+=vec[start+j*stride]*power_mod(root,i*j,P);
			sum%=P;
		}
		dest[start+i*stride]=sum;
	}
	#ifdef DEBUG
	printf("\nOutput\n  ");
	for (i=0; i*stride+start<N; i++)
		printf("%ld*x%ld + ",dest[start+i*stride],i);
	printf("\n");
	#endif
}

/*
* Computes forward NTT
* in Z[X]/p
*    /(X^N-1)
* Flags adjust parallelism
*/
void ntt(long* vec, long* dest, long N, long P, long root){
	long i,j,sum; 
    for (i=0; i<N; i++){
        sum=0;
        for (j=0; j<N; j++){
            sum+=vec[j] * power_mod(root,i*j,P);
            sum%=P;
        }
        dest[i]=sum;
    }
}

/*
* Computes reverse NTT
*/
void inv_ntt(long* vec, long* dest, long N, long P, long root){
    ntt(vec,dest,N,P,inverse(root,P));
    long inv_N = inverse(N,P);
    for (int i=0; i<N; i++)
        dest[i]=(dest[i]*inv_N)%P;
}

/*
* Multiplies two elements in Z[X]/P
*                            /(X^N-1)
* root has order N in F_P
*/
void convolution(long* vec1, long* vec2, long* dest, long* temp, long N, long P, long root){
    ntt(vec1,dest,N,P,root);
    ntt(vec2,temp,N,P,root);
	//printf("F(v)\n");
	//printPoly(dest,N);
	//printf("F(v2)\n");
	//printPoly(temp,N);
    long i;
    for (i=0; i<N; i++){
        temp[i]=(temp[i]*dest[i])%P;
    }
	//printf("F(v1)*F(v2)\n");
	//printPoly(temp,N);
    inv_ntt(temp,dest,N,P,root);
}

// Helper Methods for NTT

/*
* Slow convolution algorithm to check results of NTT
*/
void check_conv(long* vec1, long* vec2, long* dest, long N, long P){
	long i,j,di;
	for (i=0; i<N; i++)
		dest[i]=0;
	for (i=0; i<N; i++){
		for (j=0; j<N; j++){
			di=(i+j)%N;
			dest[di]+=(vec1[i]*vec2[j]);
			dest[di]%=P;
		}
	}
}
void printPoly(long* vec, long N){
    printf("	%ldX^%d ",vec[0],0);
    for (long i=1; i<N; i++){
        printf("+ %ldX^%ld",vec[i],i);
    }
    printf("\n");
}

/*
* checks if g is generator for F_P
* O(P*log(P))
*/
int is_generator(long g,long P){
	for (long i=1; i<=P/2; i++){
		if (i|(P-1) && power_mod(g,i,P)==1)
			return 0;
	}
	return 1;
}
/*
* Finds generator for F_P
* (P must be prime)
* O(P^2*log(P))
*/
long generator(long P){
	for (long i=1; i<P; i++){
		if (is_generator(i,P)){
			printf("%ld is generator for %ld\n",i,P);
			return i;
		}
	}
	printf("Unable to find generator for %ld\n",P);
	return -1;
}

/*
* finds element of F_P which has order N
* this produces the w for NTT (Equivalent to N-th root of 1 in DFT)
*/
long get_root(long N, long P){
	long gen = generator(P);
	long k = (P-1)/N;
	return power_mod(gen,k,P);
}

//implements fast power algorithm O(log(exp))
long power_mod(long base, long exp, long mod){
    long res=1;
    long cur=base;
    while (exp>0){
        if (exp%2==1)
            res=(res*cur)%mod;
        cur=(cur*cur)%mod;
        exp/=2;
    }
    return res;
}
// returns positive lift in mod
long true_mod(long val, long mod){
    long ret=val%mod;
    if (ret<0)
        return ret+mod;
    return ret;
}

/*
 Finds inverse of val in F_mod, i.e. i : i*val=1 (mod)
 Takes O(log(mod)) time,
 Uses extended euclidian algorithm (although it does not look like it)
*/
long inverse(long val, long mod){
	long x,y,a,b,temp;
	x=mod;
	y=val;
	a=0;
	b=1;
	while (y!=0){
	  	temp=b;
	 	b=a-(x/y)*b;
	 	a=temp;
	  	temp=y;
	  	y=x%y;
	  	x=temp;
	}
	if (x==1)
		return true_mod(a,mod);
	else
		printf("Unable to find inverse for %ld in F%ld\n",val,mod);
}

int is_prime(long P){
	if (P==1)
		return 0;
	for (long i=2; i*i<=P; i++){
		if (P%i==0)
			return 0;
	}
	return 1;
}

/* 
* finds a usable mod for NTT algorithm
* i.e. P : P>midMod && N|P-1
*/
long usable_mod(long N, long minMod){
	long cur=(minMod/N)*N;
	while (!is_prime(cur+1))
		cur+=N;
	return cur+1;
}

//Some test cases
void test_ntt(){
	long vec[3]={2,3,3};
	long trans[3]={1,6,6};
	long dest[3];
	ntt(vec,dest,3,7,2);
	printf("Forward transform result\n");
	printPoly(dest,3);
	printf("Correct\n");
	printPoly(trans,3);
	printf("Inverse transform result\n");
	inv_ntt(trans,dest,3,7,2);
	printPoly(dest,3);
	printf("Correct\n");
	printPoly(vec,3);
}
void test_IsPrime(){
	int primes[16]={1,1,0,1,0,1,0,0,0,1,0,1,0,0,0,1};
	for (int i=2; i<=17; i++){
		if (is_prime(i)!=primes[i-2])
			printf("is prime failed for %d\n",i);
	}
}
void test_inverse(){
    long mods[5]={3,5,7,11,13};
    long val=2;
    long invs[5]={2,3,4,6,7};
    for (int i=0; i<5; i++){
        long res=inverse(val,mods[i]);
        if (true_mod(res,mods[i])!=invs[i])
            printf("Inverse Error: 2*%ld!=1 mod %ld\n",true_mod(res,mods[i]),mods[i]);
    }
}
long power(long l1, long l2){
	long res=1;
	while(l2-->0)
		res*=l1;
	return res;
}
void test_power_mod(){
    long mods[5]={3,5,7,11,13};
	for (int i=0; i<5; i++){
		for (int j=1; j<6; j++){
			long res = power_mod(2,j,mods[i]);
			long cor = power(2,j)%mods[i];
			if (res!=cor)
				printf("2^%ld!=%ld mod %ld\n",res,cor,mods[i]);
		}
	}
}

/*Requirements
*  kN+1=P (k in Z)
*  
*/
int main(int argc, char** argv){
	//test_ntt();
	test_IsPrime();
	test_inverse();
	test_power_mod();
	printf("Tests Completed\n");
	if (argc<3){
        printf("Recieved %d args. Usage: ntt <N> <P>\n",argc);
        exit(0);
    }
    long N,P,root,i;
    N = atoi(argv[1]);
    P = atoi(argv[2]);
    P = usable_mod(N,P);
	root = get_root(N,P);
	printf("root: %ld\n",root);

    long* vec1 = (long*)malloc(sizeof(long)*N);
    long* vec2= (long*)malloc(sizeof(long)*N);
    long* res = (long*)malloc(sizeof(long)*N);
    long* temp = (long*)malloc(sizeof(long)*N);

    srand(time(NULL));
    for (i=0; i<N; i++){
        vec1[i] = rand()%P;
        vec2[i] = rand()%P;
    }

    #ifdef CONV_TEST
	printf("Polynomial 1\n");
    printPoly(vec1,N);
    printf("Polynomial 2\n");
    printPoly(vec2,N);

    convolution(vec1,vec2,res,temp,N,P,root);

    printf("NTT Result\n");
    printPoly(res,N);
	
	check_conv(vec1,vec2,temp,N,P);

	printf("Correct Result\n");
	printPoly(temp,N);
	#endif
	
	#ifdef FNTT_TEST
	
	fntt(vec1,res,temp,N,P,root,1,0,3);

	ntt(vec1,temp,N,P,root);

	printf("Polynomial\n");
	printPoly(vec1,N);

	printf("FNTT Result\n");
    printPoly(res,N);
	
	printf("Correct Result\n");
	printPoly(temp,N);
	#endif
	
	//check results
	int success=1;
	for (i=0; i<N; i++){
		if (true_mod(res[i],P)!=true_mod(temp[i],P))
			success=0;
	}
	if (success)
		printf("Test Successful\n");
	else
		printf("Test is Failure.\n");
	
	free(vec1);
	free(vec2);
	free(res);
	free(temp);
}
