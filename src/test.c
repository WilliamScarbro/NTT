#include <stdio.h>
#include "ntt.h"

long true_mod(long val, long mod){
	long ret=val%mod;
	if (ret<0)
		return ret+mod;
	return ret;
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
void test_power_mod(){
	long res=power_mod(3,5,4);
	if (res!=1)
		printf("Power_mod error: 3^5!=%ld mod 4\n",res);
}

int main(){
	test_inverse();
	test_power_mod();
}


