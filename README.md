# RSA_GMPLIB
 Maximum support for encryption length 2048 bit

First, you should install gmp library in you system and then complie and run this rsa program.
Please follow gmplib.org step to install gmp lib: https://gmplib.org/manual/Introduction-to-GMP.html#Introduction-to-GMP 

complie command:
$ gcc -g miller_rabin_primality_test.c random_number.c rsa_algorithm.c select_parameter.c main.c -lgmp -o rsa

run command:
$ ./rsa
