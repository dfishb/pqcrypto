#!/bin/sh

gcc  -c gfp2.c -o gfp2.o -lm -w -lpthread -lrt -marm -lgmp -lcrypto 

gcc -c arm-768.s -o arm-768.o

gcc arm-768.o gfp2.o -o gfp2 -lm -w -lpthread -lrt -marm -lgmp -lcrypto


# francisco flags not working -flto -mfloat-abi=softfp

# -ftree-vectorizer-verbose=1 -march=armv7-a -Ofast --save-temps -O3  -fprefetch-loop-arrays -funroll-all-loops -ftree-vectorize -mthumb-interwork -mvectorize-with-neon-quad -mtune=cortex-a15 -lgmp -lcrypto -mfpu=neon -Wa,-mimplicit-it=thumb
#-L/usr/lib/arm-linux-gnueabihf#-flto -ftree-vectorizer-verbose=1 -march=armv7-a -Ofast --save-temps -O3 -marm -fprefetch-loop-arrays -funroll-all-loops -ftree-vectorize -mthumb-interwork -mvectorize-with-neon-quad -mtune=cortex-a15 -lgmp -lcrypto -mfpu=neon -Wa,-mimplicit-it=thumb  -lm -w -lpthread -lrt