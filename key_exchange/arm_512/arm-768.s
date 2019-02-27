
.align	2
.section .text

PRIME_OLD:
    .word  0xffffffff @first 12 limbs (12 least significant limbs) in the prime are 0xffffffff
    .word  0x5ce89647
    .word  0xf007669a
    .word  0x484504f9
    .word  0xade00d91
    .word  0xc24486e3
    .word  0x0979d570
    .word  0x9a4c7025
    .word  0x8bbae367
    .word  0x9f6808b4
    .word  0xa06a805a
    .word  0x87fabdfa
    .word  0xe69ebefa
    .word  0x5


PRIME:
    .word 0xffffffff @first 8 limbs (8 least significant limbs) in the prime are 0xffffffff
    .word 0x14a0a4b7
    .word 0xa8a062c9
    .word 0x21dad274
    .word 0xdbbaf789
    .word 0x5bb8d489
    .word 0x45957692
    .word 0x4ef2a188
    .word 0x213fbe56
    .word 0x1a5


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ 25 word modular addition
@ Usage:
@    void addm(int *res, int *a, int *b)
@ Operation:
@    res = b + a mod p
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
.align	2
.global	addm
.type	addm, %function
addm:


    push {r4-r12}

    LDM r1!, {r3-r7}
    LDM r2!, {r8-r12}
    ADDS r3,r8
    ADCS r4,r9
    ADCS r5,r10
    ADCS r6,r11
    ADCS r7,r12
    STM r0!, {r3-r7}

    LDM r1!, {r3-r7}
    LDM r2!, {r8-r12}
    ADCS r3,r8
    ADCS r4,r9
    ADCS r5,r10
    ADCS r6,r11
    ADCS r7,r12
    STM r0!, {r3-r7}

    LDM r1!, {r3-r7}
    LDM r2!, {r8-r12}
    ADCS r3,r8
    ADCS r4,r9
    ADCS r5,r10
    ADCS r6,r11
    ADCS r7,r12
    STM r0!, {r3-r7}

    LDM r1!, {r3-r4}
    LDM r2!, {r8-r9}
    ADCS r3,r8
    ADCS r4,r9
    STM r0!, {r3-r4}


    @d//decrementing r0 to the beginning, then subtracting the prime and palcing the result in r11
    STR r3, [r0], #-68
    STR r3, [r13], #-68

    LDM r0!, {r4-r12}

    LDR r3, PRIME
    SUBS r4, r3
    SBCS r5, r3
    SBCS r6, r3
    SBCS r7, r3
    SBCS r8, r3
    SBCS r9, r3
    SBCS r10, r3
    SBCS r11, r3

    LDR r3, PRIME+4
    SBCS r12, r3

    STM r13!, {r4-r12} 


    LDM r0!, {r4-r11}

    LDR r3, PRIME+8
    SBCS r4, r3

    LDR r3, PRIME+12
    SBCS r5, r3

    LDR r3, PRIME+16
    SBCS r6, r3

    LDR r3, PRIME+20
    SBCS r7, r3

    LDR r3, PRIME+24
    SBCS r8, r3

    LDR r3, PRIME+28
    SBCS r9, r3

    LDR r3, PRIME+32
    SBCS r10, r3

    LDR r3, PRIME+36
    SBCS r11, r3

    STM r13!, {r4-r11}



    @//if there was not a carry, move the difference into the output register

    STRCS r3, [r0], #-68
    STRCS r3, [r13], #-68

    LDMCS r13!, {r3-r12}
    STMCS r0!, {r3-r12}

    LDMCS r13!, {r3-r9}
    STMCS r0!, {r3-r9}


    pop {r4-r12}
    bx lr


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ 25 word modular addition
@ Usage:
@    void addm_old(int *res, int *a, int *b)
@ Operation:
@    res = b + a mod p
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
.align	2
.global	addm_old
.type	addm_old, %function
addm_old:


push {r4-r12}

LDM r1!, {r3-r7}
LDM r2!, {r8-r12}
ADDS r3,r8
ADCS r4,r9
ADCS r5,r10
ADCS r6,r11
ADCS r7,r12
STM r0!, {r3-r7}

LDM r1!, {r3-r7}
LDM r2!, {r8-r12}
ADCS r3,r8
ADCS r4,r9
ADCS r5,r10
ADCS r6,r11
ADCS r7,r12
STM r0!, {r3-r7}

LDM r1!, {r3-r7}
LDM r2!, {r8-r12}
ADCS r3,r8
ADCS r4,r9
ADCS r5,r10
ADCS r6,r11
ADCS r7,r12
STM r0!, {r3-r7}

LDM r1!, {r3-r7}
LDM r2!, {r8-r12}
ADCS r3,r8
ADCS r4,r9
ADCS r5,r10
ADCS r6,r11
ADCS r7,r12
STM r0!, {r3-r7}

LDM r1!, {r3-r7}
LDM r2!, {r8-r12}
ADCS r3,r8
ADCS r4,r9
ADCS r5,r10
ADCS r6,r11
ADC r7,r12
STM r0!, {r3-r7}

@d//decrementing r0 to the beginning, then subtracting the prime and palcing the result in r11
STR r3, [r0], #-100
STR r3, [r13], #-100

LDM r0!, {r4-r12}
LDR r3, PRIME
SUBS r4, r3
SBCS r5, r3
SBCS r6, r3
SBCS r7, r3
SBCS r8, r3
SBCS r9, r3
SBCS r10, r3
SBCS r11, r3
SBCS r12, r3
STM r13!, {r4-r12}

LDM r0!, {r4-r12}
SBCS r4, r3
SBCS r5, r3
SBCS r6, r3


LDR r3, PRIME+4
SBCS r7, r3

LDR r3, PRIME+8
SBCS r8, r3

LDR r3, PRIME+12
SBCS r9, r3

LDR r3, PRIME+16
SBCS r10, r3

LDR r3, PRIME+20
SBCS r11, r3

LDR r3, PRIME+24
SBCS r12, r3

STM r13!, {r4-r12}
LDM r0!, {r4-r10}

LDR r3, PRIME+28
SBCS r4, r3

LDR r3, PRIME+32
SBCS r5, r3

LDR r3, PRIME+36
SBCS r6, r3

LDR r3, PRIME+40
SBCS r7, r3

LDR r3, PRIME+44
SBCS r8, r3

LDR r3, PRIME+48
SBCS r9, r3

LDR r3, PRIME+52
SBCS r10, r3

STM r13!, {r4-r10}


@//if there was not a carry, move the difference into the output register

STRCS r3, [r0], #-100
STRCS r3, [r13], #-100

LDMCS r13!, {r3-r12}
STMCS r0!, {r3-r12}

LDMCS r13!, {r3-r12}
STMCS r0!, {r3-r12}

LDMCS r13!, {r3-r7}
STMCS r0!, {r3-r7}

pop {r4-r12}
bx lr







