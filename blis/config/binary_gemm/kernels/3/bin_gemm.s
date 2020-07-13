# mark_description "Intel(R) C Intel(R) 64 Compiler XE for applications running on Intel(R) 64, Version 14.0.3.174 Build 2014042";
# mark_description "2";
# mark_description "-I/root/blis/include/blis/ -S -msse4.2 -std=c99";
	.file "bin_gemm.c"
	.text
..TXTST0:
# -- Begin  bin_gemm_2x2
# mark_begin;
       .align    16,0x90
	.globl bin_gemm_2x2
bin_gemm_2x2:
# parameter 1: %rdi
# parameter 2: %rsi
# parameter 3: %rdx
# parameter 4: %rcx
# parameter 5: %r8
# parameter 6: %r9
# parameter 7: 176 + %rsp
# parameter 8: 184 + %rsp
# parameter 9: 192 + %rsp
..B1.1:                         # Preds ..B1.0
..___tag_value_bin_gemm_2x2.1:                                  #13.1
        pushq     %r13                                          #13.1
..___tag_value_bin_gemm_2x2.3:                                  #
        pushq     %r14                                          #13.1
..___tag_value_bin_gemm_2x2.5:                                  #
        pushq     %r15                                          #13.1
..___tag_value_bin_gemm_2x2.7:                                  #
        pushq     %rbx                                          #13.1
..___tag_value_bin_gemm_2x2.9:                                  #
        pushq     %rbp                                          #13.1
..___tag_value_bin_gemm_2x2.11:                                 #
        subq      $128, %rsp                                    #13.1
..___tag_value_bin_gemm_2x2.13:                                 #
        movq      %rdi, %rax                                    #14.25
        shrq      $63, %rax                                     #14.25
        xorl      %r8d, %r8d                                    #24.2
        addq      %rdi, %rax                                    #14.25
        xorl      %ebx, %ebx                                    #24.11
        sarq      $1, %rax                                      #14.25
        movq      %rsi, %r10                                    #13.1
        xorl      %esi, %esi                                    #25.2
        xorl      %ebp, %ebp                                    #25.11
        lea       (%rax,%rax), %r11                             #15.23
        subq      %r11, %rdi                                    #15.19
        xorl      %r11d, %r11d                                  #27.13
        testq     %rax, %rax                                    #27.23
        je        ..B1.5        # Prob 10%                      #27.23
                                # LOE rax rdx rcx rbx rbp rsi rdi r8 r9 r10 r11 r12 r13 r14 r15
..B1.2:                         # Preds ..B1.1
        movq      %rbp, 48(%rsp)                                #
        movq      %rsi, 56(%rsp)                                #
        movq      %rbx, 40(%rsp)                                #
        movq      %r8, 64(%rsp)                                 #
        movq      %rdi, (%rsp)                                  #
        movq      %rax, 32(%rsp)                                #
        movq      %r10, 8(%rsp)                                 #
        movq      %r9, 16(%rsp)                                 #
        movq      %r12, 24(%rsp)                                #
..___tag_value_bin_gemm_2x2.14:                                 #
                                # LOE rdx rcx r11
..B1.3:                         # Preds ..B1.3 ..B1.2
        movq      40(%rdx), %rdi                                #62.11
        movq      48(%rdx), %rsi                                #69.11
        movq      %rsi, %r12                                    #73.11
        movq      %rdi, 120(%rsp)                               #62.11
        movq      48(%rcx), %rdi                                #70.11
        andq      %rdi, %r12                                    #73.11
        movq      %r11, 80(%rsp)                                #
        movq      16(%rdx), %rbx                                #41.11
        movq      (%rcx), %r11                                  #29.10
        movq      16(%rcx), %r15                                #42.11
        popcnt    %r12, %r14                                    #73.11
        movq      %rbx, %r12                                    #45.11
        movq      (%rdx), %rax                                  #28.10
        andq      %r15, %r12                                    #45.11
        movq      %r11, 88(%rsp)                                #29.10
        movq      %r15, 104(%rsp)                               #42.11
        movq      %rax, %r15                                    #32.11
        andq      88(%rsp), %r15                                #32.11
        movq      32(%rdx), %rbp                                #55.11
        movq      %rbp, %r13                                    #59.11
        popcnt    %r15, %r15                                    #32.11
        movq      8(%rdx), %r10                                 #35.11
        movq      24(%rdx), %r9                                 #48.11
        movq      32(%rcx), %r8                                 #56.11
        andq      %r8, %r13                                     #59.11
        popcnt    %r12, %r12                                    #45.11
        addq      64(%rsp), %r15                                #32.4
        movq      %r10, 96(%rsp)                                #35.11
        addq      %r15, %r12                                    #45.4
        movq      %r9, 112(%rsp)                                #48.11
        movq      %rcx, 72(%rsp)                                #
        movq      8(%rcx), %r11                                 #30.11
        andq      %r11, %rax                                    #33.11
        movq      24(%rcx), %r10                                #43.11
        andq      %r10, %rbx                                    #46.11
        movq      40(%rcx), %r9                                 #57.11
        andq      %r9, %rbp                                     #60.11
        movq      56(%rcx), %rcx                                #71.11
        andq      %rcx, %rsi                                    #74.11
        popcnt    %r13, %r13                                    #59.11
        addq      %r12, %r13                                    #59.4
        popcnt    %rsi, %r12                                    #74.11
        popcnt    %rbp, %rsi                                    #60.11
        popcnt    %rbx, %rbp                                    #46.11
        popcnt    %rax, %rbx                                    #33.11
        addq      %r13, %r14                                    #73.4
        addq      40(%rsp), %rbx                                #33.4
        addq      %rbx, %rbp                                    #46.4
        addq      %rbp, %rsi                                    #60.4
        movq      56(%rdx), %rax                                #76.11
        andq      %rax, %rcx                                    #79.11
        movq      %r14, 64(%rsp)                                #73.4
        andq      %rax, %rdi                                    #78.11
        movq      120(%rsp), %r14                               #64.11
        andq      %r14, %r8                                     #64.11
        movq      88(%rsp), %rbx                                #38.11
        lea       (%rsi,%r12), %r13                             #74.4
        movq      96(%rsp), %r12                                #38.11
        andq      %r12, %r11                                    #39.11
        movq      %r13, 40(%rsp)                                #74.4
        andq      %r12, %rbx                                    #38.11
        movq      112(%rsp), %r13                               #50.11
        andq      %r13, %r10                                    #51.11
        popcnt    %rcx, %rax                                    #79.11
        andq      %r14, %r9                                     #65.11
        popcnt    %r11, %rcx                                    #39.11
        addq      $64, %rdx                                     #81.4
        popcnt    %rdi, %r15                                    #78.11
        popcnt    %r8, %rdi                                     #64.11
        movq      104(%rsp), %r8                                #50.11
        popcnt    %rbx, %rbp                                    #38.11
        andq      %r13, %r8                                     #50.11
        addq      48(%rsp), %rcx                                #39.4
        popcnt    %r10, %r10                                    #51.11
        addq      56(%rsp), %rbp                                #38.4
        addq      %rcx, %r10                                    #51.4
        popcnt    %r8, %rsi                                     #50.11
        popcnt    %r9, %r9                                      #65.11
        addq      %rbp, %rsi                                    #50.4
        addq      %r10, %r9                                     #65.4
        addq      %rsi, %rdi                                    #64.4
        movq      80(%rsp), %r11                                #27.33
        incq      %r11                                          #27.33
        lea       (%r9,%rax), %rcx                              #79.4
        movq      %rcx, 48(%rsp)                                #79.4
        movq      72(%rsp), %rcx                                #82.4
        lea       (%rdi,%r15), %r8                              #78.4
        addq      $64, %rcx                                     #82.4
        movq      %r8, 56(%rsp)                                 #78.4
        cmpq      32(%rsp), %r11                                #27.23
        jne       ..B1.3        # Prob 82%                      #27.23
                                # LOE rdx rcx r11
..B1.4:                         # Preds ..B1.3
        movq      48(%rsp), %rbp                                #
        movq      56(%rsp), %rsi                                #
        movq      40(%rsp), %rbx                                #
        movq      64(%rsp), %r8                                 #
        movq      (%rsp), %rdi                                  #
        movq      8(%rsp), %r10                                 #
        movq      16(%rsp), %r9                                 #
        movq      24(%rsp), %r12                                #
..___tag_value_bin_gemm_2x2.15:                                 #
                                # LOE rdx rcx rbx rbp rsi rdi r8 r9 r10 r12 r13 r14 r15
..B1.5:                         # Preds ..B1.4 ..B1.1
        xorl      %eax, %eax                                    #86.13
        testq     %rdi, %rdi                                    #86.23
        je        ..B1.9        # Prob 10%                      #86.23
        .align    16,0x90
                                # LOE rax rdx rcx rbx rbp rsi rdi r8 r9 r10 r12 r13 r14 r15
..B1.7:                         # Preds ..B1.5 ..B1.7
        movq      (%rdx), %r13                                  #87.10
        movq      %r13, %r15                                    #91.11
        movq      8(%rcx), %r11                                 #89.11
        andq      %r11, %r13                                    #92.11
        popcnt    %r13, %r13                                    #92.11
        incq      %rax                                          #86.33
        movq      (%rcx), %r14                                  #88.10
        addq      %r13, %rbx                                    #92.4
        movq      8(%rdx), %r13                                 #94.11
        andq      %r14, %r15                                    #91.11
        andq      %r13, %r14                                    #97.11
        andq      %r13, %r11                                    #98.11
        popcnt    %r15, %r15                                    #91.11
        addq      $8, %rdx                                      #100.4
        popcnt    %r14, %r14                                    #97.11
        addq      %r15, %r8                                     #91.4
        popcnt    %r11, %r11                                    #98.11
        addq      %r14, %rsi                                    #97.4
        addq      %r11, %rbp                                    #98.4
        addq      $8, %rcx                                      #101.4
        cmpq      %rdi, %rax                                    #86.23
        jne       ..B1.7        # Prob 82%                      #86.23
                                # LOE rax rdx rcx rbx rbp rsi rdi r8 r9 r10 r12
..B1.9:                         # Preds ..B1.7 ..B1.5
        movsd     (%r10), %xmm7                                 #109.13
        movq      %r8, %rdx                                     #109.22
        movq      %r8, %r10                                     #109.22
        andq      $1, %rdx                                      #109.22
        shrq      $1, %r10                                      #109.22
        xorps     %xmm1, %xmm1                                  #109.22
        orq       %r10, %rdx                                    #109.22
        xorps     %xmm0, %xmm0                                  #109.22
        cvtsi2sdq %rdx, %xmm1                                   #109.22
        cvtsi2sdq %r8, %xmm0                                    #109.22
        addsd     %xmm1, %xmm1                                  #109.22
        movq      %rbx, %r11                                    #110.31
        movq      %rbx, %rdi                                    #110.31
        shrq      $1, %rdi                                      #110.31
        andq      $1, %r11                                      #110.31
        orq       %rdi, %r11                                    #110.31
        xorps     %xmm3, %xmm3                                  #110.31
        testq     %r8, %r8                                      #109.22
        xorps     %xmm2, %xmm2                                  #110.31
        cvtsi2sdq %r11, %xmm3                                   #110.31
        cvtsi2sdq %rbx, %xmm2                                   #110.31
        addsd     %xmm3, %xmm3                                  #110.31
        jl        ..L16         # Prob 50%                      #109.22
        movaps    %xmm0, %xmm1                                  #109.22
..L16:                                                          #
        testq     %rbx, %rbx                                    #110.31
        movq      %rsi, %rdx                                    #111.31
        movq      %rsi, %rbx                                    #111.31
        xorps     %xmm5, %xmm5                                  #111.31
        mulsd     %xmm7, %xmm1                                  #109.22
        xorps     %xmm4, %xmm4                                  #111.31
        cvtsi2sdq %rsi, %xmm4                                   #111.31
        addsd     (%r9), %xmm1                                  #109.3
        jl        ..L17         # Prob 50%                      #110.31
        movaps    %xmm2, %xmm3                                  #110.31
..L17:                                                          #
        shrq      $1, %rbx                                      #111.31
        andq      $1, %rdx                                      #111.31
        orq       %rbx, %rdx                                    #111.31
        testq     %rsi, %rsi                                    #111.31
        cvtsi2sdq %rdx, %xmm5                                   #111.31
        mulsd     %xmm7, %xmm3                                  #110.31
        addsd     %xmm5, %xmm5                                  #111.31
        movq      176(%rsp), %rcx                               #13.1
        movq      %rbp, %rsi                                    #112.38
        movsd     %xmm1, (%r9)                                  #109.3
        xorps     %xmm8, %xmm8                                  #112.38
        jl        ..L18         # Prob 50%                      #111.31
        movaps    %xmm4, %xmm5                                  #111.31
..L18:                                                          #
        mulsd     %xmm7, %xmm5                                  #111.31
        lea       (%r9,%rcx,8), %r8                             #110.5
        addsd     (%r8), %xmm3                                  #110.5
        movq      184(%rsp), %rax                               #13.1
        xorps     %xmm6, %xmm6                                  #112.38
        movsd     %xmm3, (%r8)                                  #110.5
        shrq      $1, %rsi                                      #112.38
        cvtsi2sdq %rbp, %xmm6                                   #112.38
        addsd     (%r9,%rax,8), %xmm5                           #111.5
        movsd     %xmm5, (%r9,%rax,8)                           #111.5
        movq      %rbp, %r9                                     #112.38
        andq      $1, %r9                                       #112.38
        orq       %rsi, %r9                                     #112.38
        testq     %rbp, %rbp                                    #112.38
        cvtsi2sdq %r9, %xmm8                                    #112.38
        addsd     %xmm8, %xmm8                                  #112.38
        jl        ..L19         # Prob 50%                      #112.38
        movaps    %xmm6, %xmm8                                  #112.38
..L19:                                                          #
        mulsd     %xmm7, %xmm8                                  #112.38
        addsd     (%r8,%rax,8), %xmm8                           #112.5
        movsd     %xmm8, (%r8,%rax,8)                           #112.5
        addq      $128, %rsp                                    #113.1
..___tag_value_bin_gemm_2x2.20:                                 #
        popq      %rbp                                          #113.1
..___tag_value_bin_gemm_2x2.22:                                 #
        popq      %rbx                                          #113.1
..___tag_value_bin_gemm_2x2.24:                                 #
        popq      %r15                                          #113.1
..___tag_value_bin_gemm_2x2.26:                                 #
        popq      %r14                                          #113.1
..___tag_value_bin_gemm_2x2.28:                                 #
        popq      %r13                                          #113.1
..___tag_value_bin_gemm_2x2.30:                                 #
        ret                                                     #113.1
        .align    16,0x90
..___tag_value_bin_gemm_2x2.31:                                 #
                                # LOE
# mark_end;
	.type	bin_gemm_2x2,@function
	.size	bin_gemm_2x2,.-bin_gemm_2x2
	.data
# -- End  bin_gemm_2x2
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
	.4byte 0x00000014
	.8byte 0x7801000100000000
	.8byte 0x0000019008070c10
	.4byte 0x00000000
	.4byte 0x0000008c
	.4byte 0x0000001c
	.8byte ..___tag_value_bin_gemm_2x2.1
	.8byte ..___tag_value_bin_gemm_2x2.31-..___tag_value_bin_gemm_2x2.1
	.byte 0x04
	.4byte ..___tag_value_bin_gemm_2x2.3-..___tag_value_bin_gemm_2x2.1
	.4byte 0x028d100e
	.byte 0x04
	.4byte ..___tag_value_bin_gemm_2x2.5-..___tag_value_bin_gemm_2x2.3
	.4byte 0x038e180e
	.byte 0x04
	.4byte ..___tag_value_bin_gemm_2x2.7-..___tag_value_bin_gemm_2x2.5
	.4byte 0x048f200e
	.byte 0x04
	.4byte ..___tag_value_bin_gemm_2x2.9-..___tag_value_bin_gemm_2x2.7
	.4byte 0x0583280e
	.byte 0x04
	.4byte ..___tag_value_bin_gemm_2x2.11-..___tag_value_bin_gemm_2x2.9
	.4byte 0x0686300e
	.byte 0x04
	.4byte ..___tag_value_bin_gemm_2x2.13-..___tag_value_bin_gemm_2x2.11
	.4byte 0x0401b00e
	.4byte ..___tag_value_bin_gemm_2x2.14-..___tag_value_bin_gemm_2x2.13
	.2byte 0x138c
	.byte 0x04
	.4byte ..___tag_value_bin_gemm_2x2.15-..___tag_value_bin_gemm_2x2.14
	.2byte 0x04cc
	.4byte ..___tag_value_bin_gemm_2x2.20-..___tag_value_bin_gemm_2x2.15
	.4byte 0x04c6300e
	.4byte ..___tag_value_bin_gemm_2x2.22-..___tag_value_bin_gemm_2x2.20
	.4byte 0x04c3280e
	.4byte ..___tag_value_bin_gemm_2x2.24-..___tag_value_bin_gemm_2x2.22
	.4byte 0x04cf200e
	.4byte ..___tag_value_bin_gemm_2x2.26-..___tag_value_bin_gemm_2x2.24
	.4byte 0x04ce180e
	.4byte ..___tag_value_bin_gemm_2x2.28-..___tag_value_bin_gemm_2x2.26
	.4byte 0x04cd100e
	.4byte ..___tag_value_bin_gemm_2x2.30-..___tag_value_bin_gemm_2x2.28
	.8byte 0x000000000000080e
	.byte 0x00
# End
