V34 :0x4 cufft_m
21 ../common/cufft_m.cuf S624 0
07/27/2019  19:55:15
use pgi_acc_common private
use iso_c_binding private
enduse
B 525 iso_c_binding c_loc
B 526 iso_c_binding c_funloc
B 527 iso_c_binding c_associated
B 528 iso_c_binding c_f_pointer
B 529 iso_c_binding c_f_procpointer
B 608 iso_c_binding c_sizeof
D 58 20 27
D 60 26 663 8 662 7
D 69 26 666 8 665 7
D 78 20 87
D 80 26 663 8 662 7
D 101 26 746 8 745 7
D 3786 23 13 1 3 0 0 0 0 1 0
 0 0 3 3 0 6768
D 3789 23 13 1 3 0 0 0 0 1 0
 0 0 3 3 0 6769
D 3792 23 14 1 3 0 0 0 0 1 0
 0 0 3 3 0 6770
D 3795 23 14 1 3 0 0 0 0 1 0
 0 0 3 3 0 6771
D 3798 23 9 1 3 0 0 0 0 1 0
 0 0 3 3 0 6772
D 3801 23 13 1 3 0 0 0 0 1 0
 0 0 3 3 0 6773
D 3804 23 10 1 3 0 0 0 0 1 0
 0 0 3 3 0 6774
D 3807 23 14 1 3 0 0 0 0 1 0
 0 0 3 3 0 6775
D 3810 23 9 1 3 0 0 0 0 1 0
 0 0 3 3 0 6776
D 3813 23 9 1 3 0 0 0 0 1 0
 0 0 3 3 0 6777
D 3816 23 10 1 3 0 0 0 0 1 0
 0 0 3 3 0 6778
D 3819 23 10 1 3 0 0 0 0 1 0
 0 0 3 3 0 6779
D 3822 20 5745
D 3824 20 125
D 3826 20 5745
S 624 24 0 0 0 9 1 0 5015 10005 0 A 0 0 0 0 B 0 18 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 18 0 0 0 0 0 0 cufft_m
S 625 6 4 0 0 6 627 624 5023 80000c 0 A 0 0 0 0 B 0 20 0 0 0 0 0 0 0 0 0 0 11223 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 cufft_forward
S 627 6 4 0 0 6 628 624 5037 80000c 0 A 0 0 0 0 B 0 21 0 0 0 4 0 0 0 0 0 0 11223 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 cufft_inverse
S 628 6 4 0 0 6 630 624 5051 80000c 0 A 0 0 0 0 B 0 22 0 0 0 8 0 0 0 0 0 0 11223 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 cufft_r2c
S 630 6 4 0 0 6 632 624 5061 80000c 0 A 0 0 0 0 B 0 23 0 0 0 12 0 0 0 0 0 0 11223 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 cufft_c2r
S 632 6 4 0 0 6 634 624 5071 80000c 0 A 0 0 0 0 B 0 24 0 0 0 16 0 0 0 0 0 0 11223 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 cufft_c2c
S 634 6 4 0 0 6 636 624 5081 80000c 0 A 0 0 0 0 B 0 25 0 0 0 20 0 0 0 0 0 0 11223 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 cufft_d2z
S 636 6 4 0 0 6 638 624 5091 80000c 0 A 0 0 0 0 B 0 26 0 0 0 24 0 0 0 0 0 0 11223 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 cufft_z2d
S 638 6 4 0 0 6 1 624 5101 80000c 0 A 0 0 0 0 B 0 27 0 0 0 28 0 0 0 0 0 0 11223 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 cufft_z2z
S 640 19 0 0 0 9 1 624 5111 4000 0 A 0 0 0 0 B 0 30 0 0 0 0 0 0 0 643 0 0 0 0 0 0 9 1 0 0 0 0 0 624 0 0 0 0 cufftdestroy
O 640 1 643
S 641 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 642 3 0 0 0 58 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 5124 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 12 63 75 66 66 74 44 65 73 74 72 6f 79
S 643 14 5 0 0 0 1 624 5111 0 18000 A 1000000 0 0 0 B 0 0 0 0 0 0 0 1 1 640 0 0 0 0 0 0 0 0 0 0 0 31 0 624 0 0 642 0 cufftdestroy
F 643 1 644
S 644 1 3 0 0 60 1 643 5137 2004 6000 A 0 0 0 0 B 0 31 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
R 662 25 6 iso_c_binding c_ptr
R 663 5 7 iso_c_binding val c_ptr
R 665 25 9 iso_c_binding c_funptr
R 666 5 10 iso_c_binding val c_funptr
R 700 6 44 iso_c_binding c_null_ptr$ac
R 702 6 46 iso_c_binding c_null_funptr$ac
R 703 26 47 iso_c_binding ==
R 705 26 49 iso_c_binding !=
S 732 19 0 0 0 9 1 624 5861 4000 0 A 0 0 0 0 B 0 37 0 0 0 0 0 0 0 735 0 0 0 0 0 0 987 1 0 0 0 0 0 624 0 0 0 0 cufftsetstream
O 732 1 735
S 733 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 14 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 734 3 0 0 0 78 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 5876 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 14 63 75 66 66 74 53 65 74 53 74 72 65 61 6d
S 735 14 5 0 0 0 1 624 5861 0 18000 A 1000000 0 0 0 B 0 0 0 0 0 0 0 14 2 732 0 0 0 0 0 0 0 0 0 0 0 38 0 624 0 0 734 0 cufftsetstream
F 735 2 736 737
S 736 1 3 0 0 60 1 735 5137 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
S 737 1 3 0 0 7 1 735 5891 2004 6000 A 0 0 0 0 B 0 38 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 stream
R 745 25 6 pgi_acc_common c_devptr
R 746 5 7 pgi_acc_common cptr c_devptr
R 752 6 13 pgi_acc_common c_null_devptr$ac
R 756 26 17 pgi_acc_common =
S 810 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 5607 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 11143 19 0 0 0 9 1 624 75607 4000 0 A 0 0 0 0 B 0 46 0 0 0 0 0 0 0 0 0 0 0 0 0 0 993 6 0 0 0 0 0 624 0 0 0 0 cufftexec
O 11143 6 11187 11180 11173 11165 11157 11145
S 11144 3 0 0 0 58 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 75617 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 12 63 75 66 66 74 45 78 65 63 43 32 43
S 11145 14 5 0 0 0 1 624 75630 0 18000 A 1000000 0 0 0 B 0 48 0 0 0 0 0 4994 4 0 0 0 0 0 0 0 0 0 0 0 0 48 0 624 0 0 11144 0 cufftexecc2c
F 11145 4 11146 11147 11148 11149
S 11146 1 3 0 0 60 1 11145 5137 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
S 11147 7 3 0 0 3786 1 11145 75643 80a104 2000 A 0 0 0 0 B 0 48 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 idata
S 11148 7 3 0 0 3789 1 11145 75649 80a104 2000 A 0 0 0 0 B 0 48 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 odata
S 11149 1 3 0 0 6 1 11145 75655 2004 6000 A 0 0 0 0 B 0 48 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 direction
S 11154 6 1 0 0 6 1 11145 9404 40800006 2000 A 0 0 0 0 B 0 55 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0
S 11155 6 1 0 0 6 1 11145 6233 40800006 2000 A 0 0 0 0 B 0 55 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_1
S 11156 3 0 0 0 58 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 75717 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 12 63 75 66 66 74 45 78 65 63 5a 32 5a
S 11157 14 5 0 0 0 1 624 75730 0 18000 A 1000000 0 0 0 B 0 58 0 0 0 0 0 4998 4 0 0 0 0 0 0 0 0 0 0 0 0 58 0 624 0 0 11156 0 cufftexecz2z
F 11157 4 11158 11159 11160 11161
S 11158 1 3 0 0 60 1 11157 5137 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
S 11159 7 3 0 0 3792 1 11157 75643 80a104 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 idata
S 11160 7 3 0 0 3795 1 11157 75649 80a104 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 odata
S 11161 1 3 0 0 6 1 11157 75655 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 direction
S 11162 6 1 0 0 6 1 11157 6239 40800006 2000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2
S 11163 6 1 0 0 6 1 11157 6245 40800006 2000 A 0 0 0 0 B 0 65 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3
S 11164 3 0 0 0 58 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 75743 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 12 63 75 66 66 74 45 78 65 63 52 32 43
S 11165 14 5 0 0 0 1 624 75756 0 18000 A 1000000 0 0 0 B 0 68 0 0 0 0 0 5002 3 0 0 0 0 0 0 0 0 0 0 0 0 68 0 624 0 0 11164 0 cufftexecr2c
F 11165 3 11166 11167 11168
S 11166 1 3 0 0 60 1 11165 5137 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
S 11167 7 3 0 0 3798 1 11165 75643 80a104 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 idata
S 11168 7 3 0 0 3801 1 11165 75649 80a104 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 odata
S 11170 6 1 0 0 6 1 11165 9448 40800006 2000 A 0 0 0 0 B 0 75 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_4
S 11171 6 1 0 0 6 1 11165 9454 40800006 2000 A 0 0 0 0 B 0 76 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5
S 11172 3 0 0 0 58 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 75769 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 12 63 75 66 66 74 45 78 65 63 44 32 5a
S 11173 14 5 0 0 0 1 624 75782 0 18000 A 1000000 0 0 0 B 0 79 0 0 0 0 0 5005 3 0 0 0 0 0 0 0 0 0 0 0 0 79 0 624 0 0 11172 0 cufftexecd2z
F 11173 3 11174 11175 11176
S 11174 1 3 0 0 60 1 11173 5137 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
S 11175 7 3 0 0 3804 1 11173 75643 80a104 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 idata
S 11176 7 3 0 0 3807 1 11173 75649 80a104 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 odata
S 11178 6 1 0 0 6 1 11173 9460 40800006 2000 A 0 0 0 0 B 0 86 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6
S 11179 6 1 0 0 6 1 11173 9466 40800006 2000 A 0 0 0 0 B 0 87 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_7
S 11180 14 5 0 0 0 1 624 75795 0 18000 A 1000000 0 0 0 B 0 90 0 0 0 0 0 5008 3 0 0 0 0 0 0 0 0 0 0 0 0 90 0 624 0 0 11164 0 cufftexecr2cinplace
F 11180 3 11181 11182 11183
S 11181 1 3 0 0 60 1 11180 5137 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
S 11182 7 3 0 0 3810 1 11180 75643 80a104 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 idata
S 11183 7 3 0 0 3813 1 11180 75649 80a104 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 odata
S 11185 6 1 0 0 6 1 11180 9489 40800006 2000 A 0 0 0 0 B 0 97 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_8
S 11186 6 1 0 0 6 1 11180 9495 40800006 2000 A 0 0 0 0 B 0 98 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_9
S 11187 14 5 0 0 0 1 624 75815 0 18000 A 1000000 0 0 0 B 0 101 0 0 0 0 0 5011 3 0 0 0 0 0 0 0 0 0 0 0 0 101 0 624 0 0 11172 0 cufftexecd2zinplace
F 11187 3 11188 11189 11190
S 11188 1 3 0 0 60 1 11187 5137 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
S 11189 7 3 0 0 3816 1 11187 75643 80a104 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 idata
S 11190 7 3 0 0 3819 1 11187 75649 80a104 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 odata
S 11191 6 1 0 0 6 1 11187 9501 40800006 2000 A 0 0 0 0 B 0 107 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_10
S 11192 6 1 0 0 6 1 11187 9508 40800006 2000 A 0 0 0 0 B 0 108 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_11
S 11193 19 0 0 0 9 1 624 75835 4000 0 A 0 0 0 0 B 0 113 0 0 0 0 0 0 0 11195 0 0 0 0 0 0 994 1 0 0 0 0 0 624 0 0 0 0 cufftplan1d
O 11193 1 11195
S 11194 3 0 0 0 3822 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 75847 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 11 63 75 66 66 74 50 6c 61 6e 31 64
S 11195 14 5 0 0 0 1 624 75835 0 18000 A 1000000 0 0 0 B 0 0 0 0 0 0 0 5014 4 11193 0 0 0 0 0 0 0 0 0 0 0 114 0 624 0 0 11194 0 cufftplan1d
F 11195 4 11196 11197 11198 11199
S 11196 1 3 0 0 60 1 11195 5137 2004 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
S 11197 1 3 0 0 6 1 11195 75859 2004 6000 A 0 0 0 0 B 0 114 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nx
S 11198 1 3 0 0 6 1 11195 56400 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 type
S 11199 1 3 0 0 6 1 11195 75862 2004 6000 A 0 0 0 0 B 0 114 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 batch
S 11200 19 0 0 0 9 1 624 75868 4000 0 A 0 0 0 0 B 0 123 0 0 0 0 0 0 0 11202 0 0 0 0 0 0 995 1 0 0 0 0 0 624 0 0 0 0 cufftplanmany
O 11200 1 11202
S 11201 3 0 0 0 3824 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 75882 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 13 63 75 66 66 74 50 6c 61 6e 4d 61 6e 79
S 11202 14 5 0 0 0 1 624 75868 0 18000 A 1000000 0 0 0 B 0 0 0 0 0 0 0 5018 11 11200 0 0 0 0 0 0 0 0 0 0 0 124 0 624 0 0 11201 0 cufftplanmany
F 11202 11 11203 11204 11205 11206 11207 11208 11209 11210 11211 11212 11213
S 11203 1 3 0 0 60 1 11202 5137 2004 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
S 11204 1 3 0 0 6 1 11202 75896 2004 6000 A 0 0 0 0 B 0 124 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rank
S 11205 1 3 0 0 6 1 11202 75901 2004 2000 A 0 0 0 0 B 0 124 0 0 0 0 0 0 0 0 0 0 0 0 0 0 47 0 0 0 0 0 0 0 0 0 0 0 n
S 11206 1 3 0 0 6 1 11202 75903 2004 2000 A 0 0 0 0 B 0 124 0 0 0 0 0 0 0 0 0 0 0 0 0 0 47 0 0 0 0 0 0 0 0 0 0 0 inembed
S 11207 1 3 0 0 6 1 11202 75911 2004 6000 A 0 0 0 0 B 0 124 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 istride
S 11208 1 3 0 0 6 1 11202 75919 2004 6000 A 0 0 0 0 B 0 124 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 idist
S 11209 1 3 0 0 6 1 11202 75925 2004 2000 A 0 0 0 0 B 0 124 0 0 0 0 0 0 0 0 0 0 0 0 0 0 47 0 0 0 0 0 0 0 0 0 0 0 onembed
S 11210 1 3 0 0 6 1 11202 75933 2004 6000 A 0 0 0 0 B 0 124 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ostride
S 11211 1 3 0 0 6 1 11202 75941 2004 6000 A 0 0 0 0 B 0 124 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 odist
S 11212 1 3 0 0 6 1 11202 56400 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 type
S 11213 1 3 0 0 6 1 11202 75862 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 batch
S 11214 19 0 0 0 9 1 624 75947 4000 0 A 0 0 0 0 B 0 136 0 0 0 0 0 0 0 11218 0 0 0 0 0 0 997 1 0 0 0 0 0 624 0 0 0 0 cufftplan2d
O 11214 1 11215
S 11215 27 0 0 0 9 11224 624 75959 0 400000 A 0 0 0 0 B 0 137 0 0 0 0 999 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 cufftplan2dswap
Q 11215 11214 0
S 11216 19 0 0 0 9 1 624 75975 4000 0 A 0 0 0 0 B 0 140 0 0 0 0 0 0 0 0 0 0 0 0 0 0 998 1 0 0 0 0 0 624 0 0 0 0 cufftplan2dc
O 11216 1 11218
S 11217 3 0 0 0 3826 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 75988 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 11 63 75 66 66 74 50 6c 61 6e 32 64
S 11218 14 5 0 0 0 1 624 75947 0 18000 A 1000000 0 0 0 B 0 0 0 0 0 0 0 5029 4 11214 0 0 0 0 0 0 0 0 0 0 0 141 0 624 0 0 11217 0 cufftplan2d
F 11218 4 11219 11220 11221 11222
S 11219 1 3 0 0 60 1 11218 5137 2004 2000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
S 11220 1 3 0 0 6 1 11218 75859 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nx
S 11221 1 3 0 0 6 1 11218 76000 2004 6000 A 0 0 0 0 B 0 141 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ny
S 11222 1 3 0 0 6 1 11218 56400 2004 6000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 type
S 11223 11 0 0 0 9 760 624 76003 40800000 805000 A 0 0 0 0 B 0 149 0 0 0 32 0 0 625 638 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _cufft_m$8
S 11224 23 5 0 0 0 11226 624 75959 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cufftplan2dswap
S 11225 1 3 0 0 60 1 11224 5137 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 plan
S 11226 14 5 0 0 0 1 11224 75959 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 5034 4 0 0 0 0 0 0 0 0 0 0 0 0 151 0 624 0 0 0 0 cufftplan2dswap
F 11226 4 11225 11227 11228 11229
S 11227 1 3 0 0 6 1 11224 76014 4 7000 A 0 0 0 0 B 0 156 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _V_nx
S 11228 1 3 0 0 6 1 11224 76020 4 7000 A 0 0 0 0 B 0 156 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _V_ny
S 11229 1 3 0 0 6 1 11224 76026 4 7000 A 0 0 0 0 B 0 156 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _V_type
A 27 2 0 0 0 6 641 0 0 0 27 0 0 0 0 0 0 0 0 0 0 0
A 83 1 0 0 0 60 700 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 86 1 0 0 0 69 702 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 87 2 0 0 0 6 733 0 0 0 87 0 0 0 0 0 0 0 0 0 0 0
A 103 1 0 0 0 101 752 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 125 2 0 0 0 6 810 0 0 0 125 0 0 0 0 0 0 0 0 0 0 0
A 5745 2 0 0 5255 6 5607 0 0 0 5745 0 0 0 0 0 0 0 0 0 0 0
A 6768 1 0 0 6152 6 11154 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6769 1 0 0 5098 6 11155 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6770 1 0 0 6547 6 11162 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6771 1 0 0 4487 6 11163 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6772 1 0 0 5133 6 11170 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6773 1 0 0 5550 6 11171 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6774 1 0 0 6646 6 11178 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6775 1 0 0 6277 6 11179 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6776 1 0 0 6710 6 11185 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6777 1 0 0 6664 6 11186 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6778 1 0 0 5150 6 11191 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 6779 1 0 0 6167 6 11192 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 131 1 1
V 83 60 7 0
S 0 60 0 0 0
A 0 6 0 0 1 2 0
J 132 1 1
V 86 69 7 0
S 0 69 0 0 0
A 0 6 0 0 1 2 0
J 36 1 1
V 103 101 7 0
S 0 101 0 0 0
A 0 80 0 0 1 83 0
Z
