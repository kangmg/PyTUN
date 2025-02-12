# easyTUN

## 특징
* pyTUN 스크립트를 수정하여 여러 온도에 대해서 transmission coefficient와 rate constant를 계산할 수 있도록 했습니다.
* 현재는 Gaussian output 파일만을 지원하며 Wigner, Skodje, Eckart correction을 지원합니다.

## 실행 방법
```
> python easyTUN.py test.inp
```

## Test data
* 반응 종류 : ClFHP 분자의 configurational inversion
* 계산 레벨 : CCSD/def2SVP
* 프로그램  : Gaussian 16

## 파일
easyTUN : 실행파일
|File Name|Content|
|:-:|:-:|
|test.inp|example input file|
|test.out|example output file|
|calculations/*.xyz| xyz geometry|
|calculations/*_sp.out| Single Point calculation file |
|calculations/*.out| Freqeucny calculation file |


## 인풋 구조
```
# Rxn symmetry number #
RXN Symmetry Number = auto


# Temperatures #
#   Style #1  :  Temp_1   Temp_2   Temp_3   Temp_4 ...
#   Style #2  :  start_Temp   end_Temp   interval
Temperatures list [K] = 300 500 10


# Special Temperatures (optional) #
#   Style #1  :  Stemp_1 Stemp_2 Stemp_3 ... e.g. 273.15 293.15 298.15
#   Style #2  :  'auto'  will  set  273.15(STP) 293.15(NTP) 298.15(SATP)
Special Temperatures [K] = auto


# Input Files #
# ---------------------------------------------------------------------- #
Reactant Frequency file = calculations/reactant.out
Reactant Single Point file = calculations/reactant_sp.out
Reactant XYZ file = calculations/reactant.xyz

TS Frequency file = calculations/ts.out
TS Single Point file = calculations/ts_sp.out
TS XYZ file = calculations/ts.xyz

Product Frequency file = calculations/product.out
Product Single Point file = calculations/product_sp.out
Product XYZ file = calculations/product.xyz
# ---------------------------------------------------------------------- #
```

## 실행 결과
```
# Level of thoery Info. # 
Freqency calculation  : RCCSD-FC/def2SVP 
SP calculation        : RCCSD-FC/def2SVP

# RXN Info. # 
Detected Point group  : Reactant [ C1 ] / TS [ C1 ] 
RXN Symmetry Number   : 1.0 
Rxn Barrier           : 339456.1985099471  [J/mol]

TEMPERATURE      kappa_wigner     kappa_skodje     kappa_eckart     k_wiger             k_skodje              k_eckart             k_uncorrtd 
------------------------------------------------------------------------------------------------------------------------------------------------
273.15           2.5646           40.3500          15.7541          1.7827E-52           2.8049E-51           1.0951E-51           6.9515E-53
293.15           2.3584           10.0935          8.3539           4.7211E-48           2.0206E-47           1.6723E-47           2.0018E-48
298.15           2.3132           8.5469           7.4562           4.8676E-47           1.7985E-46           1.5690E-46           2.1043E-47
300.00           2.2970           8.0926           7.1719           1.1316E-46           3.9867E-46           3.5331E-46           4.9263E-47
310.00           2.2147           6.3123           5.9576           9.0916E-45           2.5913E-44           2.4457E-44           4.1051E-45
320.00           2.1400           5.2063           5.1158           5.5581E-43           1.3522E-42           1.3287E-42           2.5973E-43
330.00           2.0719           4.4554           4.5031           2.6505E-41           5.6996E-41           5.7605E-41           1.2793E-41
340.00           2.0098           3.9139           4.0396           1.0078E-39           1.9626E-39           2.0257E-39           5.0145E-40
350.00           1.9529           3.5062           3.6781           3.1153E-38           5.5931E-38           5.8673E-38           1.5952E-38
360.00           1.9007           3.1888           3.3889           7.9651E-37           1.3363E-36           1.4201E-36           4.1906E-37
370.00           1.8527           2.9354           3.1525           1.7105E-35           2.7102E-35           2.9106E-35           9.2328E-36
380.00           1.8084           2.7287           2.9558           3.1283E-34           4.7203E-34           5.1131E-34           1.7299E-34
390.00           1.7675           2.5572           2.7896           4.9326E-33           7.1366E-33           7.7850E-33           2.7907E-33
400.00           1.7296           2.4129           2.6473           6.7805E-32           9.4594E-32           1.0378E-31           3.9203E-32
410.00           1.6944           2.2899           2.5241           8.2078E-31           1.1092E-30           1.2227E-30           4.8440E-31
420.00           1.6617           2.1840           2.4164           8.8290E-30           1.1604E-29           1.2839E-29           5.3131E-30
430.00           1.6313           2.0920           2.3214           8.5092E-29           1.0912E-28           1.2109E-28           5.2161E-29
440.00           1.6030           2.0114           2.2370           7.4031E-28           9.2893E-28           1.0331E-27           4.6184E-28
450.00           1.5765           1.9402           2.1615           5.8538E-27           7.2046E-27           8.0261E-27           3.7133E-27
460.00           1.5517           1.8770           2.0935           4.2332E-26           5.1208E-26           5.7114E-26           2.7282E-26
470.00           1.5284           1.8206           2.0320           2.8157E-25           3.3538E-25           3.7434E-25           1.8422E-25
480.00           1.5067           1.7699           1.9761           1.7316E-24           2.0341E-24           2.2711E-24           1.1493E-24
490.00           1.4862           1.7241           1.9251           9.8929E-24           1.1477E-23           1.2815E-23           6.6566E-24
500.00           1.4669           1.6827           1.8784           5.2741E-23           6.0500E-23           6.7535E-23           3.5954E-23
```


## pyTUN에서의 변화점
* xyz file로부터 RXN symmetry를 자동으로 예측하도록 코드를 추가했고, input 파일을 읽는 방식으로 바꿨습니다.
* python3에서 구동되도록 코드 스타일을 조금 수정했습니다.
## Todo

- [ ] Enables the parsing of Psi4 and ORCA output files

[Oiginal README file](https://github.com/kangmg/PyTUN/blob/master/information.md)
