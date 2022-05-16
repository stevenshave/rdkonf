# rdkonf v6

Implementation of high quality RDKit conformer generator as described in: Ebejer, Jean-Paul, Garrett M. Morris, and Charlotte M. Deane. "Freely available conformer generation methods: how good are they?." Journal of chemical information and modeling 52.5 (2012): 1146-1158.

Please reference this repository in publications using this software

## Usage
----
> python rdkonf.py -h

> usage: rdkonf.py [-h] inFile outFile numConformers [removeHydrogens]
>
> rdkonf conformer generation
>
> positional arguments:
>
>  inFile           In file (SDF, or SMI)
>
>  outFile          Out file (SDF)
>
>  numConformers    The number of conformers per molecule to make
>
>  removeHydrogens  Should molecules be written using implicit hydrogens
>                   (optional, default = 1)
> 
> optional arguments:
>
>   -h, --help       show this help message and exit

Using the supplied test.smi file, we can gerate 1 low energy conformer for each as follows:

>python rdkonf.py test.smi test_3d.sdf 1

or 10:
>python rdkonf.py test.smi test_3d.sdf 10

Optionally, we may request the 10 diverse conformers  include implicit hydrogens:

>python rdkonf.py test.smi test_3d_H.sdf 10 0

If not enough diverse conformers are found, then less than requested may be returned.

Molecules have an underscore and conformer number added to their title, output for diclofenac included as the first molecule in the test.smi becomes:
>
>Diclofenac_1
>
>   RDKit          3D
>
>
>
> 19 20  0  0  0  0  0  0  0  0999 V2000
>
>   -4.1029   -1.3436   -0.9698 O   0  0  0  0  0  0  0  0  0  0  0  0
>
>   -3.1212   -1.1835   -0.1992 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>   -2.2528   -2.2392    0.0446 O   0  0  0  0  0  0  0  0  0  0  0  0
>
>   -2.8689    0.1488    0.4203 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>   -1.9358    0.9521   -0.4277 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>   -2.4673    1.7860   -1.4186 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>   -1.6205    2.5338   -2.2401 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>   -0.2364    2.4506   -2.0762 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>    0.3012    1.6164   -1.0933 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>   -0.5409    0.8531   -0.2692 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>    0.0060    0.0133    0.7574 N   0  0  0  0  0  0  0  0  0  0  0  0
>
>    1.3809   -0.3846    0.8433 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>    2.0990   -0.1307    2.0249 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>    1.3498    0.6950    3.3814 Cl  0  0  0  0  0  0  0  0  0  0  0  0
>
>    3.4385   -0.5149    2.1288 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>    4.0653   -1.1728    1.0698 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>    3.3509   -1.4574   -0.0949 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>    2.0121   -1.0736   -0.2148 C   0  0  0  0  0  0  0  0  0  0  0  0
>
>    1.1432   -1.5488   -1.6665 Cl  0  0  0  0  0  0  0  0  0  0  0  0
>
>  1  2  2  0
>
>  2  3  1  0
>
>  2  4  1  0
>
>  4  5  1  0
>
>  5  6  2  0
>
>  6  7  1  0
>
>  7  8  2  0
>
>  8  9  1  0
>
>  9 10  2  0
>
> 10 11  1  0
>
> 11 12  1  0
>
> 12 13  2  0
>
> 13 14  1  0
>
> 13 15  1  0
>
> 15 16  2  0
>
> 16 17  1  0
>
> 17 18  2  0
>
> 18 19  1  0
>
> 10  5  1  0
>
> 18 12  1  0
>
>M  END
>
>$$$$