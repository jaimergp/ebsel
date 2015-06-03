ebsel
=============================
(EMSL_Basis_Set_Exchange_Local)

Create a local copy of the famous [EMSL Basis Set Exchange](https://bse.pnl.gov/bse/portal) and use it easily through the API.

* Make a compact copy (~90 MB of Sqlite3 database files) of the EMSL Basis Set Exchange website
* __API for scripting__
* Quick local access without delay
* Only needs [Python](https://www.python.org/)

##Dependencies
* Python 2.7

###### Optional
If you plan to re-download data from the BSE - not just use the included snapshots - you need:
* The [requests](http://docs.python-requests.org/en/latest/) python module. ```pip install requests``` (Install it in a virtualenv or with sudo.)

##Installation
* Download the git repository (```cd ~ && git clone https://github.com/mattbernst/ebsel.git``` for example).
* You can now use ```EMSL_api.py``` mostly as described for [TApplencourt's original  EMSL_Basis_Set_Exchange_Local.](https://github.com/TApplencourt/EMSL_Basis_Set_Exchange_Local)
* This fork is more tailored for use as a library called by other programs.

##Usage:
_export PYTHONPATH=$HOME:$PYTHONPATH (presuming you cloned into your home directory)_
```
In [1]: from ebsel.src import EMSL_local

In [2]: el_gaussian94 = EMSL_local.EMSL_local(fmt="g94")

In [3]: print el_gaussian94.get_available_elements("6-31G*")
['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']

In [4]: print el_gaussian94.get_available_basis_sets(elements=["C", "H", "O", "N", "Se"])[:5]
[(u'3-21G', u"VDZ Valence Double Zeta: 2 Funct.'s/Valence AO"), (u'6-311G', u"VTZ Valence Triple Zeta: 3 Funct.'s/Valence AO"), (u'6-311G*', u'VTZP Valence Triple Zeta + Polarization on Heavy Atoms'), (u'6-311G**', u'VTZP Valence Triple Zeta + Polarization on All Atoms'), (u'6-311G** Polarization', u'1P Polarization functions associated with 6-311G')]

In [5]: print "\n".join(el_gaussian94.get_basis("6-311G", elements=["O", "Se"]))
****
O     0 
S   6   1.00
   8588.5000000              0.00189515       
   1297.2300000              0.0143859        
    299.2960000              0.0707320        
     87.3771000              0.2400010        
     25.6789000              0.5947970        
      3.7400400              0.2808020        
SP   3   1.00
     42.1175000              0.1138890              0.0365114        
      9.6283700              0.9208110              0.2371530        
      2.8533200             -0.00327447             0.8197020        
SP   1   1.00
      0.9056610              1.0000000              1.0000000        
SP   1   1.00
      0.2556110              1.0000000              1.0000000        
****
Se     0 
S   6   1.00
 405400.0000000              0.0008310        
  60850.0000000              0.0062331        
  13910.0000000              0.0320710        
   3989.0000000              0.1263600        
   1324.0000000              0.3889900        
    487.0000000              0.5488800        
S   3   1.00
    487.0000000              0.1796300        
    193.2000000              0.6222700        
     82.1900000              0.2506700        
S   1   1.00
     30.5200000              1.0000000        
S   1   1.00
     13.2500000              1.0000000        
S   1   1.00
      4.5100000              1.0000000        
S   1   1.00
      1.8670000              1.0000000        
S   1   1.00
      0.3190000              1.0000000        
S   1   1.00
      0.1120000              1.0000000        
P   3   1.00
   2706.0000000              0.0221460        
    638.6000000              0.1818000        
    203.8000000              0.8615500        
P   3   1.00
     75.9600000              0.3424300        
     31.0200000              0.5054100        
     13.0500000              0.2623700        
P   3   1.00
     13.0500000              0.0701630        
      6.9860000              0.3841900        
      3.2340000              0.6036800        
P   1   1.00
      1.4750000              1.0000000        
P   1   1.00
      0.7275000              1.0000000        
P   1   1.00
      0.2869000              1.0000000        
P   1   1.00
      0.0967900              1.0000000        
D   4   1.00
     99.0100000              0.0255960        
     28.4100000              0.1545900        
      9.8630000              0.4287800        
      3.5140000              0.5862000        
D   1   1.00
      1.1710000              1.0000000        
****

#could also repeat all above and get NWChem or GAMESS-US formats using
#el_nwchem=EMSL_local.EMSL_local(fmt="nwchem") or
#el_gamess=EMSL_local.EMSL_loccal(fmt="gamess-us")
```

###### Supplementing data from the Basis Set Exchange
The EMSL Basis Set Exchange hosts many but not all basis sets. Sometimes you may want to use a basis set that is not in the BSE, but keep using the same API. If you add an NWChem format basis file in the db/nwchem/ directory, with a name ending in .nwbas, you can use that basis set data from the API for __all__ supported formats. The NWChem data will be automatically extracted and reformatted for use in GAMESS-US or Gaussian 94 formats as necessary. See the test cases involving "g3mp2large" in test_local.py for examples.

Feel free to fork/pull request. 

In papers where you use the basis sets obtained from the Basis Set Exchange please cite this :
>The Role of Databases in Support of Computational Chemistry Calculations
>
>>--<cite>Feller, D.; J. Comp. Chem., 17(13), 1571-1586, 1996.</cite>

>Basis Set Exchange: A Community Database for Computational Sciences
>
>>--<cite>Schuchardt, K.L., Didier, B.T., Elsethagen, T., Sun, L., Gurumoorthi, V., Chase, J., Li, J., and Windus ; T.L.
>>J. Chem. Inf. Model., 47(3), 1045-1052, 2007, doi:10.1021/ci600510j.</cite>

And don't forget: 
>These documents may be freely distributed and used for non-commercial, scientific and educational purposes. 
>-- <cite>http://www.pnl.gov/notices.asp</cite>

