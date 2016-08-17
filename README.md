# poh5

I/O routines for poh5 (PaNDa on HDF5) file format written by C.

# What is poh5?

Poh5 file format is re-implementation of PaNDa format on HDF5.

PaNDa is an original format for NICAM datafile, which is called
'ADVANCED' in NICAM source/namelist.

Different from original PaNDa format, you can include
horizontal grid (hgrid) in each file.

# APIs

There are currentry two implementation of APIs, by C and by python,
the latter is in other repository.

# Details of poh5 format

## Global Attributes
- glevel
- rlevel
- grid\_topology
- is\_complete
- description
- note
- my_pe
- num\_of\_pe
- num\_of\_var
- rgnid

## Group Hierarchy
The root group of HDF5 ("/") has two group, "Var" and "Grd".

"Grd" group is to contain grid coordinates, currently only horizontal `hgrid` is supported.
`hgrid` consists of two arrays, "grd_x" and "grd_xt", these are stored as 4/5 dimension dataset in a group "Hgrid".

"Var" group is to contain multiple variable data, one sub group for one variable.
Each variable group has two datasets, "data" and "time".
"data" dataset is 5-dimensioned, as `var[t,l,k,j,i]`.
"time" dataset is 1-dimensioned of compound datatype, two member struct.

![](http://g.gravizo.com/g?
  digraph G {
    rankdir = LR
    root[shape=folder,label="/"]
    Grd[shape=folder]
    Hgrid[shape=folder]
    Var[shape=folder]
    var1[shape=folder,label="var1"]
    var2[shape=folder,label="var2"]
    grdx[shape=note,label="grd_x"]
    grdxt[shape=note,label="grd_xt"]
    data1[shape=note,label="data"]
    data2[shape=note,label="data"]
    time1[shape=note,label="time"]
    time2[shape=note,label="time"]
    root -> Grd
    Grd -> Hgrid
    Hgrid -> grdx
    Hgrid -> grdxt
    root -> Var
    Var -> var1
    var1 -> data1
    var1 -> time1
    Var -> var2
    var2 -> data2
    var2 -> time2
  }
)

