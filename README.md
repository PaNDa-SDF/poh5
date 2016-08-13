# poh5

I/O routines for poh5 (PaNDa on HDF5) file format written by C.

# What is poh5?

Poh5 file format is re-implementation of PaNDa format on HDF5.

PaNDa is an original format for NICAM datafile, which is called
'ADVANCED' in NICAM source/namelist.

Different from original PaNDa format, from version 94, you can include
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

__Add figure here__

