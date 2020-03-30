# ConformalBarycenter

A Mathematica package  for

 - computing conformal barycenters of spherical point clouds;
 - computing the conformal centralizations of  spherical points clouds, 
 - closing polygonal lines while preserving the edge lengths; and     
 - computing Douady-Earle extensions of closed curves into spheres

in two and three dimensions.

Tested with Mathematica 12.0.

Copyright: (c)  2020 Jason Cantarella and Henrik Schumacher



To load the package, use the `Get` command with the full path to the file _ConformalBarycenter.m_. 

Alternatively, one may install the package by copying the file _ConformalBarycenter.m_ to the folder that Mathematica return upon evaluating

    FileNameJoin[{$UserBaseDirectory, "Applications"}]

Afterwards, the package can be loaded with

    Needs["ConformalBarycenter`"]

The files

    Example_ConformalClosure.nb
    Example_DouadyEarleExtension.nb

provide a couple of usage examples. Please make sure that those files are contained in the same folder as _ConformalBarycenter.m_ (or install the package and evaluate ``Needs["ConformalBarycenter`"]``).
