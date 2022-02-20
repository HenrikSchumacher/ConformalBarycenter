(* ::Package:: *)

 (* ::Package:: *)

 (* Mathematica Package  *)


 (* ::Title:: *)
 (*ConformalBarycenter.m*)


 (* ::Text:: *)
 (*A Mathematica package for *)
 (*    - computing conformal barycenters of spherical point clouds;*)
 (*    - computing the conformal centralizations of  spherical points clouds, *)
 (*    - closing polygonal lines while preserving the edge lengths; and     *)
 (*    - computing Douady-Earle extensions of closed curves into spheres*)
 (*in two and three dimensions.*)
 (**)
 (*Tested with Mathematica 12.0.*)
 (**)
 (*Copyright: (c)  2020 Jason Cantarella and Henrik Schumacher*)


 BeginPackage["ConformalBarycenter`"];


 (* ::Chapter:: *)
 (*Public*)


 ConformalBarycenter::usage=
 "ConformalBarycenter[pts] computes the conformal barycenter
 of a point cloud on \!\(\*SuperscriptBox[\(S\), \(1\)]\) or \!\(\*SuperscriptBox[\(S\), \(2\)]\).
 The input points are expected to be an array of dimensions n\[Times]2, n\[Times]3, m\[Times]n\[Times]2 or m\[Times]n\[Times]3.
      
 ConformalBarycenter[pts,weights] computes the conformal barycenter
 of a weighted point cloud on \!\(\*SuperscriptBox[\(S\), \(1\)]\) or \!\(\*SuperscriptBox[\(S\), \(2\)]\). The weights are expected
 to be an m\[Times]n array. They must be positive, but do not have to
 sum to 1.";

 ConformalClosure::usage=
 "ConformalClosure[Line[{\!\(\*
 StyleBox[\"pts\",\nFontSlant->\"Italic\"]\)}]] computes the conformal barycenter closure of the
 polygonal line with vertices given by \!\(\*
 StyleBox[\"pts\",\nFontSlant->\"Italic\"]\). The result is a Line with the same number of
 vertices. The output has the same edgelengths and the same Euclidean barycenter as the input.
 ConformalClosure[Tube[{\!\(\*
 StyleBox[\"pts\",\nFontSlant->\"Italic\"]\)},\!\(\*
 StyleBox[\"thickness\",\nFontSlant->\"Italic\"]\)]] does the same when the vertex set of is wrapped by a Tube. The result is a Tube with the same number of
 vertices and the same thickness.
 ConformalClosure[\!\(\*
 StyleBox[\"x\",\nFontSlant->\"Italic\"]\)] interprets \!\(\*
 StyleBox[\"x\",\nFontSlant->\"Italic\"]\) as the list of (not neccessarily unit) edge vectors. It returns the edge vectors \!\(\*
 StyleBox[\"y\",\nFontSlant->\"Italic\"]\) of the conformal barycenter closure. Corresponding edge vectors of \!\(\*
 StyleBox[\"x\",\nFontSlant->\"Italic\"]\) and \!\(\*
 StyleBox[\"y\",\nFontSlant->\"Italic\"]\) have the same length.
 ConformalClosure[{ \!\(\*
 StyleBox[\"x\",\nFontSlant->\"Italic\"]\), \!\(\*
 StyleBox[\"\[Omega]\",\nFontSlant->\"Italic\"]\) }] computes the conformal barycenter closure of the
 polygonal line whose unit edge vectors are given by \!\(\*
 StyleBox[\"x\",\nFontSlant->\"Italic\"]\) and whose edglengths are given by \!\(\*
 StyleBox[\"\[Omega]\",\nFontSlant->\"Italic\"]\). The result is the pair {\!\(\*
 StyleBox[\"y\",\nFontSlant->\"Italic\"]\),\!\(\*
 StyleBox[\"\[Omega]\",\nFontSlant->\"Italic\"]\)} where \!\(\*
 StyleBox[\"y\",\nFontSlant->\"Italic\"]\) is the array of new  unit edge vectors.
 ";

 $DouadyEarleData::usage="Returns the results of the last call to the backend function ConformalBaycenter`Private`cConformalBarycenter3D for later inspection when it was called from ConformalBarycenter or DouadyEarleExtension.";
 $DouadyEarleData=Association[];

 $DouadyEarleReport::usage="Returns the short summary of the results of the backend function ConformalBaycenter`Private`cConformalBarycenter3D for later inspection when it was called from ConformalBarycenter or DouadyEarleExtension.";
 $DouadyEarleReport=Association[];

 DouadyEarleDisk::usage="DouadyEarleDisk[] generates a discrete disk that can be used as first argument of DouadyEarleExtension.";

 DouadyEarleGrid::usage="DouadyEarleGrid[] generates a discrete meshlines that can be used as first argument of DouadyEarleExtension.";

 DouadyEarleExtension::usage="DouadyEarleExtension[\!\(\*
 StyleBox[\"diskpts\",\nFontSlant->\"Italic\"]\),\!\(\*
 StyleBox[\"curvepts\",\nFontSlant->\"Italic\"]\)] returns the Douady-Earle of the piecewise-linear interpolaton of the list of points \!\(\*
 StyleBox[\"curvepts\",\nFontSlant->\"Italic\"]\), evaluated at each point in the list \!\(\*
 StyleBox[\"diskpts\",\nFontSlant->\"Italic\"]\).
 DouadyEarleExtension[\!\(\*
 StyleBox[\"meshregion\",\nFontSlant->\"Italic\"]\),\!\(\*
 StyleBox[\"curvepts\",\nFontSlant->\"Italic\"]\)] does the same for all vertex coordinates in the MeshRegion \!\(\*
 StyleBox[\"meshregion\",\nFontSlant->\"Italic\"]\) and generates the output surface as a MeshRegion";


 VonMisesFisherDistribution::usage="VonMisesFisherDistribution[\[Mu],\[Kappa]] represents the spherical con Mises-Fisher distribution on the 2-sphere with mean direction \[Mu] and concentration parameter \[Kappa]."


 (* ::Chapter:: *)
 (*Private*)


 Begin["ConformalBarycenter`Private`"];


 Edgelengths[verts_?MatrixQ] := Norm /@ Differences[verts];
 Edgelengths[lineverts_Line] := Edgelengths[List @@ lineverts];
 Edgelengths[tubeverts_Tube] := Edgelengths[tubeverts[[1]]];


 (* ::Subsection:: *)
 (*ConformalBarycenter*)


 (* ::Subsubsection::Closed:: *)
 (*ConformalBarycenter*)


 Unprotect[$GiveUpThreshold];
 $GiveUpThreshold=1. 10^-14;
 Protect[$GiveUpThreshold];

 ConformalBarycenter::dim="Error: Dimensions of point cloud must be 2 or 3.";
 ConformalBarycenter::wdim="Error: Dimensions of weights is `1` does not match the dimension of the points. Dimension of weights should be  `2` or `3`.";
 ConformalBarycenter::unknwn="Warning: Unknown method `1`. Using default method.";
 ConformalBarycenter::nconv="Warning: Error estimator not reduced below tolerance `1` after `2` iterations.";

 ConformalBarycenter::nconv="Warning: Tolerance `1` not met after `2` iterations. Maximal error estimator is `3`.";

 ConformalBarycenter::noNK="Warning: Kantorovich condition not met at final configuration. Try to increase the value of the option \"MaxIterations\" or specify a more efficient \"Method\".";

 ConformalBarycenter::illcon="Warning: Smallest eigenvalue detected at final configuration is `1` while the respective residual is `2`. The problem is probably too instable to be treated in machine precision.";

 ConformalBarycenter::hugeres="Warning: Final maximal residual `1` is not below tolerance `2`. Try to increase the value of the option \"MaxIterations\" or specify a more efficient \"Method\".";

 ConformalBarycenter::tinyres="Warning: Final maximal residual `1` is below give-up tolerance `2`. Requested tolerance might be too small to be treated in machine precision.";

 Options[ConformalBarycenter]={
     "StartingVector"->Automatic
     ,"Tolerance"->1. 10^-12
     ,"Method"->{"Newton"->True,"ArmijoConstant"->0.01,"CheckKantorovichCondition"->True,"Regularization"->1.,"HyperbolicUpdate"->True}
     ,"MaxIterations"->1000
     ,"Verbose"->False
     ,"ReturnMeasures"->False
     ,"CreateReport"->True
     ,"PrintWarnings"->True
 };

 ConformalBarycenter[
     pts_?MatrixQ,
     Optional[weights_?VectorQ,Automatic],
     opts:OptionsPattern[]
 ]:=If[ListQ[#],First[#],#]&/@ConformalBarycenter[{pts},weights,opts];

 ConformalBarycenter[
     pts_?(ArrayQ[#,_,NumericQ]&&ArrayDepth[#]==3&),
     Optional[weights_?(MatrixQ[#,NumericQ]||VectorQ[#,NumericQ]||Automatic===#&),Automatic],
     opts:OptionsPattern[]
 ]:= Module[{m,n,d,\[Omega],initialguess,TOL,maxiter,method,defaultmethod,timing,data,cf,\[Epsilon]data,\[Rho]data,qdata,\[Lambda]data,q,\[Lambda],\[Rho],\[Epsilon]},
     {m,n,d}=Dimensions[pts];
     
     Switch[d,
         2,cf=cConformalBarycenter2D;,
         3,cf=cConformalBarycenter3D;,
         _,Message[ConformalBarycenter::dim];Return[$Failed];
     ];

     (*Option handling.*)
     initialguess=OptionValue["StartingVector"];
     If[initialguess===Automatic||Not[VectorQ[initialguess]&&Length[initialguess]==d],
         initialguess=ConstantArray[1.,d]
     ];
     If[weights===Automatic,
         \[Omega]=ConstantArray[1./n,n];
         ,
         If[!Dimensions[weights]=={m,n}&&!Dimensions[weights]=={n},
             Message[ConformalBarycenter::wdim,Dimensions[weights],{m,n},{n}];
             Return[$Failed];
             ,
             \[Omega]=weights/Total[weights,{-1}];
         ];
     ];

     maxiter=OptionValue["MaxIterations"];
     If[Not[NumericQ[maxiter]&&TrueQ[maxiter>=0]]&&Not[VectorQ[maxiter,NumericQ]&&Length[maxiter]!=m],
         TOL="MaxIterations"/.Options[ConformalBarycenter]
     ];
     TOL=OptionValue["Tolerance"];
     If[Not[NumericQ[TOL]&&TrueQ[TOL>=0.]],TOL="Tolerance"/.Options[ConformalBarycenter]];
     method=OptionValue["Method"];
     defaultmethod=Association["Method"/.Options[ConformalBarycenter]];
     If[StringQ[method],
         method=Switch[method,
             "Default",
             defaultmethod,
             "Newton",
             <|"Newton"->True,"ArmijoConstant"->0.01,"Regularization"->0.,"CheckKantorovichCondition"->True,"HyperbolicUpdate"->True|>,
             "RegularizedNewton",
             defaultmethod,
             "MilnorAbikoffYe"|"Milnor"|"AbikoffYe",
             <|"Newton"->False,"ArmijoConstant"->0.,"Regularization"->0.,"CheckKantorovichCondition"->True,"HyperbolicUpdate"->False|>,
             "Gradient",
             <|"Newton"->False,"ArmijoConstant"->0.125,"Regularization"->0.,"CheckKantorovichCondition"->False,"HyperbolicUpdate"->True|>,
             _,
             (Message[DouadyEarleExtension::unknwn,StringJoin["\"",method,"\""]];defaultmethod)
         ];
         ,
         If[method===Automatic,
             method=defaultmethod,
             method=Merge[{Association@method,defaultmethod},First]
         ];
     ];
     If[!NumericQ[method["Regularization"]],method["Regularization"]=defaultmethod["Regularization"]];
     If[!NumericQ[method["ArmijoConstant"]],method["ArmijoConstant"]=defaultmethod["ArmijoConstant"]];
     If[TrueQ[OptionValue["Verbose"]],Print["Running with the following options:\n",method];];

     (*The 10 percent of code that actually count.*)
     timing=AbsoluteTiming[
         data=cf[
             Sequence@@Transpose[#,RotateLeft[Range@TensorRank[#]]]&[pts],\[Omega],initialguess,
             maxiter,TOL,method["Regularization"],method["ArmijoConstant"],
             !TrueQ[!method["Newton"]],
             !TrueQ[!method["CheckKantorovichCondition"]],
             !TrueQ[!method["HyperbolicUpdate"]],
             TrueQ[OptionValue["ReturnMeasures"]]
         ];
     ][[1]];

     (*Message handling.*)
     {\[Epsilon]data,\[Rho]data,qdata,\[Lambda]data}=Transpose[data[[All,d+1;;d+4]]];
     If[OptionValue["PrintWarnings"],
         \[Epsilon]=Max[\[Epsilon]data];
         If[\[Epsilon]>TOL,
             q=Max[qdata];
             \[Lambda]=Min[\[Lambda]data];
             \[Rho]=Max[Extract[\[Rho]data,Position[\[Lambda]data,\[Lambda]]]];
             If[\[Rho]<TOL,
                 If[q<1.,
                     Message[ConformalBarycenter::nconv,TOL,maxiter,\[Epsilon]];
                     If[\[Rho]<$GiveUpThreshold,Message[ConformalBarycenter::illcon,\[Lambda],\[Rho]]];
                     ,
                     Message[ConformalBarycenter::noNK,q];
                 ],
                 If[\[Rho]<$GiveUpThreshold,
                     Message[ConformalBarycenter::tinyres,\[Rho],$GiveUpThreshold];,
                     Message[ConformalBarycenter::hugeres,\[Rho],TOL];
                 ];
             ];
         ];
     ];

     (*Constructing the return value.*)
     If[TrueQ[OptionValue["Verbose"]],Print["Computations finished."]];
     Association@@{
         "Barycenters"-> data[[All,1;;d]],
         "Timing"->timing,
         "ErrorEstimators"->\[Epsilon]data,
         "Residuals"->\[Rho]data,
         "q"->qdata,
         "\[Lambda]"->\[Lambda]data,
         "Iterations"->Round@data[[All,d+5]],
         "Backtrackings"->Round@data[[All,d+6]],
         "Eigen evaluations"->Round@data[[All,d+7]],
         "Method"->method,
         If[TrueQ[OptionValue["ReturnMeasures"]],
             "CentralizedMeasures"->ArrayReshape[data[[All,d+8;;]],{m,n,d}],
             Nothing
         ]
     }
 ];


 (* ::Subsubsection::Closed:: *)
 (*cConformalBarycenter3D*)


 Block[{t,getSmallestEigenvalue3D},

 With[{
         id=N@IdentityMatrix[3],
         eps=10. $MachineEpsilon
     },
     (* Function to compute the smallest eigenvalue of a symmetric positive-semidefinite 3x3 matrix
     A = {{a11,a12,a13},{a12,a22,a23},{a13,a23,a33}}. *)
     getSmallestEigenvalue3D=Compile[{{a11,_Real},{a12,_Real},{a13,_Real},{a22,_Real},{a23,_Real},{a33,_Real}},
         Block[{\[Lambda]1,\[Lambda]2,\[Lambda]3,p1,p2,q,p,pinv,\[Phi],r,b11,b12,b13,b22,b23,b33,\[Lambda]min},
             p1=a12^2+a13^2+a23^2;

             If[Sqrt[p1]<1. eps Sqrt[a11^2+a22^2+a33^2],
                 (* A is diagonal. *)
                 (*{\[Lambda]3,\[Lambda]2,\[Lambda]1}=Sort[{a11,a22,a33}];*)
                 \[Lambda]min=Min[a11,a22,a33];
                 ,
                 q=(a11+a22+a33)/3.;
                 p2=(a11-q)^2+(a22-q)^2+(a33-q)^2+2. p1;
                 p=Sqrt[p2/6.];
                 pinv=1./p;
                 b11=(a11-q)pinv;
                 b22=(a22-q)pinv;
                 b33=(a33-q)pinv;
                 b12=a12 pinv;
                 b13=a13 pinv;
                 b23=a23 pinv;

                 r=0.5(-b13^2 b22+2.b12 b23 b13-b11 b23^2-b12^2 b33+b11 b22 b33);

                 (* In exact arithmetic for a symmetric matrix-1\[LessEqual]r\[LessEqual]1 but computation error can leave it slightly outside this range.*)
                 If[ (r<=-1.),
                     \[Phi]=Pi/3.,
                     If[(r>=1.),
                         \[Phi]=0.,
                         \[Phi]=ArcCos[r]/3.
                     ];
                 ];
                 (* the eigenvalues satisfy eig3\[LessEqual]eig2\[LessEqual]eig1*)
                 (*\[Lambda]1=q+2. p Cos[\[Phi]];
                 \[Lambda]3=q+2.p Cos[\[Phi]+2.Pi/3.];
                 \[Lambda]2=3.q-\[Lambda]1-\[Lambda]3; *)
                 \[Lambda]min=q+2.p Cos[\[Phi]+2.Pi/3.];
             ];
             \[Lambda]min
         ]
     ];
 ];

 With[{
         eps=$MachineEpsilon,
         giveupTOL=$GiveUpThreshold,
         smallone=1.-16$MachineEpsilon,
         bigone=1.+16$MachineEpsilon,
         infty=$MaxMachineNumber,
         gfactor=4.,
         Tanhsmall=N[HornerForm/@PadeApproximant[Tanh[t]/t,{t,0,{8,9}}]],
         get\[Lambda]min=getSmallestEigenvalue3D,
         normalizationthreshold=(0.99^2+10$MachineEpsilon)
     },

     cConformalBarycenter3D=Compile[{
             {xinit1,_Real,1},{xinit2,_Real,1},{xinit3,_Real,1},{\[Omega],_Real,1},{winit,_Real,1},
             {maxiter,_Integer},
             {TOL,_Real},{regularization,_Real},{ArmijoC,_Real},
             {NewtonQ,True|False},{KantorovichQ,True|False},
             {HyperbolicUpdateQ,True|False},{ReturnMeasureQ,True|False}
         },
         Block[{x1,x2,x3,z1,z2,z3,zz,zx2,wz2,w1,w2,w3,ww,wx2,w1new,w2new,w3new,
                 u1,u2,u3,F1,F2,F3,V1,V2,V3,\[Gamma],maxbiter,biter,metadata,invnormx,linesearchQ,
                 a,b,c,aa,bb,cc,dd,iter,reg,residual,
                 scale,normu,\[Tau],t,A11,A12,A13,A22,A23,A33,detAinv,
                 \[Phi]\[Tau],D\[Phi]0,res2,q,\[Lambda]min,errorestimator,backtrackingcounter,eigcounter,succeeded,continue
             },
             q=1.;
             normu=0.;
             \[Lambda]min=eps;
             linesearchQ=ArmijoC>0;
             continue=True;
             succeeded=False;
             errorestimator=infty;
             residual=infty;

             eigcounter=0;

             (*Initialization of matrix parameters.*)
             detAinv=1.;
             A11=A22=A33=1.;
             A12=A13=A23=0.;

             (*Initialization of line search parameters.*)
             \[Gamma]=0.1;
             maxbiter=12;
             backtrackingcounter=0;

             (*Initialization of current iterate w.*)
             w1=Compile`GetElement[winit,1];
             w2=Compile`GetElement[winit,2];
             w3=Compile`GetElement[winit,3];
             (*Use Euclidean barycenter as initial guess if the supplied initial guess does not make sense.*)
             If[w1^2+w2^2+w3^2>=smallone,
                 w1=\[Omega] . xinit1;
                 w2=\[Omega] . xinit2;
                 w3=\[Omega] . xinit3;
             ];

             (*Shift measure along w.*)
             wx2=(2.w1)xinit1+(2.w2)xinit2+(2.w3)xinit3;
             ww=w1^2+w2^2+w3^2;
             cc=1./Subtract[bigone+ww,wx2];
             aa=bigone-ww;
             bb=wx2-2.;
             x1=Times[aa xinit1+bb w1,cc];
             x2=Times[aa xinit2+bb w2,cc];
             x3=Times[aa xinit3+bb w3,cc];
             (*If w lies close to the boundary of the ball, then normalizing x is a good idea.*)
             If[ww>normalizationthreshold,
                 invnormx=1./Sqrt[x1^2+x2^2+x3^2];
                 x1*=invnormx;
                 x2*=invnormx;
                 x3*=invnormx;
             ];

             (*Start optimization loop.*)
             iter=0;
             While[(
                 V1=0.5x1;V2=0.5x2;V3=0.5x3;
                 F1=\[Omega] . V1;F2=\[Omega] . V2;F3=\[Omega] . V3;
                 res2=gfactor (F1^2+F2^2+F3^2);
                 residual=Sqrt[res2];
                 If[NewtonQ,
                     (*-\[Del]F(0) with respect to measure y*)
                     A11=-gfactor \[Omega] . (V1 V1)+1.;
                     A22=-gfactor \[Omega] . (V2 V2)+1.;
                     A33=-gfactor \[Omega] . (V3 V3)+1.;
                     A12=-gfactor \[Omega] . (V1 V2);
                     A13=-gfactor \[Omega] . (V1 V3);
                     A23=-gfactor \[Omega] . (V2 V3);
                 ];
                 If[KantorovichQ,
                     If[residual<=TOL||NewtonQ,
                         (*Performing Kantorovich condition check only when have to... or if we really want it.*)
                         If[!NewtonQ,
                             (*-\[Del]F(0) with respect to measure y*)
                             A11=-gfactor \[Omega] . (V1 V1)+1.;
                             A22=-gfactor \[Omega] . (V2 V2)+1.;
                             A33=-gfactor \[Omega] . (V3 V3)+1.;
                             A12=-gfactor \[Omega] . (V1 V2);
                             A13=-gfactor \[Omega] . (V1 V3);
                             A23=-gfactor \[Omega] . (V2 V3);
                         ];
                         ++eigcounter;
                         \[Lambda]min=get\[Lambda]min[A11,A12,A13,A22,A23,A33];
                         q=4.residual/\[Lambda]min^2;
                         If[(q<1),
                             (*"Kantorovich condition satisfied; this allows to compute an error estimator."*)
                             errorestimator=0.5\[Lambda]min q;
                             (*"And we should deactivate line search. Otherwise, we may run into precision issues."*)
                             linesearchQ=False;
                             continue=(errorestimator >TOL);
                             succeeded=Not[continue]
                         ,
                             errorestimator=infty;
                             linesearchQ=ArmijoC>0;
                             (*There is no way to reduce the residual below machine epsilon.*)
                             (*If the algorithm reaches here, the problem is probably too ill-conditioned to be solved in machine precision.*)
                             continue=residual>giveupTOL;
                         ];
                     ,
                         (*Skipped Kantorovich condition check because residual was too high.*)
                         linesearchQ=ArmijoC>0;
                         q=1.;
                         \[Lambda]min=eps;
                         errorestimator=infty;
                         continue=(residual>Max[giveupTOL,TOL])
                     ];
                 ,
                     (*"Kantorovich condition check deactivated."*)
                     linesearchQ=ArmijoC>0;
                     q=1.;
                     \[Lambda]min=eps;
                     errorestimator=infty;
                     continue=(residual>Max[giveupTOL,TOL])
                 ];

                 (*Condition for While*)
                 continue&&(iter<maxiter)
             ),
                 (*Executing body of While*)
                 iter++;
                 (*Compute update direction u.*)
                 If[NewtonQ,
                     (*Compute (regularized) Newton search direction.*)
                     reg=regularization res2;
                     (*reg=regularization residual;*)
                     A11+=reg;A22+=reg;A33+=reg;
                     (*a very crude way of solving A.u = F*)
                     detAinv=1./(-A13 A13 A22+2. A12 A13 A23-A11 A23 A23-A12 A12 A33+A11 A22 A33);
                     u1=detAinv(-A23 A23 F1+A22 A33 F1+A13 A23 F2-A12 A33 F2-A13 A22 F3+A12 A23 F3);u2=detAinv(A13 A23 F1-A12 A33 F1-A13 A13 F2+A11 A33 F2+A12 A13 F3-A11 A23 F3);u3=detAinv(-A13 A22 F1+A12 A23 F1+A12 A13 F2-A11 A23 F2-A12 A12 F3+A11 A22 F3);
                 ,
                     (*Compute gradient descent direction.*)
                     If[HyperbolicUpdateQ,
                         u1=F1/gfactor;
                         u2=F2/gfactor;
                         u3=F3/gfactor;
                     ,
                         (*The factor 2 is here to reproduce the Abikoff-Ye algorithm (in the absence of linesearch).*)
                         u1=2.F1;
                         u2=2.F2;
                         u3=2.F3;
                     ];
                 ];
                 \[Tau]=1.;
                 normu=Sqrt[u1^2+u2^2+u3^2];
                 scale=If[HyperbolicUpdateQ,
                     (*Riemannian exponential at the origin.*)
                     t=\[Tau] normu;
                     \[Tau] If[t<=1.,Tanhsmall,If[t<=50.,Tanh[t]/t,1./t]]
                 ,
                     (*Euclidean exponential at the origin.*)
                     \[Tau]
                 ];
                 z1=scale u1;z2=scale u2;z3=scale u3;
                 If[linesearchQ,
                     (*Linesearch with potential as merit function.*)
                     D\[Phi]0=-gfactor (F1 u1+F2 u2+F3 u3);
                     biter=0;
                     While[(
                         (*Compute potential and check Armijo condition.*)
                         zz=z1^2+z2^2+z3^2;
                         \[Phi]\[Tau]=\[Omega] . (Log[Abs[Subtract[bigone+zz,(2.z1)x1+(2.z2)x2+(2.z3)x3]/(bigone-zz)]]);
                         (*\[Phi]0=0.;*)
                         (\[Phi]\[Tau](*-\[Phi]0*)-ArmijoC \[Tau] D\[Phi]0)>0&&biter<maxbiter
                     ),
                         ++biter;
                         (*Utilizing quadratic interpolation of merit function,*)
                         \[Tau]=Max[\[Gamma] \[Tau],-0.5 D\[Phi]0 \[Tau]^2/(\[Phi]\[Tau](*-\[Phi]0*)-D\[Phi]0 \[Tau]) ];
                         (*Update*)
                         scale=If[HyperbolicUpdateQ,
                             (*Riemannian exponential at the origin.*)
                             t=\[Tau] normu;
                             \[Tau] If[t<=1.,Tanhsmall,If[t<=50.,Tanh[t]/t,1./t]]
                         ,
                             (*Euclidean exponential at the origin.*)
                             \[Tau]
                         ];
                         z1=scale u1;z2=scale u2;z3=scale u3;
                     ];
                     backtrackingcounter+=biter;
                 ];
                 (*Shift the point z along -w to get new updated point w.*)
                 wz2=(2.w1)z1+(2.w2)z2+(2.w3)z3;
                 ww=w1^2+w2^2+w3^2;
                 zz=z1^2+z2^2+z3^2;
                 c=1./(bigone+ww zz+wz2);
                 a=1.-ww;
                 b=1.+zz+wz2;
                 w1=Times[a z1+b w1,c];
                 w2=Times[a z2+b w2,c];
                 w3=Times[a z3+b w3,c];

                 (*Shift the input measure along w to 0 to simplify gradient, Hessian, and update computation.*)
                 wx2=(2.w1)xinit1+(2.w2)xinit2+(2.w3)xinit3;
                 ww=w1^2+w2^2+w3^2;
                 cc=1./Subtract[bigone+ww,wx2];
                 aa=bigone-ww;
                 bb=wx2-2.;
                 x1=Times[aa xinit1+bb w1,cc];
                 x2=Times[aa xinit2+bb w2,cc];
                 x3=Times[aa xinit3+bb w3,cc];
                 (*If w lies close to the boundary of the ball, then normalizing x is a good idea.*)
                 If[ww>normalizationthreshold,
                     invnormx=1./Sqrt[x1^2+x2^2+x3^2];
                     x1*=invnormx;
                     x2*=invnormx;
                     x3*=invnormx;
                 ];
             ];
             (*Optimization loop terminated.*)

             (*Compute error estimator in the case that the algorithm has not succeeded...*)
             If[Not[succeeded]||iter>=maxiter,
                 ++eigcounter;
                 (*-\[Del]F(0) with respect to measure y*)
                 V1=0.5x1;V2=0.5x2;V3=0.5x3;
                 F1=\[Omega] . V1;F2=\[Omega] . V2;F3=\[Omega] . V3;
                 res2=gfactor (F1^2+F2^2+F3^2);
                 residual=Sqrt[res2];
                 A11=-gfactor \[Omega] . (V1 V1)+1.;
                 A22=-gfactor \[Omega] . (V2 V2)+1.;
                 A33=-gfactor \[Omega] . (V3 V3)+1.;
                 A12=-gfactor \[Omega] . (V1 V2);
                 A13=-gfactor \[Omega] . (V1 V3);
                 A23=-gfactor \[Omega] . (V2 V3);
                 \[Lambda]min=get\[Lambda]min[A11,A12,A13,A22,A23,A33];
                 q=4.residual/\[Lambda]min^2;
                 errorestimator=If[(q<1),0.5\[Lambda]min q,infty];
             ];

             metadata={
                 w1,w2,w3,errorestimator,residual,q,\[Lambda]min,N[iter],N[backtrackingcounter],N[eigcounter]
             };
             If[!ReturnMeasureQ,
                 (*Just return the conformal barycenter.*)
                 metadata
             ,
                 (*Return also centralized measure.*)
                 Join[metadata,Flatten[Transpose[{x1,x2,x3}]]]
             ]
         ],
         CompilationTarget->"C",
         RuntimeAttributes->{Listable},
         Parallelization->True,
         RuntimeOptions->"Speed",
         CompilationOptions->{"InlineCompiledFunctions" -> True}
     ];
 ]
 ]



 (* ::Subsubsection::Closed:: *)
 (*cConformalBarycenter2D*)


 Block[{t},

 With[{
     eps=$MachineEpsilon,
     giveupTOL=$GiveUpThreshold,
     smallone=1.-16$MachineEpsilon,
     bigone=1.+16$MachineEpsilon,
     infty=$MaxMachineNumber,
     gfactor=4.,
     Tanhsmall=N[HornerForm/@PadeApproximant[Tanh[t]/t,{t,0,{8,9}}]],
     normalizationthreshold=(0.99^2+10$MachineEpsilon)
 },

 cConformalBarycenter2D=Compile[{
         {xinit1,_Real,1},{xinit2,_Real,1},{\[Omega],_Real,1},{winit,_Real,1},
         {maxiter,_Integer},
         {TOL,_Real},{regularization,_Real},{ArmijoC,_Real},
         {NewtonQ,True|False},{KantorovichQ,True|False},
         {HyperbolicUpdateQ,True|False},{ReturnMeasureQ,True|False}
     },
     Block[{x1,x2,z1,z2,zz,zx2,wz2,w1,w2,ww,wx2,w1new,w2new,
             u1,u2,F1,F2,V1,V2,\[Gamma],maxbiter,biter,metadata,invnormx,linesearchQ,
             a,b,c,aa,bb,cc,dd,iter,reg,residual,
             scale,normu,\[Tau],t,A11,A12,A22,detAinv,
             \[Phi]\[Tau],D\[Phi]0,res2,q,\[Lambda]min,errorestimator,backtrackingcounter,eigcounter,succeeded,continue
         },
         q=1.;
         normu=0.;
         \[Lambda]min=eps;
         linesearchQ=ArmijoC>0;
         continue=True;
         succeeded=False;
         errorestimator=infty;
         residual=infty;

         eigcounter=0;

         (*Initialization of matrix parameters.*)
         detAinv=1.;
         A11=A22=1.;
         A12=0.;

         (*Initialization of line search parameters.*)
         \[Gamma]=0.1;
         maxbiter=12;
         backtrackingcounter=0;

         (*Initialize the currect iterate w.*)
         w1=Compile`GetElement[winit,1];
         w2=Compile`GetElement[winit,2];
         (*Use Euclidean barycenter as initial guess if the supplied initial guess does not make sense.*)
         If[w1^2+w2^2>=smallone,
             w1=\[Omega] . xinit1;
             w2=\[Omega] . xinit2;
         ];

         (*Shift measure along w.*)
         wx2=(2.w1)xinit1+(2.w2)xinit2;
         ww=w1^2+w2^2;
         cc=1./Subtract[bigone+ww,wx2];
         aa=bigone-ww;
         bb=wx2-2.;
         x1=Times[aa xinit1+bb w1,cc];
         x2=Times[aa xinit2+bb w2,cc];
         (*If w lies close to the boundary of the ball, then normalizing x is a good idea.*)
         If[ww>normalizationthreshold,
             invnormx=1./Sqrt[x1^2+x2^2];
             x1*=invnormx;
             x2*=invnormx;
         ];

         (*Start optimization loop.*)
         iter=0;
         While[(
             V1=0.5x1;V2=0.5x2;
             F1=\[Omega] . V1;F2=\[Omega] . V2;
             res2=gfactor (F1^2+F2^2);
             residual=Sqrt[res2];
             If[NewtonQ,
                 (*-\[Del]F(0) with respect to measure y*)
                 A11=-gfactor \[Omega] . (V1 V1)+bigone;
                 A22=-gfactor \[Omega] . (V2 V2)+bigone;
                 A12=-gfactor \[Omega] . (V1 V2);
             ];
             If[KantorovichQ,
                 If[residual<=TOL||NewtonQ,
                     (*Performing Kantorovich condition check only if we have to... or if we really want it.*)
                     If[!NewtonQ,
                         (*-\[Del]F(0) with respect to measure y*)
                         A11=-gfactor \[Omega] . (V1 V1)+bigone;
                         A22=-gfactor \[Omega] . (V2 V2)+bigone;
                         A12=-gfactor \[Omega] . (V1 V2);
                     ];
                     ++eigcounter;
                     \[Lambda]min=0.5 (A11+A22-Sqrt[Abs[A11^2+4. A12^2 - 2. A11 A22+A22^2]]);
                     q=4.residual/\[Lambda]min^2;
                     If[(q<1),
                         (*"Kantorovich condition satisfied; this allows to compute an error estimator."*)
                         errorestimator=0.5\[Lambda]min q;
                         (*"And we should deactivate line search. Otherwise, we may run into precision issues."*)
                         linesearchQ=False;
                         continue=(errorestimator >TOL);
                         succeeded=Not[continue]
                     ,
                         linesearchQ=ArmijoC>0;
                         errorestimator=infty;
                         (*There is no way to reduce the residual below machine epsilon.*)
                         (*If the algorithm reaches here, the problem is probably too ill-conditioned to be solved in machine precision.*)
                         continue=residual>giveupTOL;
                     ];
                 ,
                     (*Skipped Kantorovich condition check because residual was too high.*)
                     linesearchQ=ArmijoC>0;
                     q=1.;
                     \[Lambda]min=eps;
                     errorestimator=infty;
                     continue=(residual>Max[giveupTOL,TOL])
                 ];
             ,
                 (*"Kantorovich condition check deactivated."*)
                 linesearchQ=ArmijoC>0;
                 q=1.;
                 \[Lambda]min=eps;
                 errorestimator=infty;
                 continue=(residual>Max[giveupTOL,TOL])
             ];

             (*Condition for While*)
             continue&&(iter<maxiter)
         ),
             (*Executing body of While*)
             iter++;
             (*Compute update direction u.*)
             If[NewtonQ,
                 (*Compute (regularized) Newton search direction.*)
                 reg=regularization res2;
                 (*reg=regularization residual;*)
                 A11+=reg;A22+=reg;
                 (*a very crude way of solving A.u = F*)
                 detAinv=1./(A11 A22-A12^2);
                 u1=detAinv(A22 F1-A12 F2);
                 u2=detAinv(-A12 F1+A11 F2);
             ,
                 (*Compute gradient descent direction.*)
                 If[HyperbolicUpdateQ,
                     u1=F1/gfactor;
                     u2=F2/gfactor;
                 ,
                     (*The factor 2 is here to reproduce the Abikoff-Ye algorithm (in the absence of linesearch).*)
                     u1=2.F1;
                     u2=2.F2;
                 ];
             ];
             \[Tau]=1.;
             normu=Sqrt[u1^2+u2^2];
             scale=If[HyperbolicUpdateQ,
                 (*Riemannian exponential at the origin.*)
                 t=\[Tau] normu;
                 \[Tau] If[t<=1.,Tanhsmall,If[t<=50.,Tanh[t]/t,1./t]]
             ,
                 (*Euclidean exponential at the origin.*)
                 \[Tau]
             ];
             z1=scale u1;z2=scale u2;
             If[linesearchQ,
                 (*Linesearch with potential as merit function.*)
                 D\[Phi]0=-gfactor (F1 u1+F2 u2);
                 biter=0;
                 While[(
                     (*Compute potential and check Armijo condition.*)
                     zz=z1^2+z2^2;
                     \[Phi]\[Tau]=\[Omega] . (Log[Abs[Subtract[bigone+zz,(2.z1)x1+(2.z2)x2]/(bigone-zz)]]);
                     \[Phi]\[Tau]>(*-\[Phi]0*)-ArmijoC \[Tau] D\[Phi]0&&biter<maxbiter
                 ),
                     ++biter;
                     (*Utilizing quadratic interpolation of merit function,*)
                     \[Tau]=Max[\[Gamma] \[Tau],-0.5 D\[Phi]0 \[Tau]^2/(\[Phi]\[Tau](*-\[Phi]0*)-D\[Phi]0 \[Tau]) ];
                     (*Update*)
                     scale=If[HyperbolicUpdateQ,
                         (*Riemannian exponential at the origin.*)
                         t=\[Tau] normu;
                         \[Tau] If[t<=1.,Tanhsmall,If[t<=50.,Tanh[t]/t,1./t]]
                     ,
                         (*Euclidean exponential at the origin.*)
                         \[Tau]
                     ];
                     z1=scale u1;z2=scale u2;
                 ];
                 backtrackingcounter+=biter;
             ];

             (*Shift the point z along -w to get new updated point w.*)
             wz2=(2.w1)z1+(2.w2)z2;
             ww=w1^2+w2^2;
             zz=z1^2+z2^2;
             c=1./(bigone+ww zz+wz2);
             a=1.-ww;
             b=1.+zz+wz2;
             w1=Times[a z1+b w1,c];
             w2=Times[a z2+b w2,c];

             (*Shift the input measure along w to 0 to simplify gradient, Hessian, and update computation.*)
             wx2=(2.w1)xinit1+(2.w2)xinit2;
             ww=w1^2+w2^2;
             cc=1./Subtract[bigone+ww,wx2];
             aa=bigone-ww;
             bb=wx2-2.;
             x1=Times[aa xinit1+bb w1,cc];
             x2=Times[aa xinit2+bb w2,cc];
             (*If w lies close to the boundary of the ball, then normalizing x is a good idea.*)
             If[ww>normalizationthreshold,
                 invnormx=1./Sqrt[x1^2+x2^2];
                 x1*=invnormx;
                 x2*=invnormx;
             ];
         ];
         (*Optimization loop terminated.*)
         (*Compute error estimator in the case that the algorithm has not succeeded...*)
         If[Not[succeeded]||iter>=maxiter,
             ++eigcounter;
             (*-\[Del]F(0) with respect to measure y*)
             V1=0.5x1;V2=0.5x2;
             F1=\[Omega] . V1;F2=\[Omega] . V2;
             res2=gfactor (F1^2+F2^2);
             residual=Sqrt[res2];
             A11=-gfactor \[Omega] . (V1 V1)+1.;
             A22=-gfactor \[Omega] . (V2 V2)+1.;
             A12=-gfactor \[Omega] . (V1 V2);
             \[Lambda]min=0.5 (A11+A22-Sqrt[Abs[A11^2+4. A12^2 - 2. A11 A22+A22^2]]);
             q=4.residual/\[Lambda]min^2;
             errorestimator=If[(q<1),0.5\[Lambda]min q,infty];
         ];

         metadata={
             w1,w2,errorestimator,residual,q,\[Lambda]min,N[iter],N[backtrackingcounter],N[eigcounter]
         };
         If[!ReturnMeasureQ,
             (*Just return the conformal barycenter.*)
             metadata
         ,
             (*Return also centralized measure.*)
             Join[metadata,Flatten[Transpose[{x1,x2}]]]
         ]
     ],
     CompilationTarget->"C",
     RuntimeAttributes->{Listable},
     Parallelization->True,
     RuntimeOptions->"Speed",
     CompilationOptions->{"InlineCompiledFunctions" -> True}
     ]
 ];
 ]


 (* ::Subsection:: *)
 (*ConformalClosure*)


 (* ::Subsubsection::Closed:: *)
 (*ConformalClosure*)


 Options[ConformalClosure]:=Options[ConformalBarycenter];

 ConformalClosure[Line[v_?(MatrixQ[#]||(ArrayQ[#]&&TensorRank[#]===3)&)],opts:OptionsPattern[]]:=Module[{len,x,y},
     x=cUnitEdgeVectors[v];
     len=cEdgeLengths[v];
     y=ConformalBarycenter[x,len(*/Total[len,{-1}]*),"ReturnMeasures"->True,opts]["CentralizedMeasures"];
     Line[cFromEdges[cMean[v,len],y,len,True]]
 ];

 ConformalClosure[Tube[v_?(MatrixQ[#]||(ArrayQ[#]&&TensorRank[#]===3)&),\[Theta]_],opts:OptionsPattern[]]:=Module[{len,x,y},
     x=cUnitEdgeVectors[v];
     len=cEdgeLengths[v];
     y=ConformalBarycenter[x,len(*/Total[len,{-1}]*),"ReturnMeasures"->True,opts]["CentralizedMeasures"];
     Tube[cFromEdges[cMean[v,len],y,len,True],\[Theta]]
 ];

 ConformalClosure[{x_?(MatrixQ[#]||(ArrayQ[#]&&TensorRank[#]===3)&),len_?ArrayQ},opts:OptionsPattern[]]:={
     ConformalBarycenter[Divide[x,Sqrt[Dot[x^2,ConstantArray[1.,Dimensions[x][[-1]]]]]],len(*/Total[len,{-1}]*),"ReturnMeasures"->True,opts]["CentralizedMeasures"],len
 };

 ConformalClosure[x_?(MatrixQ[#]||(ArrayQ[#]&&TensorRank[#]===3)&),opts:OptionsPattern[]]:=Module[{len},
     len=Sqrt[Dot[x^2,ConstantArray[1.,Dimensions[x][[-1]]]]];
     ConformalBarycenter[Divide[x,len],len(*/Total[len,{-1}]*),"ReturnMeasures"->True,opts]["CentralizedMeasures"]len
 ];


 (* ::Subsubsection::Closed:: *)
 (*Vertex-Edge-Conversion*)


 cFromEdges := cFromEdges=Compile[{{p0,_Real,1},{x,_Real,2},{len,_Real,1},{centeredQ,True|False}},
     Block[{p,c,n,\[Eta],Linv},
         c=p0;
         n=Length[len];
         Linv=1./Total[len];
         p=Table[0.,{n+1},{Length[p0]}];
         If[centeredQ,
             \[Eta]=Table[0.,n];
             \[Eta][[;;-1]]+=0.5len;
             \[Eta][[1;;-2]]+=0.5len[[2;;-1]];
             Do[\[Eta][[i]]+=\[Eta][[i+1]],{i,n-1,1,-1}];
             c-=((\[Eta] len) . x) Linv^2;
         ];
         p[[1]]=c;
         Do[
             c=c+(Compile`GetElement[len,i] Linv) Compile`GetElement[x,i];
             p[[i+1]]=c,
             {i,1,n}
         ];
         p
     ],
     CompilationTarget->"C",
     RuntimeAttributes->{Listable},
     Parallelization->True,
     RuntimeOptions->"Speed"
 ];

 cUnitEdgeVectors := cUnitEdgeVectors = Compile[{{v,_Real,2}},
     # Power[Dot[#^2,Table[1.,{Dimensions[v][[-1]]}]],-0.5]&[Differences[v]],
     CompilationTarget->"C",
     RuntimeAttributes->{Listable},
     Parallelization->True,
     RuntimeOptions->"Speed"
 ];

 cEdgeLengths := cEdgeLengths = Compile[{{v,_Real,2}},
     Power[Dot[#^2,Table[1.,{Dimensions[v][[-1]]}]],0.5]&[Differences[v]],
     CompilationTarget->"C",
     RuntimeAttributes->{Listable},
     Parallelization->True,
     RuntimeOptions->"Speed"
 ];

 cMean := cMean = Compile[{{v,_Real,2},{len,_Real,1}},
     Block[{m,Linv},
         Linv=1./Total[len];
         m=Table[0.,{Dimensions[v][[2]]}];
         Do[m+=(Compile`GetElement[v,i]+Compile`GetElement[v,i+1])(0.5Compile`GetElement[len,i]),{i,1,Length[len]}];
         m Linv
     ],
     CompilationTarget->"C",
     RuntimeAttributes->{Listable},
     Parallelization->True,
     RuntimeOptions->"Speed"
 ];


 (* ::Subsection:: *)
 (*DouadyEarleExtension*)


 (* ::Subsubsection::Closed:: *)
 (*DouadyEarleExtension*)


 DouadyEarleExtension::dim="Error: Dimensions of curve points must be 2 or 3.";
 DouadyEarleExtension::unknwn="Warning: Unknown method `1`. Using default method.";
 DouadyEarleExtension::nconv="Warning: Error estimator not reduced below tolerance `1` after `2` iterations.";

 Options[DouadyEarleExtension]={
     "StartingVector"->Automatic
     ,"QuadraturePoints"->720
     ,"Tolerance"->1. 10^-8
     ,"Method"->{"Newton"->True,"ArmijoConstant"->0.25,"CheckKantorovichCondition"->True,"Regularization"->1.,"HyperbolicUpdate"->True}
     ,"MaxIterations"->1000
     ,"Verbose"->False
     ,"StoreData"->True
     ,"CreateReport"->True
 };

 DouadyEarleExtension[M_MeshRegion,curvepts_?MatrixQ,opts:OptionsPattern[]]:= MeshRegion[
     DouadyEarleExtension[MeshCoordinates[M],curvepts,opts],
     MeshCells[M,All,"Multicells"->True],
     Sequence@@Options[M]
 ];

 DouadyEarleExtension[diskpts_?MatrixQ,curvepts_?MatrixQ,OptionsPattern[]]:=
 Module[{xt,d,bndplist,intplist,initialguess,maxiter,TOL,regularization,data,pts,m,method,defaultmethod,timing,\[Omega],cf,c\[Gamma],r2,boole},
     d=Dimensions[curvepts][[2]];
     Switch[d,
         2,c\[Gamma]=cDouadyEarleBoundaryCurve2D;cf=cDouadyEarleExtension2D;,
         3,c\[Gamma]=cDouadyEarleBoundaryCurve3D;cf=cDouadyEarleExtension3D;,
         _,(Message[DouadyEarleExtension::dim];Return[$Failed];)
     ];

     m=OptionValue["QuadraturePoints"];
     If[Not[IntegerQ[m]&&TrueQ[m>0]],m="QuadraturePoints"/.Options[DouadyEarleExtension]];
     initialguess=OptionValue["StartingVector"];
     If[initialguess===Automatic||Not[VectorQ[initialguess]&&Length[initialguess]==d],
         initialguess=ConstantArray[1.,d]
     ];

     maxiter=OptionValue["MaxIterations"];
     If[Not[NumericQ[maxiter]&&TrueQ[maxiter>=0]],TOL="MaxIterations"/.Options[DouadyEarleExtension]];
     TOL=OptionValue["Tolerance"];
     If[Not[NumericQ[TOL]&&TrueQ[TOL>=0.]],TOL="Tolerance"/.Options[DouadyEarleExtension]];
     method=OptionValue["Method"];
     defaultmethod=Association["Method"/.Options[DouadyEarleExtension]];
     If[StringQ[method],
         method=Switch[method,
             "Newton",
             <|"Newton"->True,"ArmijoConstant"->0.25,"Regularization"->0.,"CheckKantorovichCondition"->True,"HyperbolicUpdate"->True|>,
             "RegularziedNewton",defaultmethod,
             "MilnorAbikoffYe"|"Milnor"|"AbikoffYe",
             <|"Newton"->False,"ArmijoConstant"->0.,"Regularization"->0.,"CheckKantorovichCondition"->True,"HyperbolicUpdate"->False|>,
             "Gradient",
             <|"Newton"->False,"ArmijoConstant"->0.25,"Regularization"->0.,"CheckKantorovichCondition"->False,"HyperbolicUpdate"->True|>,
             _,(Message[DouadyEarleExtension::unknwn];defaultmethod)
         ];
     ,
         If[method===Automatic,
             method=defaultmethod
         ,
             method=Merge[{Association@method,defaultmethod},First]
         ];
     ];
     regularization=If[NumericQ[method["Regularization"]],
         method["Regularization"],
         defaultmethod["Regularization"]
     ];
     If[TrueQ[OptionValue["Verbose"]],Print["Running DouadyEarleExtension with the following options:\n",method];];

     pts=ConstantArray[0.,{Length[diskpts],d}];
     r2=Total[diskpts^2,{2}];
     boole=UnitStep[r2-(1.-100$MachineEpsilon)];
     bndplist=Random`Private`PositionsOf[boole,1];
     intplist=Random`Private`PositionsOf[boole,0];

     (*Scramble the points so that the workload is better distributed over the compute kernels.*)
     intplist=intplist[[RandomSample[1;;Length[intplist]]]];

     (*Apply the boundary curve to all boundary points.*)
     If[TrueQ[OptionValue["Verbose"]],Print["Computing boundary curve..."];];
     If[
         Length[bndplist]>0,
         pts[[bndplist]]=c\[Gamma][
             Sequence@@Transpose[curvepts],
             ArcTan@@Transpose[diskpts][[1;;2,bndplist]]
         ];
     ];

     (*Apply the actual Douady-Earle extension to all interior points.*)
     If[TrueQ[OptionValue["Verbose"]],Print["Computing interior surface..."];];

     timing=AbsoluteTiming[
         If[
             Length[intplist]>0,
             $ConformalBarycenterData=data=cf[
                 Sequence@@Transpose[curvepts],
                 diskpts[[intplist]],
                 Sequence@@Transpose[ReIm[Exp[I Most@Subdivide[0.,2.Pi,m]]]],
                 ConstantArray[1./m,m],
                 initialguess,
                 maxiter,
                 TOL,regularization,method["ArmijoConstant"],
                 method["Newton"],method["CheckKantorovichCondition"],method["HyperbolicUpdate"]
             ];
             pts[[intplist]]=data[[All,1;;d]];
         ];
     ][[1]];
     $DouadyEarleReport=With[{
         f=x|->(Through[<|
             "Max"->Max[#]&,
             "Mean"->Mean[N[#]]&,
             "Median"->Median[#/.\[Infinity]->$MaxMachineNumber]&|>[x]
         ]),
         nsucceeded=(Length[data]-Total[UnitStep[data[[All,d+1]]-TOL]])
         }
     ,
         If[
             nsucceeded<Length[data],
             Message[DouadyEarleExtension::nconv,TOL,maxiter]
         ];
         Grid[{{
             Dataset@Association[
                 "Succeeded"->nsucceeded>=Length[data],
                 "Time elapsed"->timing,
                 "Points processed"->Length[data],
                 "Succeeded points"->nsucceeded,
                 "Success rate"->N[nsucceeded]/Length[data]
             ]
         ,
             Dataset[Association[
                 "Iterations"->f[Round@data[[All,d+5]]],
                 "Error estimator"->f[data[[All,d+1]]/.$MaxMachineNumber->\[Infinity]],
                 "Residual"->f[data[[All,d+2]]],
                 "q"->f[data[[All,d+3]]],
                 "\[Lambda]min"->f[1/data[[All,d+4]]],
                 "Backtrackings"->f[Round@data[[All,d+5]]],
                 "Eigen evaluations"->f[Round@data[[All,d+6]]]
             ]]
         }}]
     ];
     If[TrueQ[OptionValue["Verbose"]],Print["Computations finished. Printing report:"]];
     If[TrueQ[OptionValue["Verbose"]],Print[$DouadyEarleReport]];
     Developer`ToPackedArray[pts]
 ];


 (* ::Subsubsection::Closed:: *)
 (*cDoudaryEarleExtension2D*)



 cDouadyEarleExtension2D := cDouadyEarleExtension2D = With[{
     bigone=1.+16$MachineEpsilon,
     twoPi=2.Pi,
     getcenter=cConformalBarycenter2D
     },
     Compile[{
         {c1,_Real,1},{c2,_Real,1},
         {w,_Real,1},
         {z01,_Real,1},{z02,_Real,1},{\[Omega],_Real,1},
         {winit,_Real,1},
         {maxiter,_Integer},
         {TOL,_Real},{regularization,_Real},{ArmijoC,_Real},
         {NewtonQ,True|False},{KantorovichQ,True|False},{HyperbolicUpdateQ,True|False}
         },
         Block[{m,w1,w2,ww,wz2,aa,bb,cc,\[Theta],y1,y2,\[Lambda],hinv,T,idx},
             m=Min[Length[c1],Length[c2]];
             (*Transform the boundary points z with the inverse of the shift that moves w to the origin.*)
             w1=-Compile`GetElement[w,1];
             w2=-Compile`GetElement[w,2];
             wz2=(2.w1)z01+(2.w2)z02;
             ww=w1^2+w2^2;
             cc=1./Subtract[bigone+ww,wz2];
             aa=bigone-ww;
             bb=wz2-2.;
             (*Transform to angular coordinates.*)
             \[Theta]=Mod[ArcTan[Times[aa z01+bb w1,cc],Times[aa z02+bb w2,cc]],twoPi];

             (*employing piecewise-linear interpolation*)
             hinv=m/twoPi;
             T=Mod[\[Theta],twoPi]hinv;
             idx=Mod[Floor[T]+1,m,1];

             \[Lambda]=Mod[T,1.];
             y1=Subtract[1.,\[Lambda]]c1[[idx]]+\[Lambda] c1[[Mod[idx+1,m,1]]];
             y2=Subtract[1.,\[Lambda]]c2[[idx]]+\[Lambda] c2[[Mod[idx+1,m,1]]];
             \[Lambda]=1./Sqrt[y1^2+y2^2];
             y1*=\[Lambda];
             y2*=\[Lambda];

             getcenter[
                 y1,y2,\[Omega],winit,
                 maxiter,
                 TOL,regularization,ArmijoC,
                 NewtonQ,KantorovichQ,HyperbolicUpdateQ,False
             ]
         ],
         CompilationTarget->"C",
         RuntimeAttributes->{Listable},
         Parallelization->True,
         RuntimeOptions->"Speed",
         CompilationOptions->{"InlineCompiledFunctions"->True}
     ]
 ];


 (* ::Subsubsection::Closed:: *)
 (*cDouadyEarleBoundaryCurve2D*)


 cDouadyEarleBoundaryCurve2D := cDouadyEarleBoundaryCurve2D = Compile[{{c1,_Real,1},{c2,_Real,1},{\[Theta],_Real}},
     Block[{\[Lambda],idx,T,x1,x2,m,hinv},
         m=Min[Length[c1],Length[c2]];
         hinv=m/(2.Pi);
         T=Mod[\[Theta],2.Pi]hinv;
         idx=Mod[Floor[T]+1,m,1];
         \[Lambda]=Mod[T,1.];
         x1=(1.-\[Lambda]) Compile`GetElement[c1,idx]+\[Lambda] Compile`GetElement[c1,Mod[idx+1,m,1]];
         x2=(1.-\[Lambda]) Compile`GetElement[c2,idx]+\[Lambda] Compile`GetElement[c2,Mod[idx+1,m,1]];
         \[Lambda]=1./Sqrt[x1^2+x2^2];
         {\[Lambda] x1,\[Lambda] x2}
     ],
     CompilationTarget->"C",
     RuntimeAttributes->{Listable},
     Parallelization->True,
     RuntimeOptions->"Speed"
 ];


 (* ::Subsubsection::Closed:: *)
 (*cDouadyEarleExtension3D*)



 cDouadyEarleExtension3D := cDouadyEarleExtension3D = With[{
     bigone=1.+16$MachineEpsilon,
     twoPi=2.Pi,
     getcenter=cConformalBarycenter3D
     },
     Compile[{
         {c1,_Real,1},{c2,_Real,1},{c3,_Real,1},
         {w,_Real,1},
         {z01,_Real,1},{z02,_Real,1},{\[Omega],_Real,1},
         {winit,_Real,1},
         {maxiter,_Integer},
         {TOL,_Real},{regularization,_Real},{ArmijoC,_Real},
         {NewtonQ,True|False},{KantorovichQ,True|False},{HyperbolicUpdateQ,True|False}
         },
         Block[{m,w1,w2,ww,wz2,aa,bb,cc,\[Theta],y1,y2,y3,\[Lambda],hinv,T,idx},
             m=Min[Length[c1],Length[c2],Length[c3]];
             (*Transform the boundary points z with the inverse of the shift that moves w to the origin.*)
             w1=-Compile`GetElement[w,1];
             w2=-Compile`GetElement[w,2];
             wz2=(2.w1)z01+(2.w2)z02;
             ww=w1^2+w2^2;
             cc=1./Subtract[bigone+ww,wz2];
             aa=bigone-ww;
             bb=wz2-2.;
             (*Transform to angular coordinates.*)
             \[Theta]=Mod[ArcTan[Times[aa z01+bb w1,cc],Times[aa z02+bb w2,cc]],twoPi];

             (*employing piecewise-linear interpolation*)
             hinv=m/twoPi;
             T=Mod[\[Theta],twoPi]hinv;
             idx=Mod[Floor[T]+1,m,1];

             \[Lambda]=Mod[T,1.];
             y1=Subtract[1.,\[Lambda]]c1[[idx]]+\[Lambda] c1[[Mod[idx+1,m,1]]];
             y2=Subtract[1.,\[Lambda]]c2[[idx]]+\[Lambda] c2[[Mod[idx+1,m,1]]];
             y3=Subtract[1.,\[Lambda]]c3[[idx]]+\[Lambda] c3[[Mod[idx+1,m,1]]];
             \[Lambda]=1./Sqrt[y1^2+y2^2+y3^2];
             y1*=\[Lambda];
             y2*=\[Lambda];
             y3*=\[Lambda];

             getcenter[
                 y1,y2,y3,\[Omega],winit,
                 maxiter,
                 TOL,regularization,ArmijoC,
                 NewtonQ,KantorovichQ,HyperbolicUpdateQ,False
             ]
         ],
         CompilationTarget->"C",
         RuntimeAttributes->{Listable},
         Parallelization->True,
         RuntimeOptions->"Speed",
         CompilationOptions->{"InlineCompiledFunctions"->True}
     ]
 ];


 (* ::Subsubsection::Closed:: *)
 (*cDouadyEarleBoundaryCurve3D*)


 cDouadyEarleBoundaryCurve3D := cDouadyEarleBoundaryCurve3D = Compile[{{c1,_Real,1},{c2,_Real,1},{c3,_Real,1},{\[Theta],_Real}},
     Block[{\[Lambda],idx,T,x1,x2,x3,m,hinv},
         m=Min[Length[c1],Length[c2],Length[c3]];
         hinv=m/(2.Pi);
         T=Mod[\[Theta],2.Pi]hinv;
         idx=Mod[Floor[T]+1,m,1];
         \[Lambda]=Mod[T,1.];
         x1=(1.-\[Lambda])Compile`GetElement[c1,idx]+\[Lambda] Compile`GetElement[c1,Mod[idx+1,m,1]];
         x2=(1.-\[Lambda])Compile`GetElement[c2,idx]+\[Lambda] Compile`GetElement[c2,Mod[idx+1,m,1]];
         x3=(1.-\[Lambda])Compile`GetElement[c3,idx]+\[Lambda] Compile`GetElement[c3,Mod[idx+1,m,1]];
         \[Lambda]=1./Sqrt[x1^2+x2^2+x3^2];
         {\[Lambda] x1,\[Lambda] x2,\[Lambda] x3}
     ],
     CompilationTarget->"C",
     RuntimeAttributes->{Listable},
     Parallelization->True,
     RuntimeOptions->"Speed"
 ];


 (* ::Subsubsection::Closed:: *)
 (*DouadyEarleGrid*)


 Options[DouadyEarleGrid]={
     "Mesh"->{24,36},
     "Subdivisions"->{12,6}
 };

 DouadyEarleGrid[OptionsPattern[]]:=Module[{ncircs,nrays,kcircs,krays,pts,edges},
     {ncircs,nrays}=OptionValue["Mesh"];
     {kcircs,krays}=OptionValue["Subdivisions"];
     pts=Join[
         Flatten[Subdivide[0.,1.,ncircs][[2;;]] ConstantArray[CirclePoints[{1.,0.},kcircs nrays],ncircs],1],
         Flatten[Compile[{{x,_Real,1},{y,_Real,1}},
             Partition[x,1] . {y},
             RuntimeAttributes->{Listable}
         ][Subdivide[0.,1.,krays ncircs],CirclePoints[{1.,0.},nrays]],1]
         ];
     edges=Join[
         Join@@(Partition[#,2,1,1]&/@Partition[Range[ncircs kcircs nrays],{kcircs nrays}]),
         Join@@(Partition[#,2,1]&/@Partition[Range[ncircs kcircs nrays+1,ncircs kcircs nrays+nrays (krays ncircs+1)],{krays ncircs+1}])
     ];
     MeshRegion[pts,Line[edges],
         MeshCellShapeFunction->{1->(Line[#]&),0->({}&)},
         MeshCellStyle->{{1,All}->Black}
     ]
 ];


 (* ::Subsubsection::Closed:: *)
 (*DouadyEarleDisk*)


 Options[DouadyEarleDisk]={"Mesh"->720};

 DouadyEarleDisk[OptionsPattern[]]:=MeshRegion[
     DiscretizeRegion[
         Disk[],
         MaxCellMeasure->(1->4 Pi/OptionValue["Mesh"])
     ],
     MeshCellShapeFunction->{0->None,1->None},
     BaseStyle->{EdgeForm[None]},
     MeshCellStyle->{1->None}
 ];


 (* ::Subsection::Closed:: *)
 (*SampleSphericalVonMisesDistribution*)


 cSphericalVonMises:=cSphericalVonMises=Compile[{{A,_Real,2},{\[Kappa],_Real},{s,_Real},{\[Phi],_Real}},
     Block[{z},
         z=Log[1.+(Exp[2. \[Kappa]]-1.) s]/\[Kappa]-1.;
         A . {Sqrt[Abs[1.-z^2]]Cos[\[Phi]],Sqrt[1.-z^2]Sin[\[Phi]],z}
     ],
     CompilationTarget->"C",
     RuntimeAttributes->{Listable},
     Parallelization->True,
     RuntimeOptions->"Speed"
 ];

 VonMisesFisherDistribution/: Random`DistributionVector[VonMisesFisherDistribution[\[Mu]_,\[Kappa]_], dims_,prec_]:=Module[{A,\[Phi],s},
     A=RotationMatrix[{{0,0,1},\[Mu]}];
     \[Phi]=RandomReal[{-Pi,Pi},dims];
     s=RandomReal[{0,1},dims];
     cSphericalVonMises[A,\[Kappa],s,\[Phi]]
 ];


 (* ::Subsection::Closed:: *)
 (*cSphericalGaussProcess*)


 (* ::Subsubsection::Closed:: *)
 (*cSphericalGaussProces2D*)


 cSphericalGaussProcess2D:=cSphericalGaussProcess2D=Compile[{{\[Nu]0,_Real,1},{X,_Real,1}},
     Block[{\[Nu],u,v,w,i,j,\[Alpha]},
         \[Nu]=\[Nu]0;
         Table[
             If[k==0,
                 \[Nu]
             ,
                 w=Compile`GetElement[X,k]{-Compile`GetElement[\[Nu],2], Compile`GetElement[\[Nu],1]};
                 \[Alpha]=Sqrt[w . w];
                 \[Nu]=\[Nu] Cos[\[Alpha]]+Sinc[\[Alpha]]w
             ],
             {k,0,Length[X]}
         ]
     ],
     CompilationTarget->"C",
     RuntimeAttributes->{Listable},
     Parallelization->True
 ];


 (* ::Subsubsection::Closed:: *)
 (*cSphericalGaussProcess3D*)


 cSphericalGaussProcess3D := cSphericalGaussProcess3D = With[{id=IdentityMatrix[3,WorkingPrecision->MachinePrecision]},
     Compile[{{\[Nu]0,_Real,1},{X,_Real,2}},
         Block[{\[Nu],u,v,w,i,j,\[Alpha]},
             \[Nu]=\[Nu]0;
             Table[
                 If[k==0,
                     \[Nu]
                 ,
                     {i,j}=Ordering[Abs[\[Nu]]][[1;;2]];
                     u=id[[i]];
                     u=u-\[Nu] \[Nu] . u;
                     u/=Sqrt[u . u];
                     v=id[[j]];
                     v=v-\[Nu] \[Nu] . v;
                     v/=Sqrt[v . v];
                     w=u Compile`GetElement[X,k,1]+v Compile`GetElement[X,k,2];
                     \[Alpha]=Sqrt[w . w];
                     \[Nu]=\[Nu] Cos[\[Alpha]]+Sinc[\[Alpha]]w
                 ]
             ,
                 {k,0,Length[X]}
             ]
         ],
         CompilationTarget->"C",
         RuntimeAttributes->{Listable},
         Parallelization->True
     ]
 ];


 (* ::Chapter::Closed:: *)
 (*End*)


 End[];


 EndPackage[];
