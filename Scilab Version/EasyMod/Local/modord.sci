function [p] = modord(n,ni)

// ------------------   This file is part of EasyMod   ----------------------------
//  Internal function
//
//  This function calculates the modal order p as a function of:
//    - the number of modes (n),
//    - the number of outputs (ni).
//
// Copyright (C) 2012 David WATTIAUX, Georges KOUROUSSIS, Delphine LUPANT


// Building of data matrix
N = [1:40] ; 
N = N' ; 
Ni1 = 2*N ; 
Ni2 = N ; 
Ni3 = [1 2 2 3 4 4 5 6 6 7 8 8 9 10 10 11 12 12 13 14 14 15 16 16 17 18 18 19 20 20 21 22 22 23 24 24 25 26 26 27] ; 
Ni3 = Ni3' ; 
Ni4 = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20] ; 
Ni4 = Ni4' ; 
Ni5 = [1 1 2 2 2 3 3 4 4 4 5 5 6 6 6 7 7 8 8 8 9 9 10 10 10 11 11 12 12 12 13 13 14 14 14 15 15 16 16 16] ; 
Ni5 = Ni5' ; 
Ni6 = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12 13 13 13 14] ; 
Ni6 = Ni6' ; 
Ni7 = [1 1 1 2 2 2 2 3 3 3 4 4 4 4 5 5 5 6 6 6 6 7 7 7 8 8 8 8 9 9 9 10 10 10 10 11 11 11 12 12] ; 
Ni7 = Ni7' ; 
Ni8 = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8 9 9 9 9 10 10 10 10] ; 
Ni8 = Ni8' ; 
Ni9 = [1 1 1 1 2 2 2 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 6 6 6 6 6 7 7 7 7 8 8 8 8 8 9 9 9 9] ; 
Ni9 = Ni9' ; 
Ni10 = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8] ; 
Ni10 = Ni10' ; 
Ni11 = [1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 6 6 6 6 6 6 7 7 7 7 7 8 8] ; 
Ni11 = Ni11' ; 
P = [N Ni1 Ni2 Ni3 Ni4 Ni5 Ni6 Ni7 Ni8 Ni9 Ni10 Ni11] ; 

// Value searching
p = P(n,ni+1) ; 

endfunction