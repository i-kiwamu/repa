{-# OPTIONS -fno-warn-incomplete-patterns #-}
{-# LANGUAGE PackageImports #-}
{-# LANGUAGE FlexibleContexts #-}

-- | Algorithms operating on matrices.
-- 
--   These functions should give performance comparable with nested loop C
--   implementations. 
-- 
--   If you care deeply about runtime performance then you
--   may be better off using a binding to LAPACK, such as hvector.
--
module Data.Array.Repa.Algorithms.Matrix
        ( --  * Projections
          row
        , col

          -- * Matrix Multiplication.
        , mmultP,      mmultS,      mmult

          -- * inner product between Vector
        , dotP,        dotS

          -- * Matrix-Vector Multiplication.
        , mvmultP,     mvmultS,     mvmult

          -- * outer product between Vector
        , outerP,      outerS,      outer

          -- * Transposition.
        , transpose2P, transpose2S

          -- * Trace.
        , trace2P,     trace2S

          -- * Determinant.
        , det2P,       det2S

          -- * LU decomposition.
        , luP,         luS,         lu

          -- * diagonal elements of matrix.
        , mdiagP,      mdiagS,      mdiag)


where
import Data.Array.Repa                  as R
import Data.Array.Repa.Eval             as R
import Data.Array.Repa.Unsafe           as R
import Control.Monad
import Control.Monad.ST.Strict
import Data.Vector.Unboxed.Base (Unbox)
import Data.List                        as L


-- Projections ----------------------------------------------------------------
-- | Take the row number of a rank-2 index.
row :: DIM2 -> Int
row (Z :. r :. _) = r
{-# INLINE row #-}


-- | Take the column number of a rank-2 index.
col :: DIM2 -> Int
col (Z :. _ :. c) = c
{-# INLINE col #-}


-- MMult ----------------------------------------------------------------------
-- | Matrix matrix multiply, in parallel.
mmultP  :: (Num a, Elt a, Unbox a, Source r a, Monad m)
        => Array r DIM2 a
        -> Array r DIM2 a
        -> m (Array U DIM2 a)

mmultP arr brr 
 = [arr, brr] `deepSeqArrays` 
   do   trr      <- transpose2P brr
        let (Z :. h1  :. _)  = extent arr
        let (Z :. _   :. w2) = extent brr
        computeP 
         $ fromFunction (Z :. h1 :. w2)
         $ \ix   -> R.sumAllS 
                  $ R.zipWith (*)
                        (unsafeSlice arr (Any :. (row ix) :. All))
                        (unsafeSlice trr (Any :. (col ix) :. All))
{-# NOINLINE mmultP #-}


-- | Matrix matrix multiply, sequentially.
mmultS  :: (Num a, Elt a, Unbox a, Source r a)
        => Array r DIM2 a 
        -> Array r DIM2 a
        -> Array U DIM2 a

mmultS arr brr
 = [arr, brr]  `deepSeqArrays` (runST $
   do   trr     <- R.now $ transpose2S brr
        let (Z :. h1  :. _)  = extent arr
        let (Z :. _   :. w2) = extent brr
        return $ computeS 
         $ fromFunction (Z :. h1 :. w2)
         $ \ix   -> R.sumAllS 
                  $ R.zipWith (*)
                        (unsafeSlice arr (Any :. (row ix) :. All))
                        (unsafeSlice trr (Any :. (col ix) :. All)))
{-# NOINLINE mmultS #-}


-- | Matrix matrix multiply, delayed.
mmult :: (Num a, Elt a, Unbox a, Source r1 a, Source r2 a)
      => Array r1 DIM2 a
      -> Array r2 DIM2 a
      -> Array D  DIM2 a
mmult arr brr = fromFunction (Z :. h1 :. w2)
                $ \ix -> R.sumAllS
                         $ R.zipWith (*)
                           (unsafeSlice arr (Any :. (row ix) :. All))
                           (unsafeSlice trr (Any :. (col ix) :. All))
  where trr = R.transpose brr
        (Z :. h1 :. _ ) = extent arr
        (Z :. _  :. w2) = extent brr
{-# NOINLINE mmult #-}


-- inner product --------------------------------------------------------------
-- | inner product in parallel
dotP :: (Num a, Elt a, Unbox a, Monad m, Source r1 a, Source r2 a)
     => Array r1 DIM1 a
     -> Array r2 DIM1 a
     -> m a
dotP av bv = R.sumAllP $ R.zipWith (*) av bv
{-# NOINLINE dotP #-}


-- | inner product sequentially
dotS :: (Num a1, Elt a1, Unbox a1, Source r1 a1,
         Num a2, Elt a2, Unbox a2, Source r2 a2,
         Num a3, Elt a3, Unbox a3)
     => Array r1 DIM1 a1
     -> Array r2 DIM1 a2
     -> a3
dotS av bv = R.sumAllS $ R.zipWith (*) av bv
{-# NOINLINE dotS #-}


-- Matrix vector multiply -----------------------------------------------------
-- | matrix vector multiply in parallel
mvmultP :: (Num a1, Num a2, Unbox a2, Source r1 a1, Source r2 a2, Monad m)
        => Array r1 DIM2 a1
        -> Array r2 DIM1 a2
        -> m (Array U DIM1 a2)
mvmultP arr vrr
 = arr `deepSeqArray` vrr `deepSeqArray` computeP $ mvmult arr vrr
{-# NOINLINE mvmultP #-}


-- | matrix vector multiply sequentially
mvmultS :: (Num a1, Num a2, Unbox a2, Source r1 a1, Source r2 a2)
        => Array r1 DIM2 a1
        -> Array r2 DIM1 a2
        -> Array U  DIM1 a2
mvmultS arr vrr
 = arr `deepSeqArray` vrr `deepSeqArray` (runST $
   return $ computeS $ mvmult arr vrr)
{-# NOINLINE mvmultS #-}


-- | matrix vector multiply delayed
mvmult :: (Num a1, Elt a1, Unbox a1, Source r1 a1, Source r2 a2)
       => Array r1 DIM2 a1
       -> Array r2 DIM1 a2
       -> Array D  DIM1 a3
mvmult arr vrr
 = R.map doti rs
 where nr = row $ extent arr
       rs = fromListUnboxed (Z :. nr) [0..(nr-1)]
       doti :: Int -> a2
       doti r = dotS (slice arr (Any :. r :. All)) vrr
       {-# NOINLINE doti #-}
{-# NOINLINE mvmult #-}


-- outer product --------------------------------------------------------------
-- | outer product in parallel
outerP :: (Num a, Unbox a, Source r1 a, Source r2 a, Monad m)
       => Array r1 DIM1 a
       -> Array r2 DIM1 a
       -> m (Array U DIM2 a)
outerP va vb = va `deepSeqArray` vb `deepSeqArray` computeP $ outer va vb
{-# NOINLINE outerP #-}


-- | outer product sequentially
outerS :: (Num a, Unbox a, Source r1 a, Source r2 a)
       => Array r1 DIM1 a
       -> Array r2 DIM1 a
       -> Array U DIM2 a
outerS va vb = va `deepSeqArray` vb `deepSeqArray` computeS $ outer va vb
{-# NOINLINE outerS #-}


-- | outer product delayed
outer :: (Num a, Source r1 a, Source r2 a)
       => Array r1 DIM1 a
       -> Array r2 DIM1 a
       -> Array D DIM2 a
outer va vb = ma *^ mb
 where na = size $ extent va
       nb = size $ extent vb
       ma = extend (Z :. All :. nb) va
       mb = extend (Z :. na :. All) vb
{-# NOINLINE outer #-}


-- Transpose ------------------------------------------------------------------
-- | Transpose a 2D matrix, in parallel.
transpose2P :: (Unbox a, Source r a, Monad m)
    => Array r DIM2 a 
    -> m (Array U DIM2 a)

transpose2P arr
 = arr `deepSeqArray` computeUnboxedP $
   unsafeBackpermute new_extent swap arr
 where  swap (Z :. i :. j)      = Z :. j :. i
        new_extent              = swap (extent arr)
{-# NOINLINE transpose2P #-}


-- | Transpose a 2D matrix, sequentially.
transpose2S :: (Unbox a, Source r a)
            => Array r DIM2 a 
            -> Array U DIM2 a

transpose2S arr
 = arr `deepSeqArray` computeUnboxedS $
   unsafeBackpermute new_extent swap arr
 where  swap (Z :. i :. j)      = Z :. j :. i
        new_extent              = swap (extent arr)
{-# NOINLINE transpose2S #-}



-- Trace ----------------------------------------------------------------------
-- | Get the trace of a (square) 2D matrix, in parallel.
trace2P :: (Num a, Elt a, Unbox a, Source r a, Monad m)
        => Array r DIM2 a
        -> m a
trace2P x 
 = liftM (safeHead . toList) $ sumP $ slice y (Z :. (0 :: Int) :. All)
 where
    safeHead []     = error "repa-algorithms: trace2P empty list"
    safeHead (x':_) = x'
    y               = unsafeBackpermute (extent x) f x
    f (Z :. i :. j) = Z :. (i - j) `mod` nRows:. j
    Z :. nRows :. _nCols = extent x
{-# NOINLINE trace2P #-}


-- | Get the trace of a (square) 2D matrix, sequentially.
trace2S :: (Num a, Elt a, Unbox a, Source r a)
        => Array r DIM2 a
        -> a
trace2S x 
 = safeHead $ toList $ sumS $ slice y (Z :. (0 :: Int) :. All)
 where
    safeHead []     = error "repa-algorithms: trace2S empty list"
    safeHead (x':_) = x'
    y               =  unsafeBackpermute (extent x) f x
    f (Z :. i :. j) = Z :. (i - j) `mod` nRows:. j
    Z :. nRows :. _nCols = extent x
{-# NOINLINE trace2S #-}


-- Determinant ----------------------------------------------------------------
-- | Get the determinant of a (square) 2D matrix, in parallel.
det2P :: (Num a, Elt a, Fractional a, Unbox a, Source r a, Monad m)
      => Array r DIM2 a
      -> m a
det2P arr = arr `deepSeqArray` do
      luarr <- luP arr
      ds <- mdiagP luarr
      x <- R.foldP (*) 1 ds
      return $ x ! Z
{-# NOINLINE det2P #-}


-- | Get the determinant of a (square) 2D matrix sequentially
det2S :: (Num a, Elt a, Fractional a, Unbox a, Source r a)
      => Array r DIM2 a
      -> a
det2S arr = arr `deepSeqArray` (runST $
      do luarr <- R.now $ luS arr
         ds <- R.now $ mdiagS luarr
         x <- R.now $ R.foldS (*) 1 ds
         return $ x ! Z)
{-# NOINLINE det2S #-}


-- LU decomposition -----------------------------------------------------------
-- | LU decomposition by Crout algorithms (L(i,i) = 1) in parallel
luP :: (Num a, Fractional a, Unbox a, Source r a, Monad m)
    => Array r DIM2 a
    -> m (Array U DIM2 a)
luP arr = arr `deepSeqArray` computeP $ lu arr
{-# NOINLINE luP #-}


-- | LU decomposition by Crout algorithms (L(i,i) = 1) sequentially
luS :: (Num a, Fractional a, Unbox a, Source r a)
    => Array r DIM2 a
    -> Array U DIM2 a
luS arr = arr `deepSeqArray` computeS $ lu arr
{-# NOINLINE luS #-}


-- | LU decomposition by Crout algorithms (L(i,i) = 1) delayed
lu :: (Num a, Fractional a, Source r a)
   => Array r DIM2 a
   -> Array D DIM2 a
lu arr = arr'
  where arr' = R.fromFunction e (\(Z :. i :. j) -> f i j)
        e = R.extent arr
        f i j = if i > j
                -- lower triangles
                then (aij - sum [a'i_ k * a'_j k | k <- [0..(j-1)]]) / a'jj
                -- upper triangles
                else  aij - sum [a'i_ k * a'_j k | k <- [0..(i-1)]]
          where
            aij    = arr  R.! (Z :. i :. j)
            a'i_ k = arr' R.! (Z :. i :. k)
            a'_j k = arr' R.! (Z :. k :. j)
            a'jj   = arr' R.! (Z :. j :. j)
{-# NOINLINE lu #-}


-- Swap row -----------------------------------------------------------
-- | Swap row delayed
swapRow2D :: (Source r Double)
          => Int
          -> Int
          -> Array r DIM2 Double
          -> Array D DIM2 Double
swapRow2D i j arr = undefined


-- Diagonal -------------------------------------------------------------------
-- | Diagonal elements in parallel
mdiagP :: (Unbox a, Source r a, Monad m)
       => Array r DIM2 a
       -> m (Array U DIM1 a)
mdiagP arr = arr `deepSeqArray` computeP $ mdiag arr
{-# NOINLINE mdiagP #-}


-- | Diagonal elements sequentially
mdiagS :: (Unbox a, Source r a)
       => Array r DIM2 a
       -> Array U DIM1 a
mdiagS arr = arr `deepSeqArray` computeS $ mdiag arr
{-# NOINLINE mdiagS #-}


-- | Diagonal elements delayed
mdiag :: (Source r a)
      => Array r DIM2 a
      -> Array D DIM1 a
mdiag arr = slice md (Any :. (0::Int))
  where md = R.traverse arr f g
        f :: DIM2 -> DIM2
        f (Z :. r :. c) = Z :. (min r c) :. (1::Int)
        g :: (DIM2 -> a) -> DIM2 -> a
        g get (Z :. r :. _) = get (Z :. r :. r)
{-# NOINLINE mdiag #-}


-- Solve linear equations by LU decomposition ------------------------------------
-- | Solve in parallel
solveLUP :: (Num a1, Num a2, Source r1 a1, Source r2 a2, Monad m)
         => Array r1 DIM2 a1
         -> Array r2 DIM1 a2
         -> m (Array U DIM1 a2)
solveLUP arr vec = arr `deepSeqArray` vec `deepSeqArray`
         do luarr <- luP arr
            let larr = R.fromFunction (R.extent arr) (extractL luarr)
                uarr = R.fromFunction (R.extent arr) (extractU luarr)
                y = R.fromFunction e (\(Z :. i) -> fy larr i)
            computeP
             $ R.fromFunction e (\(Z :. i) -> fx uarr y (n-i))
  where e = R.extent vec
        n = size e - 1
        extractL :: (Num a, Fractional a, Unbox a)
                 => Array U DIM2 a -> DIM2 -> a
        extractL mat (Z :. r :. c)
            | r == c    = 1.0
            | r > c     = mat R.! (Z :. r :. c)
            | otherwise = 0.0
        extractU :: (Num a, Fractional a, Unbox a)
                 => Array U DIM2 a -> DIM2 -> a
        extractU mat (Z :. r :. c)
            | r <= c    = mat R.! (Z :. r :. c)
            | otherwise = 0.0
        fy :: (Num a) => Array D DIM2 a -> Int -> a
        fy _ 0    = vec ! (Z :. 0 :: DIM1)
        fy larr i = vec ! (Z :. i) - dotS ys ls
          where ys = R.fromFunction (Z :. i) (\(Z :. k) -> fy larr k)
                ls = slice (extract (Z :. i :. 0 :: DIM2)
                                    (Z :. 1 :. i :: DIM2) larr)
                           (Any :. (0::Int) :. All)
        fx :: (Num a, Fractional a, Elt a, Unbox a)
           => Array D DIM2 a -> Array D DIM1 a -> Int -> a
        fx uarr y 0 = y ! (Z :. n) / uarr ! (Z :. n :. n)
        fx uarr y i = (y ! (Z :. n-i) - dotS xs us) / uarr ! (Z :. n-i :. n-i)
          where xs = R.fromListUnboxed (Z :. i)
                     $ L.map (fx uarr y) $ reverse [0..(i-1)]
                us = slice (extract (Z :. n-i :. n-i+1) (Z :. 1 :. i) uarr)
                           (Any :. (0::Int) :. All)
