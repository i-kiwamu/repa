{-# OPTIONS -fno-warn-incomplete-patterns #-}
{-# LANGUAGE PackageImports #-}

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
        , mmultP,      mmultS

          -- * inner product between Vector
        , dotP,        dotS

          -- * Matrix-Vector Multiplication.
        , mvmultP,     mvmultS

          -- * outer product between Vector
        , outerP,        outerS

          -- * Transposition.
        , t2P, t2S

          -- * Trace.
        , trace2P, trace2S

          -- * Transposition.
        , det2P)


where
import Data.Array.Repa                  as R
import Data.Array.Repa.Eval             as R
import Data.Array.Repa.Unsafe           as R
import Control.Monad
import Control.Monad.ST.Strict
import Data.Vector.Unboxed.Base (Unbox)
import Data.Ix (range)


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
mmultP  :: Monad m
        => Array U DIM2 Double 
        -> Array U DIM2 Double 
        -> m (Array U DIM2 Double)

mmultP arr brr 
 = [arr, brr] `deepSeqArrays` 
   do   trr      <- t2P brr
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
mmultS  :: Array U DIM2 Double 
        -> Array U DIM2 Double 
        -> Array U DIM2 Double

mmultS arr brr
 = [arr, brr]  `deepSeqArrays` (runST $
   do   trr     <- R.now $ t2S brr
        let (Z :. h1  :. _)  = extent arr
        let (Z :. _   :. w2) = extent brr
        return $ computeS 
         $ fromFunction (Z :. h1 :. w2)
         $ \ix   -> R.sumAllS 
                  $ R.zipWith (*)
                        (unsafeSlice arr (Any :. (row ix) :. All))
                        (unsafeSlice trr (Any :. (col ix) :. All)))
{-# NOINLINE mmultS #-}


dotP :: (Monad m, Source r1 a, Source r2 a, Num a, Elt a, Unbox a)
         => Array r1 DIM1 a
         -> Array r2 DIM1 a
         -> m a
dotP av bv = R.sumAllP $ R.zipWith (*) av bv
{-# NOINLINE dotP #-}


dotS :: (Source r1 a, Source r2 a, Num a, Elt a, Unbox a)
         => Array r1 DIM1 a
         -> Array r2 DIM1 a
         -> a
dotS av bv = R.sumAllS $ R.zipWith (*) av bv
{-# NOINLINE dotS #-}


mvmultP :: Monad m
           => Array U DIM2 Double
           -> Array U DIM1 Double
           -> m (Array U DIM1 Double)
mvmultP arr vrr
 = arr `deepSeqArray` vrr `deepSeqArray` (computeP $ R.map doti rs)
 where nr = row $ extent arr
       rs = fromListUnboxed (Z :. nr) [0..(nr-1)]
       doti :: Int -> Double
       doti r = dotS (slice arr (Any :. r :. All)) vrr
       {-# NOINLINE doti #-}
{-# NOINLINE mvmultP #-}


mvmultS :: Array U DIM2 Double
           -> Array U DIM1 Double
           -> Array U DIM1 Double
mvmultS arr vrr
 = arr `deepSeqArray` vrr `deepSeqArray` (runST $
   return $ computeS $ R.map doti rs)
 where nr = row $ extent arr
       rs = fromListUnboxed (Z :. nr) [0..(nr-1)] 
       doti :: Int -> Double
       doti r = dotS (slice arr (Any :. r :. All)) vrr
       {-# NOINLINE doti #-}
{-# NOINLINE mvmultS #-}


outerP :: Monad m
          => Array U DIM1 Double
          -> Array U DIM1 Double
          -> m (Array U DIM2 Double)
outerP va vb
 = [va, vb] `deepSeqArrays` (computeP $ ma *^ mb)
 where na = size $ extent va
       nb = size $ extent vb
       ma = extend (Z :. All :. nb) va
       mb = extend (Z :. na :. All) vb
{-# NOINLINE outerP #-}


outerS :: Array U DIM1 Double
          -> Array U DIM1 Double
          -> Array U DIM2 Double
outerS va vb
 = [va, vb] `deepSeqArrays` (computeS $ ma *^ mb)
 where na = size $ extent va
       nb = size $ extent vb
       ma = extend (Z :. All :. nb) va
       mb = extend (Z :. na :. All) vb
{-# NOINLINE outerS #-}


-- Transpose ------------------------------------------------------------------
-- | Transpose a 2D matrix, in parallel.
t2P :: Monad m 
       => Array U DIM2 Double 
       -> m (Array U DIM2 Double)

t2P arr
 = arr `deepSeqArray` computeUnboxedP $
   unsafeBackpermute new_extent swap arr
 where  swap (Z :. i :. j)      = Z :. j :. i
        new_extent              = swap (extent arr)
{-# NOINLINE t2P #-}


-- | Transpose a 2D matrix, sequentially.
t2S :: Array U DIM2 Double 
       -> Array U DIM2 Double

t2S arr
 = arr `deepSeqArray` computeUnboxedS $
   unsafeBackpermute new_extent swap arr
 where  swap (Z :. i :. j)      = Z :. j :. i
        new_extent              = swap (extent arr)
{-# NOINLINE t2S #-}


-- Trace ------------------------------------------------------------------------
-- | Get the trace of a (square) 2D matrix, in parallel.
trace2P :: Monad m => Array U DIM2 Double -> m Double
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
trace2S :: Array U DIM2 Double -> Double
trace2S x 
 = safeHead $ toList $ sumS $ slice y (Z :. (0 :: Int) :. All)
 where
    safeHead []     = error "repa-algorithms: trace2S empty list"
    safeHead (x':_) = x'
    y               =  unsafeBackpermute (extent x) f x
    f (Z :. i :. j) = Z :. (i - j) `mod` nRows:. j
    Z :. nRows :. _nCols = extent x
{-# NOINLINE trace2S #-}


-- Determinant -------------------------------------------------------------------
-- | Get the determinant of a (square) 2D matrix, in parallel.
det2P :: Monad m => Array U DIM2 Double -> m Double
det2P arr = undefined


-- LU decomposition --------------------------------------------------------------
-- | LU decomposition (L(i,i) = 1) in parallel
luP :: Monad m
       => Array U DIM2 Double 
       -> m (Array U DIM2 Double)
luP arr = arr `deepSeqArray` computeP arr'
  where arr' = R.fromFunction e (\(Z :. i :. j) -> f i j)
        e = R.extent arr
        n = row e - 1
        f i j = if i > j
                then (aij - sum [a'i_ k * a'_j k | k <- [1..(j-1)]]) / a'jj
                else  aij - sum [a'i_ k * a'_j k | k <- [1..(i-1)]]
          where
            aij    = arr  R.! (Z :. i :. j)
            a'i_ k = arr' R.! (Z :. i :. k)
            a'_j k = arr' R.! (Z :. k :. j)
            a'jj   = arr' R.! (Z :. j :. j)
{-# NOINLINE luP #-}


-- | LU decomposition (L(i,i) = 1) sequentially
luS :: Array U DIM2 Double 
       -> Array U DIM2 Double
luS arr = arr `deepSeqArray` computeS arr'
  where arr' = R.fromFunction e (\(Z :. i :. j) -> f i j)
        e = R.extent arr
        n = row e - 1
        f i j = if i > j
                then (aij - sum [a'i_ k * a'_j k | k <- [1..(j-1)]]) / a'jj
                else  aij - sum [a'i_ k * a'_j k | k <- [1..(i-1)]]
          where
            aij    = arr  R.! (Z :. i :. j)
            a'i_ k = arr' R.! (Z :. i :. k)
            a'_j k = arr' R.! (Z :. k :. j)
            a'jj   = arr' R.! (Z :. j :. j)
{-# NOINLINE luS #-}
