{-# OPTIONS -fno-warn-incomplete-patterns #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}

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
        , dotP,        dotS,        dot

          -- * Matrix-Vector Multiplication.
        , mvmultP,     mvmultS,     mvmult

          -- * outer product between Vector
        , outerP,      outerS,      outer

          -- * Transposition.
        , transpose2P, transpose2S, transpose2

          -- * Append vertically
        , vappend, (-+-)

          -- * Trace.
        , trace2P,     trace2S,     trace2

          -- * Determinant.
        , det2P,       det2S

          -- * LU factorization.
        , luP,         luS,         lu

          -- * diagonal elements of matrix.
        , mdiagP,      mdiagS,      mdiag

          -- * diagonal matrix from vector
        , vdiagP,      vdiagS,      vdiag

          -- * solve linear system by LU factorization
        , solveLUP,    solveLUS,    solveLU)


where
import Data.Array.Repa                  as R
import Data.Array.Repa.Eval             as R
import Data.Array.Repa.Unsafe           as R
import Control.Monad
import Control.Monad.ST.Strict
import qualified Data.Vector.Unboxed    as V
import Data.List                        as L
import Data.Ord (comparing)
-- import Debug.Trace


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
mmultP arr brr = [arr, brr] `R.deepSeqArrays` do
       trr <- R.computeP $ R.transpose brr
       R.computeP $ R.fromFunction (Z :. h1 :. w2) (fun arr trr)
  where (Z :. h1 :. _ ) = R.extent arr
        (Z :. _  :. w2) = R.extent brr
        fun :: Array U DIM2 Double
            -> Array U DIM2 Double
            -> DIM2
            -> Double
        fun ma mb (Z :. r :. c) = [ma, mb] `R.deepSeqArrays`
            let va = slicedUnbox ma r
            in V.sum $ V.zipWith (*) va $ slicedUnbox mb c
{-# NOINLINE mmultP #-}


-- | Matrix matrix multiply, sequentially.
mmultS :: Array U DIM2 Double
       -> Array U DIM2 Double
       -> Array U DIM2 Double
mmultS arr brr = [arr, brr] `R.deepSeqArrays` runST (do
       trr <- R.now $ R.computeS $ R.transpose brr
       return $ R.fromUnboxed (Z :. h1 :. w2)
              $ V.concat $ L.map (fun arr trr) [0..(h1-1)])
  where (Z :. h1 :. _ ) = R.extent arr
        (Z :. _  :. w2) = R.extent brr
        fun :: Array U DIM2 Double
            -> Array U DIM2 Double
            -> Int
            -> V.Vector Double
        fun ma mb h = [ma, mb] `R.deepSeqArrays`
                      let va = slicedUnbox ma h
                          vw = V.enumFromN 0 (V.length va)
                          f = V.sum . V.zipWith (*) va
                      in V.map (f . slicedUnbox mb) vw
{-# NOINLINE mmultS #-}


slicedUnbox :: Array U DIM2 Double
            -> Int
            -> V.Vector Double
slicedUnbox mat r = if r >= nr || r < 0
                       then V.empty
                       else V.take nc $ V.drop (r*nc) vm
  where (Z :. nr :. nc) = extent mat
        vm = toUnboxed mat
{-# INLINE slicedUnbox #-}


-- | Matrix matrix multiply, delayed.
mmult :: (Source r1 Double, Source r2 Double)
      => Array r1 DIM2 Double
      -> Array r2 DIM2 Double
      -> Array D  DIM2 Double
mmult arr brr =
    let trr = R.transpose brr
    in trr `R.deepSeqArray` R.fromFunction (Z :. h1 :. w2) (fun arr trr)
  where (Z :. h1 :. _ ) = R.extent arr
        (Z :. _  :. w2) = R.extent brr
        fun :: (Source r1 Double, Source r2 Double)
            => Array r1 DIM2 Double
            -> Array r2 DIM2 Double
            -> DIM2
            -> Double
        fun ma mb (Z :. r :. c) =
            let va = unsafeSlice ma (Z :. r :. All)
                vb = unsafeSlice mb (Z :. c :. All)
            in dot va vb
{-# NOINLINE mmult #-}


-- Inner prodduct -------------------------------------------------------------
-- | inner product, in parallel.
dotP :: (Monad m)
     => Array U DIM1 Double
     -> Array U DIM1 Double
     -> m Double
dotP av bv = av `deepSeqArray` bv `deepSeqArray`
             R.sumAllP $ R.zipWith (*) av bv
{-# NOINLINE dotP #-}


-- | inner product, sequentially.
dotS :: Array U DIM1 Double
     -> Array U DIM1 Double
     -> Double
dotS av bv = [av, bv] `deepSeqArrays`
             R.sumAllS $ R.zipWith (*) av bv
{-# NOINLINE dotS #-}


-- | inner product, delayed.
dot :: (Source r1 Double, Source r2 Double)
    => Array r1 DIM1 Double
    -> Array r2 DIM1 Double
    -> Double
dot av bv = R.sumAllS $ R.zipWith (*) av bv
{-# NOINLINE dot #-}


-- Matrix vector mutiply ------------------------------------------------------
-- | Matrix vector multiply, in parallel.
mvmultP :: (Monad m)
       => Array U DIM2 Double
       -> Array U DIM1 Double
       -> m (Array U DIM1 Double)
mvmultP arr vec = arr `deepSeqArray` vec `deepSeqArray`
        R.computeP $ mvmult arr vec
{-# NOINLINE mvmultP #-}


-- | Matrix vector multiply, sequentially.
mvmultS :: Array U DIM2 Double
        -> Array U DIM1 Double
        -> Array U DIM1 Double
mvmultS arr vec = arr `deepSeqArray` vec `deepSeqArray`
                  R.fromUnboxed (Z :. h1) $ V.map (fun arr vb) vh
  where h1 = row $ R.extent arr
        vh = V.enumFromN 0 h1
        vb = R.toUnboxed vec
        fun :: Array U DIM2 Double
            -> V.Vector Double
            -> Int
            -> Double
        fun m v h = let va = slicedUnbox m h
                    in V.sum $ V.zipWith (*) va v
{-# NOINLINE mvmultS #-}


-- | Matrix vector multiply, delayed.
mvmult :: (Source r1 Double, Source r2 Double)
       => Array r1 DIM2 Double
       -> Array r2 DIM1 Double
       -> Array D  DIM1 Double
mvmult arr vec = R.fromFunction e
                 $ \(Z :. r) ->
                       R.sumAllS 
                       $ R.zipWith (*) (unsafeSlice arr (Any :. r :. All))
                                       vec
  where (e :. _) = R.extent arr
{-# NOINLINE mvmult #-}


-- Outer product --------------------------------------------------------------
-- | outer product, in parallel.
outerP :: (Monad m)
       => Array U DIM1 Double
       -> Array U DIM1 Double
       -> m (Array U DIM2 Double)
outerP va vb = [va, vb] `deepSeqArrays` computeP $ outer va vb
{-# NOINLINE outerP #-}


-- | outer product, sequentially.
outerS :: Array U DIM1 Double
       -> Array U DIM1 Double
       -> Array U DIM2 Double
outerS va vb = [va, vb] `deepSeqArrays` computeS $ outer va vb
{-# NOINLINE outerS #-}


-- | outer product, delayed.
outer :: (Source r1 Double, Source r2 Double)
      => Array r1 DIM1 Double
      -> Array r2 DIM1 Double
      -> Array D DIM2 Double
outer va vb = ma *^ mb
  where na = size $ extent va
        nb = size $ extent vb
        ma = extend (Z :. All :. nb) va
        mb = extend (Z :. na :. All) vb
{-# NOINLINE outer #-}



-- Transpose ------------------------------------------------------------------
-- | Transpose a 2D matrix, in parallel.
transpose2P :: (Monad m)
            => Array U DIM2 Double
            -> m (Array U DIM2 Double)
transpose2P arr
 = arr `deepSeqArray` computeUnboxedP $ transpose2 arr
{-# NOINLINE transpose2P #-}


-- | Transpose a 2D matrix, sequentially.
transpose2S :: Array U DIM2 Double
            -> Array U DIM2 Double
transpose2S arr
 = arr `deepSeqArray` computeUnboxedS
   $ unsafeBackpermute new_extent swap arr
 where  swap (Z :. i :. j)      = Z :. j :. i
        new_extent              = swap (extent arr)
{-# NOINLINE transpose2S #-}


-- | Transpose a 2D matrix, delayed
transpose2 :: (Source r a)
           => Array r DIM2 a
           -> Array D DIM2 a
transpose2 arr = unsafeBackpermute new_extent swap arr
 where  swap (Z :. i :. j)      = Z :. j :. i
        new_extent              = swap (extent arr)
{-# NOINLINE transpose2 #-}


-- append two matrix vertically -----------------------------------------------
vappend, (-+-) :: (Source r1 a, Source r2 a)
        => Array r1 DIM2 a
        -> Array r2 DIM2 a
        -> Array D  DIM2 a
vappend arr1 arr2 = transpose2 $ R.append ta1 ta2
  where ta1 = transpose2 arr1
        ta2 = transpose2 arr2
{-# NOINLINE vappend #-}

(-+-) = vappend
{-# NOINLINE (-+-) #-}


-- Trace ----------------------------------------------------------------------
-- | Get the trace of a (square) 2D matrix, in parallel.
trace2P :: (Monad m)
        => Array U DIM2 Double
        -> m Double
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
trace2S :: Array U DIM2 Double
        -> Double
trace2S x 
 = safeHead $ toList $ sumS $ slice y (Z :. (0 :: Int) :. All)
 where
    safeHead []     = error "repa-algorithms: trace2S empty list"
    safeHead (x':_) = x'
    y               =  unsafeBackpermute (extent x) f x
    f (Z :. i :. j) = Z :. (i - j) `mod` nRows:. j
    Z :. nRows :. _nCols = extent x
{-# NOINLINE trace2S #-}


-- | Get the trace of a (square) 2D matrix, sequentially.
trace2 :: (Source r Double)
       => Array r DIM2 Double
       -> Double
trace2 x 
 = safeHead $ toList $ sumS $ slice y (Z :. (0 :: Int) :. All)
 where
    safeHead []     = error "repa-algorithms: trace2S empty list"
    safeHead (x':_) = x'
    y               =  unsafeBackpermute (extent x) f x
    f (Z :. i :. j) = Z :. (i - j) `mod` nRows:. j
    Z :. nRows :. _nCols = extent x
{-# NOINLINE trace2 #-}


-- Determinant ----------------------------------------------------------------
-- | Get the determinant of a (square) 2D matrix, in parallel.
det2P :: (Monad m)
      => Array U DIM2 Double
      -> m Double
det2P arr = arr `deepSeqArray` do
      let luarr = fst $ lu arr
          ds = mdiag luarr
      x <- R.foldP (*) 1 ds
      return $ x ! Z
{-# NOINLINE det2P #-}


-- | Get the determinant of a (square) 2D matrix sequentially
det2S :: Array U DIM2 Double
      -> Double
det2S arr = arr `deepSeqArray` runST (do 
         let luarr = fst $ lu arr
             ds = mdiag luarr
             x = R.foldS (*) 1 ds
         return $ x ! Z)
{-# NOINLINE det2S #-}


-- LU factorization by Crout algorithms (L(i,i) = 1) --------------------------
-- | LU factorization, in parallel
luP :: (Monad m)
    => Array U DIM2 Double
    -> m (Array U DIM2 Double,
          Array U DIM1 (Int, Int))
luP arr = arr `deepSeqArray` do
    let (mat, p) = luS arr
    return (mat, p)
{-# NOINLINE luP #-}


-- | LU factorization, sequentially
luS :: Array U DIM2 Double
    -> (Array U DIM2 Double,
        Array U DIM1 (Int, Int))
luS arr = arr `deepSeqArray` (computeS mat, computeS p)
  where (mat, p) = lu arr
{-# NOINLINE luS #-}


lu :: (Source r Double)
   => Array r DIM2 Double
   -> (Array D DIM2 Double,
       Array D DIM1 (Int, Int))
lu arr = iter (nr-1) ma0 md0 mu0 ml0 [(0,p0)]
  where nr = row (R.extent arr)
        l_0 = R.slice arr (Any :. (0::Int))
        (pivot0, p0) = V.maximumBy (comparing (abs . fst))
                       $ V.zip (R.toUnboxed $ R.computeS l_0)
                               (V.fromList [0..(nr-1)])
        ma0 = swapRow (0,p0) arr
        md0 = R.extract (Z :. 0 :. 0) (Z :. 1 :. 1) ma0
        ml0 = R.map (/ pivot0)
              $ R.extract (Z :. 1 :. 0) (Z :. nr-1 :. 1) ma0
        mu0 = R.extract (Z :. 0 :. 1) (Z :. 1 :. nr-1) ma0
        iter :: (Source r1 Double, Source r2 Double, Source r3 Double)
             => Int
             -> Array r1 DIM2 Double
             -> Array D  DIM2 Double
             -> Array r2 DIM2 Double
             -> Array r3 DIM2 Double
             -> [(Int, Int)]
             -> (Array D DIM2 Double,
                 Array D DIM1 (Int, Int))
        iter 0 _  md _  _  swaps = (md,
                                    R.delay
                                    $ R.fromListUnboxed (Z :. nr) swaps)
        iter q ma md mu ml swaps = iter (q-1) ma' md' mu'' ml'' ((k,pj+k):swaps)
          where k = nr - q  -- k = 1..(nr-1) when q = (nr-1)..1
                mu_k = R.slice mu (Any :. (0::Int))
                mu' = R.extract (Z :. 0 :. 1) (Z :. k :. nr-k-1) mu
                l_k = R.fromFunction (Z :. nr-k)
                      $ \(Z :. i) -> ma R.! (Z :. k+i :. k)
                                     - dot mu_k (R.slice ml (Any :. i :. All))
                (pivot, pj) = V.maximumBy (comparing (abs . fst))
                              $ V.zip (R.toUnboxed $ R.computeS l_k)
                                      (V.fromList [0..(nr-k-1)])
                ma' = swapRow (k,pj+k) ma
                mlS = swapRow (0,pj) ml
                ml' = R.extract (Z :. 1 :. 0) (Z :. nr-k-1 :. k) mlS
                mlk_ = R.slice mlS (Any :. (0::Int) :. All)
                l_k' = R.extract (Z :. 1 :. 0) (Z :. nr-k-1 :. 1)
                       $ swapRow (0,pj)
                       $ R.reshape (Z :. nr-k :. 1 :: DIM2)
                       $ R.map (/ pivot) l_k
                ml'' = ml' R.++ l_k'
                uk_ = R.fromFunction (Z :. nr-k-1)
                      $ \(Z :. j) -> ma' R.! (Z :. k :. k+j+1)
                                     - dot mlk_ (R.slice mu' (Any :. j))
                mu'' = mu' -+- (R.reshape (Z :. 1 :. nr-k-1 :: DIM2) uk_)
                md' = (md R.++ (R.reshape (Z :. k :. 1 :: DIM2) mu_k))
                      -+- (R.reshape (Z :. 1 :. k+1 :: DIM2)
                           $ mlk_ R.++ (R.fromListUnboxed (Z :. 1 :: DIM1) [pivot]))
{-# NOINLINE lu #-}


swapRow :: (Source r Double)
        => (Int, Int)
        -> Array r DIM2 Double
        -> Array D DIM2 Double
swapRow (r1, r2) mat
    | r1 < 0 || r1 >= nr || r2 < 0 || r2 >= nr || r1 == r2 = R.delay mat
    | otherwise = R.traverse mat id
                  $ \_ (Z :. r :. c) ->
                       case () of _
                                    | r == b1   -> mat R.! (Z :. b2 :. c)
                                    | r == b2   -> mat R.! (Z :. b1 :. c)
                                    | otherwise -> mat R.! (Z :. r  :. c)
  where b1 = min r1 r2
        b2 = max r1 r2
        nr = row (R.extent mat)
{-# NOINLINE swapRow #-}


-- Diagonal vector ------------------------------------------------------------
-- | Diagonal elements, in parallel
mdiagP :: (Monad m)
       => Array U DIM2 Double
       -> m (Array U DIM1 Double)
mdiagP arr = arr `deepSeqArray` computeP $ mdiag arr
{-# NOINLINE mdiagP #-}


-- | Diagonal elements, sequentially
mdiagS :: Array U DIM2 Double
       -> Array U DIM1 Double
mdiagS arr = arr `deepSeqArray` computeS $ mdiag arr
{-# NOINLINE mdiagS #-}


-- | Diagonal elements, delayed
mdiag :: (Source r Double)
      => Array r DIM2 Double
      -> Array D DIM1 Double
mdiag arr = slice md (Any :. (0::Int))
  where md = R.traverse arr f g
        f :: DIM2 -> DIM2
        f (Z :. r :. c) = Z :. (min r c) :. (1::Int)
        g :: (DIM2 -> a) -> DIM2 -> a
        g get (Z :. r :. _) = get (Z :. r :. r)
{-# NOINLINE mdiag #-}


-- Diagonal matrix ------------------------------------------------------------
-- | Diagonal matrix, in parallel
vdiagP :: (Monad m)
       => Array U DIM1 Double
       -> m(Array U DIM2 Double)
vdiagP vec = vec `deepSeqArray` computeP $ vdiag vec
{-# NOINLINE vdiagP #-}


-- | Diagonal matrix, sequentially
vdiagS :: Array U DIM1 Double
       -> Array U DIM2 Double
vdiagS vec = vec `deepSeqArray` computeS $ vdiag vec
{-# NOINLINE vdiagS #-}


-- | Diagonal matrix, delayed
vdiag :: (Source r Double)
      => Array r DIM1 Double
      -> Array D DIM2 Double
vdiag vec = R.fromFunction (Z :. n :. n)
            (\(Z :. r :. c) -> if r == c
                               then vec R.! (Z :. r)
                               else 0.0)
  where n = R.size $ R.extent vec
{-# NOINLINE vdiag #-}


-- Solve linear equations by LU factorization ---------------------------------
-- | Solve linear equation, in parallel
solveLUP :: (Monad m)
         => Array U DIM2 Double
         -> Array U DIM1 Double
         -> m (Array U DIM1 Double)
solveLUP arr vec = arr `deepSeqArray` vec `deepSeqArray`
                   computeP $ solveLU arr vec
{-# NOINLINE solveLUP #-}


-- | Solve linear equation, sequentially
solveLUS :: Array U DIM2 Double
         -> Array U DIM1 Double
         -> Array U DIM1 Double
solveLUS arr vec = arr `deepSeqArray` vec `deepSeqArray`
                   computeS $ solveLU arr vec
{-# NOINLINE solveLUS #-}


-- | Solve linear equation, delayed
solveLU :: forall r1 r2 .
           (Source r1 Double, Source r2 Double)
        => Array r1 DIM2 Double
        -> Array r2 DIM1 Double
        -> Array D  DIM1 Double
solveLU arr vec = xvec
  where (luarr, parr) = lu arr
        larr = R.fromFunction (R.extent luarr)
               $ \(Z :. r :. c) -> if r > c
                                   then luarr R.! (Z :. r :. c)
                                   else if r == c then 1.0 else 0.0
        uarr = R.fromFunction (R.extent luarr)
               $ \(Z :. r :. c) -> if r <= c
                                   then luarr R.! (Z :. r :. c)
                                   else 0.0
        pvec = R.toUnboxed $ R.computeS parr
        vec' = R.reshape e
               $ V.foldr swapRow (R.reshape (e :. 1) vec) pvec
        e = R.extent vec
        n = size e - 1
        fy :: Array D DIM2 Double
           -> Int
           -> Double
        fy _ 0    = vec' ! (Z :. 0)
        fy la i = vec' ! (Z :. i) - dot ys ls
          where ys = R.fromFunction (Z :. i) (\(Z :. k) -> fy la k)
                ls = R.slice (R.extract (Z :. i :. 0) (Z :. 1 :. i) la)
                             (Any :. (0::Int) :. All)
        fx :: Array D DIM2 Double
           -> Array D DIM1 Double
           -> Int
           -> Double
        fx ua y 0 = y ! (Z :. n) / ua ! (Z :. n :. n)
        fx ua y i = (y ! (Z :. n-i) - dot xs us)
                         / ua ! (Z :. n-i :. n-i)
          where xs = R.fromListUnboxed (Z :. i)
                     $ L.map (fx ua y) $ reverse [0..(i-1)]
                us = R.slice (R.extract (Z :. n-i :. n-i+1) (Z :. 1 :. i) ua)
                             (Any :. (0::Int) :. All)
        yvec = R.fromFunction e (\(Z :. i) -> fy larr i)
        xvec = R.fromFunction e (\(Z :. i) -> fx uarr yvec (n-i))
{-# NOINLINE solveLU #-}
