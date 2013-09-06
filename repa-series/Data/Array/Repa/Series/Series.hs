
module Data.Array.Repa.Series.Series
        ( Series (..)
        , index
        , indexFloatX4
        , indexDoubleX2
        , length
        , toVector
        , rateOfSeries
        , down4
        , tail4
        , unsafeFromVector)
where
import Data.Array.Repa.Series.Rate
import Data.Array.Repa.Series.Vector            (Vector)
import Data.Vector.Primitive                    (Prim)
import Data.Primitive.ByteArray
import System.IO.Unsafe
import GHC.Exts
import GHC.Prim
import qualified Data.Array.Repa.Series.Vector  as V
import qualified Data.Vector.Primitive          as P
import qualified Data.Vector.Primitive.Mutable  as PM
import Prelude hiding (length)
import Debug.Trace

-- | A `Series` is a source of element data that is tagged by rate variable,
--   which is a type level version of its length.
--
--   Although the manifest representation of a series supports random-access
--   indexing, all fusable series process must consume their series sequentially.
--
--   The rate parameter @k@ represents the abstract length of the series.
--
data Series k a
        = Series 
        { seriesStart           :: Word#
        , seriesLength          :: Word#
        , seriesMBA_start       :: Word#
        , seriesMBA             :: MutableByteArray# RealWorld
        , seriesVector          :: !(P.Vector a) }  

-- | Take the length of a series.
length :: Series k a -> Word#
length (Series start len _ _ d) = len
{-# INLINE [1] length #-}


-- | Get the Rate / Length of a series.
rateOfSeries :: Series k a -> RateNat k
rateOfSeries s
 = RateNat (seriesLength s)
{-# INLINE [1] rateOfSeries #-}


-- | Window a series to the initial range of 4 elements.
down4 :: forall k a. RateNat (Down4 k) -> Series k a -> Series (Down4 k) a
down4 r (Series start len mba_start mba vec)        
       = Series start len mba_start mba vec
{-# INLINE [1] down4 #-}


-- | Window a series to the ending elements.
tail4 :: forall k a. RateNat (Tail4 k) -> Series k a -> Series (Tail4 k) a
tail4 r (Series start len mba_start mba vec)        
        = Series (quotWord# len (int2Word# 4#) `timesWord#` (int2Word# 4#))
                 len mba_start mba vec
{-# INLINE [1] tail4 #-}


-- | Index into a series.
index :: Prim a => Series k a -> Word# -> a
index s ix
 = P.unsafeIndex (seriesVector s) 
                 (I# (word2Int# (ix `plusWord#` seriesStart s)))
{-# INLINE [1] index #-}


-- | Retrieve a packed FloatX4 from a `Series`.
indexFloatX4  :: Series (Down4 k) Float -> Word# -> FloatX4#
indexFloatX4 s ix
 = let  !mba            = seriesMBA s
        !offset         = word2Int#
                        ( plusWord# (seriesMBA_start s)
                        ( plusWord# (seriesStart s)
                                    (timesWord# ix (int2Word# 4#))))

        (# _, f4 #)     = readFloatX4Array# mba offset realWorld#
   in   f4
{-# INLINE [1] indexFloatX4 #-}


-- | Retrieve a packed DoubleX2 from a `Series`.
indexDoubleX2 :: Series (Down2 k) Float -> Word# -> DoubleX2#
indexDoubleX2 s ix
        = doubleToDoubleX2# (int2Double# 5#)                    -- TODO: fixme
{-# INLINE [1] indexDoubleX2 #-}


-- | Convert a series to a vector, discarding the rate information.
toVector :: Prim a => Series k a -> Vector a
toVector s
 = unsafePerformIO
 $ do   V.fromPrimitive (seriesVector s)
{-# INLINE [1] toVector #-}


-- | Unsafely convert a vector to a series of an arbitrary rate.
--   
--   The rate variable in the result is arbitrary,
--   so a series created this way may not have the same length as others
--   of the same rate.
unsafeFromVector :: Prim a => Vector a -> IO (Series k a)
unsafeFromVector (V.Vector len mv)
 = do   let !pv@(P.MVector (I# start) (I# _) (MutableByteArray mba))
                = mv
        v       <- P.unsafeFreeze mv
        return $ Series (int2Word# 0#) len (int2Word# start) mba v
{-# NOINLINE unsafeFromVector #-}


