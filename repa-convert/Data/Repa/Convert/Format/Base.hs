
module Data.Repa.Convert.Format.Base
        ( Format   (..)
        , Packable (..)
        , Packables(..)

        -- * Conversion
        , packToList,   unpackFromList
        , packToString, unpackFromString

        -- * Constraints
        , forFormat
        , listFormat

        -- * Strict products
        , (:*:)(..))
where
import Data.Repa.Product
import Data.Char
import Data.Word
import Control.Monad
import System.IO.Unsafe
import qualified Foreign.Storable               as S
import qualified Foreign.Marshal.Alloc          as S
import qualified Foreign.Ptr                    as S


---------------------------------------------------------------------------------------------------
-- | Relates a storage format to the Haskell type of the value
--   that is stored in that format.
class Format f where

 -- | Get the type of a value with this format.
 type Value f  

 -- | Yield the minumum number of bytes that a value of this
 --   format will take up. 
 -- 
 --   Packing a value into this format
 --   is guaranteed to use at least this many bytes.
 --   This is exact for fixed-size formats.
 minSize    :: f -> Int


 -- | For fixed size formats, yield their size (length) in bytes.
 --
 --   Yields `Nothing` if this is not a fixed size format.
 --
 fixedSize  :: f -> Maybe Int


 -- | Yield the size of a value in the given format.
 --
 --   Yields `Nothing` when a collection of values is to be packed into a
 --   fixed length format, but the size of the collection does not match
 --   the format.
 --
 --   If `fixedSize` returns a size then `packedSize` returns the same size.
 --
 packedSize :: f -> Value f -> Maybe Int


 -- | Yield the number of separate fields in this format.
 fieldCount :: f -> Maybe Int

  
---------------------------------------------------------------------------------------------------
-- | Class of storage formats that can have values packed and unpacked
--   from foreign bufferes. 
-- 
--   The methods are written using continuations to make it easier for
--   GHC to optimise its core code when packing/unpacking many fields.
--
class Format   format 
   => Packable format where


 -- | Pack a value into a buffer using the given format.
 -- 
 --   The buffer must be at least as long as the size returned by
 --   `fixedSize` / `packedSize`. 
 -- 
 --   If the format contains fixed width fields and the corresponding
 --   value has too many elements, then this function returns `False`, 
 --   otherwise `True`.
 --
 pack   :: S.Ptr Word8                  -- ^ Target Buffer.
        -> format                       -- ^ Storage format.
        -> Value format                 -- ^ Value to pack.
        -> (Int -> IO (Maybe a))        -- ^ Continue, given the number of bytes written.
        -> IO (Maybe a)


 -- | Unpack a value from a buffer using the given format.
 --
 --   This is the inverse of `pack` above.
 -- 
 --   PRECONDITION: The length of the buffer must be at least the
 --   minimum size required of the format (minSize). This allows
 --   us to avoid repeatedly checking for buffer overrun when
 --   unpacking fixed size formats.
 --
 unpack :: S.Ptr Word8                  -- ^ Source buffer.
        -> Int                          -- ^ Length of buffer.
        -> format                       -- ^ Format of buffer.
        -> ((Value format, Int) -> IO (Maybe a)) 
                                        -- ^ Continue, given the unpacked value and the 
                                        --   number of bytes read. 
        -> IO (Maybe a)


-- | Class of field containers, eg comma or pipe-separated fields.
class Format format 
   => Packables container format where

 packs  :: S.Ptr Word8                  -- ^ Target Buffer.
        -> container                    -- ^ Field container.
        -> format                       -- ^ Storage format.
        -> Value format                 -- ^ Value to pack.
        -> (Int -> IO (Maybe a))        -- ^ Continue, given the number of bytes written.
        -> IO (Maybe a)

 unpacks
        :: S.Ptr Word8                  -- ^ Target Buffer.
        -> Int                          -- ^ Length of buffer.
        -> container                    -- ^ Field container.
        -> format                       -- ^ Storage format.
        -> ((Value format, Int) -> IO (Maybe a))        
                                        -- ^ Continue, given the number of bytes written.
        -> IO (Maybe a)


---------------------------------------------------------------------------------------------------
-- | Pack a value into a list of `Word8`.
packToList 
        :: Packable format
        => format -> Value format -> Maybe [Word8]
packToList f x
 | Just lenMax  <- packedSize f x
 = unsafePerformIO
 $ do   buf     <- S.mallocBytes lenMax
        mResult <- pack buf f x (\o -> return (Just o))
        case mResult of
         Nothing      -> return Nothing
         Just lenUsed -> do
                xs      <- mapM (S.peekByteOff buf) [0 .. lenUsed - 1]
                S.free buf
                return $ Just xs

 | otherwise    = Nothing


-- | Unpack a value from a list of `Word8`.
unpackFromList
        :: Packable format
        => format -> [Word8] -> Maybe (Value format)

unpackFromList f xs
 = unsafePerformIO
 $ do   let len = length xs
        buf     <- S.mallocBytes len
        mapM_ (\(o, x) -> S.pokeByteOff buf o x)
                $ zip [0 .. len - 1] xs
        unpack buf len f $ \(v, _) -> return (Just v)


-- | Pack a value to a String.
packToString
        :: Packable format
        => format -> Value format -> Maybe String
packToString f v
        = liftM (map (chr . fromIntegral)) $ packToList f v


-- | Unpack a value from a String.
unpackFromString 
        :: Packable format
        => format -> String -> Maybe (Value format)
unpackFromString f s
        = unpackFromList f $ map (fromIntegral . ord) s


---------------------------------------------------------------------------------------------------
-- | Constrain the type of a value to match the given format.
-- 
--   The value itself is not used.
--
forFormat :: format -> Value format  -> Value format
forFormat _ v = v
{-# INLINE forFormat #-}


-- | Constrain the type of some values to match the given format.
--
--   The value itself is not used.
--
listFormat :: format -> [Value format] -> [Value format]
listFormat _ v = v
{-# INLINE listFormat #-}

