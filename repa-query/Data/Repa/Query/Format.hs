
-- | Conversion between the open type-index representation of formats,
--   and closed data-type view. 
module Data.Repa.Query.Format
        ( Row   (..)
        , Delim (..)
        , Field (..)

        , FieldBox (..)
        , flattens
        , flattensBox

        , showField
        , readField)
where
import Data.Char
import Data.Word
import Data.Int
import Data.Repa.Bits.Date32                    as Date32
import qualified Data.Map                       as Map
import qualified Data.Repa.Product              as P


-- | Row format.
data Row aa
        = Row Delim (Field aa)
        deriving (Eq, Show)


-- | How the rows and fields are delimited.
data Delim
        -- | Format with fixed-length rows.
        --   The fields are all concatenated together, with no delimitors.
        = Fixed

        -- | Format with a single field on each lines.
        | Lines

        -- | Format with multiple fields on each line,
        --   where the fields are separated by a special character.
        | LinesSep Char
        deriving (Eq, Show)


-- | Field formats supported by Repa tables.
--
--   The repa-convert library defines several data formats using a singleton
--   type to represent each format. We enumerate all the formats here in a
--   closed data type to make it easier to write parsers and pretty pritners
--   for queries that mention them.
--
--   The names of the constructors have the same names as the ones in 
--   repa-convert, so when we read/show them they match.
--
data Field a where

        -- Compound field.
        (:*:)           :: Field x -> Field y -> Field (x P.:*: y)

        -- Big-endian 8-bit unsigned word.
        Word8be         :: Field Word8

        -- Big-endian 8-bit signed integer.
        Int8be          :: Field Int8

        -- Big-endian 16-bit unsigned word.
        Word16be        :: Field Word16

        -- Big-endian 16-bit signed integer.
        Int16be         :: Field Int16

        -- Big-endian 32-bit unsigned word.
        Word32be        :: Field Word32

        -- Big-endian 32-bit signed integer.
        Int32be         :: Field Int32

        -- Big-endian 64-bit unsigned word.
        Word64be        :: Field Word64

        -- Big-endian 64-bit signed integer.
        Int64be         :: Field Int64

        -- Big-endian 32-bit IEEE 754 float.
        Float32be       :: Field Float

        -- Big-endian 64-bit IEEE 754 float.
        Float64be       :: Field Double

        -- Date in ASCII YYYYsMMsDD format.
        YYYYsMMsDD      :: Char -> Field Date32

        -- Date in ASCII DDsMMsYYYY format.
        DDsMMsYYYY      :: Char -> Field Date32

        -- Human readable ASCII integer.
        IntAsc          :: Field Int  

        -- Human readable ASCII double.
        DoubleAsc       :: Field Double

        -- Fixed length ASCII string.
        FixAsc          :: Int -> Field String

        -- Variable length ASCII string.
        VarAsc          :: Field String


deriving instance Show (Field a)
deriving instance Eq   (Field a)

infixr :*:


---------------------------------------------------------------------------------------------------
-- | Existential container for field formats,
--   and dictionaries to work with them.
data FieldBox
        =  forall a
        .  FieldBox (Field a)

instance Show FieldBox where
 show (FieldBox f)      = show f


-- | Flatten compound fields into their parts and box up the components.
flattens :: Field a -> [FieldBox]
flattens ff
 = case ff of
        (:*:) f1 f2     -> flattens f1 ++ flattens f2
        _               -> [FieldBox ff]


-- | Like `flattens`, but start with a boxed field format.
flattensBox :: FieldBox -> [FieldBox]
flattensBox (FieldBox f)
        = flattens f


---------------------------------------------------------------------------------------------------
-- | Show a field format.
showField :: Field a -> String
showField ff = show ff


-- | Parse a field format.
readField :: String -> Maybe FieldBox
readField ss
 | Just f       <- Map.lookup ss atomic
 = Just f

 | ["YYYYsMMsDD", sep]  <- words ss
 , ['\'', c, '\'']      <- sep
 = Just $ FieldBox $ YYYYsMMsDD c

 | ["DDsMMsYYYY", sep]  <- words ss
 , ['\'', c, '\'']      <- sep
 = Just $ FieldBox $ DDsMMsYYYY c

 | ["FixAsc", sLen]     <- words ss
 , all isDigit sLen
 = Just $ FieldBox $ FixAsc (read sLen)

 | otherwise
 = Nothing

 where
        atomic
         = Map.fromList
         [ ("Word8be",    FieldBox   Word8be)
         , ("Int8be",     FieldBox    Int8be) 
         , ("Word16be",   FieldBox  Word16be)
         , ("Int16be",    FieldBox   Int16be)
         , ("Word32be",   FieldBox  Word32be)
         , ("Int32be",    FieldBox   Int32be)
         , ("Word64be",   FieldBox  Word64be)
         , ("Int64be",    FieldBox   Int64be)
         , ("Float32be",  FieldBox Float32be)
         , ("Float64be",  FieldBox Float64be)
         , ("IntAsc",     FieldBox    IntAsc)
         , ("DoubleAsc",  FieldBox DoubleAsc)
         , ("VarAsc",     FieldBox    VarAsc) ]


