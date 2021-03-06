{-# LANGUAGE NoMonomorphismRestriction #-}

-- | Primitive functions used by produced Repa code when compiling the queries.
--   All names used by the emitted query code are defined here so we can keep
--   track of what is being used.
--
module Data.Repa.Query.Runtime.Primitive
        ( -- * From Prelude
          (>>=), (=<<), return
        , error
        , negate, abs, signum
        , add, sub, mul, div
        , eq,  neq
        , gt,  ge,  lt,  le

          -- * From Data.Repa.Flow.Auto
        , map_i
        , folds_i
        , groupsBy_i

          -- * From Data.Repa.Flow.Auto.IO
        , fromFiles
        , sourceLinesFormat
        , sourceFixedFormat

          -- * From Data.Repa.Flow.Auto.Format
        , concatPackFormat_i
        , unlinesPackFormat_i

          -- * From Data.Repa.Convert.Format
        , pattern App
        , pattern Sep
        , pattern Word8be,      pattern Int8be
        , pattern Word16be,     pattern Int16be
        , pattern Word32be,     pattern Int32be
        , pattern Word64be,     pattern Int64be
        , pattern Float32be
        , pattern Float64be
        , pattern YYYYsMMsDD
        , pattern DDsMMsYYYY
        , pattern IntAsc
        , pattern DoubleAsc
        , pattern FixAsc
        , pattern VarAsc)
where
import qualified Prelude                        as P
import qualified Data.Repa.Flow.Auto            as F
import qualified Data.Repa.Flow.Auto.IO         as F
import qualified Data.Repa.Flow.Auto.Format     as F
import qualified Data.Repa.Convert.Format       as C


-- Prelude
(>>=)                   = (\x y -> x P.>>= y)
(=<<)                   = (\x y -> x P.=<< y)
return                  = P.return
error                   = P.error

negate                  = (\x   -> P.negate x)
abs                     = (\x   -> P.abs    x)
signum                  = (\x   -> P.signum x)
add                     = (\x y -> x P.+  y)
sub                     = (\x y -> x P.-  y)
mul                     = (\x y -> x P.*  y)
div                     = (\x y -> x P./  y)
eq                      = (\x y -> x P.== y)
neq                     = (\x y -> x P./= y)
gt                      = (\x y -> x P.>  y)
ge                      = (\x y -> x P.>= y)
lt                      = (\x y -> x P.<  y)
le                      = (\x y -> x P.<= y)


-- Data.Repa.Flow.Auto
map_i                   = F.map_i
folds_i                 = F.folds_i
groupsBy_i              = F.groupsBy_i


-- Data.Repa.Flow.Auto.IO
fromFiles               = F.fromFiles
sourceLinesFormat       = F.sourceLinesFormat
sourceFixedFormat       = F.sourceFixedFormat


-- Data.Repa.Flow.Auto.Format
concatPackFormat_i      = F.concatPackFormat_i
unlinesPackFormat_i     = F.unlinesPackFormat_i


-- Data.Repa.Convert.Format
pattern App   fs        = C.App fs
pattern Sep s fs        = C.Sep s fs

pattern Word8be         = C.Word8be
pattern  Int8be         =  C.Int8be

pattern Word16be        = C.Word16be
pattern  Int16be        =  C.Int16be

pattern Word32be        = C.Word32be
pattern  Int32be        =  C.Int32be

pattern Word64be        = C.Word64be
pattern  Int64be        =  C.Int64be

pattern Float32be       = C.Float32be
pattern Float64be       = C.Float64be

pattern YYYYsMMsDD c    = C.YYYYsMMsDD c
pattern DDsMMsYYYY c    = C.DDsMMsYYYY c

pattern IntAsc          = C.IntAsc
pattern DoubleAsc       = C.DoubleAsc

pattern FixAsc len      = C.FixAsc len
pattern VarAsc          = C.VarAsc

