
-- | Drivers used by queries compiled via Repa.
--
--   This code is imported by the generated query code and used by the 
--   running queries, rather than being used at query compile time.
--
module Data.Repa.Query.Runtime.Driver
        (streamSourcesToStdout)
where
import Data.Repa.Flow                                   as F
import Data.Word
import System.IO
import qualified Data.Repa.Flow.Generic                 as FG
import qualified Data.Repa.Array.Generic                as AG
import qualified Data.Repa.Array.Material.Foreign       as AF
import qualified Data.Repa.Array.Material.Auto          as AA
import qualified Foreign.Ptr                            as Foreign
import qualified Foreign.ForeignPtr                     as Foreign


-- | Read data from a bundle of sources and write it to stdout.
--  
--   This function only works for sources bundles containing a single stream.
--   If this is not true then the function returns False, and no data is
--   written to stdout.
--
streamSourcesToStdout
        :: Sources Word8 -> IO Bool

streamSourcesToStdout (FG.Sources 1 pullX)
 = do   go
        return True

 where  go 
         = pullX 0 eat_streamSource eject_streamSource

        eat_streamSource (chunk :: AG.Array AA.A Word8)
         = do   let (start, len, fptr :: Foreign.ForeignPtr Word8) 
                        = AF.toForeignPtr $ AG.convert AF.F chunk

                Foreign.withForeignPtr fptr $ \ptr 
                 -> hPutBuf stdout (ptr `Foreign.plusPtr` start) len

                hFlush stdout
                go

        eject_streamSource
         = do   hClose stdout
                return ()

streamSourcesToStdout _
 = return False
