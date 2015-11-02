module Paths_repa_algorithms (
    version,
    getBinDir, getLibDir, getDataDir, getLibexecDir,
    getDataFileName, getSysconfDir
  ) where

import qualified Control.Exception as Exception
import Data.Version (Version(..))
import System.Environment (getEnv)
import Prelude

catchIO :: IO a -> (Exception.IOException -> IO a) -> IO a
catchIO = Exception.catch

version :: Version
version = Version [3,4,0,2] []
bindir, libdir, datadir, libexecdir, sysconfdir :: FilePath

bindir     = "/Users/kiwamu/.cabal/bin"
libdir     = "/Users/kiwamu/.cabal/lib/x86_64-osx-ghc-7.10.2/repa-algorithms-3.4.0.2-65pcTUu7L9n9HROwhdY1Qj"
datadir    = "/Users/kiwamu/.cabal/share/x86_64-osx-ghc-7.10.2/repa-algorithms-3.4.0.2"
libexecdir = "/Users/kiwamu/.cabal/libexec"
sysconfdir = "/Users/kiwamu/.cabal/etc"

getBinDir, getLibDir, getDataDir, getLibexecDir, getSysconfDir :: IO FilePath
getBinDir = catchIO (getEnv "repa_algorithms_bindir") (\_ -> return bindir)
getLibDir = catchIO (getEnv "repa_algorithms_libdir") (\_ -> return libdir)
getDataDir = catchIO (getEnv "repa_algorithms_datadir") (\_ -> return datadir)
getLibexecDir = catchIO (getEnv "repa_algorithms_libexecdir") (\_ -> return libexecdir)
getSysconfDir = catchIO (getEnv "repa_algorithms_sysconfdir") (\_ -> return sysconfdir)

getDataFileName :: FilePath -> IO FilePath
getDataFileName name = do
  dir <- getDataDir
  return (dir ++ "/" ++ name)
