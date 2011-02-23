{-# OPTIONS_HADDOCK hide #-}
{-# LANGUAGE TypeOperators, ExplicitForAll, FlexibleContexts #-}

module Data.Array.Repa.Operators.IndexSpace
	( reshape
	, append, (+:+)
	, transpose
	, replicate
	, slice
	, backpermute
	, backpermuteDft)
where
import Data.Array.Repa.Index
import Data.Array.Repa.Slice
import Data.Array.Repa.Internals.Elt
import Data.Array.Repa.Internals.Base
import Data.Array.Repa.Operators.Traverse
import Data.Array.Repa.Shape		as S
import Prelude				hiding (replicate)

stage	= "Data.Array.Repa.Operators.IndexSpace"

-- Index space transformations --------------------------------------------------------------------
-- | Impose a new shape on the elements of an array.
--   The new extent must be the same size as the original, else `error`.
--   Reshaping a `DelayedCursor` array makes it `Delayed`.
--   TODO: Finish this for Partitioned arrays.
reshape	:: (Shape sh, Shape sh', Elt a) 
	=> sh'
	-> Array sh a
	-> Array sh' a

{-# INLINE reshape #-}
reshape sh' arr
	| not $ S.size sh' == S.size (extent arr)
	= error $ stage ++ ".reshape: reshaped array will not match size of the original"

	| otherwise
	= case arr of
		Manifest _  uarr 
		 -> Manifest sh' uarr

		Delayed  sh f
		 -> Delayed sh' (f . fromIndex sh . toIndex sh')

		DelayedCursor sh makeCursor _ loadElem
		 -> Delayed sh' (loadElem . makeCursor . fromIndex sh . toIndex sh')

		-- TODO: we need to reshape all of the contained functions.
		Partitioned{}
		 -> error $ stage ++ ".reshape: not finished for segmented"


-- | Append two arrays.
--
--   OBLIGATION: The higher dimensions of both arrays must have the same extent.
--
--   @tail (listOfShape (shape arr1)) == tail (listOfShape (shape arr2))@
--
append, (+:+)	
	:: (Shape sh, Elt a)
	=> Array (sh :. Int) a
	-> Array (sh :. Int) a
	-> Array (sh :. Int) a

{-# INLINE append #-}
append arr1 arr2 
 = unsafeTraverse2 arr1 arr2 fnExtent fnElem
 where
 	(_ :. n) 	= extent arr1

	fnExtent (sh :. i) (_  :. j) 
		= sh :. (i + j)

	fnElem f1 f2 (sh :. i)
      		| i < n		= f1 (sh :. i)
  		| otherwise	= f2 (sh :. (i - n))

{-# INLINE (+:+) #-}
(+:+) arr1 arr2 = append arr1 arr2


-- | Transpose the lowest two dimensions of an array. 
--	Transposing an array twice yields the original.
transpose 
	:: (Shape sh, Elt a) 
	=> Array (sh :. Int :. Int) a
	-> Array (sh :. Int :. Int) a

{-# INLINE transpose #-}
transpose arr 
 = unsafeTraverse arr
	(\(sh :. m :. n) 	-> (sh :. n :.m))
	(\f -> \(sh :. i :. j) 	-> f (sh :. j :. i))


-- | Replicate an array, according to a given slice specification.
replicate
	:: ( Slice sl
	   , Shape (FullShape sl)
	   , Shape (SliceShape sl)
	   , Elt e)
	=> sl
	-> Array (SliceShape sl) e
	-> Array (FullShape sl) e

{-# INLINE replicate #-}
replicate sl arr
	= backpermute 
		(fullOfSlice sl (extent arr)) 
		(sliceOfFull sl)
		arr

-- | Take a slice from an array, according to a given specification.
slice	:: ( Slice sl
	   , Shape (FullShape sl)
	   , Shape (SliceShape sl)
	   , Elt e)
	=> Array (FullShape sl) e
	-> sl
	-> Array (SliceShape sl) e

{-# INLINE slice #-}
slice arr sl
	= backpermute 
		(sliceOfFull sl (extent arr))
		(fullOfSlice sl)
		arr


-- | Backwards permutation of an array's elements.
--	The result array has the same extent as the original.
backpermute
	:: forall sh sh' a
	.  (Shape sh, Shape sh', Elt a) 
	=> sh' 				-- ^ Extent of result array.
	-> (sh' -> sh) 			-- ^ Function mapping each index in the result array
					--	to an index of the source array.
	-> Array sh a 			-- ^ Source array.
	-> Array sh' a

{-# INLINE backpermute #-}
backpermute newExtent perm arr
	= traverse arr (const newExtent) (. perm) 
	

-- | Default backwards permutation of an array's elements.
--	If the function returns `Nothing` then the value at that index is taken
--	from the default array (@arrDft@)
backpermuteDft
	:: forall sh sh' a
	.  (Shape sh, Shape sh', Elt a) 
	=> Array sh' a			-- ^ Default values (@arrDft@)
	-> (sh' -> Maybe sh) 		-- ^ Function mapping each index in the result array
					--	to an index in the source array.
	-> Array sh  a			-- ^ Source array.
	-> Array sh' a

{-# INLINE backpermuteDft #-}
backpermuteDft arrDft fnIndex arrSrc
	= Delayed (extent arrDft) fnElem
	where	fnElem ix	
		 = case fnIndex ix of
			Just ix'	-> arrSrc ! ix'
			Nothing		-> arrDft ! ix
