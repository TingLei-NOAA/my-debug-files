I created review_profile_level_range.py for this.

Example:

python review_profile_level_range.py \
  --profiles-dir dr-normprofiles-dbz-1G-2var \
  --k 1 \
  --i-start 1 \
  --i-end 5 \
  --j-start 1229 \
  --j-end 1469 \
  --top-n 10
It prints:

how many profile files fall inside that i/j window
the minimum and maximum values at level k
the corresponding subdomain_index, i, j, and file
the smallest/largest top-n entries
So yes, the script explicitly shows how i,j are recovered: from the second line written by sample_profiles_from_fv3.py
